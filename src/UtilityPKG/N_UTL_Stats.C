//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Math.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <limits>

#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_FormatTime.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace Stats {

Stat StatTop::s_statTop(0);

namespace {

MetricsMask s_enabledMetricsMask = METRICS_LAP_COUNT | METRICS_CPU_TIME | METRICS_WALL_TIME;        ///< Bit mask of enabled metrics

template <class T>
typename MetricTraits<T>::Type
value_now() {
  if (MetricTraits<T>::METRIC & getEnabledStatMetricsMask())
    return MetricTraits<T>::value_now();
  else
    return 0;
}

std::vector<std::string> &
split(
  const std::string &           path,
  char                          separator,
  std::vector<std::string> &    path_vector)
{
  for (std::string::const_iterator it = path.begin(); ; ) {
    std::string::const_iterator it2 = std::find(it, path.end(), separator);
    path_vector.push_back(std::string(it, it2));
    if (it2 == path.end())
      break;
    it = it2 + 1;
  }

  return path_vector;
}

} // namespace <empty>


MetricsMask getEnabledStatMetricsMask()
{
  return s_enabledMetricsMask;
}

void
setEnabledStatMetricsMask(
  MetricsMask   stat_mask)
{
  s_enabledMetricsMask = stat_mask | METRICS_LAP_COUNT;
}


///
/// Class <b>StatImpl</b> is the core stat class.  The Stat class is a
/// wrapper around StatImpl so that the buried references can be constructed more easily.
///
/// Each stat has a lap counter, cpu stat, wall stat and other metrics.  Each time a stat is
/// started, the cpu start time, wall start time and other metrics, set to the process' current
/// values.  When the stat is stopped, the lap counter is incremented, and the cpu, wall, and other
/// values are accumulated with the difference between now and the start time.
///
/// Each stat may have a list of subordinate stats.  The relationship is purely
/// hierarchical in that a there is no timing relationship assumed between the stats other
/// than the grouping.  There is no relation between the starting and stopping of parent
/// and subordinate stats.
///
/// The subordinate stats are stored as pointers to a new stat on the heap, since the
/// calling function will be receiving a reference to this memory which can never change
/// location.  The subordinate stats are not sorted in the list as they should very
/// rarely be created or looked up by name, rather the calling function stores the
/// reference via the Stat class.
///
////
class StatImpl
{
  friend class Stat;
  friend class StatTop;

public:
  static void updateRootStat(StatImpl *root_stat);

  static Stat createRootStat(const std::string &name, const StatSet &stat_set);

  static void deleteRootStat(StatImpl *root_stat);

  static std::vector<Stat> &findStats(StatImpl *root_stat, const std::string &path_tail, std::vector<Stat> &found_stats);

  static void findStat(StatImpl *stat, std::vector<std::string> &path_tail_vector, std::vector<Stat> &found_stats);

private:
  ///
  /// Static function <b>reg</b> returns a reference to an existing stat or newly
  /// created stat of the specified <b>name</b> which is subordinate to the
  /// <b>parent</b> stat.
  ///
  /// @return      a <b>StatImpl</b> reference to the stat with the
  ///        specified name that is subordinate to the
  ///        <b>parent</b> stat.
  ///
  static StatImpl *reg(const std::string &name, StatMask stat_mask, StatImpl *parent_stat, const StatSet &stat_set) {
    return parent_stat->addSubstat(name, stat_mask, stat_set);
  }

  ///
  /// Creates a new <b>Stat</b> instance.
  ///
  /// @param name    a <b>std::string</b> const reference to the name of
  ///        the stat.
  ///
  ///
  StatImpl(const std::string &name, StatMask stat_mask, StatImpl *parent_stat, const StatSet &stat_set);

  ///
  /// Destroys a <b>StatImpl</b> instance.
  ///
  ///
  ~StatImpl();

  StatImpl(const StatImpl &StatImpl);
  StatImpl &operator=(const StatImpl &StatImpl);

  ///
  /// Class <b>finder</b> is a binary predicate for finding a subordinate stat.
  ///
  /// Note that the subordinate stat is an unsorted list as there are very few stats
  /// created and should rarely be looked up by name.
  ///
  class finder : private std::unary_function<Stat, bool>
  {
  public:
    explicit finder(const std::string &name)
      : m_name(name)
    {}

    bool operator()(Stat stat) const {
      return equal_nocase(stat.getName(), m_name);
    }

  private:
    std::string    m_name;
  };

public:
  ///
  /// Member function <b>getName</b> returns the name of the stat.
  ///
  /// @return      a <b>std::string</b> const reference to the stat's
  ///        name.
  ///
  const std::string &getName() const {
    return m_name;
  }

  ///
  /// Member function <b>getStatMask</b> returns the stat mask of the stat.
  ///
  /// @return      a <b>StatMask</b> value to the stat mask.
  ///
  StatMask getStatMask() const {
    return m_statMask;
  }

  ///
  /// Member function <b>getStatSet</b> returns the stat set of the stat.
  ///
  /// @return      a <b>StatSet</b> const reference to the stat set.
  ///
  const StatSet &getStatSet() const {
    return m_statSet;
  }

  ///
  /// Member function <b>shouldRecord</b> returns true if any of the specified stat
  /// bit masks are set in the enable stat bit mask.
  ///
  /// @param stat_mask    a <b>StatMask</b> value to test the enable stat
  ///        bit mask against.
  ///
  ///
  bool shouldRecord() const {
    return m_statSet.shouldRecord(m_statMask) && s_enabledMetricsMask;
  }

  ///
  /// Member function <b>getSubstatLapCount</b> returns the substat lap counter.
  ///
  /// @return      a <b>Counter</b> value of the substat lap
  ///        counter.
  ///
  double getSubstatLapCount() const {
    return m_substatLapCount;
  }

  void setSubstatLapCount(double value) {
    m_substatLapCount = value;
  }

  ///
  /// Member function <b>getLapCount</b> returns the lap counter metric.  The lap
  /// count metric is the number of times the stop function has been executed.
  ///
  /// @return      a <b>CounterMetric</b> const reference of the lap counter
  ///        metric.
  ///
  template <class T>
  const Stat::Metric<T> &getMetric() const;

  ///
  /// Member function <b>getStatList</b> returns the substats associated with
  /// this stat.
  ///
  /// @return      a <b>StatList</b> const reference to the sub
  ///        time list.
  ///
  const StatList &getStatList() const {
    return m_substatList;
  }

  StatList::iterator begin() {
    return m_substatList.begin();
  }

  StatList::const_iterator begin() const {
    return m_substatList.begin();
  }

  StatList::iterator end() {
    return m_substatList.end();
  }

  StatList::const_iterator end() const {
    return m_substatList.end();
  }

  ///
  /// Member function <b>reset</b> resets the accumulated time and lap times.
  ///
  ///
  void reset();

  ///
  /// Member function <b>checkpoint</b> checkpoints the stat and all substats.
  ///
  ///
  void checkpoint() const;

  ///
  /// Member function <b>start</b> sets the start stat.
  ///
  /// @return      a <b>StatImpl</b> reference to the stat.
  ///
  StatImpl &start();

  ///
  /// Member function <b>lap</b> sets the stop stat.
  ///
  /// @return      a <b>StatImpl</b> reference to the stat.
  ///
  StatImpl &lap();

  ///
  /// Member function <b>stop</b> sets the stop stat and sums the just completed lap
  /// time to the stat.
  ///
  /// @return      a <b>StatImpl</b> reference to the stat.
  ///
  StatImpl &stop();

  ///
  /// Member function <b>accumulateSubstatLapCounts</b> sums the lap counter of all
  /// subordinate stats.  This is used to determin which stats have been activated at all.
  ///
  /// @return      an <b>int</b> value of the number of subordinate
  ///        stat laps.
  ///
  double accumulateSubstatLapCounts() const;

  Stat getSubstat(const std::string &name);

private:
  ///
  /// Member function <b>addSubstat</b> returns a reference to an existing or new
  /// substat with the specified name.
  ///
  /// @param name    a <b>std::string</b> value of the stat's name.
  ///
  /// @param stat_mask    a <b>StatMask</b> value of the class of the stat.
  ///
  /// @return      a <b>StatImpl</b> reference to the stat with
  ///        specified name.
  ///
  StatImpl *addSubstat(const std::string &name, StatMask stat_mask, const StatSet &stat_set);

private:
  std::string                   m_name;               ///< Name of the stat
  StatMask                      m_statMask;           ///< Bit mask to enable stat
  StatImpl *                    m_parentStat;         ///< Parent stat
  mutable double                m_substatLapCount;    ///< Sum of substat lap counts and m_lapCount
  unsigned                      m_lapStartCount;      ///< Number of pending lap stops

  StatList                      m_substatList;        ///< List of subordinate stats

  const StatSet &               m_statSet;            ///< Stat enabled mask
  Stat::Metric<LapCount>        m_lapCount;           ///< Number of laps accumulated
  Stat::Metric<CPUTime>         m_cpuTime;            ///< CPU time
  Stat::Metric<WallTime>        m_wallTime;           ///< Wall time
  Stat::Metric<MPICount>        m_MPICount;           ///< MPI call count
  Stat::Metric<MPIByteCount>    m_MPIByteCount;       ///< MPI byte count
  Stat::Metric<HeapAlloc>       m_heapAlloc;          ///< Heap allocated
  std::string                   m_text;               ///< Simple text data (not accumulated)
};

void
updateRootStat(
  Stat                  root_stat)
{
  StatImpl::updateRootStat(root_stat.m_statImpl);
}


Stat
createRootStat(
  const std::string &   name,
  const StatSet &       stat_set)
{
  return StatImpl::createRootStat(name, stat_set);
}


void
deleteRootStat(
  Stat                  stat)
{
  StatImpl::deleteRootStat(stat.m_statImpl);
  stat.m_statImpl = 0;
}


std::vector<Stat> &
findStats(Stat root_stat, const std::string &path_tail, std::vector<Stat> &found_stats) {
  StatImpl::findStats(root_stat.m_statImpl, path_tail, found_stats);
  return found_stats;
}


StatImpl::StatImpl(
  const std::string &   name,
  StatMask              stat_mask,
  StatImpl *            parent_stat,
  const StatSet &       stat_set)
  : m_name(name),
    m_statMask(stat_mask),
    m_parentStat(parent_stat),
    m_substatLapCount(0.0),
    m_lapStartCount(0),
    m_substatList(),
    m_statSet(stat_set)
{}


StatImpl::~StatImpl()
{
  try {
    for (StatList::iterator it = m_substatList.begin(); it != m_substatList.end(); ++it)
      delete (*it).m_statImpl;
  }
  catch (std::exception &) {
  }
}


template<>
const Stat::Metric<LapCount> &
StatImpl::getMetric<LapCount>() const {
  return m_lapCount;
}


template<>
const Stat::Metric<CPUTime> &
StatImpl::getMetric<CPUTime>() const {
  return m_cpuTime;
}


template<>
const Stat::Metric<WallTime> &
StatImpl::getMetric<WallTime>() const {
  return m_wallTime;
}


template<>
const Stat::Metric<MPICount> &
StatImpl::getMetric<MPICount>() const {
  return m_MPICount;
}


template<>
const Stat::Metric<MPIByteCount> &
StatImpl::getMetric<MPIByteCount>() const {
  return m_MPIByteCount;
}


template<>
const Stat::Metric<HeapAlloc> &
StatImpl::getMetric<HeapAlloc>() const {
  return m_heapAlloc;
}


void
StatImpl::reset()
{
  m_lapStartCount = 0;

  m_lapCount.reset();
  m_cpuTime.reset();
  m_wallTime.reset();
  m_MPICount.reset();
  m_MPIByteCount.reset();
  m_heapAlloc.reset();
}


Stat
StatImpl::getSubstat(
  const std::string &  name)
{
  StatList::iterator it = std::find_if(m_substatList.begin(), m_substatList.end(), finder(name));

  if (it == m_substatList.end())
    throw std::runtime_error("Stat not found");
  else
    return *it;
}


StatImpl *
StatImpl::addSubstat(
  const std::string &  name,
  StatMask          stat_mask,
  const StatSet &      stat_set)
{
  StatList::iterator it = std::find_if(m_substatList.begin(), m_substatList.end(), finder(name));

  if (it == m_substatList.end()) {
    StatImpl *stat_impl = new StatImpl(name, stat_mask, this, stat_set);
    m_substatList.push_back(Stat(stat_impl));
    return stat_impl;
  }
  else
    return (*it).m_statImpl;
}


StatImpl &
StatImpl::start()
{
  if (shouldRecord()) {
    if (m_lapStartCount++ == 0) {
      m_lapCount.m_lapStart = m_lapCount.m_lapStop;

      m_cpuTime.m_lapStop = m_cpuTime.m_lapStart = value_now<CPUTime>();
      m_wallTime.m_lapStop = m_wallTime.m_lapStart = value_now<WallTime>();
      m_MPICount.m_lapStop = m_MPICount.m_lapStart = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = m_MPIByteCount.m_lapStart = value_now<MPIByteCount>();
      m_heapAlloc.m_lapStop = m_heapAlloc.m_lapStart = value_now<HeapAlloc>();
    }
  }

  return *this;
}


StatImpl &
StatImpl::lap()
{
  if (shouldRecord()) {
    if (m_lapStartCount > 0) {
      m_cpuTime.m_lapStop = value_now<CPUTime>();
      m_wallTime.m_lapStop = value_now<WallTime>();
      m_MPICount.m_lapStop = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();
    }
  }

  return *this;
}


StatImpl &
StatImpl::stop()
{
  if (shouldRecord()) {
    if (--m_lapStartCount <= 0) {
      m_lapStartCount = 0;
      m_lapCount.m_lapStop++;

      m_cpuTime.m_lapStop = value_now<CPUTime>();
      m_wallTime.m_lapStop = value_now<WallTime>();
      m_MPICount.m_lapStop = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();

      m_lapCount.addLap();
      m_cpuTime.addLap();
      m_wallTime.addLap();
      m_MPICount.addLap();
      m_MPIByteCount.addLap();
      m_heapAlloc.addLap();
    }
  }

  return *this;
}


double
StatImpl::accumulateSubstatLapCounts() const
{
  m_substatLapCount = m_lapCount.getAccumulatedLap(false);

  for (StatList::const_iterator it = m_substatList.begin(); it != m_substatList.end(); ++it)
    (*it).m_statImpl->accumulateSubstatLapCounts();

  for (StatList::const_iterator it = m_substatList.begin(); it != m_substatList.end(); ++it)
    m_substatLapCount += (*it).m_statImpl->m_substatLapCount;

  return m_substatLapCount;
}


void
StatImpl::checkpoint() const
{
  m_lapCount.checkpoint();
  m_cpuTime.checkpoint();
  m_wallTime.checkpoint();
  m_MPICount.checkpoint();
  m_MPIByteCount.checkpoint();
  m_heapAlloc.checkpoint();

  for (StatList::const_iterator it = m_substatList.begin(); it != m_substatList.end(); ++it)
    (*it).m_statImpl->checkpoint();
}


void
StatImpl::updateRootStat(
  StatImpl *                    root_stat)
{
  root_stat->m_lapCount.m_lapStop = value_now<LapCount>();
  root_stat->m_cpuTime.m_lapStop = value_now<CPUTime>();
  root_stat->m_wallTime.m_lapStop = value_now<WallTime>();
  root_stat->m_MPICount.m_lapStop = value_now<MPICount>();
  root_stat->m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
  root_stat->m_heapAlloc.m_lapStop = value_now<HeapAlloc>();

  root_stat->m_lapCount.m_accumulatedLap = root_stat->m_lapCount.m_lapStop - root_stat->m_lapCount.m_lapStart;
  root_stat->m_cpuTime.m_accumulatedLap = root_stat->m_cpuTime.m_lapStop - root_stat->m_cpuTime.m_lapStart;
  root_stat->m_wallTime.m_accumulatedLap = root_stat->m_wallTime.m_lapStop - root_stat->m_wallTime.m_lapStart;
  root_stat->m_MPICount.m_accumulatedLap = root_stat->m_MPICount.m_lapStop - root_stat->m_MPICount.m_lapStart;
  root_stat->m_MPIByteCount.m_accumulatedLap = root_stat->m_MPIByteCount.m_lapStop - root_stat->m_MPIByteCount.m_lapStart;
  root_stat->m_heapAlloc.m_accumulatedLap = root_stat->m_heapAlloc.m_lapStop - root_stat->m_heapAlloc.m_lapStart;
}



Stat
StatImpl::createRootStat(
  const std::string &           name,
  const StatSet &               stat_set)
{
  StatImpl *stat_impl = new StatImpl(name, 0, 0, stat_set);
  stat_impl->m_parentStat = stat_impl;
  StatTop::s_statTop.m_statImpl = stat_impl;
  return Stat(stat_impl);
}


void
StatImpl::deleteRootStat(
  StatImpl *                    root_stat)
{
  delete root_stat;
  StatTop::s_statTop.m_statImpl = 0;
}


void
StatImpl::findStat(
  StatImpl *                    stat,
  std::vector<std::string> &    path_tail_vector,
  std::vector<Stat> &           found_stats)
{
  if (stat->begin() == stat->end()) { // at leaf

  }
  else
    for (StatList::const_iterator it = stat->begin(); it != stat->end(); ++it)
      findStat((*it).m_statImpl, path_tail_vector, found_stats);
}


std::vector<Stat> &
StatImpl::findStats(
  StatImpl *                    root_stat,
  const std::string &           path_tail,
  std::vector<Stat> &           found_stats)
{
  std::vector<std::string> path_tail_vector;

  findStat(root_stat, split(path_tail, '.', path_tail_vector), found_stats);

  return found_stats;
}



Stat::Stat(const std::string &name, const Stat parent)
  : m_statImpl(StatImpl::reg(name, parent.getStatMask(), parent.m_statImpl, parent.getStatSet()))
{}

Stat::Stat(const std::string &name, StatMask stat_mask, const Stat parent)
  : m_statImpl(StatImpl::reg(name, stat_mask, parent.m_statImpl, parent.getStatSet()))
{}

Stat::Stat(const std::string &name, const Stat parent, const StatSet &stat_set)
  : m_statImpl(StatImpl::reg(name, parent.getStatMask(), parent.m_statImpl, stat_set))
{}

Stat::Stat(const std::string &name, StatMask stat_mask, const Stat parent, const StatSet &stat_set)
  : m_statImpl(StatImpl::reg(name, stat_mask, parent.m_statImpl, stat_set))
{}


const std::string &
Stat::getName() const {
  return m_statImpl->m_name;
}

StatMask
Stat::getStatMask() const {
  return m_statImpl->getStatMask();
}

const StatSet &
Stat::getStatSet() const {
  return m_statImpl->getStatSet();
}

double
Stat::getSubstatLapCount() const {
  return m_statImpl->getSubstatLapCount();
}

const StatList &
Stat::getStatList() const {
  return m_statImpl->getStatList();
}

template<class T>
const Stat::Metric<T> &
Stat::getMetric() const {
  return m_statImpl->getMetric<T>();
}

template const Stat::Metric<LapCount> &Stat::getMetric<LapCount>() const;
template const Stat::Metric<CPUTime> &Stat::getMetric<CPUTime>() const;
template const Stat::Metric<WallTime> &Stat::getMetric<WallTime>() const;
template const Stat::Metric<MPICount> &Stat::getMetric<MPICount>() const;
template const Stat::Metric<MPIByteCount> &Stat::getMetric<MPIByteCount>() const;
template const Stat::Metric<HeapAlloc> &Stat::getMetric<HeapAlloc>() const;


bool
Stat::shouldRecord() const
{
  return m_statImpl->shouldRecord();
}

StatList::iterator
Stat::begin()
{
  return m_statImpl->begin();
}

StatList::const_iterator
Stat::begin() const
{
  return m_statImpl->begin();
}

StatList::iterator
Stat::end()
{
  return m_statImpl->end();
}

StatList::const_iterator
Stat::end() const
{
  return m_statImpl->end();
}

double
Stat::accumulateSubstatLapCounts() const {
  return m_statImpl->accumulateSubstatLapCounts();
}

Stat &
Stat::start() {
  m_statImpl->start();
  return *this;
}

Stat &
Stat::lap() {
  m_statImpl->lap();
  return *this;
}

Stat &
Stat::stop() {
  m_statImpl->stop();
  return *this;
}

void
Stat::checkpoint() const {
  m_statImpl->checkpoint();
}


StatTop::StatTop(Stat stat)
{
  s_statTop.m_statImpl = stat.m_statImpl;
}

StatTop::StatTop(const std::string &name)
{
  Stat new_stat(name, s_statTop);

  s_statTop.m_statImpl = new_stat.m_statImpl;
}

StatTop::~StatTop()
{
  if (s_statTop.m_statImpl)
    s_statTop.m_statImpl = s_statTop.m_statImpl->m_parentStat;
}

Stat
StatTop::getTop()
{
  return s_statTop;
}

namespace {

size_t
s_statNameMaxWidth = DEFAULT_STAT_NAME_MAX_WIDTH;		///< Maximum width for names

} // namespace


//
// XyceRootStat member functions:
//
XyceRootStat::XyceRootStat()
  : m_xyceStat(Xyce::Stats::createRootStat("Xyce", xyceStatSet()))
{}


XyceRootStat::~XyceRootStat()
{
  Xyce::Stats::deleteRootStat(m_xyceStat);
}


Xyce::Stats::Stat & XyceRootStat::xyceStat()
{
  return m_xyceStat;
}


StatSet &
xyceStatSet()
{
  static StatSet s_xyceStatSet(STAT_ALL);

  return s_xyceStatSet;
}


XyceRootStat *xyceRootStat()
{
  static XyceRootStat *s_xyceRootStat = new XyceRootStat();

  return s_xyceRootStat;
}


Stat &
xyceStat()
{
  return xyceRootStat()->xyceStat();
}

void
setEnabledStatMask(
  StatMask                      stat_mask)
{
  xyceStatSet().setEnabledStatMask(stat_mask);
}


StatMask
getEnabledStatMask()
{
  return xyceStatSet().getEnabledStatMask();
}


void
setTimeFormat(
  int                           time_format)
{
  setStatTimeFormat(time_format);
}


void
setTimeFormatMillis()
{
  if ((getTimeFormat() & Xyce::TIMEFORMAT_STYLE_MASK ) == Xyce::TIMEFORMAT_HMS) {
    if (getXyceWallTime() > 3600.0)
      setTimeFormat(getTimeFormat() & ~Xyce::TIMEFORMAT_MILLIS);
    else
      setTimeFormat(getTimeFormat() | Xyce::TIMEFORMAT_MILLIS);
  }
  else if ((getTimeFormat() & Xyce::TIMEFORMAT_STYLE_MASK ) == Xyce::TIMEFORMAT_SECONDS) {
    if (getXyceWallTime() > 1000.0)
      setTimeFormat(getTimeFormat() & ~Xyce::TIMEFORMAT_MILLIS);
    else
      setTimeFormat(getTimeFormat() | Xyce::TIMEFORMAT_MILLIS);
  }
}


int
getTimeFormat()
{
  return getStatTimeFormat();
}


void
setStatNameMaxWidth(
  size_t                        width)
{
  s_statNameMaxWidth = width;
}


size_t
getStatNameMaxWidth()
{
  return s_statNameMaxWidth;
}


MetricTraits<CPUTime>::Type
getXyceCPUTime()
{
  return xyceStat().getMetric<CPUTime>().getAccumulatedLap(false);
}


MetricTraits<WallTime>::Type
getXyceWallTime()
{
  return xyceStat().getMetric<WallTime>().getAccumulatedLap(false);
}


MetricTraits<CPUTime>::Type
getCPULapTime(Stat stat) {
  return stat.getMetric<CPUTime>().getLap();
}

MetricTraits<CPUTime>::Type
getCPUAccumulatedLapTime(Stat stat) {
  return stat.getMetric<CPUTime>().getAccumulatedLap(false);
}

} // namespace Stats
} // namespace Xyce
