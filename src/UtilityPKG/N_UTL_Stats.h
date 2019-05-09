//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Stats_h
#define Xyce_N_UTL_Stats_h

#include <iosfwd>
#include <vector>
#include <list>
#include <string>

#include <N_UTL_StatMetricTraits.h>
#include <N_PDS_ParallelMachine.h>

///
/// @addtogroup DiagStatDetail
/// @{
///

namespace Xyce {
namespace Stats {

class Stat;
class StatSet;
class StatImpl;

typedef unsigned StatMask;        ///< Stat classification mask

/**
 * Function <b>getEnabledMetricsMask</b> retruns the stat enable bit mask.
 *
 * @return      a <b>MetricsMask</b> value of the stat enable bit
 *        mask.
 */
MetricsMask getEnabledStatMetricsMask();

/**
 * Function <b>setEnabledMetricsMask</b> set the stat enable bit mask to
 * <b>stat_mask</b>.
 *
 * @param stat_mask    a <b>MetricsMask</b> value to set the stat enable bit
 *        mask to.
 *
 */
void setEnabledStatMetricsMask(MetricsMask stat_mask);

/**
 * Function <b>updateRootStat</b> updates the root stats stop and total
 * metric values with the current time.
 *
 * @param root_stat      a <b>Stat</b> reference to the root stat.
 *
 */
void updateRootStat(Stat root_stat);

/**
 * Function <b>createRootStat</b> creates a root stat.  Root stats are the root of a stat
 * hierarchy.  The stat_set specifies the stat groupings for this root stat.  The percentage of a
 * child stat is the ratio of that stat the its root.
 *
 * @param name                  a <b>std::string</b> const reference to the name of the new root
 *                              stat.
 *
 * @param stat_set           a <b>StatSet</b> const reference of the stat set of the new root
 *                              stat.
 *
 * @return      a <b>Stat</b> value of the new root stat.
 */
Stat createRootStat(const std::string &name, const StatSet &stat_set);

/**
 * Function <b>deleteRootStat</b> deletes a root stat and all of it's children stats.  All
 * children Stats are invalidated and can no longer be used.
 *
 * @param                       a <b>Stat</b> value of the root stat to delete.
 */
void deleteRootStat(Stat stat);

/**
 * @brief Member function <code>findStat</code> return a vector of stats whose tail of the dot
 * separated name from root_time to leaf matches the specified path_tail.
 *
 * @param root_stat    a <code>Stat</code> value of the root to begin search.
 *
 * @param path_tail    a <code>std::string</code> const reference to the dot separated tail
 *                              to match.
 *
 * @param found_stat    a <code>std::vector<Stat></code> reference to the vector to store
 *                              matching stats.
 *
 * @return      a <code>std::vector<Stat></code> reference to found_stat.
 */
std::vector<Stat> &findStats(Stat root_stat, const std::string &path_tail, std::vector<Stat> &found_stats);

/**
 * @brief Class <b>StatSet</b> implements a set of stat classifications.  A time classification
 * consists of a bit mask set StatMask
 *
 */
class StatSet
{
public:
  explicit StatSet(StatMask enabled_stat_mask)
    : m_enabledStatMask(enabled_stat_mask)
  {}

  StatSet(const StatSet &stat_set)
    : m_enabledStatMask(stat_set.m_enabledStatMask)
  {}

private:
  StatSet &operator=(StatSet &stat_set) {
    m_enabledStatMask = stat_set.m_enabledStatMask;

    return *this;
  }

public:
  ~StatSet()
  {}

  /**
   * Member function <b>getEnabledStatMask</b> returns the stat enable bit mask.
   *
   * @return      a <b>StatMask</b> value of the stat enable bit
   *        mask.
   */
  StatMask getEnabledStatMask() const {
    return m_enabledStatMask;
  }

  /**
   * Member function <b>setEnabledStatMask</b> set the stat enable bit mask to
   * <b>stat_mask</b>.
   *
   * @param stat_mask    a <b>StatMask</b> value to set the stat enable bit
   *        mask to.
   *
   */
  void setEnabledStatMask(StatMask stat_mask) {
    m_enabledStatMask = stat_mask;
  }

  /**
   * Member function <b>shouldRecord</b> returns true if any of the specified stat
   * bit masks are set in the enable stat bit mask.
   *
   * @param stat_mask    a <b>StatMask</b> value to test the enable stat
   *        bit mask against.
   *
   */
  bool shouldRecord(StatMask stat_mask) const {
    return (stat_mask == 0 || (m_enabledStatMask & stat_mask));
  }

private:
  StatMask    m_enabledStatMask;  ///< Bit mask of enabled stat
};


typedef std::list<Stat> StatList;    ///< A vector of subordinate stats.

/**
 * @brief Class <b>Stat</b> implements a diagnostic stat and stat container for the
 * collection and display of execution times.
 *
 */
class Stat
{
  friend class StatImpl;
  friend class StatTop;
  friend class TimeBlock;
  friend class TimeBlockSynchronized;
  friend void updateRootStat(Stat);
  friend Stat createRootStat(const std::string &, const StatSet &);
  friend void deleteRootStat(Stat);
  friend std::vector<Stat> &findStats(Stat, const std::string &, std::vector<Stat> &);

public:
  /**
   * Class <b>Metric</b> maintains the metric data for the stat or counter.  The
   * start and stop times maintain the current lap time.  When a lap completes, its
   * time/count is accumlated to the total.  The total time/count can be stored in the
   * checkpoint member variable.  The total can be retrieved as either absolute time/count
   * the diffence from the checkpoint value.
   *
   */
  template <typename T>
  struct Metric
  {
    Metric()
      : m_lapStart(0),
        m_lapStop(0),
        m_accumulatedLap(0),
        m_checkpoint(0)
    {}

    /**
     * Member function <b>reset</b> resets the metric values to zero.
     *
     */
    void reset() {
      m_lapStart = m_lapStop = m_accumulatedLap = m_checkpoint = 0;
    }

    /**
     * Member function <b>addLap</b> adds the most recently completed lap to the total.
     *
     * @return      a <b>T</b> value of the total.
     */
    typename MetricTraits<T>::Type addLap() {
      return m_accumulatedLap += m_lapStop - m_lapStart;
    }

    /**
     * Member function <b>checkpoint</b> checkpoints the metrics by storing the
     * total time in the checkpoint value.
     *
     */
    void checkpoint() const {
      m_checkpoint = m_accumulatedLap;
    }

    /**
     * Member function <b>getLap</b> returns the value of the most recently
     * lap.
     *
     * @return      a <b>T</b> value of the most recent lap.
     */
    typename MetricTraits<T>::Type getLap() const {
      return m_lapStop - m_lapStart;
    }

    /**
     * Member function <b>getStart</b> returns the start value of the most recent lap.
     *
     * @return      a <b>T</b> value of the start of the most recent lap.
     */
    typename MetricTraits<T>::Type getStart() const {
      return m_lapStart;
    }

    /**
     * Member function <b>getStop</b> returns the stop value of the most recent lap.
     *
     * @return      a <b>T</b> value of the stop of the most recent lap.
     */
    typename MetricTraits<T>::Type getStop() const {
      return m_lapStop;
    }

    /**
     * Member function <b>getAccumulatedLap</b> returns the accumulated value of the metric.
     * If the <b>checkpoint</b> parameter if true, the value returned is the
     * difference between the accumulated value and the checkpointed value.
     *
     * @param checkpoint  a <b>bool</b> value of true of the checkpointed
     *        value is to be returned.
     *
     * @return      a <b>T</b> value of the accumulated or the
     *        checkpoint difference.
     */
    typename MetricTraits<T>::Type getAccumulatedLap(bool arg_checkpoint = false) const {
      if (arg_checkpoint)
        return m_accumulatedLap - m_checkpoint;
      else
        return m_accumulatedLap;
    }

    typename MetricTraits<T>::Type    m_lapStart;    ///< Most recent start time/count
    typename MetricTraits<T>::Type    m_lapStop;    ///< Most recent stop or lap time/count
    typename MetricTraits<T>::Type    m_accumulatedLap;  ///< Accumulated time/count
    mutable typename MetricTraits<T>::Type      m_checkpoint;    ///< Checkpointed time/count
  };

  /**
   * Creates a new <b>Stat</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the stat.
   *
   * @param parent    a <b>Stat</b> value of the parent stat.
   *
   */
  Stat(const std::string &name, const Stat parent);

  /**
   * Creates a new <b>Stat</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the stat.
   *
   * @param parent    a <b>Stat</b> value of the parent stat.
   *
   * @param stat_set    a <b>StatSet</b> value of the stat set used to interpret the
   *                            StatMask's of this and children stats.
   *
   */
  Stat(const std::string &name, const Stat parent, const StatSet &stat_set);

  /**
   * Creates a new <b>Stat</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the stat.
   *
   * @param stat_mask    a <b>StatMask</b> value which enables this stat.
   *
   * @param parent    a <b>Stat</b> value of the parent stat.
   *
   */
  Stat(const std::string &name, StatMask stat_mask, const Stat parent);

  /**
   * Creates a new <b>Stat</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the stat.
   *
   * @param stat_mask    a <b>StatMask</b> value which enables this stat.
   *
   * @param parent    a <b>Stat</b> value of the parent stat.
   *
   * @param stat_set    a <b>StatSet</b> value of the stat set used to interpret the
   *                            StatMask's of this and children stats.
   *
   */
  Stat(const std::string &name, StatMask stat_mask, const Stat parent, const StatSet &stat_set);

  /**
   * Creates the root <b>Stat</b> stat instance.
   *
   */
  explicit Stat(StatImpl &stat_impl)
    : m_statImpl(&stat_impl)
  {}

  explicit Stat(StatImpl *stat_impl)
    : m_statImpl(stat_impl)
  {}

  Stat(const Stat &stat)
    : m_statImpl(stat.m_statImpl)
  {}

  Stat &operator=(const Stat &stat) {
    if (this != &stat)
      m_statImpl = stat.m_statImpl;

    return *this;
  }

  virtual ~Stat()
  {}

  const StatList &getStatList() const;

  StatList::iterator begin();
  StatList::const_iterator begin() const;
  StatList::iterator end();
  StatList::const_iterator end() const;

  /**
   * Member function <b>getName</b> returns the name of the stat.
   *
   * @return      a <b>std::string</b> const reference to the stat's
   *        name.
   */
  const std::string &getName() const;


  /**
   * Member function <b>getStatMask</b> returns the stat mask of the stat.
   *
   * @return      a <b>StatMask</b> value to the stat mask.
   */
  const StatSet &getStatSet() const;

  /**
   * Member function <b>getStatMask</b> returns the stat mask of the stat.
   *
   * @return      a <b>StatMask</b> value to the stat mask.
   */
  StatMask getStatMask() const;

  bool shouldRecord() const;

  /**
   * Member function <b>getSubstatLapCount</b> returns the substat lap counter.
   *
   * @return      a <b>Counter</b> value of the substat lap
   *        counter.
   */
  double getSubstatLapCount() const;

  /**
   * Member function <b>getLapCount</b> returns the lap counter metric.  The lap
   * count metric is the number of times the stop function has been executed.
   *
   * @return      a <b>CounterMetric</b> const reference of the lap counter
   *        metric.
   */
  template <class T>
  const Metric<T> &getMetric() const;

  /**
   * Member function <b>accumulateSubstatLapCounts</b> accumulates the substat la
   * counts.
   *
   * @return      an <b>int</b> value of the count.
   */
  double accumulateSubstatLapCounts() const;

  /**
   * Member function <b>start</b> starts the lap stat.
   *
   * @return      a <b>Stat</b> reference to this stat.
   */
  Stat &start();

  /**
   * Member function <b>lap</b> sets the lap stop time.
   *
   * @return      a <b>Stat</b> reference to the stat.
   */
  Stat &lap();

  /**
   * Member function <b>stop</b> sets the lap stop time and sums the just completed
   * lap time to the stat.
   *
   * @return      a <b>Stat</b> reference to the stat.
   */
  Stat &stop();

  /**
   * Member function <b>checkpoint</b> checkpoints the metrics by storing the
   * total time in the checkpoint value.
   *
   */
  void checkpoint() const;

private:
  StatImpl *    m_statImpl;      ///< Reference to the actual stat
};


class StatTop 
{
  friend class StatImpl;
  
public:
  StatTop(Stat stat);

  StatTop(const std::string &name);

private:
  StatTop(const StatTop &);
  StatTop &operator=(const StatTop &);

public:
  ~StatTop();

  static Stat getTop();

private:
  static Stat   s_statTop;
};

/**
 * Class <b>TimeBlock</b> is a time sentry for timing a statement block.  The
 * stat is generally started upon construction. But, the start is delayed if the second
 * argument is false.  In this case, manually start the stat by calling the start()
 * function.  This gives the safety of using a sentry, but does not force to awkwardness
 * associated with local variables crossing the timed block.
 *
 */
class TimeBlock
{
public:
  /**
   * Creates a new <b>TimeBlock</b> instance.  The newly created instance will
   * start the stat if the <b>start</b> value is true, which is the default
   * case.  If the <b>start</b> value is false, the calling function is
   * responsible for starting the stat at the appropriate time.
   *
   * @param stat    a <b>Stat</b> reference to the stat accumulate
   *        block run times.
   *
   * @param start_stat  a <b>bool</b> value to have the stat started on
   *        construction.
   *
   */
  explicit TimeBlock(Stat &stat, bool start_stat = true)
    : m_stat(stat),
      m_started(start_stat)
  {
    if (start_stat)
      m_stat.start();
  }

  /**
   * Creates a new <b>TimeBlock</b> instance.  The newly created instance will
   * start the stat if the <b>start</b> value is true, which is the default
   * case.  If the <b>start</b> value is false, the calling function is
   * responsible for starting the stat at the appropriate time.
   *
   * @param stat    a <b>Stat</b> reference to the stat accumulate
   *        block run times.
   *
   * @param start_stat  a <b>bool</b> value to have the stat started on
   *        construction.
   *
   */
  explicit TimeBlock(StatTop &stat, bool start_stat = true)
    : m_stat(stat.getTop()),
      m_started(start_stat)
  {
    if (start_stat)
      m_stat.start();
  }

private:
  TimeBlock(const TimeBlock &);
  TimeBlock &operator=(const TimeBlock &);

public:
  /**
   * Destroys a <b>TimeBlock</b> instance.  Stops the stat if is has been started.
   *
   */
  ~TimeBlock() {
    try {
      if (m_started)
        m_stat.stop();
    }
    catch (...) {
    }
  }

  /**
   * Member function <b>start</b> starts the stat associated with the time block.
   *
   */
  void start() {
    m_started = true;
    m_stat.start();
  }

  /**
   * Member function <b>lap</b> sets the stop time of the stat associated with
   * the time block.
   *
   */
  void lap() {
    m_stat.lap();
  }

  /**
   * Member function <b>stop</b> stops the stat associated with the time block.
   *
   */
  void stop() {
    m_started = false;
    m_stat.stop();
  }

private:
  Stat          m_stat;  ///< Stat to accumulate block run times.
  bool          m_started;  ///< Stat has been started
};

/**
 * Class <b>TimeBlockSynchronized</b> is a time sentry for timing a statement
 * block.  The stat is generally started upon construction. But, the start is delayed
 * if the second argument is false.  In this case, manually start the stat by calling
 * the start() function.  This gives the safety of using a sentry, but does not force to
 * awkwardness associated with local variables crossing the timed block.
 *
 * Prior to starting the stat, an MPI synchronization barrier is set so that the
 * timing of routines which require MPI communication will all be at a known location
 * prior to executing.
 *
 */
class TimeBlockSynchronized
{
public:
  /**
   * Creates a new <b>TimeBlockSynchronized</b> instance.  If
   * <b>start_stat</b> is true, then the stat is started using the
   * <b>start()</b>.  An <b>MPI_Barrier</b> is called to synchronize the
   * start of the stat. The destructor will always stop a started stat.
   *
   * @param stat    a <b>Stat</b> reference to the stat to start.
   *
   * @param mpi_comm    a <b>MPI_Comm</b> value of the mpi communicator.

   * @param start_stat  a <b>bool</b> value to start the stat on construction.
   *
   */
  TimeBlockSynchronized(Stat &stat, Parallel::Machine mpi_comm, bool start_stat = true);

  /**
   * Destroys a <b>TimeBlockSynchronized</b> instance.  Stops the stat if it has
   * been started.
   *
   */
  ~TimeBlockSynchronized();

  /**
   * Member function <b>start</b> starts the stat associated with the time block.
   * An <b>MPI_Barrier</b> is executed prior to starting the stat.
   *
   */
  void start();

  /**
   * Member function <b>stop</b> stops the stat associated with the time block.
   *
   */
  void stop();

private:
  Stat &      m_stat;  ////< Stat to accumulate block run times.
  Parallel::Machine    m_mpiComm;  ////< MPI comm to synchronize across
  bool      m_started;  ////< Stat has been started
};


/**
 * @brief Enumeration <b><unnnamed></b> defines the bit mask values for the diagnostic
 * stat's in the <b>Diag</b> namespace.
 *
 */
enum StatSetMask {
  STAT_XYCE		= 0x00000001,		///< Enable Xyce stats
  STAT_ANALYSIS		= 0x00000002,		///< Enable analysis
  STAT_NONLINEAR	= 0x00000004,		///< Enable nonlinear solver
  STAT_SOLVER   	= 0x00000008,		///< Enable linear solver
  STAT_TIME_INTEGRATOR	= 0x00000010,		///< Enable time integrator

  STAT_PROFILE_1	= 0x00001000,		///< Enable profile 1 stats
  STAT_PROFILE_2	= 0x00002000,		///< Enable profile 2 stats
  STAT_PROFILE_3	= 0x00004000,		///< Enable profile 3 stats
  STAT_PROFILE_4	= 0x00008000,		///< Enable profile 4 stats
  STAT_APP_1		= 0x00010000,		///< Enable application defined 1
  STAT_APP_2		= 0x00020000,		///< Enable application defined 2
  STAT_APP_3		= 0x00040000,		///< Enable application defined 3
  STAT_APP_4		= 0x00080000,		///< Enable application defined 4
  STAT_ALL		= 0x000FFFFF,		///< Enable all stats
  STAT_NONE		= 0x00000000,		///< Enable no stats

  STAT_FORCE		= 0x00000000		///< Force stat to be active
};


StatSet &xyceStatSet();

Stat &xyceStat();

void xyceStatDestroy();

/**
 * @brief Enumeration <b><unnamed></b> defines some constants used for the
 * indented display of the diagnostic stats.
 *
 */
enum {
  DEFAULT_STAT_NAME_MAX_WIDTH = 40			///< Width to truncate the name
};

/**
 * Function <b>setEnabledStatMask</b> set the stat enable bit mask to
 * <b>stat_mask</b>.
 *
 * @param stat_mask		a <b>StatMask</b> value to set the stat enable bit
 *				mask to.
 *
 */
void setEnabledStatMask(StatMask stat_mask);

/**
 * Function <b>getEnabledStatMask</b> retruns the stat enable bit mask.
 *
 * @return			a <b>StatMask</b> value of the stat enable bit
 *				mask.
 */
StatMask getEnabledStatMask();

void setTimeFormat(int time_format);

void setTimeFormatMillis();

int getTimeFormat();

/**
 * @brief Member function <b>setStatNameMaxWidth</b> sets the maximum width for names
 * displayed in the stat output table.
 *
 * @param width		a <b>size_t</b> value to set for the maximum width for
 *				names displayed in the stat output table.
 *
 */
void setStatNameMaxWidth(size_t width);

/**
 * @brief Member function <b>getTimeNameMaxWidth</b> returns the width to use for the
 * name cell of the table display.
 *
 * @return			a <b>size_t</b> value of the width to use for the name
 *				cell of the table display.
 */
size_t getStatNameMaxWidth();

MetricTraits<CPUTime>::Type getCPULapTime(Stat stat);
MetricTraits<CPUTime>::Type getCPUAccumulatedLapTime(Stat stat);
MetricTraits<CPUTime>::Type getXyceCPUTime();
MetricTraits<CPUTime>::Type getXyceWallTime();


class XyceRootStat
{
  public:
    XyceRootStat();
    virtual ~XyceRootStat();
    Stat &xyceStat();

  private:
    Stat                        m_xyceStat;
};

} // namespace Stats
} // namespace Xyce

///
/// @}
///

#endif // Xyce_N_UTL_Stats_h

