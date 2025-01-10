//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <typeinfo>
#include <utility>
#include <algorithm>
#include <limits>

#include <N_UTL_PrintStats.h>
#include <N_UTL_PrintTable.h>

#include <N_UTL_NoCase.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Marshal.h>

namespace Xyce {
namespace Stats {

namespace {
struct ParallelStats;
}
}

namespace Util {

template <class T>
Marshal &operator<<(Marshal &mout, const Stats::Stat::Metric<T> &t);

Marshal &operator<<(Marshal &mout, const Stats::Stat &t);

Marshal &operator>>(Marshal &min, Stats::ParallelStats &t);
}

namespace Stats {

namespace {

///
/// Class <b>Percent</b> is a functor which display the percentage of the numerator
/// to the denominator.  The value is displayed at (xx.xx%) for values in the range 0.01
/// to 99.99, (0.00%) if zero, (<0.01%) is less than 0.01, and (100.0%) for 100 percent.
///
///
struct Percent
{
  Percent(double numerator, double denominator)
    : m_numerator(numerator),
      m_denominator(denominator)
  {}

  ///
  /// Member function <b>operator()</b> writes the percentage as a string to the
  /// output stream.
  ///
  /// @param os    a <b>std::ostream</b> reference to the output stream
  ///        to write to.
  ///
  /// @return      a <b>std::ostream</b> reference to the output stream
  ///        written to.
  ///
  std::ostream &operator()(std::ostream &os) const;

private:
  double    m_numerator;
  double    m_denominator;
};


std::ostream &
Percent::operator()(
  std::ostream &  os) const
{
  std::ostringstream strout;

  if (m_numerator == 0.0)
    strout << "( 0.00%)";
  else if (m_denominator == 0.0)
    strout << "(  NaN )";
  else {
    double ratio = m_numerator/m_denominator*100.0;
    if (ratio < 0.01)
      strout << "(<0.01%)";
    else if (ratio >= 100.0)
      strout << "(" << std::setw(5) << std::setprecision(1) << std::fixed << ratio << "%)";
    else
      strout << "(" << std::setw(5) << std::setprecision(2) << std::fixed << ratio << "%)";
  }

  return os << strout.str();
}


///
/// Member function <b>operator&lt;&lt;</b> ...
///
/// @param os      a <b>std::ostream</b> variable ...
///
/// @param p      a <b>StatsImpl::Percent</b> variable ...
///
/// @return      a <b>std::ostream</b> ...
///
inline std::ostream &operator<<(std::ostream &os, const Percent &p) {
  return p(os);
}

struct ParallelStats
{
  template <typename T>
  struct Metric
  {
    Metric()
      : m_value(0),
        m_sum(0.0),
        m_min(std::numeric_limits<double>::max()),
        m_max(0.0)
    {}

    typename MetricTraits<T>::Type  m_value;  ///< Metric value
    typename MetricTraits<T>::Type  m_checkpoint;  ///< Metric checkpointed value
    double                          m_sum;    ///< Reduction sum
    double                          m_min;    ///< Reduction min
    double                          m_max;          ///< Reduction max

    void accumulate(const Metric<T> &metric, bool checkpoint) {
      double value = static_cast<double>(metric.m_value);
      if (checkpoint)
        value -= static_cast<double>(metric.m_checkpoint);

      m_sum += value;
      m_min = std::min(m_min, value);
      m_max = std::max(m_max, value);
    }
  };

  ParallelStats()
    : m_name(),
      m_timerMask(0),
      m_subtimerLapCount(0),
      m_lapCount(),
      m_cpuTime(),
      m_wallTime(),
      m_MPICount(),
      m_MPIByteCount(),
      m_heapAlloc(),
      m_subtimerList()
  {}

  ParallelStats(const ParallelStats &parallel_timer)
    : m_name(parallel_timer.m_name),
      m_timerMask(parallel_timer.m_timerMask),
      m_subtimerLapCount(parallel_timer.m_subtimerLapCount),
      m_lapCount(parallel_timer.m_lapCount),
      m_cpuTime(parallel_timer.m_cpuTime),
      m_wallTime(parallel_timer.m_wallTime),
      m_MPICount(parallel_timer.m_MPICount),
      m_MPIByteCount(parallel_timer.m_MPIByteCount),
      m_heapAlloc(parallel_timer.m_heapAlloc),
      m_subtimerList(parallel_timer.m_subtimerList)
  {}

  ParallelStats &operator=(const ParallelStats &parallel_timer) {
    if (this != &parallel_timer)
    {
      m_name = parallel_timer.m_name;
      m_timerMask = parallel_timer.m_timerMask;
      m_subtimerLapCount = parallel_timer.m_subtimerLapCount;
      m_lapCount = parallel_timer.m_lapCount;
      m_cpuTime = parallel_timer.m_cpuTime;
      m_wallTime = parallel_timer.m_wallTime;
      m_MPICount = parallel_timer.m_MPICount;
      m_heapAlloc = parallel_timer.m_heapAlloc;
      m_subtimerList = parallel_timer.m_subtimerList;
    }

    return *this;
  }

  template <class T>
  const Metric<T> &getMetric() const;

  std::string                   m_name;                 ///< Name of the timer
  Stats::StatMask               m_timerMask;
  double                        m_subtimerLapCount;     ///< Sum of subtimer lap counts and m_lapCount

  Metric<LapCount>              m_lapCount;             ///< Number of laps accumulated
  Metric<CPUTime>               m_cpuTime;              ///< CPU time
  Metric<WallTime>              m_wallTime;             ///< Wall time
  Metric<MPICount>              m_MPICount;             ///< MPI call count
  Metric<MPIByteCount>          m_MPIByteCount;	        ///< MPI byte count
  Metric<HeapAlloc>             m_heapAlloc;            ///< MPI byte count

  std::list<ParallelStats>      m_subtimerList;         ///< Sub timers
};

template<>
const ParallelStats::Metric<LapCount> &
ParallelStats::getMetric<LapCount>() const {
  return m_lapCount;
}


template<>
const ParallelStats::Metric<CPUTime> &
ParallelStats::getMetric<CPUTime>() const {
  return m_cpuTime;
}


template<>
const ParallelStats::Metric<WallTime> &
ParallelStats::getMetric<WallTime>() const {
  return m_wallTime;
}


template<>
const ParallelStats::Metric<MPICount> &
ParallelStats::getMetric<MPICount>() const {
  return m_MPICount;
}


template<>
const ParallelStats::Metric<MPIByteCount> &
ParallelStats::getMetric<MPIByteCount>() const {
  return m_MPIByteCount;
}


template<>
const ParallelStats::Metric<HeapAlloc> &
ParallelStats::getMetric<HeapAlloc>() const {
  return m_heapAlloc;
}

class finder
{
  using result_type = bool;
  using first_argument_type = ParallelStats;

public:
  finder(const std::string &name)
    : m_name(name)
  {}

  bool operator()(const ParallelStats &parallel_timer) const {
    return equal_nocase(parallel_timer.m_name, m_name);
  }

private:
  std::string           m_name;
};


void
merge_parallel_timer(
  ParallelStats &       p0,
  const ParallelStats & p1,
  bool                  checkpoint)
{
  p0.m_timerMask = p1.m_timerMask;
  p0.m_subtimerLapCount += p1.m_subtimerLapCount;
  p0.m_lapCount.accumulate(p1.m_lapCount, checkpoint);
  p0.m_cpuTime.accumulate(p1.m_cpuTime, checkpoint);
  p0.m_wallTime.accumulate(p1.m_wallTime, checkpoint);
  p0.m_MPICount.accumulate(p1.m_MPICount, checkpoint);
  p0.m_MPIByteCount.accumulate(p1.m_MPIByteCount, checkpoint);
  p0.m_heapAlloc.accumulate(p1.m_heapAlloc, checkpoint);


  for (std::list<ParallelStats>::const_iterator p1_it = p1.m_subtimerList.begin(); p1_it != p1.m_subtimerList.end(); ++p1_it) {
    std::list<ParallelStats>::iterator p0_it = std::find_if(p0.m_subtimerList.begin(), p0.m_subtimerList.end(), finder((*p1_it).m_name));
    if (p0_it == p0.m_subtimerList.end()) {
      p0.m_subtimerList.push_back((*p1_it));
      p0_it = --p0.m_subtimerList.end();
      merge_parallel_timer(*p0_it, *p1_it, checkpoint);
    }
    else
      merge_parallel_timer(*p0_it, *p1_it, checkpoint);
  }
}


void
collect_timers(
  Stat &                root_timer,
  ParallelStats &       parallel_timer,
  bool                  checkpoint,
  Parallel::Machine     comm)
{
#ifdef Xyce_PARALLEL_MPI

  Util::Marshal mout;
  mout << root_timer;

  const int parallel_root = 0 ;
  const int parallel_size = Parallel::size(comm);
  const int parallel_rank = Parallel::rank(comm);

  // Gather the send counts on root processor
  std::string send_string(mout.str());

  ParallelStats root_parallel_timer;

  // We need to gather the timer data in a number of 'cycles' where we only receive from a portion of the other
  // processors each cycle.  This is because buffer allocation-failures have been observed for runs on very large
  // numbers of processors if the 'root' processor tries to allocate a buffer large enough to hold timing data from all
  // other procesors.
  int num_cycles = parallel_size < 1024 ? 1 : 16;

  std::vector<char> buffer;

  for (int k = 0; k < num_cycles; ++k) {
    std::vector<int> recv_count(parallel_size, 0);
    int * const recv_count_ptr = &recv_count[0] ;

    // Send_count is the amount of data this processor needs to send.
    int send_count = send_string.size();

    // Should this processor send on the current cycle ? If not, set send_count to 0.
    if ((parallel_rank + k)%num_cycles != 0) {
      send_count = 0;
    }

    int result = MPI_Gather(&send_count, 1, MPI_INT,
                            recv_count_ptr, 1, MPI_INT,
                            parallel_root, comm);
    if (MPI_SUCCESS != result) {
      std::ostringstream message ;
      message << "stk::diag::collect_timers FAILED: MPI_Gather = " << result ;
      throw std::runtime_error(message.str());
    }

    // Receive counts are only non-zero on the root processor:
    std::vector<int> recv_displ(parallel_size + 1, 0);

    for (int i = 0 ; i < parallel_size ; ++i) {
      recv_displ[i + 1] = recv_displ[i] + recv_count[i] ;
    }

    const int recv_size = recv_displ[parallel_size] ;

    buffer.assign(recv_size, 0);

    {
      const char * const send_ptr = send_string.data();
      char * const recv_ptr = recv_size ? & buffer[0] : 0;
      int * const recv_displ_ptr = & recv_displ[0] ;

      result = MPI_Gatherv((void *) send_ptr, send_count, MPI_CHAR,
                           recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                           parallel_root, comm);
      if (MPI_SUCCESS != result) {
        std::ostringstream message ;
        message << "stk::diag::collect_timers FAILED: MPI_Gatherv = " << result ;
        throw std::runtime_error(message.str());
      }

      std::vector<ParallelStats> parallel_timer_vector;
      parallel_timer_vector.reserve(parallel_size);

      if (parallel_rank == parallel_root) {
        for (int j = 0; j < parallel_size; ++j) {
          int received_count = recv_displ[j+1] - recv_displ[j];
          if (received_count > 0) {
            //grow parallel_timer_vector by 1:
            parallel_timer_vector.resize(parallel_timer_vector.size()+1);
            Util::Marshal min(std::string(recv_ptr + recv_displ[j], recv_ptr + recv_displ[j + 1]));
            //put this data into the last entry of parallel_timer_vector:
            min >> parallel_timer_vector[parallel_timer_vector.size()-1];
          }
        }

        if (parallel_rank == parallel_root && send_count > 0)
          root_parallel_timer = parallel_timer_vector[0];

        for (size_t j = 0; j < parallel_timer_vector.size(); ++j)
          merge_parallel_timer(root_parallel_timer, parallel_timer_vector[j], checkpoint);
      }
    }
  }
  parallel_timer = root_parallel_timer;
#endif
}

PrintTable &
printSubtable(
  PrintTable &  table,
  const Stat &  root_timer,
  const Stat &  timer,
  MetricsMask   metrics_mask,
  int           depth,
  bool          timer_checkpoint)
{
  if (timer.getSubstatLapCount() != 0.0) {
    if (timer.shouldRecord()) {
      if (timer.getStatMask() == 0 || timer.getMetric<LapCount>().getAccumulatedLap(timer_checkpoint) > 0) {
        table << justify(PrintTable::Cell::LEFT) << indent(depth) << timer.getName() << end_col
              << justify(PrintTable::Cell::RIGHT) << timer.getMetric<LapCount>().getAccumulatedLap(timer_checkpoint) << end_col;

        if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<CPUTime>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<WallTime>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPICount>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<HeapAlloc>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<HeapAlloc>::format(timer.getMetric<HeapAlloc>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<HeapAlloc>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<HeapAlloc>().getAccumulatedLap(timer_checkpoint)) << end_col;
      }
      else
        table << justify(PrintTable::Cell::LEFT) << indent(depth) << span << timer.getName() << end_col;

      table << end_row;
      depth++;
    }

    for (StatList::const_iterator it = timer.begin(); it != timer.end(); ++it)
      printSubtable(table, root_timer, *it, metrics_mask, depth, timer_checkpoint);
  }

  return table;
}


PrintTable &
printSubtable(
  PrintTable &      table,
  const ParallelStats &         root_timer,
  const ParallelStats &         timer,
  MetricsMask      metrics_mask,
  int        depth,
  bool        timer_checkpoint)
{
  if (timer.m_subtimerLapCount != 0.0) {
    if (timer.m_timerMask == 0 || timer.getMetric<LapCount>().m_sum > 0) {
      table << justify(PrintTable::Cell::LEFT) << indent(depth) << timer.m_name << end_col
            << justify(PrintTable::Cell::RIGHT) << timer.getMetric<LapCount>().m_sum << end_col;

      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_sum, root_timer.getMetric<CPUTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_min, root_timer.getMetric<CPUTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_max, root_timer.getMetric<CPUTime>().m_sum) << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_sum, root_timer.getMetric<WallTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_min, root_timer.getMetric<WallTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_max, root_timer.getMetric<WallTime>().m_sum) << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_sum, root_timer.getMetric<MPICount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_min, root_timer.getMetric<MPICount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_max, root_timer.getMetric<MPICount>().m_sum) << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_sum, root_timer.getMetric<MPIByteCount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_min, root_timer.getMetric<MPIByteCount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_max, root_timer.getMetric<MPIByteCount>().m_sum) << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<HeapAlloc>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<HeapAlloc>::format(timer.getMetric<HeapAlloc>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<HeapAlloc>().m_sum, root_timer.getMetric<HeapAlloc>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<HeapAlloc>::format(timer.getMetric<HeapAlloc>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<HeapAlloc>().m_min, root_timer.getMetric<HeapAlloc>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<HeapAlloc>::format(timer.getMetric<HeapAlloc>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<HeapAlloc>().m_max, root_timer.getMetric<HeapAlloc>().m_sum) << end_col;
    }
    else
      table << justify(PrintTable::Cell::LEFT) << indent(depth) << span << timer.m_name << end_col;

    table << end_row;
    depth++;
  }

  for (std::list<ParallelStats>::const_iterator it = timer.m_subtimerList.begin(); it != timer.m_subtimerList.end(); ++it)
    printSubtable(table, root_timer, *it, metrics_mask, depth, timer_checkpoint);

  return table;
}


PrintTable &
printTable(
  PrintTable &          table,
  Stat &                root_timer,
  MetricsMask           metrics_mask,
  size_t                name_width,
  bool                  timer_checkpoint)
{
  updateRootStat(root_timer);

  root_timer.accumulateSubstatLapCounts();

  if (metrics_mask & getEnabledStatMetricsMask()) {
    table.setAutoEndCol(false);

    table << cell_width(name_width) << justify(PrintTable::Cell::CENTER) << "Stats" << (timer_checkpoint ? " (delta time)" : "") << end_col
          << justify(PrintTable::Cell::CENTER) << "Count"  << end_col;

    if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<CPUTime>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col;
    if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<WallTime>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col;
    if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPICount>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col;
    if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col;
    if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<HeapAlloc>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<HeapAlloc>::table_header() << end_col;

    table << end_header;

    printSubtable(table, root_timer, root_timer, metrics_mask, 0, timer_checkpoint);

    if (timer_checkpoint)
      root_timer.checkpoint();
  }

  return table;
}


PrintTable &
printTable(
  PrintTable &          table,
  Stat &                root_timer,
  MetricsMask           metrics_mask,
  size_t                name_width,
  bool                  timer_checkpoint,
  Parallel::Machine     parallel_machine)
{
  updateRootStat(root_timer);

  root_timer.accumulateSubstatLapCounts();

  ParallelStats parallel_timer;

  Xyce::Stats::collect_timers(root_timer, parallel_timer, timer_checkpoint, parallel_machine);

  int parallel_rank = Parallel::rank(parallel_machine);
  if (parallel_rank == 0) {
    if (metrics_mask & getEnabledStatMetricsMask()) {
      table.setAutoEndCol(false);

      table << end_col << end_col;

      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<HeapAlloc>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<HeapAlloc>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<HeapAlloc>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<HeapAlloc>::table_header() << end_col;

      table << end_header;
      table << cell_width(name_width) << justify(PrintTable::Cell::CENTER) << "Stats" << (timer_checkpoint ? " (delta time)" : "") << end_col
            << justify(PrintTable::Cell::CENTER) << "Count"  << end_col;

      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledStatMetricsMask() & MetricTraits<HeapAlloc>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;

      table << end_header;

      printSubtable(table, parallel_timer, parallel_timer, metrics_mask, 0, timer_checkpoint);
    }

    if (timer_checkpoint)
      root_timer.checkpoint();
  }

  return table;
}

} // namespace <empty>


std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint)
{
  Xyce::PrintTable print_table;

  printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint);

  os << print_table;

  return os;
}


std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint, Parallel::Machine parallel_machine)
{
  Xyce::PrintTable print_table;

  int parallel_size = Parallel::size(parallel_machine);
  if (parallel_size == 1)
    printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint);
  else
    printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint, parallel_machine);

  os << print_table;

  return os;
}


// std::ostream &printXML(std::ostream &os, MPI_Comm mpi_comm, MetricsMask metrics_mask) const;
std::ostream &printXML(std::ostream &os, MetricsMask metrics_mask, bool timer_checkpoint);

std::ostream &printSubXML(std::ostream &os, MetricsMask metrics_mask, int depth, bool timer_checkpoint);

} // namespace Stats

namespace Util {

Marshal &operator<<(Marshal &mout, const Stats::Stat &t);

template <class T>
Marshal &operator<<(Marshal &mout, const Stats::Stat::Metric<T> &t) {
  mout << t.getAccumulatedLap(false) << t.getAccumulatedLap(true);

  return mout;
}

Marshal &operator<<(Marshal &mout, const Stats::Stat &t) {
  mout << t.getName() << t.getStatMask() << t.getSubstatLapCount()
       << t.getMetric<Stats::LapCount>() << t.getMetric<Stats::CPUTime>() << t.getMetric<Stats::WallTime>()
       << t.getMetric<Stats::MPICount>() << t.getMetric<Stats::MPIByteCount>() << t.getMetric<Stats::HeapAlloc>();

  mout << t.getStatList();

  return mout;
}

Marshal &operator>>(Marshal &min, Stats::ParallelStats &t) {
  min >> t.m_name >> t.m_timerMask >> t.m_subtimerLapCount
      >> t.m_lapCount.m_value
      >> t.m_lapCount.m_checkpoint
      >> t.m_cpuTime.m_value
      >> t.m_cpuTime.m_checkpoint
      >> t.m_wallTime.m_value
      >> t.m_wallTime.m_checkpoint
      >> t.m_MPICount.m_value
      >> t.m_MPICount.m_checkpoint
      >> t.m_MPIByteCount.m_value
      >> t.m_MPIByteCount.m_checkpoint
      >> t.m_heapAlloc.m_value
      >> t.m_heapAlloc.m_checkpoint;

  min >> t.m_subtimerList;

  return min;
}

} // namespace Util
} // namespace Xyce
