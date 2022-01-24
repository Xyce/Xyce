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

#ifndef Xyce_N_UTL_StatMetricTraits_h
#define Xyce_N_UTL_StatMetricTraits_h

#include <string>

namespace Xyce {
namespace Stats {

typedef unsigned long MetricsMask;

/**
 * @brief Member function <code>setStatTimeFormat</code> sets the display format of time in the
 * output tables
 *
 * @param time_format    a <code>TimeFormat</code> variable...
 */
void setStatTimeFormat(int time_format);

/**
 * @brief Member function <code>getStatTimeFormat</code>
 *
 */
int getStatTimeFormat();

/**
 * @brief Enumeration <code>Metrics</code> assigns a unique but for each type of stat.  The
 * METRICS_FORCE type allows a stat to force itself to be active, even is not enabled.
 *
 */
enum Metrics {
  METRICS_LAP_COUNT      = 0x0001,              ///< Count of stat starts
  METRICS_CPU_TIME       = 0x0002,              ///< CPU runtime 
  METRICS_WALL_TIME      = 0x0004,              ///< Wall clock time
  METRICS_MPI_COUNT      = 0x0008,              ///< MPI class count
  METRICS_MPI_BYTE_COUNT = 0x0010,              ///< MPI byte count
  METRICS_HEAP_ALLOC     = 0x0020,              ///< Heap allocation
  METRICS_ALL            = 0x7FFF,

  METRICS_FORCE          = 0x8000               ///< Force metrics to be acquired
};


struct LapCount {};                             ///< Lap counter metric tag
struct CPUTime {};                              ///< CPU runtime metric tag
struct WallTime {};                             ///< Wall clock metric tag
struct MPICount {};                             ///< MPI call count metric tag
struct MPIByteCount {};                         ///< MPI byte count metric tag
struct HeapAlloc {};                            ///< Heap allocation metric tag

template <class T>
struct MetricTraits;

template<>
struct MetricTraits<LapCount>
{
  typedef unsigned Type;
  enum {METRIC = METRICS_LAP_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<CPUTime>
{
  typedef double Type;
  enum {METRIC = METRICS_CPU_TIME};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<WallTime>
{
  typedef double Type;
  enum {METRIC = METRICS_WALL_TIME};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type time);
};

template<>
struct MetricTraits<MPICount>
{
  typedef double Type;
  enum {METRIC = METRICS_MPI_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};

template<>
struct MetricTraits<MPIByteCount>
{
  typedef double Type;
  enum {METRIC = METRICS_MPI_BYTE_COUNT};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};


template<>
struct MetricTraits<HeapAlloc>
{
  typedef double Type;
  enum {METRIC = METRICS_HEAP_ALLOC};
  static Type value_now();
  static std::string table_header();
  static std::string format(Type count);
};

template <class T>
typename MetricTraits<T>::Type now() {
  return MetricTraits<T>::value_now();
}

} // namespace Stats
} // namespace Xyce

#endif // Xyce_N_UTL_StatMetricTraits_h
