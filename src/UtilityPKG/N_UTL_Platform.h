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

#ifndef Xyce_N_UTL_Platform_h
#define Xyce_N_UTL_Platform_h

#include <iosfwd>

namespace Xyce {

///
/// @addtogroup EnvDetail
/// @{
///

//-----------------------------------------------------------------------------
// Function      : hostname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:33:24 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>hostname</b> returns the hostname of the host running the
/// application.
///
/// @return			a <b>String</b> value of the host name obtained from
/// 			the operating system.
/// 
std::string hostname();


//-----------------------------------------------------------------------------
// Function      : domainname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:33:50 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>domainname</b> returns the domainname of the domain running the
/// application.
///
/// @return			a <b>String</b> value of the domain name obtained from
/// 			the operating system.
/// 
std::string domainname();


//-----------------------------------------------------------------------------
// Function      : username
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:33:54 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>username</b> returns the username of the user running the
/// application.
///
/// @return			a <b>String</b> value of the username obtained from
/// 			the operating system.
/// 
std::string username();


//-----------------------------------------------------------------------------
// Function      : hardware
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:00 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>hardware</b> returns the hardware type of the host running the
/// application.
///
/// @return			a <b>String</b> value of the <b>machine</b>
/// 			field of the <b>uname</b> system call or equivalent
/// 			obtained from the operating system.
/// 
std::string hardware();


//-----------------------------------------------------------------------------
// Function      : osname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:04 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>osname</b> returns the operating system nameof the host running the
/// application.
///
/// @return			a <b>String</b> value of the <b>sysname</b>
/// 			field of the <b>uname</b> system call or equivalent
/// 			obtained from the operating system.
/// 
std::string osname();


//-----------------------------------------------------------------------------
// Function      : osversion
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:07 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>osversion</b> returns the hardware type of the host running the
/// application.
///
/// @return			a <b>String</b> value of the <b>release</b>
/// 			field of the <b>uname</b> system call or equivalent
/// 			obtained from the operating system.
/// 
std::string osversion();


//-----------------------------------------------------------------------------
// Function      : pid
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:12 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>pid</b> returns the process id of the process running the
/// application.
///
/// @return			a <b>int</b> value of the process id obtained from
/// 			the operating system.
/// 
int pid();


//-----------------------------------------------------------------------------
// Function      : pgrp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:18 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>pgrp</b> returns the process group id of the process running
/// the application.
///
/// @return			a <b>int</b> value of the process group id obtained from
/// 			the operating system.
/// 
int pgrp();


//-----------------------------------------------------------------------------
// Function      : get_heap_info
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:22 2015
//-----------------------------------------------------------------------------
///
/// @brief Member function <b>get_heap_info</b> returns the amount of heap
/// memory used in bytes and the largest free block of memory in bytes.
///
/// @param heap_size		a <b>size_t</b> returns the amount of heap
/// memory used in bytes.
///
/// @param largest_free		a <b>size_t</b> returns the largest free block
/// of memory.
///
/// 
void get_heap_info(size_t &heap_size, size_t &largest_free);


//-----------------------------------------------------------------------------
// Function      : get_heap_usage
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:27 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>get_heap_usage</b> returns the number of bytes used by the heap.
///
/// @return			a <b>size_t</b> value of the number of bytes used by
/// 			the heap.
/// 
inline size_t get_heap_usage() {
  size_t heap_size;
  size_t largest_free;
  get_heap_info(heap_size, largest_free);

  return heap_size;
}

//-----------------------------------------------------------------------------
// Function      : get_available_memory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:31 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>get_available_memory</b> returns an estimation of the amount of memory available to the process.
///
/// @return			a <b>size_t</b> value of the number of bytes available to the process.
/// 
size_t get_available_memory();


//-----------------------------------------------------------------------------
// Function      : get_memory_info
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:34:36 2015
//-----------------------------------------------------------------------------
///
/// @ingroup EnvRuntimeInformationDetail
/// @brief Function <b>get_memory_info</b> returns the total memory usage of the
/// process and the number of page faults accumulated by the process.
///
/// @param memory_usage		a <b>size_t</b> reference to receive the number of
/// 			bytes currently used by the process.
///
/// @param faults		a <b>size_t</b> reference to treceive the number of
/// 			page faults incurred by the process.
///
/// 
void get_memory_info(size_t &memory_usage, size_t &faults);

///
/// @}
///

} // namespace Xyce

#endif // Xyce_N_UTL_Platform_h
