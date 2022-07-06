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

} // namespace Xyce

#endif // Xyce_N_UTL_Platform_h
