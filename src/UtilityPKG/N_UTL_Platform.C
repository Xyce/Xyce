//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        :  Provide platform-specific information for metrics
//                   reporting and other uses
//                  
//                  
//
// Special Notes  : 
//                  
//
// Creator        : David Baur
//
// Creation Date  : 25 August 2014
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <string.h>
#include <sstream>
#include <cstddef>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#include <LmCons.h>
#endif

#if defined(HAVE_MALLOC_H) && defined(HAVE_MALLINFO)
#include <malloc.h>
#endif

#if defined(HAVE_PWD_H) && defined(HAVE_GETPWUID)
#include <sys/types.h>
#include <pwd.h>
#endif

#if defined(HAVE_UNISTD_H) && defined(HAVE_GETHOSTNAME)
#include <unistd.h>
#endif

#if defined(HAVE_UNISTD_H) && defined(HAVE_GETDOMAINNAME)
#include <unistd.h>
#endif

#if defined(HAVE_UNISTD_H) && defined(HAVE_SYSCONF)
#include <unistd.h>
#endif

#if defined(HAVE_SYS_UTSNAME_H) && defined(HAVE_UNAME)
#include <sys/utsname.h>
#endif

#if defined(HAVE__PROC_SELF_STAT)
#include <fstream>
#endif

namespace Xyce {

#if defined(HAVE_WINDOWS_H)
std::string windowsVersionName();
#endif

//-----------------------------------------------------------------------------
// Function      : get_heap_info
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:12 2015
//-----------------------------------------------------------------------------
void
get_heap_info(
   std::size_t &		heap_size,
   std::size_t &		largest_free)
{
  heap_size = 0;
  largest_free = 0;

#if defined(HAVE_MALLOC_H) && defined(HAVE_MALLINFO)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = (size_t) minfo.uordblks + (size_t) minfo.hblkhd;
  largest_free = (size_t) minfo.fordblks;
#endif
}


//-----------------------------------------------------------------------------
// Function      : get_available_memory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:06 2015
//-----------------------------------------------------------------------------
size_t get_available_memory()
{
#if defined(HAVE_SYSCONF) && defined(_SC_ACPHYS_PAGES)
  static size_t pagesize = getpagesize();
  size_t avail = sysconf(_SC_AVPHYS_PAGES);
  return avail * pagesize;

#else
  return 0;
#endif
}

//-----------------------------------------------------------------------------
// Function      : get_memory_info
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:09 2015
//-----------------------------------------------------------------------------
void
get_memory_info(
  size_t &		memory_usage,
  size_t &		faults)
{
  memory_usage = 0;
  faults = 0;

#if defined(HAVE__PROC_SELF_STAT)
  std::ifstream proc("/proc/self/stat", std::ios_base::in|std::ios_base::binary);
  if (proc) {

    std::string s;
    int i;
    for (i = 0; i < 11; ++i)
      proc >> s;

    proc >> faults;
    ++i;

    for (; i < 22; ++i)
      proc >> s;

    proc >> memory_usage;
    ++i;
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : hostname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:15 2015
//-----------------------------------------------------------------------------
std::string
hostname()
{
#if defined(HAVE_GETHOSTNAME)
  char buf[255];

  if (::gethostname(buf, sizeof(buf)) == 0)
    return std::string(buf);

#elif defined(HAVE_WINDOWS_H)
  TCHAR buf[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size;

  if (::GetComputerName(buf, &size) != 0)
    return std::string(buf);
#endif

  return "";
}


//-----------------------------------------------------------------------------
// Function      : domainname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:18 2015
//-----------------------------------------------------------------------------
std::string
domainname()
{
#if defined(HAVE_GETDOMAINNAME)
  char buf[255];

  if (::getdomainname(buf, sizeof(buf)) == 0)
    return std::string(buf);

#elif defined(HAVE_WINDOWS_H)
  TCHAR buf[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size = sizeof(buf);

  if (::GetComputerNameEx(ComputerNameDnsDomain, buf, &size) != 0)
    return std::string(buf);
#endif

  return "";
}


//-----------------------------------------------------------------------------
// Function      : username
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:23 2015
//-----------------------------------------------------------------------------
std::string
username()
{
#if defined(HAVE_GETPWUID)
  struct passwd *user_info = ::getpwuid(::geteuid());

  if (user_info)
    return user_info->pw_name;

#elif defined(HAVE_WINDOWS_H)
  TCHAR buf[UNLEN + 1];
  DWORD size = sizeof(buf);

  if (::GetUserName(buf, &size) != 0)
    return std::string(buf);

#endif
  return "unknown";
}


//-----------------------------------------------------------------------------
// Function      : hardware
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:40:27 2015
//-----------------------------------------------------------------------------
std::string
hardware()
{
#if defined(HAVE_UNAME)
  struct utsname	uts_name;

  ::uname(&uts_name);

  return uts_name.machine;

#elif defined(HAVE_WINDOWS_H)

  return ::getenv("PROCESSOR_ARCHITECTURE");

#else
  return "";
#endif
}


//-----------------------------------------------------------------------------
// Function      : osname
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:44 2015
//-----------------------------------------------------------------------------
std::string
osname()
{
#if defined(HAVE_UNAME)
  struct utsname	uts_name;

  if (::uname(&uts_name) != -1)
    return uts_name.sysname;
  else
    return "";

#elif defined(HAVE_WINDOWS_H)
  return "Windows";

#else
  return "";
#endif
}


//-----------------------------------------------------------------------------
// Function      : osversion
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:37 2015
//-----------------------------------------------------------------------------
std::string
osversion()
{
#if defined(HAVE_UNAME)
  struct utsname	uts_name;

  ::uname(&uts_name);

  return uts_name.release;

#elif defined(HAVE_WINDOWS_H)
  return windowsVersionName();

#else
  return "";
#endif
}


//-----------------------------------------------------------------------------
// Function      : pid
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:27 2015
//-----------------------------------------------------------------------------
int
pid()
{
#if defined(HAVE_GETPID)
  return ::getpid();
#else
  return 0;
#endif
}


//-----------------------------------------------------------------------------
// Function      : pgrp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:23 2015
//-----------------------------------------------------------------------------
int
pgrp()
{
#if defined(HAVE_GETPGRP)
  return ::getpgrp();
#else
  return 0;
#endif
}

#if defined(HAVE_WINDOWS_H)
typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);
typedef BOOL (WINAPI *PGPI)(DWORD, DWORD, DWORD, DWORD, PDWORD);
#define PRODUCT_PROFESSIONAL	0x00000030
#define VER_SUITE_WH_SERVER	0x00008000


//-----------------------------------------------------------------------------
// Function      : windowsVersionName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:39:18 2015
//-----------------------------------------------------------------------------
std::string
windowsVersionName()
{
  OSVERSIONINFOEX osvi;
  SYSTEM_INFO si;
  BOOL bOsVersionInfoEx;
  DWORD dwType;
  ZeroMemory(&si, sizeof(SYSTEM_INFO));
  ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
  bOsVersionInfoEx = GetVersionEx((OSVERSIONINFO*) &osvi);
  if (!bOsVersionInfoEx)
    return "";

  // Call GetNativeSystemInfo if supported or GetSystemInfo otherwise.
  PGNSI pGNSI = (PGNSI) GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")), "GetNativeSystemInfo");
  if(NULL != pGNSI)
    pGNSI(&si);
  else
    GetSystemInfo(&si); // Check for unsupported OS

  if (VER_PLATFORM_WIN32_NT != osvi.dwPlatformId || osvi.dwMajorVersion <= 4 ) {
    return "";
  }

  std::stringstream os;
  os << "Microsoft "; // Test for the specific product. if ( osvi.dwMajorVersion == 6 )
  {
    if (osvi.dwMinorVersion == 0 )
    {
      if (osvi.wProductType == VER_NT_WORKSTATION )
        os << "Windows Vista ";
      else
        os << "Windows Server 2008 ";
    }
    if (osvi.dwMinorVersion == 1 )
    {
      if (osvi.wProductType == VER_NT_WORKSTATION )
        os << "Windows 7 ";
      else
        os << "Windows Server 2008 R2 ";
    }

    PGPI pGPI = (PGPI) GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")), "GetProductInfo");
    pGPI(osvi.dwMajorVersion, osvi.dwMinorVersion, 0, 0, &dwType);
    switch (dwType)
    {
      case PRODUCT_ULTIMATE:
        os << "Ultimate Edition";
        break;
      case PRODUCT_PROFESSIONAL:
        os << "Professional";
        break;
      case PRODUCT_HOME_PREMIUM:
        os << "Home Premium Edition";
        break;
      case PRODUCT_HOME_BASIC:
        os << "Home Basic Edition";
        break;
      case PRODUCT_ENTERPRISE:
        os << "Enterprise Edition";
        break;
      case PRODUCT_BUSINESS:
        os << "Business Edition";
        break;
      case PRODUCT_STARTER:
        os << "Starter Edition";
        break;
      case PRODUCT_CLUSTER_SERVER:
        os << "Cluster Server Edition";
        break;
      case PRODUCT_DATACENTER_SERVER:
        os << "Datacenter Edition";
        break;
      case PRODUCT_DATACENTER_SERVER_CORE:
        os << "Datacenter Edition (core installation)";
        break;
      case PRODUCT_ENTERPRISE_SERVER:
        os << "Enterprise Edition";
        break;
      case PRODUCT_ENTERPRISE_SERVER_CORE:
        os << "Enterprise Edition (core installation)";
        break;
      case PRODUCT_ENTERPRISE_SERVER_IA64:
        os << "Enterprise Edition for Itanium-based Systems";
        break;
      case PRODUCT_SMALLBUSINESS_SERVER:
        os << "Small Business Server";
        break;
      case PRODUCT_SMALLBUSINESS_SERVER_PREMIUM:
        os << "Small Business Server Premium Edition";
        break;
      case PRODUCT_STANDARD_SERVER:
        os << "Standard Edition";
        break;
      case PRODUCT_STANDARD_SERVER_CORE:
        os << "Standard Edition (core installation)";
        break;
      case PRODUCT_WEB_SERVER:
        os << "Web Server Edition";
        break;
    }
  }

  if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 2)
  {
    if( GetSystemMetrics(SM_SERVERR2) )
      os <<  "Windows Server 2003 R2, ";
    else if ( osvi.wSuiteMask & VER_SUITE_STORAGE_SERVER )
      os <<  "Windows Storage Server 2003";
    else if ( osvi.wSuiteMask & VER_SUITE_WH_SERVER )
      os <<  "Windows Home Server";
    else if( osvi.wProductType == VER_NT_WORKSTATION && si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64)
      os <<  "Windows XP Professional x64 Edition";
    else
      os << "Windows Server 2003, ";  // Test for the server type.

    if (osvi.wProductType != VER_NT_WORKSTATION)
    {
      if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_IA64)
      {
        if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
          os <<  "Datacenter Edition for Itanium-based Systems";
        else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
          os <<  "Enterprise Edition for Itanium-based Systems";
      }
      else if  (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64)
      {
        if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
          os <<  "Datacenter x64 Edition";
        else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
          os <<  "Enterprise x64 Edition";
        else
          os <<  "Standard x64 Edition";
      }
      else
      {
        if (osvi.wSuiteMask & VER_SUITE_COMPUTE_SERVER )
          os <<  "Compute Cluster Edition";
        else if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
          os <<  "Datacenter Edition";
        else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
          os <<  "Enterprise Edition";
        else if ( osvi.wSuiteMask & VER_SUITE_BLADE )
          os <<  "Web Edition";
        else
          os <<  "Standard Edition";
      }
    }
  }
  if ( osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1 )
  {
    os << "Windows XP ";
    if( osvi.wSuiteMask & VER_SUITE_PERSONAL )
      os <<  "Home Edition";
    else os <<  "Professional";
  }
  if ( osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0 )
  {
    os << "Windows 2000 ";
    if ( osvi.wProductType == VER_NT_WORKSTATION )
      os <<  "Professional";
    else
    {
      if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
        os <<  "Datacenter Server";
      else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
        os <<  "Advanced Server";
      else
        os <<  "Server";
    }
  } // Include service pack (if any) and build number. if(wcslen(osvi.szCSDVersion) > 0) {

  os << " " << osvi.szCSDVersion;

  os << " (build " << osvi.dwBuildNumber << L")"; if ( osvi.dwMajorVersion >= 6 ) {
    if ( si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64 )
      os <<  ", 64-bit";
    else if (si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_INTEL )
      os << ", 32-bit";
  }

  return os.str();
}
#endif

} // namespace Xyce

