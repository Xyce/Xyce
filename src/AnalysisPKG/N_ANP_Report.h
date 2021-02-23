//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Report_H
#define Xyce_N_ANP_Report_H

#include <N_ANP_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_LogStream.h>

namespace Xyce {
namespace Analysis {

struct UserWarning : public Report::UserWarning
{
  UserWarning(const AnalysisBase &analysis_base);
};

struct UserWarning0 : public Report::UserWarning0
{
  UserWarning0(const AnalysisBase &analysis_base);
};

struct UserFatal : public Report::UserFatal
{
  UserFatal(const AnalysisBase &analysis_base);
};

struct UserFatal0 : public Report::UserFatal0
{
  UserFatal0(const AnalysisBase &analysis_base);
};

struct DevelFatal : public Report::DevelFatal
{
  DevelFatal(const AnalysisBase &analysis_base);
};

struct DevelFatal0 : public Report::DevelFatal0
{
  DevelFatal0(const AnalysisBase &analysis_base);
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Report_H
