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
// Purpose        : This file contains the sensitivity class.   It mostly
//                  manages the calculations of direct (and possibly later,
//                  adjoint) sensitivities.
//
// Special Notes  : The main reason that this class is derived from
//                  N_NLS_NonLinearSolver is that this class needs to
//                  do a series of linear solves, using the jacobian
//                  matrix.  This seemed similar enough to the requirements
//                  of a nonlinear solver to have one derived off the other.
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/15/2015
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_SensitivityResiduals_h
#define Xyce_N_NLS_SensitivityResiduals_h

#include<vector>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_AnalysisManager.h>
#include <N_TIA_DataStore.h>
#include <N_LOA_NonlinearEquationLoader.h>


namespace Xyce {
namespace Nonlinear {

enum sensDiffMode
{
  SENS_FWD,
  SENS_REV,
  SENS_CNT,
  NUM_DIFF_MODES
};

//-----------------------------------------------------------------------------
// Function      : getSetParamName
//
// Purpose       : This function strips curly braces off the parameter name.
//
// Special Notes : The reason this is necessary:  This was driven by the use of 
//                 global parameters for sensitivity calculations.  It is common
//                 for global_params to be specified with curly braces surrounding 
//                 them.  However, the parameter processing functions in the device
//                 package don't fully support this.  Some functions understand them
//                 and some do not.  The "correct" way to fix this is to modify those
//                 functions, but those functions have many complex layers and
//                 this fix was quicker.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 10/8/2018
//-----------------------------------------------------------------------------
inline void getSetParamName(const std::string & origName, std::string & newName)
{
  int size=origName.size(); 
  newName = origName;
  if (size > 2) 
  { 
    if (origName[0] == '{' && origName[size-1] == '}') 
    {
      newName.resize(size-2); // copy only works into properly allocated objects
      std::string::const_iterator iter0=origName.begin(); iter0++; 
      std::string::const_iterator iter1=origName.end(); iter1--;
      std::copy(iter0,iter1,newName.begin());
    } 
  }
}

//-----------------------------------------------------------------------------
bool slowNumericalDerivatives(int iparam,
  std::vector<std::string>::const_iterator & iterParam, 
  int difference, 
  double sqrtEta_, std::string & netlistFilename_,
  TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_);

//-----------------------------------------------------------------------------
bool loadSensitivityResiduals(int difference, 
    bool forceFD_, bool forceDeviceFD_, bool forceAnalytic_, 
    double sqrtEta_,  std::string & netlistFilename_,
  TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    );

//-----------------------------------------------------------------------------
bool setupOriginalParams ( TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    );

//-----------------------------------------------------------------------------
bool testForAnalyticDerivatives ( 
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    );

//-----------------------------------------------------------------------------
bool computeSparseIndices ( const int iparam,
  TimeIntg::DataStore & ds,
  const std::vector<int>& FindicesVec,
  const std::vector<int>& QindicesVec,
  const std::vector<int>& BindicesVec
    );

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_SensitivityResiduals_h

