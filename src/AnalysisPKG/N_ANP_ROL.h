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

//-----------------------------------------------------------------------------
//
// Purpose        : ROL analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_ROL_h
#define Xyce_N_ANP_ROL_h

#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_AnalysisManager.h> 
#include <N_ANP_RegisterAnalysis.h>
#include <N_IO_OptionBlock.h>

namespace Xyce {
namespace Nonlinear {

template <typename ScalarT>
class objectiveFunctionData;
}
}

namespace Xyce {
namespace Analysis {

class ROL_Objective;

typedef double RealT;
//-------------------------------------------------------------------------
// Class         : ROL
// Purpose       : ROL analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
// Revised: Timur Takhtaganov, 06/03/2015
// Revised: Heidi Thornquist, 02/22/2022
//
class ROL : public AnalysisBase
{
public:
  ROL(
      AnalysisManager &analysis_manager, 
      Nonlinear::Manager &nonlinear_manager,
      Loader::Loader &loader, 
      Linear::System & linear_system,
      Topo::Topology & topology,
      IO::InitialConditionsManager & initial_conditions_manager); 
   
  virtual ~ROL();

  void setTIAParams(const TimeIntg::TIAParams &tia_params)
  {
    tiaParams_ = tia_params;
  } 

  const TimeIntg::TIAParams &getTIAParams() const
  {
    return tiaParams_;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    return tiaParams_;
  }

  // Method to set ROL options
  bool setROLOptions(const Util::OptionBlock & option_block);
  bool setROLObjectives(const std::vector<Util::OptionBlock>& option_block);
  bool setOptParams();

  // Method to set ROL DC description
  bool setROLDCSweep(const std::vector<Util::OptionBlock>& option_block);
  bool setROLDataOptionBlock(const std::vector<Util::OptionBlock>& option_block);

  // Method to set non-HB linear solver / preconditioning options (needed for .STEP)
  bool setLinSol(const Util::OptionBlock & option_block);

  // Method to set time integrator options (needed for initial condition / startup periods)
  bool setTimeInt(const Util::OptionBlock & option_block);

  bool getDCOPFlag() const; 

protected:
  void finalExpressionBasedSetup() {}
  bool doRun();
  bool doInit();
  bool doLoopProcess();
  bool doProcessSuccessfulStep() { return false; }
  bool doProcessFailedStep() { return false; }
  bool doHandlePredictor() { return true; }
  bool doFinish() { return true; }

public:
  // Two Level specific
  bool twoLevelStep(); 

private:
  AnalysisManager &                     analysisManager_;
  Nonlinear::Manager &                  nonlinearManager_; // TT
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  Linear::System &                      linearSystem_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  TimeIntg::TIAParams                   tiaParams_;
  int                                   stepLoopSize_;
  bool                                  sensFlag_;
  std::vector<std::string>              paramNameVec_; // TT: vector of optimization parameters
  int                                   numParams_; // TT: number of optimization parameters
  int                                   numSensParams_;  // number of sensitivity parameters returned from enableSensitivity function call
  std::string                           paramFile_;  // Name of file with parameters and bounds
  std::string                           rolParamFile_;  // Name of file with parameters and bounds
  std::string                           outputFile_;  // Name of file containing ROL output

  AnalysisBase *                        currentAnalysisObject_;

  std::vector<ROL_Objective>            rolDCObjVec_;

  std::vector<RealT>                    zInitValue_;
  std::vector<RealT>                    zLowerBoundVector_;
  std::vector<RealT>                    zUpperBoundVector_;

  std::vector<RealT>                    objectiveVec_;
  std::vector<RealT>                    dOdpVec_;
  std::vector<RealT>                    dOdpAdjVec_;
  std::vector<RealT>                    scaled_dOdpVec_;
  std::vector<RealT>                    scaled_dOdpAdjVec_;

  Util::OptionBlock                     saved_lsOB_;  // Linear solver options
  Util::OptionBlock                     saved_timeIntOB_;  // Time integrator options
  std::vector<Util::OptionBlock>        saved_sweepOB_;  // DCSweep options
  std::vector<Util::OptionBlock>        saved_dataOB_;  // DCSweep data
  std::vector<Util::OptionBlock>        saved_rolObjOB_;  // ROL objectives
};


//-------------------------------------------------------------------------
// Class         : ROL_Objective
// Purpose       : Describe ROL objective
//-------------------------------------------------------------------------
class ROL_Objective
{
  public:
    std::string objType_;                // Internal objective type, not supported by SENS (ex. data fitting)
    std::string objTag_;                 // ROL objective tag, useful for combining objectives
    std::vector<std::string> objArgs_;   // Objective arguments
 
  ROL_Objective() {}

  virtual ~ROL_Objective() {}
};

template<class Real>
class ROL_Objective_Arg
{
  public:
    int objIndex;  // either the solution vector index or data column index
    std::vector< std::vector< double > > *  objDataPtr;  // point to input data, don't delete!
   
  ROL_Objective_Arg()
  : objIndex(-1),
    objDataPtr(0)
  {}
    
  virtual ~ROL_Objective_Arg() {}
};

bool registerROLFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_ROL_h
