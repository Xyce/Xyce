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
#include <N_ANP_SweepParam.h>

namespace Xyce {
namespace Analysis {

typedef double RealT;
//-------------------------------------------------------------------------
// Class         : ROL
// Purpose       : ROL analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
// Revised: Timur Takhtaganov, 06/03/2015
//
class ROL : public AnalysisBase
{
  template <class RealT> 
  friend class EqualityConstraint_ROL_DC;
  template <class RealT> 
  friend class EqualityConstraint_ROL_DC_UQ; 
  
public:
  ROL(
      AnalysisManager &analysis_manager, 
      Nonlinear::Manager &nonlinear_manager,
      Loader::Loader &loader, 
      Linear::System & linear_system,
      Topo::Topology & topology,
      IO::InitialConditionsManager & initial_conditions_manager); 
   
  virtual ~ROL();

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setTimeIntegratorOptions(const Util::OptionBlock &option_block);

  void setTIAParams(const TimeIntg::TIAParams &tia_params)
  {
    tiaParams_ = tia_params;
  } 
  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  bool getDCOPFlag() const // override
  {
    return true;
  }



protected:
  void finalExpressionBasedSetup() {};
  bool doRun();
  bool doInit();
  bool doLoopProcess();
  bool runROLAnalysis(); 
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doHandlePredictor();
  bool doFinish();

  // virtual bool doProcessSuccessfulDCOPStep() { // override
  //   return true;
  // }

  // virtual bool doProcessFailedDCOPStep() { // override
  //   return true;
  // }

  bool doAllocations(int nc, int nz);
  bool doFree();
  std::vector<Linear::Vector *> solutionPtrVector_;
  std::vector<Linear::Vector *> statePtrVector_;
  std::vector<Linear::Vector *> constraintPtrVector_;
  std::vector<Linear::Vector *> jvecPtrVector_;
  std::vector<Linear::Vector *> testPtrVector_;
  std::vector<Linear::Vector *> mydfdpPtrVector_;
  std::vector<Linear::Vector *> mydqdpPtrVector_;
  std::vector<Linear::Vector *> mydbdpPtrVector_;
  std::vector<Linear::Vector *> mysensRHSPtrVector_;

public:
  // Two Level specific
  bool twoLevelStep(); 
  void setSweepValue(int step); 

private:
  void initializeSolution_();
  void takeStep_();

  std::vector<int>      rolSweepFailures_; // TT
  
private:
  AnalysisManager &                     analysisManager_;
  Nonlinear::Manager &                  nonlinearManager_; // TT
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  Linear::System &                      linearSystem_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  TimeIntg::TIAParams                   tiaParams_;
  SweepVector                           stepSweepVector_;
  int                                   stepLoopSize_;
  bool                                  rolLoopInitialized_; // TT
  std::vector<std::string>              paramNameVec_; // TT: vector of optimization parameters
  int                                   numParams_; // TT: number of optimization parameters
  std::vector<std::string>              uncertainParams_; // TT: vector of parameters with uncertainty
  int                                   numSensParams_;  // number of sensitivity parameters returned from enableSensitivity function call

};


bool registerROLFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_ROL_h
