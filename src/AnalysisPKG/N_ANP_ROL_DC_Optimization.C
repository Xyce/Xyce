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
// Purpose        : ROL DC analysis classes
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 01/24/08
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_ROL_DC_Optimization.h>
#include <N_ANP_ROL.h>
#include <N_ANP_SweepParamFreeFunctions.h>

#include <N_NLS_ObjectiveFunctions.h>
#include <N_IO_CmdParse.h>
#include <N_UTL_ExpressionData.h>

namespace Xyce {
namespace Analysis {

#ifdef Xyce_ROL

//-----------------------------------------------------------------------------
// TT: this function gets called by the value() function in EqualityConstraint class;
//     it updates sweep parameter so that correct source term is loaded into RHS vector
//-----------------------------------------------------------------------------
void ROL_DC::setSweepValue(int step)
{ 
  bool reset = updateSweepParams(loader_, step, dcSweepVector_.begin(), dcSweepVector_.end(), false);
}

bool ROL_DC::doInit()
{ 
  bool ret = DCSweep::doInit();
 
  analysisManager_.setAnalysisMode(ANP_MODE_DC_SWEEP); 
  //nonlinearManager_.resetAll(Nonlinear::DC_SWEEP);
  
  return ret;
}

bool ROL_DC::doProcessSuccessfulStep()
{ 
  bool ret = DCSweep::doProcessSuccessfulStep();

  int currentStep = outputManagerAdapter_.getDCAnalysisStepNumber();

  TimeIntg::DataStore & ds = *(analysisManager_.getDataStore());
  *(solutionPtrVector_[currentStep]) = *(ds.currSolutionPtr);
  double norm = 0.0;
  //ds.currSolutionPtr->lpNorm(2,&norm);
  //std::cout << "Solution norm " << currentStep << " : " << norm << std::endl;    

  return ret;
}

bool ROL_DC::doAllocations(int nc, int nz)
{
  stepLoopSize_ = nc;
  numParams_ = nz;

  // Allocate space for solution and simulation space vectors
  solutionPtrVector_.resize(nc);
  statePtrVector_.resize(nc);
  constraintPtrVector_.resize(nc);

  for (int i=0;i<nc;i++)
  {
    solutionPtrVector_[i]   = linearSystem_.builder().createVector();
    statePtrVector_[i]      = linearSystem_.builder().createVector();
    constraintPtrVector_[i] = linearSystem_.builder().createVector();
  }
  // Allocate space for sensitivity vectors
  mydfdpPtrVector_.resize(nz);
  mydqdpPtrVector_.resize(nz);
  mydbdpPtrVector_.resize(nz);
  mysensRHSPtrVector_.resize(nz);
  for (int i=0;i<nz;i++)
  {
    mydfdpPtrVector_[i] = linearSystem_.builder().createVector();
    mydqdpPtrVector_[i] = linearSystem_.builder().createVector();
    mydbdpPtrVector_[i] = linearSystem_.builder().createVector();
    mysensRHSPtrVector_[i] = linearSystem_.builder().createVector();
  }

  return true;
}

bool ROL_DC::setAnalysisParams(const std::vector<Util::OptionBlock>& OB)
{
  std::vector<Util::OptionBlock>::const_iterator it = OB.begin();
  std::vector<Util::OptionBlock>::const_iterator end = OB.end();
  for( ; it != end; ++it )
    DCSweep::setAnalysisParams( *it );

  return true;
}

bool ROL_DC::createObjectives(const std::vector<ROL_Objective>& objVec)
{
  // Save the objectives from the netlist
  objVec_ = objVec;

  // Create objective function data vector
  // NOTE:  Determine if objective argument is provided by input data instead of solution data
  std::string objType_;
  std::vector<ROL_Objective_Arg<RealT> > objArgVec;

  for ( int iObj = 0; iObj < objVec_.size(); ++iObj )
  {
    objType_ = objVec_[iObj].objType_;
    objArgVec.resize( objVec_[iObj].objArgs_.size() );
    //std::cout << "objType_ = " << objType_ << std::endl;

    for (int i=0; i<objVec_[iObj].objArgs_.size(); ++i)
    { 
      std::string objArg = objVec_[0].objArgs_[i];

      // Remove brackets if it is an expression
      if( (objArg[0]=='{') && (objArg[objArg.size()-1]=='}') )
      {
        objArg = std::string( objArg.begin()+1, objArg.end()-1 );
      }

      //std::cout << "Argument: " << objArg << std::endl;

      // Check if this argument is provided by the data tables
      bool isFound = false;
      std::map< std::string, std::vector< std::vector<double> > >& dataTableMap = getDataTablesMap();
      std::map< std::string, std::vector<std::string> >& dataNames = getDataNamesMap();
      std::map< std::string, std::vector<std::string> >::iterator it_dN = dataNames.begin();
      for ( ; (it_dN != dataNames.end()) && !isFound; ++it_dN)
      {
        std::vector< std::string >::iterator it_dN2 = it_dN->second.begin();
        for ( int index=0; (it_dN2 != it_dN->second.end()) && !isFound; ++it_dN2, ++index )
        {
          if ( *it_dN2 == objArg )
          {
            isFound = true;
            objArgVec[i].objIndex = index;
            objArgVec[i].objDataPtr = &dataTableMap[ it_dN->first ];
            //std::cout << "Found objective argument " << objArg << " index " << index << " in .DATA statement!" << std::endl;
          }
        }
      }

      // If the argument is not provided by the data tables, then create an objectiveFunctionData object.
      if (!isFound)
      {
        std::vector<Xyce::Nonlinear::objectiveFunctionData<RealT> *> objFuncDataVec;  

        // Create a new objective function
        Xyce::Nonlinear::objectiveFunctionData<double> ofData;
        ofData.objFuncString = objVec_[0].objArgs_[i];
        objFuncDataVec.push_back(&ofData);
  
        // Resolve expression to get GID
        IO::OutputMgr & output_manager = outputManagerAdapter_.getOutputManager();
        Parallel::Machine comm =  analysisManager_.getPDSManager()->getPDSComm()->comm();

        Xyce::Nonlinear::setupObjectiveFunctions(
          analysisManager_.getExpressionGroup(),
            objFuncDataVec, output_manager, linearSystem_, analysisManager_.getCommandLine() );

        Xyce::Nonlinear::setupObjectiveFuncGIDs ( objFuncDataVec, comm, topology_, output_manager );

        for ( int index=0; index<ofData.expVarGIDs.size(); ++index )
        {
          objArgVec[i].objIndex = ofData.expVarGIDs[index];
          //std::cout << "objFuncDataVec[i].expVarGIDs = " << ofData.expVarGIDs[index] << std::endl;
        }
      }
    }
  }

  // Create a full space objective (as well as a penalty when solving the
  // amplifier problem).
  if (objType_ == "ERROR")
    obj_ = ::ROL::makePtr<Objective_DC_L2Norm<RealT>>(1.0e-4, stepLoopSize_, numParams_, objArgVec);
  else if (objType_ == "AMP")
    obj_ = ::ROL::makePtr<Objective_DC_AMP<RealT>>(stepLoopSize_, numParams_);
  else
    Report::UserError0() << "ROL: Objective type " << objType_ << " is not recognized";

  if (obj_ != Teuchos::null)
    return true;
  else
    return false;
}


bool ROL_DC::doFree()
{
  for (int i=0;i<stepLoopSize_;i++)
  {
    delete solutionPtrVector_[i];
    solutionPtrVector_[i] = 0;
    delete statePtrVector_[i];
    statePtrVector_[i] = 0;
    delete constraintPtrVector_[i];
    constraintPtrVector_[i] = 0;
  }
  solutionPtrVector_.clear();
  statePtrVector_.clear();
  constraintPtrVector_.clear();

  for (int i=0;i<numParams_;i++)
  {
    delete mydfdpPtrVector_[i];
    mydfdpPtrVector_[i] = 0;
    delete mydqdpPtrVector_[i];
    mydqdpPtrVector_[i] = 0;
    delete mydbdpPtrVector_[i];
    mydbdpPtrVector_[i] = 0;
    delete mysensRHSPtrVector_[i];
    mysensRHSPtrVector_[i] = 0;
  }
  mydfdpPtrVector_.clear();
  mydqdpPtrVector_.clear();
  mydbdpPtrVector_.clear();
  mysensRHSPtrVector_.clear();

  return true;
}

#endif

} // namespace Analysis
} // namespace Xyce
