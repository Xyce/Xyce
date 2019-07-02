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
// Purpose        : Outputter class to handle delivering ".print"-like
//                  data to an external program via callbacks
//
// Special Notes  :
//
// Creator        : Tom Russo
//
// Creation Date  : 20 Feb 2018
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>
#include <N_IO_OutputterExternal.h>
#include <N_IO_OutputterLocal.h>
#include <N_IO_ExtOutWrapper.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_PDS_MPI.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_NetlistLocation.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 20 Feb 2018
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : generateOpsListFromParamsList
// Purpose       : Given a ParamList representing variables to be output,
//                 create all the ops to acces those variables.
// Special Notes :
// Scope         : file-local
// Creator       : Tom Russo
// Creation Date : 20 Feb 2018
//-----------------------------------------------------------------------------
/// Generate ops from ParamList specifying output variables
///
/// This is similar to what's done in N_IO_OutputterLocal's "fixupColumns"
/// but instead of using PrintParameters and all the cruft that goes with
/// them, we just directly use the params lists.  We never expand complex
/// types, nor do we allow any time scale factor other than 1.0.
///
/// @param[in] comm The communicator for this run
/// @param[in] theOpBuilderMgr  Reference to Op Builder manager
/// @param[in] outputParamList  ParamList containing representation of output vars
/// @param[out] theOpList       Reference to an OpList, intto which we will insrert Ops
/// @param[in] NetlistLocation  A "netlist location" to be used in error reporting.  In this outputter, this will generally be faked out, since we aren't producing output based on netlist print statements.
///

void generateOpsListFromParamsList(Parallel::Machine comm,
                                   const Util::Op::BuilderManager &theOpBuilderMgr,
                                   Util::ParamList &outputParamList,
                                   Util::Op::OpList &theOpList,
                                   NetlistLocation & theLocation)
{
  createOps(comm,theOpBuilderMgr,false,1.0,theLocation,
            outputParamList.begin(),outputParamList.end(),
            std::back_inserter(theOpList));
}

//-----------------------------------------------------------------------------
// Function      : getNamesFromOpsList
// Purpose       : Given an ops list, create a vector of strings that could
//                 be used to create a header.
// Special Notes :
// Scope         : file-local
// Creator       : Tom Russo
// Creation Date : 20 Feb 2018
//-----------------------------------------------------------------------------
void getNamesFromOpsList(Util::Op::OpList & theOpList,
                         std::vector<std::string> & theNames)
{
  for (Util::Op::OpList::const_iterator it = theOpList.begin();
       it != theOpList.end();
       it++)
  {
    theNames.push_back((*it)->getName());
  }
}

}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::OutputterExternal
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 20 Feb 2018
//-----------------------------------------------------------------------------
OutputterExternal::OutputterExternal(Parallel::Machine comm,
                                     OutputMgr &output_manager,
                                     ExternalOutputWrapper * outputWrapper)
  : outputManager_(output_manager),
    theOutputWrapper_(outputWrapper),
    currentStep_(0),
    numberOfSteps_(0),
    index_(0),
    initialized_(false)
{

  // create a fake "NetlistLocation" using wrapped object's name
  NetlistLocation fakeLocation(theOutputWrapper_->getName(),0);

  generateOpsListFromParamsList(comm,
                                outputManager_.getOpBuilderManager(),
                                theOutputWrapper_->getParamList(),
                                opList_,
                                fakeLocation);
  getNamesFromOpsList(opList_,fieldNames_);
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::~OutputterExternal
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
OutputterExternal::~OutputterExternal()
{
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doOutputTime
// Purpose       : Output the current data at a time point
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
void
OutputterExternal::doOutputTime(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVec,
  const Linear::Vector &  stateVec,
  const Linear::Vector &  storeVec,
  const Linear::Vector &  leadCurrentVec,
  const Linear::Vector &  junctionVoltageVec)
{
  if (Parallel::rank(comm) == 0 && !initialized_)
  {
    initialized_=true;
    theOutputWrapper_->outputFieldNames(fieldNames_);
  }
  std::vector<complex> complexResultList;
  getValues(comm, opList_, Util::Op::OpData(index_, &solnVec, 0,
                                            &stateVec, &storeVec, 0,
                                            &leadCurrentVec, 0,
                                            &junctionVoltageVec),
            complexResultList);

  std::vector<double> doubleResultList(complexResultList.size());
  for (int i=0; i< complexResultList.size(); i++)
    doubleResultList[i]=complexResultList[i].real();

  // initialized_ will only ever get set on proc 0, so this  assures
  // we only output once
  if (initialized_)
    theOutputWrapper_->outputReal(doubleResultList);

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doOutputFrequency
// Purpose       : Output the current data for a purely FD run (e.g. AC)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 22 Feb 2018
//-----------------------------------------------------------------------------
void
OutputterExternal::doOutputFrequency(
  Parallel::Machine     comm,
  double                frequency,
  double                fStart,
  double                fStop,
  const Linear::Vector & realSolutionVector,
  const Linear::Vector & imaginarySolutionVector,
  const Util::Op::RFparamsData & RFparams)
{
  if (Parallel::rank(comm) == 0 && !initialized_)
  {
    initialized_=true;
    theOutputWrapper_->outputFieldNames(fieldNames_);
  }
  std::vector<complex> complexResultList;
  getValues(comm, opList_, Util::Op::OpData(index_, &realSolutionVector,
                                            &imaginarySolutionVector,
                                            0,0,0),
            complexResultList);

  if (initialized_)
    theOutputWrapper_->outputComplex(complexResultList);

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doOutputHB_FD
// Purpose       : Output the current frequency domain data for an HB run
// Special Notes : We will require the user to have specified separate HB_TD
//                 and HB_FD output interfaces, and this function will output
//                 only the frequency domain data appropriate for the HB_FD
//                 interface.
//
//                 If the user specifies "HB" as the output type for their
//                 interface object, we'll never be called. Only HB_FD and
//                 HB_TD will be accepted.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 27 Mar 2018
//-----------------------------------------------------------------------------
void OutputterExternal::doOutputHB_FD(
   Parallel::Machine             comm,
   const std::vector<double> &   freqPoints,
   const Linear::BlockVector &     freqDomainSolutionVecReal,
   const Linear::BlockVector &     freqDomainSolutionVecImaginary,
   const Linear::BlockVector &     freqDomainLeadCurrentVecReal,
   const Linear::BlockVector &     freqDomainLeadCurrentVecImaginary,
   const Linear::BlockVector &     freqDomainJunctionVoltageVecReal,
   const Linear::BlockVector &     freqDomainJunctionVoltageVecImaginary)
{
  if (Parallel::rank(comm) == 0 && !initialized_)
  {
    initialized_ = true;
    theOutputWrapper_->outputFieldNames(fieldNames_);
  }

  if (theOutputWrapper_->getOutputType() == OutputType::HB_FD)
  {
    int fdblockCount = freqDomainSolutionVecReal.blockCount();
    index_ = 0;
    for (int iblock = 0; iblock < fdblockCount; ++iblock)
    {
      outputManager_.setCircuitFrequency(freqPoints[iblock]);

      Linear::Vector * real_solution_vector =
        &freqDomainSolutionVecReal.block(iblock);
      Linear::Vector * imaginary_solution_vector =
        &freqDomainSolutionVecImaginary.block(iblock);
      Linear::Vector * real_lead_current_vector =
        &freqDomainLeadCurrentVecReal.block(iblock);
      Linear::Vector * imaginary_lead_current_vector =
        &freqDomainLeadCurrentVecImaginary.block(iblock);
      Linear::Vector * real_junction_voltage_vector =
        &freqDomainJunctionVoltageVecReal.block(iblock);
      Linear::Vector * imaginary_junction_voltage_vector =
        &freqDomainJunctionVoltageVecImaginary.block(iblock);

      {
        // Fourier coefficient output
        std::vector<complex> complexResultList;
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        getValues(comm, opList_,
                  Util::Op::OpData(index_,
                                   real_solution_vector,
                                   imaginary_solution_vector,
                                   0, 0, 0,
                                   real_lead_current_vector,
                                   imaginary_lead_current_vector,
                                   real_junction_voltage_vector,
                                   imaginary_junction_voltage_vector ),
                  complexResultList);
        theOutputWrapper_->outputComplex(complexResultList);
      }
      ++index_;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doOutputHB_TD
// Purpose       : Output the current data for an HB run
// Special Notes :  We will require the user to have specified separate HB_TD
//                 and HB_FD output interfaces, and this function will output
//                 only the timedomain data appropriate for the HB_TD interface.
//
//                 If the user specifies "HB" as the output type for their
//                 interface object, we'll never be called. Only HB_FD and
//                 HB_TD will be accepted.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 27 Mar 2018
//-----------------------------------------------------------------------------
void OutputterExternal::doOutputHB_TD(
   Parallel::Machine             comm,
   const std::vector<double> &   timePoints,
   const Linear::BlockVector &     timeDomainSolutionVec,
   const Linear::BlockVector &     timeDomainLeadCurrentVec,
   const Linear::BlockVector &     timeDomainJunctionVoltageVec)
{
  if (Parallel::rank(comm) == 0 && !initialized_)
  {
    initialized_ = true;
    theOutputWrapper_->outputFieldNames(fieldNames_);
  }

  if (theOutputWrapper_->getOutputType() == OutputType::HB_TD)
  {
    int blockCount = timeDomainSolutionVec.blockCount();
    // Loop over the time points of the Linear::BlockVecor:
    for (int iblock = 0; iblock < blockCount; ++iblock)
    {
      outputManager_.setCircuitTime(timePoints[iblock]);
      Linear::Vector * time_solution_vector =
        &timeDomainSolutionVec.block(iblock);
      Linear::Vector * time_lead_current_vector =
        &timeDomainLeadCurrentVec.block(iblock);
      Linear::Vector * time_junction_voltage_vector =
        &timeDomainJunctionVoltageVec.block(iblock);

      {
        // periodic time-domain steady-state output
        std::vector<complex> complexResultList;
        getValues(comm, opList_,
                  Util::Op::OpData(index_, time_solution_vector,
                                   0, 0, 0, 0,
                                   time_lead_current_vector, 0,
                                   time_junction_voltage_vector),
                  complexResultList);
        std::vector<double> doubleResultList(complexResultList.size());
        for (int i = 0; i < complexResultList.size(); ++i)
          doubleResultList[i]=complexResultList[i].real();

        if (initialized_)
          theOutputWrapper_->outputReal(doubleResultList);
      }

      ++index_;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doFinishOutput
// Purpose       : Signal to external code that output is complete
// Special Notes : Corresponds to printing a footer in other outputters
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
void OutputterExternal::doFinishOutput()
{
  theOutputWrapper_->finishOutput();
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doStartStep
// Purpose       : Signal to external code that a new step of a .STEP loop is
//                 starting.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
void OutputterExternal::doStartStep(int current_step, int number_of_steps)
{
  currentStep_   =current_step;
  numberOfSteps_ =number_of_steps;
  theOutputWrapper_->newStepOutput(currentStep_,numberOfSteps_);
}

//-----------------------------------------------------------------------------
// Function      : OutputterExternal::doSteppingComplete
// Purpose       : Signal to external code that all stepping is complete
// Special Notes : Corresponds to printing a step footer in other outputters
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
void OutputterExternal::doSteppingComplete()
{
  theOutputWrapper_->finishedStepping();
}

}
}
}
