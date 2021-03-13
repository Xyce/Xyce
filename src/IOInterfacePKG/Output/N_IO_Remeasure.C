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
// Purpose       : Functions for remeasure classes.
// Special Notes :
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 2/15/2018
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_Remeasure.h>

namespace Xyce {
namespace IO {
namespace Measure{

//-----------------------------------------------------------------------------
// Function      : RemeasureBase::RemeasureBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureBase::RemeasureBase(
  Parallel::Communicator &pds_comm, 
  Manager &measure_manager, 
  OutputMgr &output_manager,
  Analysis::AnalysisManager &analysis_manager,
  Analysis::AnalysisCreatorRegistry &analysis_registry) 
  : pds_comm(pds_comm),
    measure_manager(measure_manager),
    output_manager(output_manager),
    analysis_manager(analysis_manager),
    analysis_registry(analysis_registry),
    analysis_mode(Xyce::Analysis::ANP_MODE_INVALID),
    index(-1) 
{
}

//-----------------------------------------------------------------------------
// Function      : RemeasureBase::~RemeasureBase
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureBase::~RemeasureBase()
{
}

//-----------------------------------------------------------------------------
// Function      : RemeasureAC::RemeasureAC
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureAC::RemeasureAC(
  Parallel::Communicator &pds_comm, 
  Manager &measure_manager, 
  OutputMgr &output_manager,
  Analysis::AnalysisManager &analysis_manager,
  Analysis::AnalysisCreatorRegistry &analysis_registry) 
  : RemeasureBase(pds_comm, measure_manager, output_manager, analysis_manager, analysis_registry),
    frequency(0)  
{
  // this enum is used to set the file suffix (.ma) and also during checking of the
  // measure mode vs. the analysis mode.
  setAnalysisMode(Xyce::Analysis::ANP_MODE_AC);

  // Note: allocateAnalysisObject(*analysisRegistry_) will likely fail for .AC.  The AC object
  // constructor segfaults on the initialization of the bVecRealPtr and bVecImgPtr member
  // variables.
  //analysis_manager.allocateAnalysisObject(analysis_registry);
}


//-----------------------------------------------------------------------------
// Function      : RemeasureAC::RemeasureAC
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureAC::~RemeasureAC()
{
}


//-----------------------------------------------------------------------------
// Function      : RemeasureAC::setIndepVarCol
// Purpose       : Set the index variable in the base object.  This is the
//               : column index of the FREQ column in the remeasured data file.
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
void RemeasureAC::setIndepVarCol(int rank, int i, std::string colName)
{
  if ( (i<3) && (colName=="FREQ") )
  {
    (rank == 0) ? index = i : index = -1;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : RemeasureAC::checkIndepVarCol
// Purpose       : Checks that a FREQ column was found in the remeasured data file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/28/19
//-----------------------------------------------------------------------------
void RemeasureAC::checkIndepVarCol(int rank, int index)
{
  if ( (rank == 0) && index < 0 )
  {
    Report::UserFatal() << "FREQ column not found in remeasured output file for AC-mode remeasure";
  }
}

//-----------------------------------------------------------------------------
// Function      : RemeasureDC::RemeasureDC
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureDC::RemeasureDC(
  Parallel::Communicator &pds_comm, 
  Manager &measure_manager, 
  OutputMgr &output_manager,
  Analysis::AnalysisManager &analysis_manager,
  Analysis::AnalysisCreatorRegistry &analysis_registry) 
  : RemeasureBase(pds_comm, measure_manager, output_manager, analysis_manager, analysis_registry),
    SweepVecPtr(0),
    dcParamsVec(0),
    DCIndex(0),
    dcStepCount(0)  
{
  // this enum is used to set the file suffix (.ms) and also during checking of the
  // measure mode vs. the analysis mode.
  setAnalysisMode(Xyce::Analysis::ANP_MODE_DC_SWEEP);

  // If we are doing remeasure for .DC, then we need to initialize the primary analysis
  // object, so that the DC Sweep Vector can be made.  The first element in the
  // DC Sweep Vector is used by the FROM-TO qualifiers for DC measures.
  analysis_manager.allocateAnalysisObject(analysis_registry);
  Analysis::DCSweep *SweepVecPtr = dynamic_cast<Analysis::DCSweep*>(&(analysis_manager.getAnalysisObject()));
  if (SweepVecPtr == 0)
  {
    Report::DevelFatal0() << "Unknown Error making DC Sweep Vector for Remeasure" << std::endl;
  }
  else
  {
    // check if the "DATA" specification was used.  If so, then a new vector of
    // SweepParams, in the "TABLE" style.  Otherwise, error out.
    if (analysis_manager.getAnalysisObject().getDataSpecification() && !SweepVecPtr->convertDataToSweepParams())
    {
      Report::DevelFatal0() << "Error making DC Sweep Vector for Remeasure" << std::endl;
    }

    // populate the values in the DC Sweep Vector
    analysis_manager.getOutputManagerAdapter().setDCSweepVector(SweepVecPtr->getDCSweepVec());
  }

  // Make a version of the DC Params Vector that is sufficient for doing remeasure
  // for .DC.  It may not have all of the fields filled out properly though.  Also, we
  // only update the values for the first element in that vector, which is the first variable
  // on the .DC line.  That first variable is used by the FROM-TO qualifiers for DC measures.
  dcParamsVec = analysis_manager.getOutputManagerAdapter().getDCSweepVector();
  if (dcParamsVec.size() == 0)
  {
    // DC sweep vector should have at least one element in it
    Report::DevelFatal0() << "Error making DC Sweep Vector for Remeasure";
  }
  setSweepLoopVals(dcParamsVec.begin(), dcParamsVec.end());
  dcParamsVec[0].updateCurrentVal(dcStepCount);
}


//-----------------------------------------------------------------------------
// Function      : RemeasureDC::~RemeasureDC
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureDC::~RemeasureDC()
{
  //clean-up pointers
  delete SweepVecPtr;
  SweepVecPtr = 0;
}


//-----------------------------------------------------------------------------
// Function      : RemeasureDC::setIndepVarCol
// Purpose       : set the index variable in the base object.  This is the
//               : column index of the INDEX column in the remeasured data file.
// Special Notes : It is a fatal error if the remeasured data file does not
//               : have an INDEX column when re-measuring DC data.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
void RemeasureDC::setIndepVarCol(int rank, int i, std::string colName)
{
  if ( (i< 2) && (colName=="Index") )
  {
    (rank == 0) ? index = i : index = -1;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : RemeasureDC::checkIndepVarCol
// Purpose       : Checks that an INDEX column was found in the remeasured data file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/28/19
//-----------------------------------------------------------------------------
void RemeasureDC::checkIndepVarCol(int rank, int index)
{
  if ( (rank == 0) && index < 0 )
  {
    // For DC mode, the comparison file must have an Index.  It will be used
    // later to "sense" when a step, caused by a .STEP line, has occurred
    // in the data.
    Report::UserFatal() << "Index column not found in remeasured output file for DC-mode remeasure";
  }
}


//-----------------------------------------------------------------------------
// Function      : RemeasureDC::resetSweepVars
// Purpose       : reset sweep variables specified on the .DC line
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
void RemeasureDC::resetSweepVars()
{
  dcStepCount=0;
  dcParamsVec[0].updateCurrentVal(dcStepCount);
  dcParamsVec[0].count=0; // for compatibility of .STEP with .DC data=table sweeps

  return;
}


//-----------------------------------------------------------------------------
// Function      : RemeasureDC::updateSweepVars
// Purpose       : update sweep variables specified on the .DC line
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
void RemeasureDC::updateSweepVars()
{
  // used to sense a "step" caused by multiple sweep variables on the .DC line
  ++dcStepCount;

  // calculate the value of first sweep variable (on the .DC line) that will
  // be used for the next iteration of this loop.  That variable is used as 
  // by the FROM-TO qualifiers for .DC measures.
  if (dcStepCount >= dcParamsVec[0].maxStep)
  {
    // Step (in .DC line) will occur on next line in the output file
    // for the 1st variable.  So, reset its value to its "start value" 
    dcStepCount=0;
    dcParamsVec[0].updateCurrentVal(dcStepCount);
  }
  else
  {
    // update value of 1st sweep variable (on .DC line) for use in 
    // next loop iteration
    dcParamsVec[0].updateCurrentVal(dcStepCount);
  }
}


//-----------------------------------------------------------------------------
// Function      : RemeasureTRAN::RemeasureTRAN
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureTRAN::RemeasureTRAN(
  Parallel::Communicator &pds_comm, 
  Manager &measure_manager, 
  OutputMgr &output_manager,
  Analysis::AnalysisManager &analysis_manager,
  Analysis::AnalysisCreatorRegistry &analysis_registry)
  : RemeasureBase(pds_comm, measure_manager, output_manager, analysis_manager, analysis_registry),
    endSimTime(0.0)
{
  // this enum is used to set the file suffix (.mt) and also during checking of the
  // measure mode vs. the analysis mode.
  setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

  analysis_manager.allocateAnalysisObject(analysis_registry);
  endSimTime = analysis_manager.getFinalTimeForRemeasure();
}


//-----------------------------------------------------------------------------
// Function      : RemeasureTRAN::~RemeasureTRAN
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
RemeasureTRAN::~RemeasureTRAN()
{
  
}


//-----------------------------------------------------------------------------
// Function      : RemeasureTRAN::setIndepVarCol
// Purpose       : set the index variable in the base object.  This is the
//               : column index of the TIME column in the remeasured data file.
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-----------------------------------------------------------------------------
void RemeasureTRAN::setIndepVarCol(int rank, int i, std::string colName)
{
  if ( (i<3) && (colName=="TIME") )
  {
    (rank == 0) ? index = i : index = -1;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : RemeasureTRAN::checkIndepVarCol
// Purpose       : Checks that a TIME column was found in the remeasured data file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/28/19
//-----------------------------------------------------------------------------
void RemeasureTRAN::checkIndepVarCol(int rank, int index)
{
  if ( (rank == 0) && index < 0 )
  {
    Report::UserFatal() << "TIME column not found in remeasured output file for TRAN-mode remeasure";
  }
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
