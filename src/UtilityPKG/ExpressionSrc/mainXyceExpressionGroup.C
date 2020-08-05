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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <unordered_map>
#include <string>
#include <random>

#include <mainXyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_UQSupport.h>
#include <N_ANP_SweepParam.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Op.h>
#include <N_DEV_Const.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_DeviceNameConverters.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::mainXyceExpressionGroup 
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
mainXyceExpressionGroup::mainXyceExpressionGroup ( 
 N_PDS_Comm & comm, Topo::Topology & top,
 Analysis::AnalysisManager &analysis_manager,
 Device::DeviceMgr & device_manager,
 IO::OutputMgr &output_manager
 ) :
 comm_(comm),
 top_(top),
 analysisManager_(analysis_manager),
 deviceManager_(device_manager),
 outputManager_(output_manager),
 time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), gmin_(0.0), dt_(0.0), alpha_(0.0),
  randomSeed_(0), 
  randomSetup_(false)
{
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::mainXyceExpressionGroup 
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
mainXyceExpressionGroup::~mainXyceExpressionGroup ()
{
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionGID_
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
int mainXyceExpressionGroup::getSolutionGID_(const std::string & nodeName)
{
  int tmpGID=-1;
  std::vector<int> svGIDList1, dummyList;
  char type1;

  std::string nodeNameUpper = nodeName;
  Xyce::Util::toUpper(nodeNameUpper);

  bool foundLocal = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_VNODE), svGIDList1, dummyList, type1);
  bool found = static_cast<int>(foundLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found, 1);

  // if looking for this as a voltage node failed, try a "device" (i.e. current) node.  I(Vsrc)
  bool foundLocal2 = false;
  if (!found)
  {
    foundLocal2 = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_DNODE), svGIDList1, dummyList, type1);
  }
  bool found2 = static_cast<int>(foundLocal2);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found2, 1);

  // Check if this is a subcircuit interface node name, which would be found in the aliasNodeMap.
  // If it is then get the GID for the corresponding "real node name". See SRN Bug 1962 for 
  // more details but, as an example, these netlist lines:
  //
  //   X1 1 2 MYSUB
  //   .SUBCKT MYSUB a c 
  //   R1   a b 0.5
  //   R2   b c 0.5
  //  .ENDS
  //
  // would produce key-value pairs of <"X1:C","2"> and <"X1:A","1"> in the aliasNodeMap 
  bool foundAliasNodeLocal = false;
  if (!found && !found2)
  {
    IO::AliasNodeMap::const_iterator alias_it = aliasNodeMap_.find( nodeNameUpper );
    if (alias_it != aliasNodeMap_.end())
    {      
      foundAliasNodeLocal = top_.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
    }
  }
  bool foundAliasNode = static_cast<int>(foundAliasNodeLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &foundAliasNode, 1);

  if(svGIDList1.size()==1)
  {
    tmpGID = svGIDList1.front();
  }

  return tmpGID;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  retval = 0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const Linear::Vector * nextSolVector = deviceManager_.getExternData().nextSolVectorPtr;
    if (nextSolVector)
    {
      retval = nextSolVector->getElementByGlobalIndex(tmpGID, 0);
    }

  }
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &retval, 1);
  return (tmpGID>=0);
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : std::complex<double> precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getSolutionVal(const std::string & nodeName, std::complex<double> & retval)
{
  double real_val=0.0;
  double imag_val=0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const Linear::Vector * nextSolVector = deviceManager_.getExternData().nextSolVectorPtr;
    if (nextSolVector)
    { 
      real_val = nextSolVector->getElementByGlobalIndex(tmpGID, 0);
    }
  }

  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &real_val, 1);

  // ERK.  To Do:
  // need to add logic to see if this is frequency domain situation or not.
  // If not, don't bother getting the imaginary part.
  //
  // Xyce currently uses a real equivalent form for everything, rather than std::complex.
  //
  // For now, however, just set imag to zero.

  retval = std::complex<double>(real_val,imag_val);
  return (tmpGID>=0);
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : double precision version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getGlobalParameterVal(const std::string &paramName, double & retval)
{
  bool success=true;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, retval);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : std::complex<double> version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval)
{
  bool success=true;
  double tmpval;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, tmpval);
  retval = std::complex<double>(tmpval,0.0);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTimeStep ()
{
  dt_ = deviceManager_.getSolverState().currTimeStep_;
  return dt_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTime
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTime() 
{ 
  // I would have preferred to use this but as the code is currently written it
  // is not safe.  The earliest call I would need to make to getTime happens before 
  // the stepErrorControl class has been allocated.  Unfortunately, the analysis
  // manager accessor returns an invalid reference in that case, which I can't 
  // really test for.
  //
  //const TimeIntg::StepErrorControl & secControl_ = (analysisManager_.getStepErrorControl());
  //time_ = secControl_.nextTime;
  
  time_ = deviceManager_.getSolverState().currTime_;
  return time_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getTemp
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getTemp() 
{ 
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>() - CONSTCtoK;
  return temp_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getVT  () 
{ 
  VT_ = (deviceManager_.getDeviceOptions().temp.getImmutableValue<double>())*CONSTKoverQ;
  return VT_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getFreq() 
{ 
  freq_ = deviceManager_.getSolverState().currFreq_;
  return freq_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getGmin
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getGmin() 
{ 
  gmin_ = deviceManager_.getGmin();
  return gmin_;
} 

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getBpTol()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getBpTol()
{
  return deviceManager_.getSolverState().bpTol_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getStartingTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getStartingTimeStep()
{
  return deviceManager_.getSolverState().startingTimeStep_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getFinalTime()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double mainXyceExpressionGroup::getFinalTime()
{
  return deviceManager_.getSolverState().finalTime_;
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getStepNumber()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/9/2020 
//-------------------------------------------------------------------------------
unsigned int mainXyceExpressionGroup::getStepNumber()
{
  //return deviceManager_.getSolverState().timeStepNumber_; // either of these should work
  return analysisManager_.getStepNumber();
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getPhaseOutputUsesRadians
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/13/2020 
//-------------------------------------------------------------------------------
bool mainXyceExpressionGroup::getPhaseOutputUsesRadians()
{
  return outputManager_.getPhaseOutputUsesRadians();
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::getRandomOpValue
//
// Purpose       : provide a value for a single random operator
//
// Special Notes : This does NOT connect to the sampling stuff (yet?)  
//                 This is written to generate the random numbers directly in
//                 the group.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/3/2020 
//-------------------------------------------------------------------------------
void mainXyceExpressionGroup::getRandomOpValue (
    Util::astRandTypes type, 
    std::vector<double> args, 
    double & value) 
{
  if (!randomSetup_)
  {
    setupRandom_();
    randomSetup_  = true;
  }

  if (type==Util::AST_AGAUSS)
  {
    double mean, stddev, n;
    int argSize = args.size();
    if (argSize < 2) { Xyce::Report::DevelFatal() << "Error.  mainXyceExpressionGroup::getRandomOpValue" <<std::endl; }
    else { mean=args[0]; stddev=args[1]; }

    if(argSize ==3) { n=args[2];  stddev /= n; }
    if (Xyce::Util::enableRandomExpression)
    {
      value = theRandomSamplesGenerator->generateNormalSample(mean, stddev);
    }
    else
    {
      value = mean;
    }
  }
  else if (type==Util::AST_GAUSS)
  {
    double mean, stddev, n;
    int argSize = args.size();
    if (argSize < 2) { Xyce::Report::DevelFatal() << "Error.  mainXyceExpressionGroup::getRandomOpValue" <<std::endl; }
    else { mean=args[0]; stddev=args[1]*mean; }

    if(argSize ==3) { n=args[2];  stddev /= n; }
    if (Xyce::Util::enableRandomExpression)
    {
      value = theRandomSamplesGenerator->generateNormalSample(mean, stddev);
    }
    else
    {
      value = mean;
    }
  }
  else if (type==Util::AST_AUNIF)
  {
    double mean, variation, n;
    int argSize = args.size();
    if (argSize < 2) { Xyce::Report::DevelFatal() << "Error.  mainXyceExpressionGroup::getRandomOpValue" <<std::endl; }
    else { mean=args[0]; variation=std::abs(std::real(args[1])); }

    if(argSize ==3) { n=args[2];  variation /= n; }

    double min   = std::real(mean) - std::real(variation);
    double max   = std::real(mean) + std::real(variation);
    if (Xyce::Util::enableRandomExpression)
    {
      value = theRandomSamplesGenerator->generateUniformSample(min, max);
    }
    else
    {
      value = mean;
    }
  }
  else if (type==Util::AST_UNIF)
  {
    double mean, variation, n;
    int argSize = args.size();
    if (argSize < 2) { Xyce::Report::DevelFatal() << "Error.  mainXyceExpressionGroup::getRandomOpValue" <<std::endl; }
    else { mean=args[0]; variation=std::abs(std::real(args[1]*mean)); }

    if(argSize ==3) { n=args[2];  variation /= n; }

    double min   = std::real(mean) - std::real(variation);
    double max   = std::real(mean) + std::real(variation);
    if (Xyce::Util::enableRandomExpression)
    {
      value = theRandomSamplesGenerator->generateUniformSample(min, max);
    }
    else
    {
      value = mean;
    }
  }
  else if (type==Util::AST_RAND)
  {
    if (Xyce::Util::enableRandomExpression)
    {
      value = theRandomSamplesGenerator->generateUniformSample(0.0,1.0);
    }
    else
    {
      value = 0.5;
    }
  }
  else if (type==Util::AST_LIMIT)
  {
  }
  else
  {
    Xyce::Report::DevelFatal() << "Error.  mainXyceExpressionGroup::getRandomOpValue" <<std::endl;
  }

  return; 
}

//-------------------------------------------------------------------------------
// Function      : mainXyceExpressionGroup::setupRandom_
// Purpose       : provide a value for a single random operator
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/3/2020 
//-------------------------------------------------------------------------------
void mainXyceExpressionGroup::setupRandom_ ()
{
  if (theRandomSamplesGenerator == 0) 
  { 
    bool userSeedGiven=false;
    long userSeed=0;
    randomSeed_ = Analysis::UQ::getTheSeed( analysisManager_.getComm(), analysisManager_.getCommandLine(), userSeed, userSeedGiven);

    theRandomSamplesGenerator = new randomSamplesGenerator(randomSeed_); 
  }
  else
  {
    randomSeed_ = theRandomSamplesGenerator->randomSeed;
  }
}

} // Util
} // Xyce
