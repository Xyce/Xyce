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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <N_DEV_DeviceMgr.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_ExtOutInterface.h>
#include <N_IO_ExtOutWrapper.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_Op.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputFileBase.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterHomotopy.h>
#include <N_IO_OutputterDC.h>
#include <N_IO_OutputterAC.h>
#include <N_IO_OutputterSParam.h>
#include <N_IO_OutputterNoise.h>
#include <N_IO_OutputterEmbeddedSampling.h>
#include <N_IO_OutputterPCE.h>
#include <N_IO_OutputterHB.h>
#include <N_IO_OutputterMPDE.h>
#include <N_IO_OutputterTransient.h>
#include <N_IO_OutputterSensitivity.h>
#include <N_IO_OutputterRawOverride.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_mmio.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>

#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_DeviceNameConverters.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_Math.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {

namespace SensitivityOptions {

enum {
  DIRECT   = 0x01,
  ADJOINT  = 0x02,
  SCALED   = 0x04,
  UNSCALED = 0x08
};

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::OutputMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::OutputMgr(
  const CmdParse &              command_line,
  Util::Op::BuilderManager &    op_builder_manager,
  const Topo::Topology &        topology)
  : title_(command_line.getArgumentValue("netlist")),
    netlistFilename_(command_line.getArgumentValue("netlist")),
    filenameSuffix_(),
    opBuilderManager_(op_builder_manager),
    topology_(topology),
    dotOpSpecified_(false),
    enableEmbeddedSamplingFlag_(false),
    enablePCEFlag_(false),
    enableHomotopyFlag_(false),
    enableSparCalcFlag_(false),
    enableSensitivityFlag_(false),
    adjointSensitivityFlag_(false),
    sensitivityOptions_(0),
    pts_per_summary_(-1),
    pts_per_summary_Given(false),
    PRINTdcstart_(0.0),
    PRINTdcstop_(0.0),
    PRINTdcvalue_(0.0),
    PRINTdcname_(""),
    initialOutputInterval_(0.0),
    printHeader_(true),
    printFooter_(true),
    printStepNumCol_(false),
    outputVersionInRawFile_(false),
    phaseOutputUsesRadians_(true),
    outputCalledBefore_(false),
    dcLoopNumber_(0),
    maxDCSteps_(0),
    hdf5FileNameGiven_(false),
    hdf5HeaderWritten_(false),
    hdf5IndexValue_(0)
{
  if (command_line.getArgumentValue("-delim") == "TAB")
  {
    defaultPrintParameters_.delimiter_ = "\t";
  }
  else if (command_line.getArgumentValue("-delim") == "COMMA")
  {
    defaultPrintParameters_.delimiter_ = ",";
  }
  else
  {
    defaultPrintParameters_.delimiter_ = command_line.getArgumentValue("-delim");
  }

  if (command_line.argExists("-a"))
  {
    defaultPrintParameters_.asciiRaw_ = true;
  }

  // -r takes precedence over -o if both are specified on the command line.
  // For more info on the behavior of -o, and its precedence over FILE= on 
  // a .PRINT line, see SON Bug 911.
  if (command_line.argExists("-r")) 
  {
    defaultPrintParameters_.overrideRaw_ = true;
    defaultPrintParameters_.overrideRawFilename_ = command_line.getArgumentValue("-r");
    defaultPrintParameters_.format_ = 
      defaultPrintParameters_.asciiRaw_ ? Format::RAW_ASCII : Format::RAW;
  }
  else if (command_line.argExists("-o"))
  {
    defaultPrintParameters_.dashoRequested_ = true;
    defaultPrintParameters_.dashoFilename_ = command_line.getArgumentValue("-o");
  }

  sensitivityOptions_ |= SensitivityOptions::ADJOINT;
  sensitivityOptions_ |= SensitivityOptions::UNSCALED;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::~OutputMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::~OutputMgr()
{
  OutputterMap::iterator it = outputterMap_.begin();
  for ( ; it != outputterMap_.end(); ++it)
  {
    OutputterMap::mapped_type::iterator it2 = (*it).second.begin();
    for ( ; it2 != (*it).second.end(); ++it2)
    {
      delete (*it2);
    }
  }

  OpenPathStreamMap::iterator it1 = openPathStreamMap_.begin();
  for ( ; it1 != openPathStreamMap_.end(); ++it1)
  {
    delete (*it1).second.second;
  }

  for (ExternalOutputWrapperMap::iterator itW = externalOutputWrapperMap_.begin();
         itW != externalOutputWrapperMap_.end(); ++itW)
  {
    for (std::vector<ExternalOutputWrapper *>::iterator itW2=(*itW).second.begin();
         itW2 != (*itW).second.end(); ++itW2)
    {
      ExternalOutputWrapper * theExtOutWrapper=(*itW2);
      delete theExtOutWrapper;
    }
    (*itW).second.clear();
  }
  
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getVarTypes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const std::vector<char> & OutputMgr::getVarTypes() const
{
  return topology_.getVarTypes();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getSolutionNodeMap
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getSolutionNodeMap() const
{
  return topology_.getSolutionNodeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getStateNodeMap
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getStateNodeMap() const
{
  return topology_.getStateNodeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getStoreNodeMap
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getStoreNodeMap() const
{
  return topology_.getStoreNodeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getExternalNodeMap
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getExternalNodeMap() const
{
  return topology_.getExternalNodeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getBranchVarsNodeMap
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getBranchVarsNodeMap() const
{
  return topology_.getBranchVarsNodeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getNoiseDeviceNameMap
// Purpose       : This gets the indexes into the noiseDataVec_ for the 
//                 deviceName[i]->deviceName entries that were put into 
//                 the symbol table during the construction of the NOISE
//                 analysis object.  These are the names of devices that 
//                 have noise sources. 
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 11/22/2017
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getNoiseDeviceNameMap() const
{
  return topology_.getNoiseDeviceNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getNoiseTypeNameMap
// Purpose       : This gets the indexes into the noiseDataVec_ for the 
//                 noiseDataVec_[i]->noiseNames[j] entries that were put into 
//                 the symbol table during the construction of the NOISE
//                 analysis object.  These are the noise types names, where
//                 where each device may have multiple noise types. 
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 6/01/2018
//-----------------------------------------------------------------------------
const NodeNameMap & OutputMgr::getNoiseTypeNameMap() const
{
  return topology_.getNoiseTypeNameMap();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::notify
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgr::notify(
  const Analysis::StepEvent &   step_event)
{
  switch (step_event.state_) 
  {
    case Analysis::StepEvent::INITIALIZE:
      outputState_.stepMaxCount_ = step_event.count_;
      break;

    case Analysis::StepEvent::STEP_STARTED:
      outputState_.stepLoopNumber_ = step_event.count_;
      startStep(outputState_.stepLoopNumber_, outputState_.stepMaxCount_);
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      break;

    case Analysis::StepEvent::FINISH:
      steppingComplete();
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file in given mode, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream * OutputMgr::openFile(
  const std::string &           path,
  std::ios_base::openmode       mode)
{
  OpenPathStreamMap::iterator it= openPathStreamMap_.find(path);

  if (path == "CONSOLE")
  {
    return &Xyce::dout();
  }
  else if (it != openPathStreamMap_.end()) 
  {
    ++(*it).second.first;
    return (*it).second.second;
  }
  else 
  {
    std::ostream *os = new std::ofstream(path.c_str(), mode);
    openPathStreamMap_[path] = std::pair<int, std::ostream *>(1, os);

    if (!os->good())
    {
      Report::UserFatal() << "Failure opening " << path;
    }

    return os;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file for output only, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream * OutputMgr::openFile(const std::string & path)
{
  return openFile(path, std::ios_base::out);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openBinaryFile
// Purpose       : open named file in binary mode for output only, create
//                 stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream * OutputMgr::openBinaryFile(const std::string & path)
{
  return openFile(path, std::ios_base::out | std::ios_base::binary);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::closeFile
// Purpose       : Close given stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
int OutputMgr::closeFile(std::ostream * os)
{
  if (os == &Xyce::dout())
  {
    return 1;
  }

  int open_count= 0;

  for (OpenPathStreamMap::iterator it= openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it) 
  {
    if ((*it).second.second == os) 
    {
      open_count= --(*it).second.first;
      if (open_count == 0) 
      {
        delete os;
        openPathStreamMap_.erase(it);
        break;
      }
    }
  }

  return open_count;
}


namespace {

//-----------------------------------------------------------------------------
// Function      : testAndSet
// Purpose       : test if set contains element and then add the element
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jun 24 09:38:58 2014
//-----------------------------------------------------------------------------
///
/// Returns true of set contained element on entry, adds element if it was not found
///
/// @invariant set contains new element on completion
///
/// @param s Set to find element
/// @param t Element to find
///
/// @return true if set contained element before insertion
///
template<class S, class T>
inline bool testAndSet(S &s, const T &t)
{
  bool found = s.find(t) != s.end();
  if (!found)
    s.insert(t);
  return found;
}

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::prepareOutput
// Purpose       : Select a set of outputtters appropriate for analysis type
// Special Notes : 
// Scope         : public
// Creator       : David Baur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
///
/// Select appropriate set of outputters for current analysis mode
///
/// @param[in] comm   Communicator for this run
/// @param analysis_mode analysis mode for which to activate output
///
/// The outputterMap object contains all of the allocated outputters
/// that have been created for .PRINT lines.  Some of these may not be
/// needed for the current analysis.  This function is called
/// automatically when an analysis creates an "ActiveOutput" sentry
/// object.  It selects all outputters that apply to the current
/// analysis mode and makes them active.
///
/// @note The "enable*" functions that are called here are all free functions
/// of the Xyce::IO::Outputter namespace, and are defined in the various
/// "N_IO_Outputter<analysistype>.C" files.  These functions in fact
/// allocate new outputter objects for each set of print parameters
/// of the appropriate type, and add them to the output manager's outputter
/// map.
/// @author David Baur
/// @date 06/28/2013
void OutputMgr::prepareOutput(
  Parallel::Machine             comm,
  Analysis::Mode                analysis_mode)
{
  // Setup rawfile if requested
  if (defaultPrintParameters_.overrideRaw_)
  {
    // This is a temporary fix until the -r and FORMAT=RAW outputs
    // are defined for LIN and HB analyses.  Note that the enableSparCalcFlag_
    // is set via the AC object, during that object's parsing of the ACLIN
    // option block.
    if (analysis_mode == Analysis::ANP_MODE_HB || enableSparCalcFlag_ ||
        enableEmbeddedSamplingFlag_ || enablePCEFlag_)
    {
      Report::UserFatal0() << "-r and -a outputs are not supported for Embedded Sampling, .HB, .LIN or .PCE analyses";
    }

    Outputter::enableRawOverrideOutput(comm, *this, analysis_mode);

    if (activeOutputterStack_.empty())
      pushActiveOutputters();

    addActiveOutputter(PrintType::RAW_OVERRIDE, analysis_mode);

    // also make SENS and Homotopy output files, since they don't support
    // RAW output.
    if (enableHomotopyFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 100)))
        Outputter::enableHomotopyOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::HOMOTOPY, analysis_mode);
    }

    if (enableSensitivityFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 101)))
      {
        if (analysis_mode == Analysis::ANP_MODE_AC)
          Outputter::enableSensitivityACOutput(comm, *this, analysis_mode);
        else
          Outputter::enableSensitivityOutput(comm, *this, analysis_mode);
      }
      addActiveOutputter(PrintType::SENS, analysis_mode);

      // allocate an adjoint output object, but do not add to the active list.
      // Need to re-think the INVALID+102 thing.
      if ( adjointSensitivityFlag_ )
      {
        if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 102)))
          Outputter::enableAdjointSensitivityOutput(comm, *this, analysis_mode);

        // ERK.  NOT putting the transient adjoint outputters on the stack.  
        // Accessing them (for now) directly from the outputterMap.  I am punting 
        // on making the active outputter structure work for transient adjoints.
      }
    }
  }
  else if (enableEmbeddedSamplingFlag_)
  {
    if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 103)))
      Outputter::enableEmbeddedSamplingOutput(comm, *this, analysis_mode);
    addActiveOutputter(PrintType::ES, analysis_mode);

    if (!defaultPrintParameters_.dashoRequested_)
    {
      // Only .PRINT ES output is made for the -o case.  .PRINT TRAN or .PRINT DC
      // output is not made.
      if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_TRANSIENT))
        Outputter::enableTransientOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::TRAN, analysis_mode);

      if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_DC_SWEEP))
        Outputter::enableDCOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::TRAN, analysis_mode);
    }
  }
  else if (enablePCEFlag_)
  {
    if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 104)))
      Outputter::enablePCEOutput(comm, *this, analysis_mode);
    addActiveOutputter(PrintType::PCE, analysis_mode);

    if (!defaultPrintParameters_.dashoRequested_)
    {
      // Only .PRINT PCE output is made for the -o case.  .PRINT TRAN or .PRINT DC
      // output is not made.
      if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_TRANSIENT))
        Outputter::enableTransientOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::TRAN, analysis_mode);

      if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_DC_SWEEP))
        Outputter::enableDCOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::TRAN, analysis_mode);
    }
  }
  else
  {
    switch (analysis_mode)
    {
      case Analysis::ANP_MODE_INVALID:
      case Analysis::ANP_MODE_DC_OP:
      case Analysis::ANP_MODE_MOR:
        break;

      case Analysis::ANP_MODE_DC_SWEEP:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_DC_SWEEP))
          Outputter::enableDCOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::TRAN, analysis_mode);
        break;

      case Analysis::ANP_MODE_TRANSIENT:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_TRANSIENT))
          Outputter::enableTransientOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::TRAN, analysis_mode);
        break;

      case Analysis::ANP_MODE_MPDE:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_MPDE))
          Outputter::enableMPDEOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::MPDE, analysis_mode);
        addActiveOutputter(PrintType::MPDE_IC, analysis_mode);
        addActiveOutputter(PrintType::MPDE_STARTUP, analysis_mode);
        break;

      case Analysis::ANP_MODE_HB:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_HB))
          Outputter::enableHBOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::HB_FD, analysis_mode);
        addActiveOutputter(PrintType::HB_TD, analysis_mode);
        addActiveOutputter(PrintType::HB_IC, analysis_mode);
        addActiveOutputter(PrintType::HB_STARTUP, analysis_mode);
        break;

      case Analysis::ANP_MODE_AC:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_AC))
	{
          // if -o was requested then either SPARAM or AC output is made
          // depending on whether a .LIN analysis is being done, or not
          if (!(defaultPrintParameters_.dashoRequested_ && enableSparCalcFlag_))
	  {
            Outputter::enableACOutput(comm, *this, analysis_mode);
          }
          Outputter::enableSParamOutput(comm, *this, analysis_mode);
        }
        if (!(defaultPrintParameters_.dashoRequested_ && enableSparCalcFlag_))
	{
          addActiveOutputter(PrintType::AC, analysis_mode);
          addActiveOutputter(PrintType::AC_IC, analysis_mode);
        }
        addActiveOutputter(PrintType::SPARAM, analysis_mode);
        break;

      case Analysis::ANP_MODE_NOISE:
        if (!testAndSet(enabledAnalysisSet_, Analysis::ANP_MODE_NOISE))
          Outputter::enableNoiseOutput(comm, *this, analysis_mode);
        addActiveOutputter(PrintType::NOISE, analysis_mode);
        break;
    }

    if (enableHomotopyFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 100)))
        Outputter::enableHomotopyOutput(comm, *this, analysis_mode);
      addActiveOutputter(PrintType::HOMOTOPY, analysis_mode);
    }

    if (enableSensitivityFlag_)
    {
      if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 101)))
      {
        if (analysis_mode == Analysis::ANP_MODE_AC)
          Outputter::enableSensitivityACOutput(comm, *this, analysis_mode);
        else
          Outputter::enableSensitivityOutput(comm, *this, analysis_mode);
      }
      addActiveOutputter(PrintType::SENS, analysis_mode);

      // allocate an adjoint output object, but do not add to the active list.
      // Need to re-think the INVALID+102 thing.
      if ( adjointSensitivityFlag_ )
      {
        if (!testAndSet(enabledAnalysisSet_, (Analysis::Mode) (Analysis::ANP_MODE_INVALID + 102)))
          Outputter::enableAdjointSensitivityOutput(comm, *this, analysis_mode);

        // ERK.  NOT putting the transient adjoint outputters on the stack.  
        // Accessing them (for now) directly from the outputterMap.  I am punting 
        // on making the active outputter structure work for transient adjoints.
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setOPAnalysisParams
// Purpose       : set the dotOpSpecified_ flag if the netlist has 
//                 a .OP statement
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 8/20/18
//-----------------------------------------------------------------------------
bool OutputMgr::setOPAnalysisParams(
  const Util::OptionBlock &     paramsBlock)
{
  dotOpSpecified_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerNonlinearOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerNonlinearOptions(const Util::OptionBlock & OB)
{
  setEnableHomotopyFlag(true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerOutputOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/30/01
//-----------------------------------------------------------------------------
bool OutputMgr::registerOutputOptions(const Util::OptionBlock & option_block)
{

  // Need to capture a variable length of output intervals.  The use case is:
  //
  // .OPTIONS OUTPUT INITIAL_INTERVAL= <interval> [<t0> <i0> [<t1> <i1>...]]
  //
  // additional options can appear before or after the INITIAL_INTERVAL as in
  //
  //.OPTIONS OUTPUT HDF5FILE= xxxx INITIAL_INTERVAL=<interval> [<t0> <i0> [<t1> <i1>...]]
  // or
  // .OPTIONS OUTPUT INITIAL_INTERVAL= <interval> [<t0> <i0> [<t1> <i1>...]] HDF5FILE=xxxx
  //

  bool intervalSpecified = false;
  bool outputPointsSpecified = false;
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; )
  {
    if ((*it).tag() == "INITIAL_INTERVAL")
    {
      intervalSpecified = true;
      // start handling a list of intervals
      // set interval value
      initialOutputInterval_= (*it).getImmutableValue<double>();
      // look for optional time pairs.
      outputIntervalPairs_.clear();
      bool doneWithTimePairs= false;

      // need to point at next parameter
      ++it;

      while (it != end && !doneWithTimePairs)
      {
        if ((*it).tag() == "TIME")
        {
          double t = (*it).getImmutableValue<double>();
          ++it;

          double iv = 0.0;
          if (it != end) 
          {
            iv = (*it).getImmutableValue<double>();
            ++it;
          }
          else
          {
            Report::UserError0() << "Must specific intervals in pairs";
          }
          
          outputIntervalPairs_.push_back(Interval(t, iv));
        }
        else
        {
          // didn't find a time pair so bail out of this loop.
          doneWithTimePairs = true;
        }
      }
    }
    else if ((*it).tag() == "HDF5FILENAME")
    {
    // look for other option tags
      hdf5FileNameGiven_= true;
      hdf5FileName_= (*it).stringValue();
      ++it;
    }
    else if ((*it).tag()=="PRINTHEADER")
    {
      // look for flag to turn off the header
      printHeader_ = (*it).getImmutableValue<bool>();
      ++it;
    }
    else if ((*it).tag()=="PRINTFOOTER")
    {
      // look for flag to turn off "End of Xyce(TM) Simulation" line
      printFooter_= (*it).getImmutableValue<bool>();
      ++it;
    }
    else if ((*it).tag()=="ADD_STEPNUM_COL")
    {
      // look for flag to turn on STEPNUM column for FORMAT=GNUPLOT output
      printStepNumCol_= (*it).getImmutableValue<bool>();
      ++it;
    }
    else if ((*it).tag()=="OUTPUTVERSIONINRAWFILE")
    {
      // look for flag to toggle output of version in header of RAW file
     outputVersionInRawFile_ = (*it).getImmutableValue<bool>();
     ++it;
    }
    else if ((*it).tag()=="PHASE_OUTPUT_RADIANS")
    {
      // look for flag to toggle whether the VP() and IP() operators use radians
      // vs. degrees.  The default is radians now.
      phaseOutputUsesRadians_ = (*it).getImmutableValue<bool>();
      ++it;
    }

    else if ( std::string( (*it).uTag() ,0,16) == "OUTPUTTIMEPOINTS") // this is a vector
    {
      outputPointsSpecified = true;
      // do nothing
      ++it;
    }
    else
    {
      // silently ignore?
      ++it;
    }
  }

  if (outputPointsSpecified && intervalSpecified)
  {
    Report::UserError0() << "Cannot specify both .options output interval and outputtimepoints";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getVariableList
// Purpose       : Get a complete list of all variables requested on all
//                 print lines in the netlist.
// Special Notes : Used only in netlist import tool to check if any lead currents
//                 were specified, in order to set bools that determine whether
//                 to load lead current data.
// Scope         : public
// Creator       : David Baur
// Creation Date : 2014-06-23
//-----------------------------------------------------------------------------
Util::ParamList OutputMgr::getVariableList() const
{
  Util::ParamList parameter_list;

  for (OutputParameterMap::const_iterator it1 = outputParameterMap_.begin(), 
      end1 = outputParameterMap_.end(); it1 != end1; ++it1)
  {
    const OutputParameterMap::mapped_type &parameter_vector = (*it1).second;

    for (OutputParameterMap::mapped_type::const_iterator it2 = parameter_vector.begin(), 
        end2 = parameter_vector.end(); it2 != end2; ++it2)
    {
      const PrintParameters &print_parameters = (*it2);

      for (Util::ParamList::const_iterator it3 = print_parameters.variableList_.begin(), 
          end3 = print_parameters.variableList_.end(); it3 != end3; ++it3)
      {
        const Util::Param &parameter = (*it3);
        parameter_list.push_back(parameter);
      }
    }
  }
  // Now do the same thing for external output requests, which are done
  // through a wrapped interface object.
  for (ExternalOutputWrapperMap::const_iterator it = externalOutputWrapperMap_.begin();
         it != externalOutputWrapperMap_.end(); ++it)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it2=(*it).second.begin();
         it2 != (*it).second.end(); ++it2)
    {
      ExternalOutputWrapper * theExtOutWrapper=(*it2);
      Util::ParamList &theParamList = theExtOutWrapper->getParamList();
      for (Util::ParamList::const_iterator it3 = theParamList.begin();
           it3 != theParamList.end();
           ++it3)
      {
        const Util::Param &parameter = (*it3);
        parameter_list.push_back(parameter);
      }
    }
  }
  
  return parameter_list;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::checkPrintParameters
// Purpose       : Do some late stage checks of variables on print lines
// Special Notes : This is called during by "initializeLate" call of the
//                 Simulator class, after the problem is fully set up.  It
//                 performs a check on all print lines in the netlist by
//                 attempting to create Ops for every variable requested.
//                 makeOps will throw errors if it cannot form an op for any
//                 of the print variables.
//                 The ops created here are discarded after creation.
//
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 6/23/2014
//-----------------------------------------------------------------------------
void OutputMgr::checkPrintParameters(
  Parallel::Machine                     comm,
  const Util::Op::BuilderManager &      op_builder_manager)
{
  Util::Op::OpList tempOpList;

  // loop over all output parameter objects of all types, create
  // ops
  for (OutputParameterMap::const_iterator it = outputParameterMap_.begin(), 
      end = outputParameterMap_.end(); it != end; ++it)
  {
    for (std::vector<PrintParameters>::const_iterator it2 = (*it).second.begin(), 
        end2 = (*it).second.end(); it2 != end2; ++it2)
    {
      PrintParameters print_parameters = (*it2);
      fixupPrintParameters(comm, print_parameters);

      makeOps(comm, op_builder_manager, print_parameters.netlistLocation_,
              print_parameters.variableList_.begin(), 
              print_parameters.variableList_.end(), 
              std::back_inserter<Util::Op::OpList>(tempOpList));

    }
  }

  // Now do the same thing for external output requests, which are done
  // through a wrapped interface object.
  for (ExternalOutputWrapperMap::const_iterator it = externalOutputWrapperMap_.begin();
         it != externalOutputWrapperMap_.end(); ++it)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it2=(*it).second.begin();
         it2 != (*it).second.end(); ++it2)
    {
      ExternalOutputWrapper * theExtOutWrapper=(*it2);
      Util::ParamList &theParamList = theExtOutWrapper->getParamList();
      fixupOutputVariables(comm,theParamList);
      // We don't actually have a netlist line for these output requests,
      // so we fake one using the name returned by the external output
      // interface object
      makeOps(comm, op_builder_manager,
              NetlistLocation(theExtOutWrapper->getName(),0.0),
              theParamList.begin(), 
              theParamList.end(), 
              std::back_inserter<Util::Op::OpList>(tempOpList));
    }
  }

  // The sole purpose of having created the ops was to have checked the
  // validity of the output requests.  We don't really need them, so throw
  // them away.
  for (Util::Op::OpList::iterator it = tempOpList.begin(), end = tempOpList.end(); it != end; ++it)
  {
    delete *it;
  }

  if (hdf5FileNameGiven_)
  {
    prepareHDF5Output(comm);
  }
}

namespace {

// functor class for use with erase/remove_if idiom
struct IsFallback
{
  bool operator()(const PrintParameters &print_parameters)
  {
    return print_parameters.fallback_;
  }
};

// functor class for use with erase/remove_if idiom
struct IsNotFallback
{
  bool operator()(const PrintParameters &print_parameters)
  {
    return !print_parameters.fallback_;
  }
};

// Functor class for use with find_if
struct PrintNamesMatch
{
  PrintNamesMatch(const std::string & name)
    :Name(name)
  {}

  bool operator()(const PrintParameters &print_parameters)
  {
    return print_parameters.filename_ == Name;
  }
  
  const std::string &Name;
};
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addOutputPrintParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jul 23 07:21:53 2014
//-----------------------------------------------------------------------------
///
/// Adds the print parameters to the output type as follows:
///
/// First, check for an existing non-fallback print parameter for that output type.
///
/// Second, if the new print parameter is not a fallback parameter, then any
/// fallback print parameters (specific to that output type),for that output
/// type are removed.
///
/// Two types of print parameters are then inserted in the vector for the
/// specified output type:
///
///     a) non-fallback print parameters (e.g., an AC_IC print_parameters that
///        came from a .PRINT AC_IC line).
///
///     b) fallback print parameters (e.g., an AC_IC print_parameters that
///        was synthesized from a .PRINT AC line) if there are no
///        non-fallback print_parameters in the vector.
///
/// During that insertion, if a prior print parameters set is found in the vector
/// that matches the file name and file format of the given print parameter
/// object, then the variable list in the given object is appended to the existing
/// object rather than adding a new element to the vector.  In this way,
/// multiple print lines to the same file are the same as having one
/// large print line with an augmented variable list. 
///
/// @invariant
///
/// @param output_type          output type to add the print parameters too
/// @param print_parameters     print parameters to add
///
//-----------------------------------------------------------------------------
void OutputMgr::addOutputPrintParameters(
  OutputType::OutputType        output_type,
  const PrintParameters &       print_parameters)
{
  OutputParameterMap::mapped_type &print_parameters_vector = outputParameterMap_[output_type];

  // Determine if a "non-fallback" print_parameters is already in the print_parameters_vector.
  // An example is an AC_IC print_parameters that was synthesized from a .PRINT AC_IC line.
  bool nonFallbackFound=false;
  OutputParameterMap::mapped_type::iterator it=std::find_if(
       print_parameters_vector.begin(),
       print_parameters_vector.end(),
       IsNotFallback());

  if (it != print_parameters_vector.end())
  { 
    nonFallbackFound=true;
  }

  // If this a a "non-fallback" print_parameters then remove any previous "fallback"
  // parameters (e.g., an AC_IC print_parameters that was synthesized from the
  // variables on a .PRINT AC line)
  if (!print_parameters.fallback_)
      print_parameters_vector.erase(
         std::remove_if(print_parameters_vector.begin(),
                        print_parameters_vector.end(),
                        IsFallback()),
         print_parameters_vector.end());

  // Add the print_parameters if it is NOT a fallback.  Also add a fallback print_parameters
  // if there are NO fallback print_parameters in the vector for the specified output type.  
  if ( (!print_parameters.fallback_) || (print_parameters.fallback_ && !nonFallbackFound) )
  {   
    // Find if there's already a parameter list in the vector with the same
    // file name:
    OutputParameterMap::mapped_type::iterator it=std::find_if(
       print_parameters_vector.begin(),
       print_parameters_vector.end(),
       PrintNamesMatch(print_parameters.filename_));

    // If there is no such name, then add the vector
    if (it == print_parameters_vector.end())
    {
      print_parameters_vector.push_back(print_parameters);
    }
    else 
    {
      // we have found a prior parameter list with the same file name.
      // If the format matches, add this list to that pre-existing list.
      // If there is a format mismatch in a non-fallback, it is an error.
      // If the format mismatch is in a fallback, ignore both the error
      // and the attempted addition to the list.
      if ((*it).format_ != print_parameters.format_)
      {
        if (!(*it).fallback_)
        {
          
          Report::UserError() << "File " << print_parameters.filename_
                              << ": multiple print statements to this file with differing formats.";
        }
      }
      else
      {
        // Some variable lists for any print type (TRAN, DC, NOISE,
        // HOMOTOPY, SENS, AC or HB_*) may have had "STEPNUM", "INDEX",
	// "TIME" or "FREQ" pushed into the front of the
        // print_parameters.variableList_.  We do *NOT* want to add
        // them to (*it).variableList_.  (Note: the check for
        // print_parameters.variableList_.end() stops the while() loop
        // from segfaulting on a .PRINT line that doesn't have any
        // variables on it.  See SRN Bug 2073.)
        Util::ParamList::const_iterator it2=print_parameters.variableList_.begin();
        while ( (it2 != print_parameters.variableList_.end()) &&
                ( (*it2).tag() == "INDEX" ||
                  (*it2).tag() == "TIME" ||
                  (*it2).tag() == "FREQ"  ||
                  (*it2).tag() == "STEPNUM" ))
        {
          ++it2;
        }

        (*it).variableList_.insert((*it).variableList_.end(),
                                 it2,
                                 print_parameters.variableList_.end());
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addExternalOutputInterface
// Purpose       : Add an External Output Interface pointer to the list of
//                 outputters
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 26 Feb 2018
//-----------------------------------------------------------------------------
///
/// Wraps the interface in an ExternalOutputWrapper and adds it to the externalOutputWrapperMap in the appropriate slot
///
void OutputMgr::addExternalOutputInterface(ExternalOutputInterface * theOutputInterface)
{
  ExternalOutputWrapper * theOutputWrapper = new ExternalOutputWrapper(theOutputInterface);
  
  externalOutputWrapperMap_[theOutputWrapper->getOutputType()].push_back(theOutputWrapper);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::parsePRINTBlock
// Purpose       : Given an OptionBlock for a print statement, construct
//                 a PrintParameters object and add it to our set of known
//                 output requests.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
// Completely revamped by Dave Baur 06/23/2014-07/23/2014
//-----------------------------------------------------------------------------
///
/// Process an OptionBlock derived from a ".PRINT" line
///
/// @param print_block The option block to process
///
/// This function is registered with the options manager to process
/// options blocks named "PRINT".  These options blocks were produced
/// from .PRINT lines by extractPrintData and added to a circuit block.
///
/// This function's end result is to produce a PrintParameters object
/// representing the print line and add it to the output manager
/// so that it can be used to control output.
///
/// @note This function is the result of an interrupted refactor, and as
/// such is a tangled set of if statements to process all the different
/// types of .PRINT lines that might exist.  Some .PRINT lines are expected
/// to set up auxiliary output files, so while one print line results in
/// a single options block being created, it may generate more than one set
/// of print parameters.  This logic needs a redesign.
///
/// @note A lot of the complication of this function comes from the fact that
/// some output types have auxiliary outputs that may or may not be explicitly
/// requested by the user.  As an example, harmonic balance, when
/// specified as ".print hb" will output both frequency and time domain
/// solutions, each to its own file.  But one could also specify ".print hb_fd"
/// and ".print hb_td" lines separately.
/// When not explicitly requested, a "fall back" output specification for
/// auxiliary output is created.  Generally,
/// a fall-back print spec contains the same options and output variables
/// that the primary does.  If the user has specified explicit .PRINT lines
/// for these auxiliary outputs, then the "fall back" spec is discarded.
///
/// @author Dave Baur
/// @date 06/23/2014
bool OutputMgr::parsePRINTBlock(const Util::OptionBlock & print_block)
{
  PrintParameters        print_parameters = defaultPrintParameters_;
  PrintType::PrintType   print_type;
  Format::Format         format = Format::STD;
  DataFormat::DataFormat dataFormat = DataFormat::RI;
  bool                   no_index = false;

  print_parameters.netlistLocation_ = print_block.getNetlistLocation();

  // Decode each option specified in the options block and populate
  // the base print_parameters
  Util::ParamList::const_iterator iterParam = print_block.begin();
  for (; iterParam != print_block.end(); ++iterParam)
  {
    if (DEBUG_IO)
      Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;

    if (iterParam->tag() == "WIDTH")
    {
      print_parameters.streamWidth_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "TYPE")
    {
      std::string s = iterParam->stringValue();
      if (s == "AC")
        print_type = PrintType::AC;
      else if (s == "AC_IC")
        print_type = PrintType::AC_IC;
      else if (s == "SPARAM")
        print_type = PrintType::SPARAM;
      else if (s == "DC")
        print_type = PrintType::DC;
      else if (s == "ES")
        print_type = PrintType::ES;
      else if (s == "HB")
        print_type = PrintType::HB;
      else if (s == "HB_TD")
        print_type = PrintType::HB_TD;
      else if (s == "HB_FD")
        print_type = PrintType::HB_FD;
      else if (s == "HB_IC")
        print_type = PrintType::HB_IC;
      else if (s == "HB_STARTUP")
        print_type = PrintType::HB_STARTUP;
      else if (s == "HOMOTOPY")
        print_type = PrintType::HOMOTOPY;
      else if (s == "MPDE")
        print_type = PrintType::MPDE;
      else if (s == "MPDE_IC")
        print_type = PrintType::MPDE_IC;
      else if (s == "MPDE_STARTUP")
        print_type = PrintType::MPDE_STARTUP;
      else if (s == "PCE")
        print_type = PrintType::PCE;
      else if (s == "SENS")
        print_type = PrintType::SENS;
      else if (s == "TRANADJOINT")
        print_type = PrintType::TRANADJOINT;
      else if (s == "TRAN")
        print_type = PrintType::TRAN;
      else if (s == "NOISE")
        print_type = PrintType::NOISE;
      else if (s == "ROL")
        print_type = PrintType::DC; // TT
      else
      {
        Report::DevelFatal0() << "Unrecognized analysis type " << s;
      }
    }
    else if (iterParam->tag() == "PRECISION")
    {
      print_parameters.streamPrecision_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "FILTER")
    {
      print_parameters.filter_ = iterParam->getImmutableValue<double>();
    }
    else if (iterParam->tag() == "FORMAT")
    {
      std::string s = iterParam->stringValue();
      if (s == "STD")
        format = Format::STD;
      else if (s == "TECPLOT")
        format = Format::TECPLOT;
      else if (s == "DAKOTA")
        format = Format::DAKOTA;
      else if (s == "PROBE")
        format = Format::PROBE;
      else if (s == "NOINDEX")
      {
        format = Format::STD;
        no_index = true;
      }
      else if (s == "GNUPLOT")
      {
        format = Format::STD;
        print_parameters.addGnuplotSpacing_ = true;
      }
      else if (s == "SPLOT")
      {
        format = Format::STD;
        print_parameters.addSplotSpacing_ = true;
      }
      else if (s == "CSV")
      {
        format = Format::CSV;
        print_parameters.delimiter_ = ",";
      }
      else if (s == "RAW")
      {
        format = defaultPrintParameters_.asciiRaw_ ? Format::RAW_ASCII : Format::RAW;
      }
      else if (s == "TOUCHSTONE")
      {
        format = Format::TS1;
      }
      else if (s == "TOUCHSTONE2")
      {
        format = Format::TS2;
      }
      else
      {
        Report::DevelFatal0() << "Unrecognized print format " << s;
      }
      print_parameters.format_ = format;
    }
    else if (iterParam->tag() == "DATAFORMAT")
    {
      // data format for Touchstone output
      std::string s = iterParam->stringValue();
      if (s == "RI")
        dataFormat = DataFormat::RI;
      else if (s == "MA")
        dataFormat = DataFormat::MA;
      else if (s == "DB")
        dataFormat = DataFormat::DB;
       else
      {
        Report::DevelFatal0() << "Unrecognized data format " << s <<  " for Touchstone output requested on .LIN ";
      }
      print_parameters.dataFormat_ = dataFormat;
    }
    else if (iterParam->tag() == "LINTYPE")
    {
      // Parameter type in Touchstone output
      std::string s = iterParam->stringValue();
      if ( (s == "S") || (s == "Y") || (s == "Z") )
      {
        print_parameters.RFparamType_ = s;
      }
      else
      {
        Report::DevelFatal0() << "Unrecognized or unsupported parameter type " << s << " for Touchstone output requested on .LIN";
      }
    }
    else if (iterParam->tag() == "TIMEWIDTH")
    {
      print_parameters.timeWidth_ = iterParam->getImmutableValue<int>();
    }
    else if (iterParam->tag() == "TIMESCALEFACTOR")
    {
      print_parameters.outputTimeScaleFactor_ = iterParam->getImmutableValue<double>();
    }
    else if (iterParam->tag() == "FILE")
    {
      // netlistFilename_ should be the default unless FILE was set to
      // something other than "" which is its default value.
      if (iterParam->stringValue() != "")
      {
        print_parameters.filename_ = iterParam->stringValue();
      }
    }
    else if (iterParam->tag() == "DELIMITER")
    {
      if (iterParam->stringValue() == "TAB")
      {
        print_parameters.delimiter_ = "\t";
      }
      else if (iterParam->stringValue() == "COMMA")
      {
        print_parameters.delimiter_ = ",";
      }
      else if (iterParam->stringValue() != "")
      {
        Report::UserWarning0() << "Invalid value of DELIMITER in .PRINT statment, ignoring";
      }
    }
    else if (iterParam->tag() == "OUTPUTSAMPLESTATS")
    {
      print_parameters.outputPCEsampleStats_ = static_cast<bool>(iterParam->getImmutableValue<bool>());
    }
    else if (iterParam->tag() == "OUTPUTALLSAMPLES")
    {
      print_parameters.outputAllPCEsamples_ = static_cast<bool>(iterParam->getImmutableValue<bool>());
    }
    else if (iterParam->tag() == "OUTPUT_PCE_COEFFS")
    {
      print_parameters.outputPCECoeffs_ = static_cast<bool>(iterParam->getImmutableValue<bool>());
    }
    else
    {
      // This must be the first print variable.
      break;
    }
  }

  // Remaining parameters in the options block are variables,
  // copy them all into the print_parameters
  print_parameters.variableList_.assign(iterParam, print_block.end());

  // Increase width to make sure that the columns do not run together
  if (print_parameters.streamWidth_ - print_parameters.streamPrecision_  < 9)
    print_parameters.streamWidth_ = print_parameters.streamPrecision_ + 9;

  // Indicate if an index column should be added
  print_parameters.printIndexColumn_ = !no_index && print_parameters.format_ == Format::STD;

  // Indicates if a stepnum column should be added
  print_parameters.printStepNumColumn_ = printStepNumCol_ && print_parameters.format_ == Format::STD;

  // -r output is not made for these formats
  if ( (print_type == PrintType::SENS) || 
       (print_type == PrintType::HOMOTOPY) || 
       (print_type == PrintType::SPARAM) ||
       (print_type == PrintType::ES) ||
       (print_type == PrintType::PCE) )
  {
    print_parameters.formatSupportsOverrideRaw_= false; 
  }

  // Assemble the apropriate flavors of output variable lists based on the PRINT type
  // and format.  Notes are:
  // 
  //   1) Certain print types cause output of columns (e.g., INDEX) that the user
  // has not explicitly requested.  This is where we add those columns to the
  // print parameters objects.
  //
  //   2) Other print types can cause the generation of multiple output files (e.g.,
  //      frequency domain and time domain output for .PRINT AC).  This is where we 
  //      create extra print parameters objects to handle those additional outputs.
  //
  //   3) The first if clause is the simplified processing done when the -o command line
  // option is specified.  Otherwise, more complicated processing (possibly with the
  // generation of "fallback print parameters) occurs based on the PrintType.
  if (print_parameters.dashoRequested_)
  {
    // -o only produces output for .PRINT AC, .PRINT DC, .PRINT HB_FD, .PRINT NOISE and
    // .PRINT TRAN lines.  It generates a Touchstone 2 file if a .LIN analysis is done.
    // -o output defaults to Format::STD, with an INDEX column and space as the delimiter_.
    print_parameters.printIndexColumn_ = true;
    print_parameters.format_ = Format::STD;
    if (printStepNumCol_)
      print_parameters.printStepNumColumn_ = true;
    print_parameters.delimiter_ = "";

    // Xyce behavior changes request for output file <netlistName>.cir to 
    // the file <netlistName>.cir.<print_parameters.defaultExtension_>
    print_parameters.defaultExtension_=".prn";

    // Reset any FILE= specified on an individual .PRINT line back to the default (blank), 
    // so that .PRINT line concatenation still works correctly in addOutputPrintParameters()
    print_parameters.filename_ = "";  

    if (print_type == PrintType::AC)
    {
      PrintParameters freq_print_parameters = print_parameters;
      freq_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (freq_print_parameters.printStepNumColumn_)
        freq_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
      freq_print_parameters.expandComplexTypes_ = true;
      addOutputPrintParameters(OutputType::AC, freq_print_parameters); 
    }
    else if (print_type == PrintType::ES)
    {
      PrintParameters es_print_parameters = print_parameters;
      addOutputPrintParameters(OutputType::ES, es_print_parameters);
    }
    else if ( (print_type == PrintType::HB) || (print_type == PrintType::HB_FD) )
    {
      PrintParameters freq_print_parameters = print_parameters;
      freq_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (freq_print_parameters.printStepNumColumn_)
        freq_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
      freq_print_parameters.expandComplexTypes_ = true;
      addOutputPrintParameters(OutputType::HB_FD, freq_print_parameters);
    }
    else if (print_type == PrintType::NOISE)
    {  
      PrintParameters noise_print_parameters = print_parameters;
      noise_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
      noise_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (noise_print_parameters.printStepNumColumn_)
        noise_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
      noise_print_parameters.expandComplexTypes_ = true;
      std::copy(noiseVariableList_.begin(), noiseVariableList_.end(), std::back_inserter(noise_print_parameters.variableList_));    
      addOutputPrintParameters(OutputType::NOISE, noise_print_parameters);
    }
    else if (print_type == PrintType::PCE)
    {
      PrintParameters pce_print_parameters = print_parameters;
      addOutputPrintParameters(OutputType::PCE, pce_print_parameters);
    }
    else if (print_type == PrintType::SPARAM)
    {
      // Note: Both SPARAM and AC output print parameters are made.  The function
      // OutputMgr::prepareOutput() will then control which outputter is actually used,
      // for the -o case, based on whether a .LIN analysis is being done, or not.
      PrintParameters sparam_print_parameters = print_parameters;
      sparam_print_parameters.format_=Format::TS2;
      addOutputPrintParameters(OutputType::SPARAM, sparam_print_parameters);
    }
    else if (print_type == PrintType::TRAN)
    {
      // OutputType::ES takes precedence over OutputType::TRAN.  This will be
      // enforced in OutputMgr:prepareOutput()
      addOutputPrintParameters(OutputType::TRAN, print_parameters);
    }
    else if (print_type == PrintType::DC)
    {
      // OutputType::ES takes precedence over OutputType::DC.  This will be
      // enforced in OutputMgr:prepareOutput()
      addOutputPrintParameters(OutputType::DC, print_parameters);
    }
    else
    {
      Report::UserWarning0() << "-o only produces output for .PRINT AC, .PRINT DC, .PRINT ES, "
                             << ".PRINT NOISE, .PRINT PCE, .PRINT TRAN, .PRINT HB, "
                             << ".PRINT HB_FD and .LIN lines";
    }
  }
  else if (print_type == PrintType::AC)
  {
    PrintParameters freq_print_parameters = print_parameters;
    if (freq_print_parameters.format_ != Format::PROBE)
    {
      freq_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
    }

    if ( (freq_print_parameters.format_ == Format::TS1) ||
         (freq_print_parameters.format_ == Format::TS2) )
    {
      freq_print_parameters.defaultExtension_ = ".FD.prn";
      // print out the Index column, since this will be in STD format
      freq_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        freq_print_parameters.printStepNumColumn_ = true;
    }

    if (freq_print_parameters.printIndexColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (freq_print_parameters.printStepNumColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;

    addOutputPrintParameters(OutputType::AC, freq_print_parameters);

    // Only create the fallback AC_IC output if the netlist has a .OP
    // statement in it.  See SON Bug 990.
    if (dotOpSpecified_)
    {
      PrintParameters ac_ic_print_parameters = print_parameters;
      ac_ic_print_parameters.fallback_ = true;
      if (ac_ic_print_parameters.format_ != Format::PROBE)
      {
        ac_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      }

      if ( (ac_ic_print_parameters.format_ == Format::TS1) ||
           (ac_ic_print_parameters.format_ == Format::TS2) )
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.prn";
        // print out the Index column, since this will be in STD format
        ac_ic_print_parameters.printIndexColumn_ = true;
        if (printStepNumCol_)
          ac_ic_print_parameters.printStepNumColumn_ = true;
      }

      if (ac_ic_print_parameters.printIndexColumn_)
        ac_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (ac_ic_print_parameters.printStepNumColumn_)
        ac_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      addOutputPrintParameters(OutputType::AC_IC, ac_ic_print_parameters);
    }
  }
  else if (print_type == PrintType::AC_IC)
  {
    // Only create the AC_IC output if the netlist has a .OP statement in it.
    // See SON Bug 990.
    if (dotOpSpecified_)
    {
      PrintParameters ac_ic_print_parameters = print_parameters;
      if (ac_ic_print_parameters.format_ != Format::PROBE)
      {
        ac_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      }

      if ( (ac_ic_print_parameters.format_ == Format::TS1) ||
           (ac_ic_print_parameters.format_ == Format::TS2) )
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.prn";
        // print out the Index column, since this will be in STD format
        ac_ic_print_parameters.printIndexColumn_ = true;
        if (printStepNumCol_)
          ac_ic_print_parameters.printStepNumColumn_ = true;
      }

      if (ac_ic_print_parameters.printIndexColumn_)
        ac_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (ac_ic_print_parameters.printStepNumColumn_)
        ac_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      addOutputPrintParameters(OutputType::AC_IC, ac_ic_print_parameters);
    }
  }
  else if (print_type == PrintType::SPARAM)
  {
    PrintParameters sparam_print_parameters = print_parameters;
    sparam_print_parameters.defaultExtension_ = ".sp";
    addOutputPrintParameters(OutputType::SPARAM, sparam_print_parameters);
  }
  else if (print_type == PrintType::NOISE)
  {
    PrintParameters noise_print_parameters = print_parameters;
    noise_print_parameters.expandComplexTypes_ = noise_print_parameters.format_ != Format::PROBE
                                                && noise_print_parameters.format_ != Format::RAW
                                                && noise_print_parameters.format_ != Format::RAW_ASCII;

    if (noise_print_parameters.format_ == Format::STD)
    {
      noise_print_parameters.defaultExtension_ = ".NOISE.prn";
    }
    else if (noise_print_parameters.format_ == Format::CSV)
    {
      noise_print_parameters.defaultExtension_ = ".NOISE.csv";
    }
    else if (noise_print_parameters.format_ == Format::TECPLOT)
    {
      noise_print_parameters.defaultExtension_ = ".NOISE.dat";
    }
    //else if (noise_print_parameters.format_ == Format::PROBE)
    //{
    //  noise_print_parameters.defaultExtension_ = ".NOISE.csd";
    //}
    //  else if (noise_print_parameters.format_ == Format::DAKOTA)
    //{
    //  noise_print_parameters.defaultExtension_ = ".NOISE.txt";
    //}
    else if ( (noise_print_parameters.format_ == Format::RAW) ||
              (noise_print_parameters.format_ == Format::RAW_ASCII) ||
              (noise_print_parameters.format_ == Format::PROBE) ||
              (noise_print_parameters.format_ == Format::TS1) ||
              (noise_print_parameters.format_ == Format::TS2) )
    {
      noise_print_parameters.defaultExtension_ = ".NOISE.prn";
      // print out the Index column, since this will be in STD format
      noise_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        noise_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      noise_print_parameters.defaultExtension_ = ".NOISE.unknown";
    }

    // adjust which columns appear in the output file, depending on the print format
    noise_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
    if (noise_print_parameters.printIndexColumn_)
      noise_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (noise_print_parameters.printStepNumColumn_)
      noise_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

    std::copy(noiseVariableList_.begin(), noiseVariableList_.end(), std::back_inserter(noise_print_parameters.variableList_));

    addOutputPrintParameters(OutputType::NOISE, noise_print_parameters);

  }
  else if (print_type == PrintType::HOMOTOPY)
  {
    PrintParameters homotopy_print_parameters = print_parameters;
    if (homotopy_print_parameters.format_ == Format::STD)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.prn";
    }
    else if (homotopy_print_parameters.format_ == Format::CSV)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csv";
    }
    else if (homotopy_print_parameters.format_ == Format::TECPLOT)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.dat";
    }
    else if (homotopy_print_parameters.format_ == Format::PROBE)
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csd";
    }
    else if ( (homotopy_print_parameters.format_ == Format::RAW) ||
              (homotopy_print_parameters.format_ == Format::RAW_ASCII) ||
              (homotopy_print_parameters.format_ == Format::TS1) ||
              (homotopy_print_parameters.format_ == Format::TS2) )
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.prn";
      // print out the Index column
      homotopy_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        homotopy_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.unknown";
    }
    addOutputPrintParameters(OutputType::HOMOTOPY, homotopy_print_parameters);
  }
  else if (print_type == PrintType::SENS)
  {
    PrintParameters sensitivity_print_parameters = print_parameters;

    if (sensitivity_print_parameters.format_ == Format::STD)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
    }
    else if (sensitivity_print_parameters.format_ == Format::CSV)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.csv";
    }
    else if (sensitivity_print_parameters.format_ == Format::TECPLOT)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.dat";
    }
    else if (sensitivity_print_parameters.format_ == Format::DAKOTA)
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.txt";
    }
    else if ( (sensitivity_print_parameters.format_ == Format::RAW) ||
              (sensitivity_print_parameters.format_ == Format::RAW_ASCII) ||
              (sensitivity_print_parameters.format_ == Format::PROBE) ||
              (sensitivity_print_parameters.format_ == Format::TS1) ||
              (sensitivity_print_parameters.format_ == Format::TS2) )
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
      // print out the Index column, since this will be in STD format
      sensitivity_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        sensitivity_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.unknown";
    }
    std::copy(sensitivityVariableList_.begin(), sensitivityVariableList_.end(), std::back_inserter(sensitivity_print_parameters.variableList_));

    addOutputPrintParameters(OutputType::SENS, sensitivity_print_parameters);
  }
  else if (print_type == PrintType::TRANADJOINT)
  {
    PrintParameters transientAdjoint_print_parameters = print_parameters;

    if (transientAdjoint_print_parameters.format_ == Format::STD)
    {
      transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.prn";
    }
    else if (transientAdjoint_print_parameters.format_ == Format::TECPLOT)
    {
      transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.dat";
    }
    else if ( (transientAdjoint_print_parameters.format_ == Format::RAW) ||
              (transientAdjoint_print_parameters.format_ == Format::RAW_ASCII) ||
              (transientAdjoint_print_parameters.format_ == Format::PROBE) ||
              (transientAdjoint_print_parameters.format_ == Format::TS1) ||
              (transientAdjoint_print_parameters.format_ == Format::TS2) )
    {
      transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.prn";
      // print out the Index column, since this will be in STD format
      transientAdjoint_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        transientAdjoint_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.unknown";
    }      
    // ERK: for now just using the same variable list, but change later.
    std::copy(
        sensitivityVariableList_.begin(), 
        sensitivityVariableList_.end(), 
        std::back_inserter(transientAdjoint_print_parameters.variableList_));

    addOutputPrintParameters(OutputType::TRANADJOINT, transientAdjoint_print_parameters);
  }
  else if (print_type == PrintType::HB)
  {
    // create fallback print parameters for HB_FD files
    PrintParameters freq_print_parameters = print_parameters;
    freq_print_parameters.fallback_ = true;
    // set parameters based on FORMAT=<val> on .PRINT HB_FD line
    update_HB_FD_printParamsForPrintFormat(freq_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported
    //if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
    if (freq_print_parameters.printIndexColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (freq_print_parameters.printStepNumColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_FD, freq_print_parameters);

    // create fallback print parameters for HB_TD files
    PrintParameters time_print_parameters = print_parameters;
    time_print_parameters.fallback_ = true;
    // set parameters based on FORMAT=<val> on .PRINT HB_TD line
    update_HB_TD_printParamsForPrintFormat(time_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported
    //if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (time_print_parameters.printIndexColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (time_print_parameters.printStepNumColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_TD, time_print_parameters);

    // create fallback print parameters for hb_ic files
    PrintParameters hb_ic_print_parameters = print_parameters;
    hb_ic_print_parameters.fallback_ = true;

    // set parameters based on FORMAT=<val> on .PRINT HB_IC line
    update_HB_IC_printParamsForPrintFormat(hb_ic_print_parameters);
  
    // uncomment this if statement when FORMAT=PROBE is supported
    //if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_ic_print_parameters.printIndexColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (hb_ic_print_parameters.printStepNumColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_IC, hb_ic_print_parameters);

    // create fallback print parameters for startup files
    PrintParameters hb_startup_print_parameters = print_parameters;
    hb_startup_print_parameters.fallback_ = true;

    // set parameters based on FORMAT=<val> on .PRINT HB_STARTUP line
    update_HB_STARTUP_printParamsForPrintFormat(hb_startup_print_parameters);
        
    // uncomment this if statement when FORMAT=PROBE is supported  
    // if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_startup_print_parameters.printIndexColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (hb_startup_print_parameters.printStepNumColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_STARTUP, hb_startup_print_parameters);
  }
  else if (print_type == PrintType::HB_TD)
  {
    PrintParameters time_print_parameters = print_parameters;

    // set parameters based on FORMAT=<val> on .PRINT HB_TD line
    update_HB_TD_printParamsForPrintFormat(time_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported  
    // if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (time_print_parameters.printIndexColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (time_print_parameters.printStepNumColumn_)
      time_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_TD, time_print_parameters);
  }
  else if (print_type == PrintType::HB_FD)
  {
    PrintParameters freq_print_parameters = print_parameters;

    // set parameters based on FORMAT=<val> on .PRINT HB_FD line
    update_HB_FD_printParamsForPrintFormat(freq_print_parameters);
        
    // uncomment this if statement when FORMAT=PROBE is supported  
    // if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
    if (freq_print_parameters.printIndexColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (freq_print_parameters.printStepNumColumn_)
      freq_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_FD, freq_print_parameters);
  }
  else if (print_type == PrintType::HB_IC)
  {
    PrintParameters hb_ic_print_parameters = print_parameters;

    // set parameters based on FORMAT=<val> on .PRINT HB_IC line
    update_HB_IC_printParamsForPrintFormat(hb_ic_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported  
    //if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_ic_print_parameters.printIndexColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (hb_ic_print_parameters.printStepNumColumn_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_IC, hb_ic_print_parameters);
  }
  else if (print_type == PrintType::HB_STARTUP)
  {
    PrintParameters hb_startup_print_parameters = print_parameters;

    // set parameters based on FORMAT=<val> on .PRINT HB_STARTUP line
    update_HB_STARTUP_printParamsForPrintFormat(hb_startup_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported    
    //if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_startup_print_parameters.printIndexColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (hb_startup_print_parameters.printStepNumColumn_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::HB_STARTUP, hb_startup_print_parameters);
  }
  else if (print_type == PrintType::MPDE)
  {
    PrintParameters mpde_print_parameters = print_parameters;
    if (mpde_print_parameters.format_ == Format::STD)
    {
      mpde_print_parameters.defaultExtension_ = ".MPDE.prn";
    }
    else if (mpde_print_parameters.format_ == Format::CSV)
    {
      mpde_print_parameters.defaultExtension_ = ".MPDE.csv";
    }
    else if (mpde_print_parameters.format_ == Format::TECPLOT)
    {
      mpde_print_parameters.defaultExtension_ = ".MPDE.dat";
    }
    else if ( (mpde_print_parameters.format_ == Format::PROBE) ||
              (mpde_print_parameters.format_ == Format::RAW) ||
              (mpde_print_parameters.format_ == Format::RAW_ASCII) ||
              (mpde_print_parameters.format_ == Format::TS1) ||
              (mpde_print_parameters.format_ == Format::TS2) )
    {
      mpde_print_parameters.defaultExtension_ = ".MPDE.prn";
    }
    else
    {
      mpde_print_parameters.defaultExtension_ = ".MPDE.unknown";
    }
    addOutputPrintParameters(OutputType::MPDE, mpde_print_parameters);

    // add fallback for MPDE_IC output
    PrintParameters mpde_ic_print_parameters = print_parameters;
    mpde_ic_print_parameters.fallback_ = true;

    // set parameters based on FORMAT=<val> on .PRINT MPDE_IC line
    update_MPDE_IC_printParamsForPrintFormat(mpde_ic_print_parameters);

    // uncomment this if statement when FORMAT=PROBE is supported 
    //if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.printIndexColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (mpde_ic_print_parameters.printStepNumColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::MPDE_IC, mpde_ic_print_parameters);

    // add fallback for MPDE_STARTUP output
    PrintParameters mpde_startup_print_parameters = print_parameters;
    mpde_startup_print_parameters.fallback_ = true;

    // set parameters based on FORMAT=<val> on .PRINT MPDE_STARTUP line
    update_MPDE_STARTUP_printParamsForPrintFormat(mpde_startup_print_parameters);
      
    // uncomment this if statement when FORMAT=PROBE is supported 
    //if (mpde_startup_print_parameters.format_ != Format::PROBE)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_startup_print_parameters.printIndexColumn_)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (mpde_startup_print_parameters.printStepNumColumn_)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::MPDE_STARTUP, mpde_startup_print_parameters);
  }
  else if (print_type == PrintType::MPDE_IC)
  {
    PrintParameters mpde_ic_print_parameters = print_parameters;
    // set parameters based on FORMAT=<val> on .PRINT MPDE_IC line
    update_MPDE_IC_printParamsForPrintFormat(mpde_ic_print_parameters);
      
    // uncomment this if statement when FORMAT=PROBE is supported 
    //if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.printIndexColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (mpde_ic_print_parameters.printStepNumColumn_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::MPDE_IC, mpde_ic_print_parameters);
  }
  else if (print_type == PrintType::MPDE_STARTUP)
  {
    PrintParameters mpde_startup_print_parameters = print_parameters;
    // set parameters based on FORMAT=<val> on .PRINT MPDE_STARTUP line
    update_MPDE_STARTUP_printParamsForPrintFormat(mpde_startup_print_parameters);
      
    // uncomment this if statement when FORMAT=PROBE is supported 
    //if (mpde_startup_print_parameters.format_ != Format::PROBE)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_startup_print_parameters.printIndexColumn_)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    if (mpde_startup_print_parameters.printStepNumColumn_)
      mpde_startup_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));
    addOutputPrintParameters(OutputType::MPDE_STARTUP, mpde_startup_print_parameters);
  }
  else if (print_type == PrintType::TRAN)
  {
    if ( (print_parameters.format_ == Format::TS1) ||
         (print_parameters.format_ == Format::TS2) )
    {
      // print out the Index column, since this will be in STD format
      print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        print_parameters.printStepNumColumn_ = true;
    }
    addOutputPrintParameters(OutputType::TRAN, print_parameters);
  }
  else if (print_type == PrintType::DC)
  {
    PrintParameters dc_print_parameters = print_parameters;
    if (dc_print_parameters.format_ == Format::STD)
    {
      dc_print_parameters.defaultExtension_ = ".prn";
    }
    else if ((print_parameters.format_ == Format::TS1) ||
             (print_parameters.format_ == Format::TS2) )
    {
      dc_print_parameters.defaultExtension_ = ".prn";
      // print out the Index column, since this will be in STD format
      dc_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        dc_print_parameters.printStepNumColumn_ = true;
    }
    addOutputPrintParameters(OutputType::DC, dc_print_parameters);
  }
  else if (print_type == PrintType::ES)
  {
    PrintParameters es_print_parameters = print_parameters;
    if (es_print_parameters.format_ == Format::STD)
    {
      es_print_parameters.defaultExtension_ = "ES.prn";
    }
    else if (es_print_parameters.format_ == Format::CSV)
    {
      es_print_parameters.defaultExtension_ = ".ES.csv";
    }
    else if (es_print_parameters.format_ == Format::TECPLOT)
    {
      es_print_parameters.defaultExtension_ = ".ES.dat";
    }
    else if ((es_print_parameters.format_ == Format::RAW) ||
             (es_print_parameters.format_ == Format::RAW_ASCII) ||
             (es_print_parameters.format_ == Format::PROBE) ||
             (es_print_parameters.format_ == Format::TS1) ||
             (es_print_parameters.format_ == Format::TS2) )
    {
      es_print_parameters.defaultExtension_ = "ES.prn";
      // print out the Index column, since this will be in STD format
      es_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        es_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      es_print_parameters.defaultExtension_ = ".ES.unknown";
    }
    addOutputPrintParameters(OutputType::ES, es_print_parameters);
  }
  else if (print_type == PrintType::PCE)
  {
    PrintParameters pce_print_parameters = print_parameters;
    if (pce_print_parameters.format_ == Format::STD)
    {
      pce_print_parameters.defaultExtension_ = "PCE.prn";
    }
    else if (pce_print_parameters.format_ == Format::CSV)
    {
      pce_print_parameters.defaultExtension_ = ".PCE.csv";
    }
    else if (pce_print_parameters.format_ == Format::TECPLOT)
    {
      pce_print_parameters.defaultExtension_ = ".PCE.dat";
    }
    else if ((pce_print_parameters.format_ == Format::RAW) ||
             (pce_print_parameters.format_ == Format::RAW_ASCII) ||
             (pce_print_parameters.format_ == Format::PROBE) ||
             (pce_print_parameters.format_ == Format::TS1) ||
             (pce_print_parameters.format_ == Format::TS2) )
    {
      pce_print_parameters.defaultExtension_ = "PCE.prn";
      // print out the Index column, since this will be in STD format
      pce_print_parameters.printIndexColumn_ = true;
      if (printStepNumCol_)
        pce_print_parameters.printStepNumColumn_ = true;
    }
    else
    {
      pce_print_parameters.defaultExtension_ = ".PCE.unknown";
    }
    addOutputPrintParameters(OutputType::PCE, pce_print_parameters);
  }
  else
  {
    Report::UserError0() << "Unrecognized .PRINT type";
  }

  pushActiveOutputters();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_HB_FD_printParamsForPrintFormat
// Purpose       : update some of the print_parameters for .PRINT HB_FD based
//                 on FORMAT= on the .PRINT Line. This is used for both
//                 fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/17/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_HB_FD_printParamsForPrintFormat(PrintParameters &freq_print_parameters)
{
  // Print out real and imaginary part of voltages and currents for most formats.
  // The exception would be FORMAT=PROBE when/if that format is supported for HB output.
  freq_print_parameters.expandComplexTypes_ = true;

  if (freq_print_parameters.format_ == Format::STD)
  {
    freq_print_parameters.defaultExtension_ = ".HB.FD.prn";
  }
  else if (freq_print_parameters.format_ == Format::CSV)
  {
    freq_print_parameters.defaultExtension_ = ".HB.FD.csv";
  }
  else if (freq_print_parameters.format_ == Format::TECPLOT)
  {
    freq_print_parameters.defaultExtension_ = ".HB.FD.dat";
  }
  else if ( (freq_print_parameters.format_ == Format::RAW) ||
            (freq_print_parameters.format_ == Format::RAW_ASCII) ||
            (freq_print_parameters.format_ == Format::PROBE) ||
            (freq_print_parameters.format_ == Format::TS1) ||
            (freq_print_parameters.format_ == Format::TS2) )
  {
    freq_print_parameters.defaultExtension_ = ".HB.FD.prn";
    // print out the Index column, since this will be in STD format
    freq_print_parameters.printIndexColumn_ = true;

    // variable needs to be set to false, if these format are supported
    //freq_print_parameters.expandComplexTypes_ = false;
  }
  else
  {
    freq_print_parameters.defaultExtension_ = ".HB.FD.unknown";   
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_HB_TD_printParamsForPrintFormat
// Purpose       : update some of the print_parameters for .PRINT HB_TD based
//                 on FORMAT= on the .PRINT Line.  This is used for both
//                 fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/17/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_HB_TD_printParamsForPrintFormat(PrintParameters &time_print_parameters)
{
  if (time_print_parameters.format_ == Format::STD)
  {
    time_print_parameters.defaultExtension_ = ".HB.TD.prn";
  }
  else if (time_print_parameters.format_ == Format::CSV)
  {
    time_print_parameters.defaultExtension_ = ".HB.TD.csv";
  }
  else if (time_print_parameters.format_ == Format::TECPLOT)
  {
    time_print_parameters.defaultExtension_ = ".HB.TD.dat";
  }
  else if ( (time_print_parameters.format_ == Format::RAW) ||
            (time_print_parameters.format_ == Format::RAW_ASCII) ||
            (time_print_parameters.format_ == Format::PROBE) ||
            (time_print_parameters.format_ == Format::TS1) ||
            (time_print_parameters.format_ == Format::TS2) )
  {
    time_print_parameters.defaultExtension_ = ".HB.TD.prn";
    // print out the Index column, since this will be in STD format
    time_print_parameters.printIndexColumn_ = true;
  }
  else
  {
    time_print_parameters.defaultExtension_ = ".HB.TD.unknown";   
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_HB_IC_printParamsForPrintFormat
// Purpose       : update some of the print_parameters for .PRINT HB_IC based
//                 on FORMAT= on the .PRINT Line.  This is used for both
//                 fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/17/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_HB_IC_printParamsForPrintFormat(PrintParameters &hb_ic_print_parameters)
{
  if (hb_ic_print_parameters.format_ == Format::STD)
  {
    hb_ic_print_parameters.defaultExtension_ = ".hb_ic.prn";
  }
  else if (hb_ic_print_parameters.format_ == Format::CSV)
  {
    hb_ic_print_parameters.defaultExtension_ = ".hb_ic.csv";
  }
  else if (hb_ic_print_parameters.format_ == Format::TECPLOT)
  {
    hb_ic_print_parameters.defaultExtension_ = ".hb_ic.dat";
  }
  else if ( (hb_ic_print_parameters.format_ == Format::RAW) ||
            (hb_ic_print_parameters.format_ == Format::RAW_ASCII) ||
            (hb_ic_print_parameters.format_ == Format::PROBE) ||
            (hb_ic_print_parameters.format_ == Format::TS1) ||
            (hb_ic_print_parameters.format_ == Format::TS2) )
  {
    hb_ic_print_parameters.defaultExtension_ = ".hb_ic.prn";
    // print out the Index column, since this will be in STD format
    hb_ic_print_parameters.printIndexColumn_ = true;
  }
  else
  {
    hb_ic_print_parameters.defaultExtension_ = ".hb_ic.unknown";   
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_HB_STARTUP_printParamsForPrintFormat
// Purpose       : Update some of the print_parameters for .PRINT HB_STARTUP based
//                 on FORMAT= on the .PRINT Line.  This is used for both
//                 fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/17/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_HB_STARTUP_printParamsForPrintFormat(PrintParameters &hb_startup_print_parameters)
{
  if (hb_startup_print_parameters.format_ == Format::STD)
  {
    hb_startup_print_parameters.defaultExtension_ = ".startup.prn";
  }
  else if (hb_startup_print_parameters.format_ == Format::CSV)
  {
    hb_startup_print_parameters.defaultExtension_ = ".startup.csv";
  }
  else if (hb_startup_print_parameters.format_ == Format::TECPLOT)
  {
    hb_startup_print_parameters.defaultExtension_ = ".startup.dat";
  }
  else if ( (hb_startup_print_parameters.format_ == Format::RAW) ||
            (hb_startup_print_parameters.format_ == Format::RAW_ASCII) ||
            (hb_startup_print_parameters.format_ == Format::PROBE) ||
            (hb_startup_print_parameters.format_ == Format::TS1) ||
            (hb_startup_print_parameters.format_ == Format::TS2) )
  {
    hb_startup_print_parameters.defaultExtension_ = ".startup.prn";
    // print out the Index column, since this will be in STD format
    hb_startup_print_parameters.printIndexColumn_ = true;
  }
  else
  {
    hb_startup_print_parameters.defaultExtension_ = ".startup.unknown";    
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_MPDE_IC_printParamsForPrintFormat
// Purpose       : update some of the print_parameters for .PRINT MPDE_IC
//                 based on FORMAT= on the .PRINT Line.  This is used for
//                 both fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/28/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_MPDE_IC_printParamsForPrintFormat(PrintParameters &mpde_ic_print_parameters)
{
  if (mpde_ic_print_parameters.format_ == Format::STD)
  {
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";
  }
  else if (mpde_ic_print_parameters.format_ == Format::CSV)
  {
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.csv";
  }
  else if (mpde_ic_print_parameters.format_ == Format::TECPLOT)
  {
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.dat";
  }
  else if ( (mpde_ic_print_parameters.format_ == Format::PROBE) ||
            (mpde_ic_print_parameters.format_ == Format::RAW) ||
            (mpde_ic_print_parameters.format_ == Format::RAW_ASCII) ||
            (mpde_ic_print_parameters.format_ == Format::TS1) ||
            (mpde_ic_print_parameters.format_ == Format::TS2) )
  {
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";
    // print out the Index column, since this will be in STD format
    mpde_ic_print_parameters.printIndexColumn_ = true;
  }
  else
  {
    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.unknown";
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::update_MPDE_STARTUP_printParamsForPrintFormat
// Purpose       : update some of the print_parameters for .PRINT MPDE_STARTUP
//                 based on FORMAT= on the .PRINT Line.  This is used for both
//                 fallbacks and non-fallbacks.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void OutputMgr::update_MPDE_STARTUP_printParamsForPrintFormat(PrintParameters &mpde_startup_print_parameters)
{
  if (mpde_startup_print_parameters.format_ == Format::STD)
  {
    mpde_startup_print_parameters.defaultExtension_ = ".startup.prn";
  }
  else if (mpde_startup_print_parameters.format_ == Format::CSV)
  {
    mpde_startup_print_parameters.defaultExtension_ = ".startup.csv";
  }
  else if (mpde_startup_print_parameters.format_ == Format::TECPLOT)
  {
    mpde_startup_print_parameters.defaultExtension_ = ".startup.dat";
  }
  else if ( (mpde_startup_print_parameters.format_ == Format::PROBE) ||
            (mpde_startup_print_parameters.format_ == Format::RAW) ||
            (mpde_startup_print_parameters.format_ == Format::RAW_ASCII) ||
            (mpde_startup_print_parameters.format_ == Format::TS1) ||
            (mpde_startup_print_parameters.format_ == Format::TS2) )
  {
    mpde_startup_print_parameters.defaultExtension_ = ".startup.prn";
    // print out the Index column, since this will be in STD format
    mpde_startup_print_parameters.printIndexColumn_ = true;
  }
  else
  {
    mpde_startup_print_parameters.defaultExtension_ = ".startup.unknown";
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSens
// Purpose       : registers set of variables to set for .SENS.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
bool OutputMgr::registerSens(const Util::OptionBlock &option_block)
{
  bool bsuccess = true;

  std::vector<std::string> functions;
  std::vector<std::string> parameters;

  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    std::string tag = (*it).uTag();
    if ( std::string( (*it).uTag() ,0,7) == "OBJFUNC") // this is a vector
    {
      functions.push_back((*it).stringValue());
    }
    else if ( std::string( (*it).uTag() ,0,7) == "OBJVARS") // this is a vector
    {
       // do nothing for now
    }
    else if (std::string((*it).uTag(), 0, 5) == "PARAM")
    {
      parameters.push_back((*it).stringValue());
    }
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  // ERK: important to set the index properly to make multiple objectives work correctly.
  //
  // For sensitivity outputs, a bunch of new variables are automatically added to the variable
  // list.  For now, they pretend to have the keyword "SENS", as this is an acceptable keyword
  // in the parser. (I tried using something like "TMPSENS" and the parser barfed.  In addition
  // to the keyword "SENS", they also have the names of the parameters included to later flesh out
  // the sensitivity names.
  //
  // Later, these lists of variables will be properly set up in the N_IO_OpBuilders.C file, which is
  // called from the OutputMgr::checkPrintParameters function.
  int index = 0;
  int index2 = 0;
  for (std::vector<std::string>::const_iterator it1 = functions.begin(), end1 = functions.end(); it1 != end1; ++it1,++index2)
  {
    const std::string &function = (*it1);

    {
      Util::Marshal mout;
      mout << function << std::string("OBJECTIVEFUNCTION") << Util::Op::identifier<SensitivityObjFunctionOp>() << index2;
      sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
    }

    for (std::vector<std::string>::const_iterator it2 = parameters.begin(), end2 = parameters.end(); it2 != end2; ++it2, ++index)
    {
      const std::string &parameter = (*it2);

      if (sensitivityOptions_ & SensitivityOptions::DIRECT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::SCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
      }

      if (sensitivityOptions_ & SensitivityOptions::ADJOINT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::SCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }

        adjointSensitivityFlag_ = true;
      }
    }
  }

  enableSensitivityFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSensOptions
// Purpose       : registers set of variables to set for .OPTIONS SENSITIVITY
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
bool OutputMgr::registerSensOptions(const Util::OptionBlock &option_block)
{
  sensitivityOptions_ = 0;

  bool adjointGiven=false;
  bool outputUnscaledGiven=false;

  sensitivityOptions_ =0;

  Util::ParamList::const_iterator it = option_block.begin(); 
  Util::ParamList::const_iterator end = option_block.end();
  for ( ; it != end; ++it)
  {
    if ((*it).uTag() == "ADJOINT")
    {
      adjointGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::ADJOINT;
      }
    }
    else if ((*it).uTag() == "DIRECT" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::DIRECT;
    }
    else if ((*it).uTag() == "OUTPUTSCALED" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::SCALED;
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::UNSCALED;
      }
    }
  }

  // default behavior is for the code to assume adjoint sensitivites are wanted
  // even if ADJOINT=1 wasn't specified.
  if (!adjointGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::ADJOINT;
  }

  // default behavior is for the code to assume unscaled sensitivities are wanted
  // for output.
  if (!outputUnscaledGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::UNSCALED;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerNoise
// Purpose       : registers set of variables to set for .NOISE.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/26/2015
//-----------------------------------------------------------------------------
bool OutputMgr::registerNoise (const Util::OptionBlock &option_block)
{
  bool bsuccess = true;

  std::vector<std::string> parameters;
  for (Util::ParamList::const_iterator it = option_block.begin(), 
      end = option_block.end(); it != end; ++it)
  { 
    // if pts_per_summary is set, then we need to output a lot more stuff.
    // NOTE.  pts_per_summary is a spice3-ism, and probably not that useful here.
    if ((*it).uTag() == "PTS_PER_SUMMARY")
    {
      pts_per_summary_ = (*it).getImmutableValue<int>();
      pts_per_summary_Given = true;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : makeParamList
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
template<class It, class Out>
void makeParamList(const std::string name, It first, It last, Out out)
{
  for (; first != last; ++first)
  {
    *out++ = Util::Param(name, 1);
    *out++ = Util::Param(*first, 0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::removeStarVariables
// Purpose       : Process V(*) and I(*) on .print line
// Special Notes : Replaces v(*) and i(*) that appears on the .print line
//                 with a list of all v() or i() variables as appropriate.
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void removeStarVariables(
  Parallel::Machine     comm,
  Util::ParamList &     variable_list,
  const NodeNameMap &   all_nodes,
  const NodeNameMap &   external_nodes)
{
  bool vStarFound = false;
  bool iStarFound = false;

  Util::ParamList::iterator vStarPosition = variable_list.end();
  Util::ParamList::iterator iStarPosition = variable_list.end();

  // remove the v(*) and i(*) from the .print line
  Util::ParamList::iterator it = variable_list.begin();
  for ( ; it != variable_list.end(); )
  {
    // process the * entries
    if ((*it).tag() == "*")
    {
      // move to the type
      --it;

      // remember type of replacement and location of *
      if ((*it).tag() == "V")
      {
        vStarFound = true;
        vStarPosition = it;
        --vStarPosition;
      }
      else if ((*it).tag() == "I")
      {
        iStarFound = true;
        iStarPosition = it;
        --iStarPosition;
      }

      // remove the v or i
      it = variable_list.erase(it);

      // remove the *
      it = variable_list.erase(it);
    }

    // move to next item on .print line
    else
    {
      ++it;
    }
  }

  std::vector<std::string> v_list;
  std::vector<std::string> i_list;

  if (vStarFound)
  {
    NodeNameMap::const_iterator it = external_nodes.begin(); 
    for ( ; it != external_nodes.end() ; ++it)
    {
      ExtendedString tmpStr((*it).first);
      tmpStr.toUpper();

      if (tmpStr.rfind("BRANCH") == std::string::npos)
        v_list.push_back(tmpStr);
    }
    ++vStarPosition;
  }

  if (iStarFound)
  {
    NodeNameMap::const_iterator iter_a = all_nodes.begin();
    for ( ; iter_a != all_nodes.end() ; ++iter_a)
    {
      ExtendedString tmpStr((*iter_a).first);
      tmpStr.toUpper();
      char devType=tmpStr[0];

      size_t pos = tmpStr.rfind("BRANCH");
      if (pos != std::string::npos)
      {
        bool addIt=true;
        tmpStr = tmpStr.substr(0, pos - 1);
        // BUG 982:  The names of devices in the "all_nodes" map
        // are all in SPICE style because of a poor design decision.
        // Other code in the IO package assumes that the params we're passing
        // in print params lists are all in Xyce style.  So we must
        // fake things out here by converting back to Xyce style.
        tmpStr = Util::spiceDeviceNameToXyceName(tmpStr);
        // BUG 989 SON:  we are allowing Y device branches to be output by
        // I(*) now, but YMIL and YMIN (mutual inductors) have two names
        // for each branch current, and one of them is already included without
        // the unintuitive "y" device syntax.  So don't print the second
        // one.  Only add tmpStr if it is neither a YMIL or YMIN.
        if (devType == 'Y')
        {
          std::string::size_type i=tmpStr.find_last_of(':');
          i=((i == std::string::npos)?0:i+1);
          std::string basename=tmpStr.substr(i);
          if (startswith_nocase(basename,"YMIL")
              || startswith_nocase(basename,"YMIN"))
          {
            addIt=false;
          }
        }
        if (addIt)
        {
          i_list.push_back(tmpStr);
        }
      }
    }
    ++iStarPosition;
  }

  Util::Marshal mout;
  mout << v_list << i_list;

  std::vector<std::string> dest;
  Parallel::AllGatherV(comm, mout.str(), dest);

  Util::ParamList vStarList;
  Util::ParamList iStarList;

  for (int p = 0; p < Parallel::size(comm); ++p)
  {
    Util::Marshal min(dest[p]);

    std::vector<std::string> x;
    std::vector<std::string> y;
    min >> x >> y;
    makeParamList("V", x.begin(), x.end(), std::back_inserter(vStarList));
    makeParamList("I", y.begin(), y.end(), std::back_inserter(iStarList));
  }

  // append temporary lists to print block, erasing temporary lists
  variable_list.splice(vStarPosition, vStarList);
  variable_list.splice(iStarPosition, iStarList);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setSweepParameters
// Purpose       : Copy .DC or .STEP sweep parameters, and set up flags if
//                 sweeping the TEMP variable.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void OutputMgr::setStepSweepVector(const Analysis::SweepVector &sweep_vector)
{
  if (!sweep_vector.empty())
  {
    outputState_.stepSweepVector_ =  sweep_vector;
  }

  // check if one of the sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!sweep_vector.empty())
  {
    for (Analysis::SweepVector::const_iterator iterP = outputState_.stepSweepVector_.begin(), 
        end = outputState_.stepSweepVector_.end(); iterP != end; ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        outputState_.tempSweepFlag_ = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setSweepParameters
// Purpose       : Copy .DC or .STEP sweep parameters, and set up flags if
//                 sweeping the TEMP variable.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void OutputMgr::setDCSweepVector(const Analysis::SweepVector &sweep_vector)
{
  if (!sweep_vector.empty())
  {
    outputState_.dcSweepVector_ = sweep_vector;
  }

  // check if one of the DC sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!sweep_vector.empty() && !outputCalledBefore_)
  {
    Analysis::SweepVector::const_iterator iterP;
    for (iterP = outputState_.dcSweepVector_.begin();
         iterP != outputState_.dcSweepVector_.end(); ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        outputState_.tempSweepFlag_ = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::fixupPrintParameters
// Purpose       : Perform some .print line checks and munging, primarily
//                 dealing with use of V(*) and I(*)
// Special Notes : This function used to call removeStarVariables directly,
//                 but that was moved to "fixupOutputVariables."  This function
//                 now exists only to maintain a prior interface for
//                 existing outputters to call.
// Scope         : public
// Creator       : David Baur, Raytheon   (modded T. Russo, 21 Feb 2018)
// Creation Date : 6/26/2013 
//-----------------------------------------------------------------------------
void OutputMgr::fixupPrintParameters(
    Parallel::Machine comm, 
    PrintParameters &print_parameters)
{
  fixupOutputVariables(comm, print_parameters.variableList_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::fixupOutputVariables
// Purpose       : 
// Special Notes : Just a re-wrapping of removeStarVariables
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 21 Feb 2018
//-----------------------------------------------------------------------------
void OutputMgr::fixupOutputVariables(
    Parallel::Machine comm, 
    Util::ParamList &outputParamList)
{
  removeStarVariables(comm, outputParamList, 
      getSolutionNodeMap(), getExternalNodeMap());
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::output
// Purpose       : Runs specified output commands
//
// Special Notes : ERK.  I've refactored this so that it receives STEP
//                 and DC sweep parameter information(mainly name and value).
//                 This function is called from the time integrator, which has
//                 all this info.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
/// Perform all output requested by transient or DC .PRINT lines
///
/// Basically, this function loops over all outputters in the
/// activeOutputterStack_ and calls their output functions.
///
/// @param[in] skipPrintLineOutput boolean to indicate that no outputter should actually be called
///
/// @note The skipPrintLineOutput bool is something of a kludge that is used
/// when the user has set specified output intervals.  The transient analysis
/// code sets this to "true" if we've computed a solution but have not yet
/// reached the time of the next requested output.  The point of that is
/// apparently that more gets done by this output function than simply
/// doing output, and we want that extra stuff to be done even though the
/// actual output won't happen.  Perhaps this function needs not to be
/// doing more than actually outputting?
///
/// @note This function is called from the analysis manager for transient and
///  DC analyses.  Other functions are called for other analyses such as AC,
///  noise, homotopy, sensitivity, etc.
void OutputMgr::output(
  Parallel::Machine             comm,
  const double                  time,
  const double                  timeStep,
  const double                  circuit_temp, 
  const int                     stepNumber,
  const int                     maxStep,
  const Analysis::SweepVector & step_sweep_vector,
  const int                     dcNumber,
  const int                     maxDC,
  const Analysis::SweepVector & dc_sweep_vector,
  const Linear::Vector &        solnVecPtr,
  const Linear::Vector &        stateVecPtr,
  const Linear::Vector &        storeVecPtr,
  const Linear::Vector &        lead_current_vector,
  const Linear::Vector &        junction_voltage_vector,
  const Linear::Vector &        lead_current_dqdt_vector,
  const std::vector<double> &   objectiveVec,
  const std::vector<double> &   dOdpVec,
  const std::vector<double> &   dOdpAdjVec,
  const std::vector<double> &   scaled_dOdpVec,
  const std::vector<double> &   scaled_dOdpAdjVec,
  bool                          skipPrintLineOutput)
{
  // copy over time:
  outputState_.circuitTime_ = time;
  outputState_.circuitTimeStep_ = timeStep;

  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;

  // copy the new values into the locally owned vector:
  if (!step_sweep_vector.empty())
  {
    outputState_.stepSweepVector_ = step_sweep_vector;
  }

  // copy over the dc sweep information:
  dcLoopNumber_ = dcNumber;
  maxDCSteps_ = maxDC;

  if (!dc_sweep_vector.empty())
  {
    outputState_.dcSweepVector_ = dc_sweep_vector;

    // For now just have PRINTdcvalue, etc just be the first parameter.
    const Analysis::SweepParam &firstParam = dc_sweep_vector.front();
    PRINTdcname_   = firstParam.name;
    PRINTdcvalue_  = firstParam.currentVal;
    if (firstParam.type == "LIST")
    {
      PRINTdcstart_  = firstParam.valList.front();
      PRINTdcstop_   = firstParam.valList.back();
    }
    else
    {
      PRINTdcstart_  = firstParam.startVal;
      PRINTdcstop_   = firstParam.stopVal;
    }
  }

  // Check for temperature:
  outputState_.circuitTemp_ = circuit_temp;

  // Needs to pass skipPrintLineOutput
  if (!skipPrintLineOutput)
  {
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator
          it = activeOutputterStack_.back().begin();
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->output(comm, solnVecPtr, stateVecPtr, storeVecPtr, 
                      lead_current_vector, junction_voltage_vector);

        (*it)->outputSensitivity(comm, objectiveVec,
                                 dOdpVec, dOdpAdjVec, 
                                 scaled_dOdpVec, scaled_dOdpAdjVec,
                                 solnVecPtr, stateVecPtr, storeVecPtr);
      }
    }
  }

  // This will become an outputter
  if (hdf5FileNameGiven_)
  {
    updateHDF5Output(comm, solnVecPtr);
  }
  
  outputCalledBefore_ = true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputSensitivity
// Purpose       : outputs sensitivity info  
//
// Special Notes : This was set up specifically to support transient adjoint 
//                 sensitivities.  Specifically for the case of using them to 
//                 compute local sensitivities.  It is (mostly) a copy of the
//                 OutputMgr::output function, but with the various 
//                 non-sensitivity stuff stripped out.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/13/2016
//-----------------------------------------------------------------------------
void OutputMgr::outputSensitivity(
  Parallel::Machine             comm,
  const double                  time,
  const double                  timeStep,
  const double                  circuit_temp, 
  const int                     stepNumber,
  const int                     maxStep,
  const Analysis::SweepVector & step_sweep_vector,
  const int                     dcNumber,
  const int                     maxDC,
  const Analysis::SweepVector & dc_sweep_vector,
  const Linear::Vector &        solnVecPtr,
  const Linear::Vector &        stateVecPtr,
  const Linear::Vector &        storeVecPtr,
  const Linear::Vector &        lead_current_vector,
  const Linear::Vector &        junction_voltage_vector,
  const Linear::Vector &        lead_current_dqdt_vector,
  const std::vector<double> &   objectiveVec,
  const std::vector<double> &   dOdpVec,
  const std::vector<double> &   dOdpAdjVec,
  const std::vector<double> &   scaled_dOdpVec,
  const std::vector<double> &   scaled_dOdpAdjVec,
  bool                          skipPrintLineOutput)
{
  // copy over time:
  outputState_.circuitTime_ = time;
  outputState_.circuitTimeStep_ = timeStep;

  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;

  // copy the new values into the locally owned vector:
  if (!step_sweep_vector.empty())
  {
    outputState_.stepSweepVector_ = step_sweep_vector;
  }

  // Check for temperature:
  outputState_.circuitTemp_ = circuit_temp;

  // Needs to pass skipPrintLineOutput
  if (!skipPrintLineOutput)
  {
    OutputterMap::iterator find_it = outputterMap_.find( PrintType::TRANADJOINT );

    if (find_it != outputterMap_.end())
    {
      std::vector<Outputter::Interface *>::const_iterator it = (*find_it).second.begin();
      std::vector<Outputter::Interface *>::const_iterator end = (*find_it).second.end();

      for ( ; it != end; ++it)
      {
        (*it)->outputSensitivity(comm, objectiveVec,
                                 dOdpVec, dOdpAdjVec, 
                                 scaled_dOdpVec, scaled_dOdpAdjVec,
                                 solnVecPtr, stateVecPtr, storeVecPtr);
      }
    }
    else
    {
      Report::UserWarning0() << "Cannot find any transient adjoint outputters!";
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::finishSensitivityOutput
// Purpose       : Runs specified finish output commands for adjoint 
//                 sensitivities.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 5/3/2016
//-----------------------------------------------------------------------------
void OutputMgr::finishSensitivityOutput()
{
  OutputterMap::iterator find_it = outputterMap_.find( PrintType::TRANADJOINT );
  if (find_it != outputterMap_.end())
  {
    std::vector<Outputter::Interface *>::const_iterator it = (*find_it).second.begin();
    std::vector<Outputter::Interface *>::const_iterator end = (*find_it).second.end();

    for ( ; it != end; ++it)
    {
      (*it)->finishOutput();
    }
  }
  else
  {
    //Report::UserWarning0() << "Cannot find any transient adjoint outputters!";
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputAC(
  Parallel::Machine     comm,
  double                frequency,
  double                fStart,
  double                fStop,
  const Linear::Vector &  real_solution_vector,
  const Linear::Vector &  imaginary_solution_vector,
  const Util::Op::RFparamsData & RFparams)
{
  outputState_.circuitFrequency_ = frequency;

  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin(); 

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputAC(comm, frequency, fStart, fStop, 
                      real_solution_vector, imaginary_solution_vector, RFparams);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputSensitivityAC
// Purpose       : .PRINT SENS output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 4/15/19
//-----------------------------------------------------------------------------
void OutputMgr::outputSensitivityAC(
  Parallel::Machine                   comm,
  double                              frequency,
  const Linear::Vector &              real_solution_vector,
  const Linear::Vector &              imaginary_solution_vector,
  const std::vector<double> &         paramVals,
  const std::vector<std::string> &    paramNameVec,
  const std::vector<std::string> &    objFuncVars,
  const std::vector<double> &         objectiveVec,
  const std::vector<double> &         dOdpVec,
  const std::vector<double> &         dOdpAdjVec,
  const std::vector<double> &         scaled_dOdpVec,
  const std::vector<double> &         scaled_dOdpAdjVec)
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it =
      activeOutputterStack_.back().begin();

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputSensitivityAC(comm, frequency, real_solution_vector, imaginary_solution_vector,
		      paramVals, paramNameVec, objFuncVars, objectiveVec,
                      dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputSParams
// Purpose       : .LIN (S-parameter) output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 03/25/2019
//-----------------------------------------------------------------------------
void OutputMgr::outputSParams(
  Parallel::Machine     comm,
  double                frequency,
  double                numFreq,
  std::vector<double> & Z0sVec,
  const Util::Op::RFparamsData & RFparams)
{
  outputState_.circuitFrequency_ = frequency;

  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it =
      activeOutputterStack_.back().begin();

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputSParams(comm, frequency, numFreq, Z0sVec, RFparams);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputNoise
// Purpose       : .PRINT output for noise runs
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputNoise(
    Parallel::Machine     comm,
    double                frequency,
    const Linear::Vector & real_solution_vector, 
    const Linear::Vector & imaginary_solution_vector,
    double                totalOutputNoiseDens_, 
    double                totalInputNoiseDens_, 
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_)
{
  outputState_.circuitFrequency_ = frequency;

  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin();

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputNoise(comm, frequency, real_solution_vector, imaginary_solution_vector,
               totalOutputNoiseDens_, totalInputNoiseDens_, noiseDataVec_);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputEmbeddedSampling
// Purpose       : .PRINT ES for embedded sampling runs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void OutputMgr::outputEmbeddedSampling(
    Parallel::Machine comm,
    bool regressionPCEenable,
    bool projectionPCEenable,
    int  numSamples,
    const std::vector<std::string> & regressionPCEcoeffs,
    const std::vector<std::string> & projectionPCEcoeffs,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec_)
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it =
      activeOutputterStack_.back().begin();

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputEmbeddedSampling(comm, regressionPCEenable, projectionPCEenable,
	       numSamples, regressionPCEcoeffs, projectionPCEcoeffs, outFuncDataVec_);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPCE
// Purpose       : .PRINT PCE for pce runs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void OutputMgr::outputPCE(
    Parallel::Machine comm,
    int numQuadPoints,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec_)
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it =
      activeOutputterStack_.back().begin();

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputPCE(comm, numQuadPoints, outFuncDataVec_);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHomotopy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const Linear::Vector &                  solution_vector)
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator
        it = activeOutputterStack_.back().begin();
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHomotopy(comm, parameter_names, param_values, solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::finishOutput
// Purpose       : Runs specified finish output commands
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 12/13/00
//-----------------------------------------------------------------------------
void OutputMgr::finishOutput()
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator
        it = activeOutputterStack_.back().begin();
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->finishOutput();
    }
  }

  if (hdf5FileNameGiven_)
  {
    closeHDF5Output();
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHB_TD
// Purpose       : Print time-domain output for Harmonic Balance runs
// Special Notes : This function is used by .PRINT HB_TD lines.  It is not
//                 used for .PRINT HB_STARTUP or .PRINT HB_IC lines though.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/5/2013
//-----------------------------------------------------------------------------
void OutputMgr::outputHB_TD(
  Parallel::Machine                     comm,
  const int                             stepNumber,
  const int                             maxStep,
  const Analysis::SweepVector &         step_sweep_parameters,
  const std::vector<double> &           timePoints,
  const Linear::BlockVector &           timeDomainSolutionVec,
  const Linear::BlockVector &           timeDomainLeadCurrentVec,
  const Linear::BlockVector &           timeDomainJunctionVoltageVec)
{
  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;

  // copy the new values into the locally owned vector:
  if (!(step_sweep_parameters.empty()))
  {
    outputState_.stepSweepVector_ = step_sweep_parameters;
  }

  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin();

    for (; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHB_TD(comm, timePoints,
        timeDomainSolutionVec, 
        timeDomainLeadCurrentVec, 
        timeDomainJunctionVoltageVec);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHB_FD
// Purpose       : Print frequency domain output for Harmonic Balance runs.
// Special Notes : This is used by .PRINT HB_FD lines.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/5/2013
//-----------------------------------------------------------------------------
void OutputMgr::outputHB_FD(
  Parallel::Machine                     comm,
  const int                             stepNumber,
  const int                             maxStep,
  const Analysis::SweepVector &         step_sweep_parameters,
  const std::vector<double> &           freqPoints,
  const Linear::BlockVector &           freqDomainSolutionVecReal,
  const Linear::BlockVector &           freqDomainSolutionVecImaginary,
  const Linear::BlockVector &           freqDomainLeadCurrentVecReal,
  const Linear::BlockVector &           freqDomainLeadCurrentVecImaginary,
  const Linear::BlockVector &           freqDomainJunctionVoltageVecReal,
  const Linear::BlockVector &           freqDomainJunctionVoltageVecImaginary )
{
  // copy over the step sweep information:
  outputState_.stepLoopNumber_ = stepNumber;
  outputState_.stepMaxCount_ = maxStep;

  // copy the new values into the locally owned vector:
  if (!(step_sweep_parameters.empty()))
  {
    outputState_.stepSweepVector_ = step_sweep_parameters;
  }

  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin();

    for (; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHB_FD(comm, freqPoints,
        freqDomainSolutionVecReal, freqDomainSolutionVecImaginary,
        freqDomainLeadCurrentVecReal, freqDomainLeadCurrentVecImaginary,
        freqDomainJunctionVoltageVecReal, freqDomainJunctionVoltageVecImaginary);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputMPDE
// Purpose       : .PRINT output for mpde runs
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/22/03
//-----------------------------------------------------------------------------
void OutputMgr::outputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points,
  const Linear::BlockVector &   solution_vector)
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
        activeOutputterStack_.back().begin(); 

    for ( ; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputMPDE(comm, time, fast_time_points, solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::startStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::startStep(int step, int max_step)
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin();

    for (;it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->startStep(step, max_step);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::resetIndex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::resetIndex()
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it = 
      activeOutputterStack_.back().begin();

    for (; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->resetIndex();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::steppingComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::steppingComplete()
{
  if (!activeOutputterStack_.empty())
  {
    std::vector<Outputter::Interface *>::const_iterator it =
      activeOutputterStack_.back().begin();

    for (; it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->steppingComplete();
    }
  }
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : populateMetadata
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date : May 12, 2015
//-----------------------------------------------------------------------------
void populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("OUTPUT");

    parameters.insert(Util::ParamMap::value_type("INITIAL_INTERVAL", Util::Param("INITIAL_INTERVAL", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TIME", Util::Param("TIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("INTERVAL", Util::Param("INTERVAL", 0.0)));
    parameters.insert(Util::ParamMap::value_type("HDF5FILENAME", Util::Param("HDF5FILENAME", "")));
    parameters.insert(Util::ParamMap::value_type("PRINTHEADER", Util::Param("PRINTHEADER", true)));
    parameters.insert(Util::ParamMap::value_type("PRINTFOOTER", Util::Param("PRINTFOOTER", true)));
    parameters.insert(Util::ParamMap::value_type("ADD_STEPNUM_COL", Util::Param("ADD_STEPNUM_COL", true)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTVERSIONINRAWFILE", Util::Param("OUTPUTVERSIONINRAWFILE", false)));
    parameters.insert(Util::ParamMap::value_type("PHASE_OUTPUT_RADIANS", Util::Param("PHASE_OUTPUT_RADIANS", true)));

    parameters.insert(Util::ParamMap::value_type("OUTPUTTIMEPOINTS", Util::Param("OUTPUTTIMEPOINTS", "VECTOR")));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("OUTPUT-LINE");

    parameters.insert(Util::ParamMap::value_type("TIME", Util::Param("TIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("INTERVAL", Util::Param("INTERVAL", 0.0)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("PRINT");

    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "TRAN")));
    parameters.insert(Util::ParamMap::value_type("FILE", Util::Param("FILE", "")));
    parameters.insert(Util::ParamMap::value_type("FORMAT", Util::Param("FORMAT", "STD")));
    parameters.insert(Util::ParamMap::value_type("DATAFORMAT", Util::Param("DATAFORMAT", "RI")));
    parameters.insert(Util::ParamMap::value_type("LINTYPE", Util::Param("LINTYPE", "S")));
    parameters.insert(Util::ParamMap::value_type("DELIMITER", Util::Param("DELIMITER", "")));
    parameters.insert(Util::ParamMap::value_type("WIDTH", Util::Param("WIDTH", 17)));
    parameters.insert(Util::ParamMap::value_type("PRECISION", Util::Param("PRECISION", 8)));
    parameters.insert(Util::ParamMap::value_type("TIMESCALEFACTOR", Util::Param("TIMESCALEFACTOR", 1.0)));
    parameters.insert(Util::ParamMap::value_type("FILTER", Util::Param("FILTER", 0.0)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTSAMPLESTATS", Util::Param("OUTPUTSAMPLESTATS", true)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTALLSAMPLES", Util::Param("OUTPUTALLSAMPLES", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_PCE_COEFFS", Util::Param("OUTPUT_PCE_COEFFS", false)));
  }
}

//-----------------------------------------------------------------------------
// Function      : extractPrintData
// Purpose       : Extract the parameters from a netlist .PRINT line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 05/02/2002
// Completely revamped by David Baur 05/12/2015
//-----------------------------------------------------------------------------
///
/// Given a .PRINT line from the netlist, create an OptionBlock named "PRINT" for it.
///
/// @param options_manager  Reference to the options manager object.
/// @param circuit_block    Circuit block to which we will add the new options_block
/// @param netlist_filename  Name of the file that contains this .PRINT line
/// @param parsed_line       Vector of strings representing the netlist line, not so much parsed as split at spaces
///
/// The parsed_line string vector contains a .PRINT line of the form
///   .PRINT <type> [<options>] [<Variables to print>+]
/// after having had each field of the line split into separate strings.
///
/// This function is registered with the options manager through an
/// "addCommandParser" call, and is called by the options manager
/// whenever a .PRINT line is encountered.
///
/// After the option block is fully populated, it is added to the
/// circuit block.  "PRINT" option blocks are later processed by
/// parsePRINTBlock, which itself is indirectly called via a registration
/// with the options manager function addCommandProcessor.
///
/// @author David Baur
/// @date 05/12/2015
bool extractPrintData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("PRINT", 
      Util::OptionBlock::ALLOW_EXPRESSIONS, 
      netlist_filename, 
      parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  if (DEBUG_IO) 
  {
    for (int ieric=0;ieric<parsed_line.size();++ieric)
    {
      Xyce::dout() << "parsed_line["<<ieric<<"] = " 
        << parsed_line[ieric].string_ << std::endl;
    }
  }

  // Set the TYPE and add it the parameters.
  std::string tmpString = parsed_line[1].string_;
  Util::toUpper(tmpString);
  if (tmpString == "TR")
  {
    tmpString = "TRAN"; // TR is a synonym for TRAN
  }
  Util::Param typeParameter("TYPE", tmpString);

  // Add the options parameter set to the model.
  addDefaultOptionsParameters(options_manager, option_block, "PRINT");

  // Reset the default TYPE with the value found.
  Util::Param *parameterPtr = Util::findParameter(option_block.begin(), option_block.end(), typeParameter.tag());
  if( parameterPtr == NULL )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Failed to find parameter\"" << typeParameter.usVal() << "\" in optionData for .PRINT statement";
  }

  parameterPtr->setVal( typeParameter.usVal() );

  // Set no default value for the output file.
  // if the user set FILE=vale it will be picked up here.  Otherwise
  // the output manager will use the top level simulation netlist name
  // as the default.
  Util::Param fileParameter("FILE", "");
  parameterPtr = Util::findParameter(option_block.begin(), option_block.end(), fileParameter.tag());

  if( parameterPtr == NULL )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Failed to find parameter\"" << typeParameter.usVal() << "\" in optionData for .PRINT statement";
  }

  parameterPtr->setVal( fileParameter.stringValue() );

  //
  // Note: The PRINT control parameters (FORMAT, WIDTH, FILE, ...) must be
  // tagged on the print line. Tagged parameters are "better" since they can
  // be identifed and parsed based soley on metadata with generic code.
  // The code that used to allow untagged parameters has been removed.

  // Check for tagged parameters.
  int parameterStartPos  = 2;
  int position = parameterStartPos;
  if ( numFields > parameterStartPos + 1 && parsed_line[parameterStartPos + 1].string_ == "=" )
  {
     // Tagged parameters found.
     Util::Param* parameterPtr;
     while ( position+1 < (int)parsed_line.size() &&
             parsed_line[position+1].string_ == "=" )
     {
       parameterPtr = Util::findParameter(option_block.begin(), option_block.end(), parsed_line[position].string_);
       if ( parameterPtr != NULL )
       {
         if (parameterPtr->tag() != "FILE")
           parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_ ).toUpper()));
         else
           parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_)));
       }
       else
       {
         Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
           << "No PRINT parameter " << parsed_line[position].string_ << " found, parameter will be ignored.";
         return false;
       }

       position += 3;
     }
  }
  else if (DEBUG_IO) // There is no point to this warning unless we're debugging
  {
    Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "No tagged parameters found";
  }

  // Complete the PRINT line parsing.
  // Some of the remaining fields are of the
  // form I(Vname), V(node) or V(node1,node2), and these fields need
  // special treatment.
  ExtendedString field("");
  std::ostringstream msg;
  int p_err=0;
  while ( position < numFields )
  {
    if (position+1 < numFields && parsed_line[position+1].string_ == "(")
    {
      // Note: the next if statement has to support I, IR, II, IM, IP and IDB.  That is why it
      // checks for strings of size 2 and 3.  A similar comment applies to the block for
      // V below.
      if (toupper(parsed_line[position].string_[0]) == 'I' && parsed_line[position].string_.size() <= 3)
      {
        if ((position+3 < numFields && parsed_line[position+3].string_ == ")") ||
            (position+4 < numFields && parsed_line[position+4].string_ == ")"))
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,1.0));

          if (parsed_line[position+3].string_ == ")")
          {
            field = parsed_line[position+2].string_;
            field.toUpper();
            option_block.addParam( Util::Param(field, 0.0) );

            position += 4;
          }
          else
          {
            // Note:  This block is here to handle the case where a user has
            // asked for I(YSOMETHING NAME) instead of I(YSOMETHING!NAME)
            field = parsed_line[position+2].string_ + " " + parsed_line[position+3].string_;
            field.toUpper();
            option_block.addParam( Util::Param(field,0.0) );

            position += 5;
          }
        }
        else
        {
          msg << "Unrecognized current specification";
          p_err = position;
        }
      }
      else if ( toupper(parsed_line[position].string_[0]) == 'V' && parsed_line[position].string_.size() <= 3 )
      {
        // position+3 < numFields test required to prevent a core dump on 
        // something invalid like .PRINT TRAN V(1
        if (position+3 < numFields && parsed_line[position+3].string_ == ")")
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,1.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 4;
        }
        else if (position+5 < numFields && parsed_line[position+5].string_ == ")")
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,2.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          field = parsed_line[position+4].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 6;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else if ( toupper(parsed_line[position].string_[0]) == 'N'
                && parsed_line[position].string_.size() == 1 )
      {
        if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
        {
          option_block.addParam( Util::Param("N",1.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 4;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else if( (parsed_line[position].string_.size() == 3) &&
             (parsed_line[position].string_[0] == 'D' || parsed_line[position].string_[0] == 'd') )
      {
        if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,1.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 4;
        }
        else if (position+5 < numFields && parsed_line[position+5].string_ == ")")
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,2.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          field = parsed_line[position+4].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 6;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else if( parsed_line[position].string_ == "W" || parsed_line[position].string_ == "w" ||
               parsed_line[position].string_ == "P" || parsed_line[position].string_ == "p")
      {
        if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,1.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 4;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else if ( ((toupper(parsed_line[position].string_[0]) == 'S') ||
                 (toupper(parsed_line[position].string_[0]) == 'Y') ||
                 (toupper(parsed_line[position].string_[0]) == 'Z'))
                && parsed_line[position].string_.size() <= 3 )
      {
        if (position+5 < numFields && parsed_line[position+5].string_ == ")")
        {
          field = parsed_line[position].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,2.0) );

          field = parsed_line[position+2].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          field = parsed_line[position+4].string_;
          field.toUpper();
          option_block.addParam( Util::Param(field,0.0) );

          position += 6;
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else
      {
        msg << "Unrecognized parenthetical specification";
        p_err = position;
      }
    }
    else // This is either an expression, or a STEP parameter.
    {
      if (parsed_line[position].string_ == "(" || parsed_line[position].string_ == ")")
      {
        msg << "Unrecognized parenthesis";
        p_err = position;
      }
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam ( Util::Param(field,0.0) );
      ++position;
    }

    if (!msg.str().empty())
    {
      msg << " in .print near ";
      position = p_err+4;
      p_err -= 2;
      if (p_err<0)
        p_err = 0;
      if (position >= numFields)
        position = numFields;
      while (p_err < position)
      {
        msg << parsed_line[p_err].string_ + " ";
        ++p_err;
      }
      Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
        << msg.str();
    }
  }

  circuit_block.addOptions(option_block);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : extractLINData
// Purpose       : Extract the parameters from a netlist .LIN line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2019
//-----------------------------------------------------------------------------
///
/// Given a .LIN line from the netlist, create OptionBlocks named "PRINT"
/// and "ACLIN" for it.
bool extractLINData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  int numFields = parsed_line.size();

  Util::OptionBlock print_option_block("PRINT", 
      Util::OptionBlock::ALLOW_EXPRESSIONS, 
      netlist_filename, 
      parsed_line[0].lineNumber_);

  Util::OptionBlock option_block_aclin("ACLIN", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

 // add SPARAMS as the print type from for a .LIN line
  Util::Param typeParameter("TYPE", "SPARAM");

  // Add the options parameter set to the model.
  addDefaultOptionsParameters(options_manager, print_option_block, "PRINT");

  // Reset the default TYPE with the value found.
  addDefaultOptionsParameters(options_manager, option_block_aclin, "ACLIN");

  Util::Param *parameterPtr = Util::findParameter(print_option_block.begin(), print_option_block.end(), typeParameter.tag());
  if( parameterPtr == NULL )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Failed to find parameter\"" << typeParameter.usVal() << "\" in optionData for .PRINT statement";
  }

  parameterPtr->setVal( typeParameter.usVal() );

  // used to default to FORMAT=TOUCHSTONE2, for a .LIN line without any parameters on it
  bool foundFormatParam = false;

  // Check for tagged parameters.   In particular, the WIDTH, PRECISION and FILE/FILENAME parameters will
  // be allowed.
  int parameterStartPos  = 1;
  int position = parameterStartPos;
  if ( numFields > parameterStartPos + 1 && parsed_line[parameterStartPos + 1].string_ == "=" )
  {
     // Tagged parameters found.
     Util::Param* parameterPtr;
     while ( position+1 < (int)parsed_line.size() &&
             parsed_line[position+1].string_ == "=" )
     {
       // FILENAME is an allowed synonym for FILE on a .LIN line, for HSPICE compatibility
       ExtendedString paramName = parsed_line[position].string_;
       paramName.toUpper();

       if (paramName == "FILENAME") { paramName = "FILE";}
       if (paramName == "FORMAT") { foundFormatParam = true;}

       if (paramName  == "SPARCALC" )
       {
         // this parameter only goes into the aclin option blocks
         parameterPtr = Util::findParameter(option_block_aclin.begin(), option_block_aclin.end(), paramName);
         parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_ ).toUpper()));
       }
       else if (paramName == "LINTYPE")
       {
         // this parameter goes into both the aclin and print option blocks
         parameterPtr = Util::findParameter(option_block_aclin.begin(), option_block_aclin.end(), paramName);
         parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_ ).toUpper()));
         parameterPtr = Util::findParameter(print_option_block.begin(), print_option_block.end(), paramName);
         parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_ ).toUpper()));
       }
       else
       {
         // all of the other parameters only go into the print option block
         parameterPtr = Util::findParameter(print_option_block.begin(), print_option_block.end(), paramName);

         if ( parameterPtr != NULL )
         {
           if (parameterPtr->tag() == "DELIMITER")
	   {
             // DELIMITER parameter is not supported for the Touchstone formats.
             Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
               << "DELIMITER parameter not supported on .LIN line";
           }
           else if (parameterPtr->tag() != "FILE")
             parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_ ).toUpper()));
           else
             parameterPtr->setVal(std::string(ExtendedString(parsed_line[position+2].string_)));
         }
         else
         {
           Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
             << "No PRINT parameter " << parsed_line[position].string_ << " found, parameter will be ignored.";
         }
       }

       position += 3;
     }
  }
  else if (DEBUG_IO) // There is no point to this warning unless we're debugging
  {
    Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "No tagged parameters found";
  }

  // Allow for a .LIN line with no parameters on it.  In that case, we want to change the default
  // print format to TOUCHSTONE2 at this point rather than in parsePRINTBlock().  This prevents
  // a useless warning message from being generated from parsePRINTBlock().
  if ( !foundFormatParam )
  {
    parameterPtr = Util::findParameter(print_option_block.begin(), print_option_block.end(), "FORMAT");
    parameterPtr->setVal(std::string("TOUCHSTONE2"));
  }

  circuit_block.addOptions(option_block_aclin);

  circuit_block.addOptions(print_option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       : Associate specific parsing functions with each .print type
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 5/12/2015
//-----------------------------------------------------------------------------
///
/// Register various OutputManager functions with the options manager.
///
/// @param output_manager Reference to an OutputMgr object
/// @param options_manager Reference to a PkgOptionsMgr object
///
/// Registers a parser for ".PRINT" lines.
/// Registers options block processors for several other options block
/// types, not all of them associated with output.
///
/// @note Many of the functions registered here as command processors (or
///  rather, options block processors) have "command parsers" registered
///  in completely different files.  Command parsers take input netlist
///  lines and create options blocks, command processors take options blocks
///  and do something with them.  The organization of the code is somewhat
///  haphazard, and for some reason (probably an abandoned refactor)
///  the command processors and command parsers are not consistently in the
///  same files.  For example, the "registerSens" function is registered here
///  to process "SENS" options blocks, but the associated command parser
///  to create "SENS" options blocks from .SENS lines is actually in the
///  nonlinear solver package's N_NLS_Sensitivity.C file.
///  It does make this code design rather hard to follow.
///  This may be a good reason to go through a reorganization of
///  code at a later date.
bool registerPkgOptionsMgr(OutputMgr & output_manager, PkgOptionsMgr &options_manager)
{
  populateMetadata(options_manager);

  options_manager.addCommandParser(".PRINT", extractPrintData);
  options_manager.addCommandParser(".LIN", extractLINData);

  options_manager.addCommandProcessor("PRINT", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::parsePRINTBlock));

  options_manager.addCommandProcessor("OP", 
    IO::createRegistrationOptions(output_manager, &OutputMgr::setOPAnalysisParams));

  // These various registrations associate command processors for
  // options blocks whose command *parsers* live in some other file.
  // The command parser for ".SENS" is "extractSENSData" in the
  // NonlinearSolverPKG/N_NLS_Sensitivity.C file.
  options_manager.addCommandProcessor("SENS", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::registerSens));
  
  // The command parser for ".NOISE" is "extractNOISEData" in the
  // AnalysisPKG/N_ANP_NOISE.C file.
  options_manager.addCommandProcessor("NOISE", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::registerNoise));

  // These are registrations of options processors for options blocks
  // of the form ".OPTIONS <pkg>".  These are all processed by a single
  // command parser (extractOptionsData) set up in
  // UtilityPKG/N_IO_OptionsBlock.C
  options_manager.addOptionsProcessor("NONLIN", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::registerNonlinearOptions));
  options_manager.addOptionsProcessor("OUTPUT", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::registerOutputOptions));
  options_manager.addOptionsProcessor("SENSITIVITY", 
      IO::createRegistrationOptions(output_manager, &OutputMgr::registerSensOptions));

  return true;
}


} // namespace IO
} // namespace Xyce
