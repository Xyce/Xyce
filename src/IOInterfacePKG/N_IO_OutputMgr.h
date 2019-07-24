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

//-----------------------------------------------------------------------------
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMgr_h
#define Xyce_N_IO_OutputMgr_h

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
using std::unordered_map;
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#else
#error neither unordered_map or tr1/unordered_map found
#endif
#include <iterator>
#include <list>
#include <set>
#include <string>
#include <vector>

#ifdef Xyce_USE_HDF5
#include <hdf5.h>
#endif

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_TOP_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

#include <N_IO_PrintTypes.h>
#include <N_IO_OutputTypes.h>

#include <N_ANP_SweepParam.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_Outputter.h>

#include <N_UTL_NoCase.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Op.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>

class N_MPDE_Manager;

namespace Xyce {
namespace IO {

typedef std::map<PrintType::PrintType, std::vector<Outputter::Interface *> > OutputterMap;
typedef std::vector<std::vector<Outputter::Interface *> > ActiveOutputterStack;
typedef std::set<Analysis::Mode> EnabledAnalysisSet;
typedef std::vector<std::pair<std::string, std::string> > StringPairVector;

enum Domain {DOMAIN_TIME, DOMAIN_FREQUENCY};

//-----------------------------------------------------------------------------
// Class         : OutputMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
class OutputMgr : public Util::Listener<Analysis::StepEvent>
{
public:
  struct OutputterKey
  {
    OutputterKey(int analysis_mode, int output_type, int format)
      : analysisMode_(analysis_mode),
        outputType_(output_type),
        format_(format)
    {}

    const int   analysisMode_;
    const int   outputType_;
    const int   format_;
  };

  typedef Outputter::Interface *(*OutputterFactory)(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);
  typedef std::map<OutputterKey, OutputterFactory> OutputterFactoryMap;
  typedef std::map<std::string, std::pair<int, std::ostream *> > OpenPathStreamMap;

  OutputMgr(const CmdParse &command_line, Util::Op::BuilderManager &op_builder_manager, const Topo::Topology &topology);
  ~OutputMgr();

private:
  OutputMgr(const OutputMgr & rhs);
  OutputMgr &operator=(const OutputMgr &rhs);

public:
  double getInitialOutputInterval() const
  {
    return initialOutputInterval_;
  }

  const IntervalVector &getOutputIntervals() const
  {
    return outputIntervalPairs_;
  }

  // helps set the dotOpSpecified_ flag if the netlist has a .OP statement
  bool setOPAnalysisParams(const Util::OptionBlock & paramsBlock);

  bool registerNonlinearOptions(const Util::OptionBlock & option_block);
  //bool registerTran(const Util::OptionBlock & option_block);
  //bool registerTranOptions(const Util::OptionBlock & option_block);
//  bool registerMPDETranOptions(const Util::OptionBlock & option_block);
//  bool registerHBOptions(const Util::OptionBlock & option_block);
  bool registerOutputOptions(const Util::OptionBlock & option_block);
//  bool registerDeviceOptions(const Util::OptionBlock & option_block);

  bool parsePRINTBlock(const Util::OptionBlock & print_block);
  bool registerSens(const Util::OptionBlock & option_block);
  bool registerSensOptions(const Util::OptionBlock & option_block);

  bool registerNoise(const Util::OptionBlock & option_block);

  void notify(const Analysis::StepEvent &step_event);

  void fixupPrintParameters(Parallel::Machine comm, PrintParameters &print_parameters);
  void fixupOutputVariables(Parallel::Machine comm, Util::ParamList &outputParamList);

public:
  void checkPrintParameters(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);
  void prepareOutput(Parallel::Machine comm, Analysis::Mode analysis_mode);
  void setStepSweepVector(const Analysis::SweepVector &sweep_vector);
  void setDCSweepVector(const Analysis::SweepVector &sweep_vector);

  const std::vector<char> &getVarTypes() const;

  const NodeNameMap &getSolutionNodeMap() const;

  const NodeNameMap &getStateNodeMap() const;

  const NodeNameMap &getStoreNodeMap() const;

  const NodeNameMap &getExternalNodeMap() const;

  const NodeNameMap &getBranchVarsNodeMap() const;

  // used to get the devices names (that have noise sources) 
  // from the Symbol Table
  const NodeNameMap &getNoiseDeviceNameMap() const;
  
  // used to get the noise type names (where each device may have
  // multiple noise sources/types) from the Symbol Table
  const NodeNameMap &getNoiseTypeNameMap() const;

  void setOutputFilenameSuffix( std::string newSuffix )
  {
    filenameSuffix_ = newSuffix;
  };

  const std::string &getFilenameSuffix() const
  {
    return filenameSuffix_;
  }

  const std::string &getTitle() const
  {
    return title_;
  }

  void setTitle(const std::string &title)
  {
    title_ = title;
  }

  const std::string &getNetlistFilename() const
  {
    return netlistFilename_;
  }

  // Runs specified output commands
  void output(
    Parallel::Machine                   comm,
    const double                        time,
    const double                        timeStep,
    const double                        circuit_temp, 
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const int                           dcNumber,
    const int                           maxDC,
    const Analysis::SweepVector &       dcParamVec1,
    const Linear::Vector &              solnVecPtr,
    const Linear::Vector &              stateVecPtr,
    const Linear::Vector &              storeVecPtr,
    const Linear::Vector &              lead_current_vector,
    const Linear::Vector &              junction_voltage_vector,
    const Linear::Vector &              lead_current_q_vector,
    const std::vector<double> &         objectiveVec,
    const std::vector<double> &         dOdpVec,
    const std::vector<double> &         dOdpAdjVec,
    const std::vector<double> &         scaled_dOdpVec,
    const std::vector<double> &         scaled_dOdpAdjVec,
    bool                                skipPrintLineOutput = false);

  void outputAC(
    Parallel::Machine                   comm,
    double                              freq,
    double                              fStart,
    double                              fStop,
    const Linear::Vector &                freqDomainSolnVecReal,
    const Linear::Vector &                freqDomainSolnVecImaginary,
    const Util::Op::RFparamsData &      RFparams);

  void outputSensitivityAC(
     Parallel::Machine                   comm,
     double                              freq,
     const Linear::Vector &              freqDomainSolnVecReal,
     const Linear::Vector &              freqDomainSolnVecImaginary,
     const std::vector<double> &         paramVals,
     const std::vector<std::string> &    paramNameVec,
     const std::vector<std::string> &    objFuncVars,
     const std::vector<double> &         objectiveVec,
     const std::vector<double> &         dOdpVec,
     const std::vector<double> &         dOdpAdjVec,
     const std::vector<double> &         scaled_dOdpVec,
     const std::vector<double> &         scaled_dOdpAdjVec);

  void outputSParams(
    Parallel::Machine                   comm,
    double                              freq,
    double                              numFreq,
    std::vector<double> &               Z0sVec,
    const Util::Op::RFparamsData & RFparams);

  void outputNoise(
    Parallel::Machine                   comm,
    double                              freq,
    const Linear::Vector &              freqDomainSolnVecReal,
    const Linear::Vector &              freqDomainSolnVecImaginary,
    double                              totalOutputNoiseDens_, 
    double                              totalInputNoiseDens_, 
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_);

  void outputMPDE(
    Parallel::Machine                   comm,
    double                              time,
    const std::vector<double> &         fast_time_points,
    const Linear::BlockVector &           solution_vector);

  void outputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         param_values,
    const Linear::Vector &                solution_vector);

  // Used for HB time-domain output such as .PRINT HB_TD lines.  This is
  // not used for .PRINT HB_STARTUP or .PRINT HB_IC lines though.
  void outputHB_TD(
    Parallel::Machine                   comm,
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const std::vector<double> &         timePoints,
    const Linear::BlockVector &         timeDomainSolutionVec,
    const Linear::BlockVector &         timeDomainLeadCurrentVec,
    const Linear::BlockVector &         timeDomainJunctionVoltageVecReal);

  // Used for HB frequency-domain output such as .PRINT HB_FD lines.
  void outputHB_FD(
    Parallel::Machine                   comm,
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const std::vector<double> &         freqPoints,
    const Linear::BlockVector &         freqDomainSolutionVecReal,
    const Linear::BlockVector &         freqDomainSolutionVecImaginary,
    const Linear::BlockVector &         freqDomainLeadCurrentVecReal,
    const Linear::BlockVector &         freqDomainLeadCurrentVecImaginary,
    const Linear::BlockVector &         freqDomainJunctionVoltageVecReal,
    const Linear::BlockVector &         freqDomainJunctionVoltageVecImaginary );


  // Runs specified output commands
  void outputSensitivity(
    Parallel::Machine                   comm,
    const double                        time,
    const double                        timeStep,
    const double                        circuit_temp, 
    const int                           stepNumber,
    const int                           maxStep,
    const Analysis::SweepVector &       stepParamVec1,
    const int                           dcNumber,
    const int                           maxDC,
    const Analysis::SweepVector &       dcParamVec1,
    const Linear::Vector &              solnVecPtr,
    const Linear::Vector &              stateVecPtr,
    const Linear::Vector &              storeVecPtr,
    const Linear::Vector &              lead_current_vector,
    const Linear::Vector &              junction_voltage_vector,
    const Linear::Vector &              lead_current_q_vector,
    const std::vector<double> &         objectiveVec,
    const std::vector<double> &         dOdpVec,
    const std::vector<double> &         dOdpAdjVec,
    const std::vector<double> &         scaled_dOdpVec,
    const std::vector<double> &         scaled_dOdpAdjVec,
    bool                                skipPrintLineOutput = false);

  void finishOutput();
  void finishSensitivityOutput();

  void startStep(
    int                           step,
    int                           max_step);

  void resetIndex();

  void setAliasNodeMap(const AliasNodeMap &alias_node_map)
  {
    aliasNodeMap_ = alias_node_map;
  }

  const IO::AliasNodeMap &getAliasNodeMap() const
  {
    return aliasNodeMap_;
  }

  template <class It>
  void setMainContextFunctionMap(It first, It last)
  {
    mainContextFunctionMap_.insert(first, last);
  }

  const Util::ParamMap &getMainContextFunctionMap() const
  {
    return mainContextFunctionMap_;
  }

  template <class It>
  void setMainContextParamMap(It first, It last)
  {
    for (; first != last; ++first)
      mainContextParamMap_.insert(Util::ParamMap::value_type((*first).tag(), *first));
  }

  const Util::ParamMap &getMainContextParamMap() const
  {
    return mainContextParamMap_;
  }

  template <class It>
  void setMainContextGlobalParamMap(It first, It last)
  {
    for (; first != last; ++first)
      mainContextGlobalParamMap_.insert(Util::ParamMap::value_type((*first).tag(), *first));
  }

  const Util::ParamMap &getMainContextGlobalParamMap() const
  {
    return mainContextGlobalParamMap_;
  }

  const Analysis::SweepVector &getStepSweepVector() const
  {
    return outputState_.stepSweepVector_;
  }

  int getStepLoopNumber() const
  {
    return outputState_.stepLoopNumber_;
  }

  int getMaxParamSteps() const
  {
    return outputState_.stepMaxCount_;
  }

  const Analysis::SweepVector &getDCSweepVector() const
  {
    return outputState_.dcSweepVector_;
  }

  int getDCLoopNumber() const
  {
    return dcLoopNumber_;
  }

  int getMaxDCSteps() const
  {
    return maxDCSteps_;
  }

  bool getTempSweepFlag() const
  {
    return outputState_.tempSweepFlag_;
  }

  Util::ParamList getVariableList() const;

  double getPRINTDCvalue() const
  {
    return PRINTdcvalue_;
  }

  double getPRINTDCstart() const
  {
    return PRINTdcstart_;
  }

  double getPRINTDCstop() const
  {
    return PRINTdcstop_;
  }

  const std::string &getPRINTDCname() const
  {
    return PRINTdcname_;
  }

  bool getPrintHeader () const
  {
    return printHeader_;
  }

  bool getPrintFooter () const
  {
    return printFooter_;
  }

  bool getOutputVersionInRawFile() const
  {
    return outputVersionInRawFile_;
  }

  void setEnableHomotopyFlag(bool value)
  {
    enableHomotopyFlag_ = value;
  }

  // used to determine, for -r output, whether a
  // .LIN analysis is being done.
  void setEnableSparCalcFlag(bool value)
  {
    enableSparCalcFlag_ = value;
  }

  double getCircuitTime() const
  {
    return outputState_.circuitTime_;
  }

  void setCircuitTime(double time)
  {
    outputState_.circuitTime_ = time;
  }

  double getCircuitTimeStep() const
  {
    return outputState_.circuitTimeStep_;
  }

  void setCircuitTimeStep(double timeStep)
  {
    outputState_.circuitTimeStep_ = timeStep;
  }

  double getCircuitFrequency() const
  {
    return outputState_.circuitFrequency_;
  }

  void setCircuitFrequency(double frequency)
  {
    outputState_.circuitFrequency_ = frequency;
  }

  double getCircuitTemp() const
  {
    return outputState_.circuitTemp_;
  }

  double getTemperature() const
  {
    return outputState_.circuitTemp_;
  }

  double getTime() const
  {
    return outputState_.circuitTime_;
  }

  double getFrequency() const
  {
    return outputState_.circuitFrequency_;
  }

  double getStepSweep(size_t index) const
  {
    return outputState_.stepSweepVector_[index].currentVal;
  }

  double getDCSweep(size_t index) const
  {
    return outputState_.dcSweepVector_[index].currentVal;
  }

  // outputter processing functions
  void clearActiveOutputters();
  void pushActiveOutputters();
  void popActiveOutputters();
  void addActiveOutputter(PrintType::PrintType print_type, Analysis::Mode analysis_mode);
  void addOutputter(PrintType::PrintType print_type, Outputter::Interface *outputter);

  void addOutputPrintParameters(OutputType::OutputType output_type, const PrintParameters &print_parameters);
  void addExternalOutputInterface(ExternalOutputInterface * theOutputInterface);

  std::pair<OutputParameterMap::const_iterator, bool> findOutputParameter(OutputType::OutputType output_type) const
  {
    OutputParameterMap::const_iterator it = outputParameterMap_.find(output_type);
    return std::pair<OutputParameterMap::const_iterator, bool>(it, it != outputParameterMap_.end());
  }

    std::pair<ExternalOutputWrapperMap::const_iterator, bool> findExternalOutputWrapper(OutputType::OutputType output_type) const
  {
    ExternalOutputWrapperMap::const_iterator it = externalOutputWrapperMap_.find(output_type);
    return std::pair<ExternalOutputWrapperMap::const_iterator, bool>(it, it != externalOutputWrapperMap_.end());
  }

  const PrintParameters &getDefaultPrintParameters() const
  {
    return defaultPrintParameters_;
  }

  const Util::Op::BuilderManager &getOpBuilderManager() const
  {
    return opBuilderManager_;
  }

  std::ostream *openFile(const std::string &path, std::ios_base::openmode mode);
  std::ostream *openFile(const std::string &path);
  std::ostream *openBinaryFile(const std::string &path);
  int closeFile(std::ostream *os);

  bool prepareHDF5Output(Parallel::Machine comm);
  bool updateHDF5Output(Parallel::Machine comm, const Linear::Vector &solnVecPtr);
  bool closeHDF5Output();

  const OutputParameterMap &getOutputParameterMap() const
  {
    return outputParameterMap_;
  }

private:
  void steppingComplete();

  void update_HB_FD_printParamsForPrintFormat(PrintParameters &freq_print_parameters);
  void update_HB_TD_printParamsForPrintFormat(PrintParameters &time_print_parameters);
  void update_HB_IC_printParamsForPrintFormat(PrintParameters &hb_ic_print_parameters);
  void update_HB_STARTUP_printParamsForPrintFormat(PrintParameters &hb_startup_print_parameters);

  void update_MPDE_IC_printParamsForPrintFormat(PrintParameters &mpde_ic_print_parameters);
  void update_MPDE_STARTUP_printParamsForPrintFormat(PrintParameters &mpde_startup_print_parameters);

private:
  std::string                   title_;
  std::string                   netlistFilename_;
  std::string                   filenameSuffix_;

  Util::Op::BuilderManager &    opBuilderManager_;
  const Topo::Topology &        topology_;
  OutputterMap                  outputterMap_;
  OutputParameterMap            outputParameterMap_;
  ExternalOutputWrapperMap      externalOutputWrapperMap_;
  ActiveOutputterStack          activeOutputterStack_;  ///< stack of outputter sets.  
                                                        ///< Only the "back" of this stack is the "active" set. 
  EnabledAnalysisSet            enabledAnalysisSet_;

  struct OutputState
  {
    OutputState()
      : circuitTime_(0.0),
        circuitTimeStep_(0.0),
        circuitTemp_(27.0),
        circuitFrequency_(0.0),
        stepLoopNumber_(0),
        stepMaxCount_(0),
        tempSweepFlag_(false)
    {}

    double                      circuitTime_;           ///< transient circuit time:
    double                      circuitTimeStep_;       ///< transient circuit time step, dt:
    double                      circuitTemp_;           ///< circuit temperature
    double                      circuitFrequency_;      ///< ac current frequency
    int                         stepLoopNumber_;        ///< step current loop
    int                         stepMaxCount_;          ///< step max count

    bool                        tempSweepFlag_;
    Analysis::SweepVector       stepSweepVector_;
    Analysis::SweepVector       dcSweepVector_;
  };

  OutputState                   outputState_;

  // print statement vars
  bool                  dotOpSpecified_; // flag to indicate if the netlist has a .OP statement
  bool                  enableHomotopyFlag_;
  bool                  enableSparCalcFlag_;
  bool                  enableSensitivityFlag_;
  bool                  adjointSensitivityFlag_;
  int                   sensitivityOptions_;
  Util::ParamList       sensitivityVariableList_;
  Util::ParamList       transientAdjointVariableList_;

  int                   pts_per_summary_;
  bool                  pts_per_summary_Given;
  Util::ParamList       noiseVariableList_;

  PrintParameters       defaultPrintParameters_;

  double                PRINTdcstart_;
  double                PRINTdcstop_;
  double                PRINTdcvalue_;
  std::string           PRINTdcname_;

  double                initialOutputInterval_;
  IntervalVector        outputIntervalPairs_;

  bool                  STEPEnabledFlag_;

  bool                  printHeader_;               // flag to indicate if user wants the header in output file.
  bool                  printFooter_;               // flag to indicate if user wants the "End of Xyce(TM)" line in the output.
  bool                  outputVersionInRawFile_;    // flag to indicate that Version should be output in the header of a RAW file.

  bool outputCalledBefore_;

  // dc loop information
  int dcLoopNumber_;
  int maxDCSteps_;

  bool isStarredPrintLineProcessed;

  AliasNodeMap          aliasNodeMap_;
  Util::ParamMap        mainContextFunctionMap_;
  Util::ParamMap        mainContextParamMap_;
  Util::ParamMap        mainContextGlobalParamMap_;

  OutputterFactoryMap   outputterFactoryMap_;
  OpenPathStreamMap     openPathStreamMap_;

  // HDF5 file support vars
  bool hdf5FileNameGiven_;
  bool hdf5HeaderWritten_;
  std::string hdf5FileName_;
  int hdf5IndexValue_;

#ifdef Xyce_USE_HDF5
  hid_t hdf5FileId_;
  hid_t hdf5PlistId_;
#endif  // Xyce_USE_HDF5

};

//-----------------------------------------------------------------------------
// Function      : OutputMgr::clearActiveOutputters
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
inline void OutputMgr::clearActiveOutputters()
{
  activeOutputterStack_.clear();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::pushActiveOutputters
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// @brief push an empty vector onto active outputter stack
///
/// This method pushes an empty vector of outputter pointers onto the
/// activeOutputterStack_ (a vector of vectors of pointers).  The "back"
/// of this stack is the vector of currently selected ("active") outputters.
///
/// New outputters added to the active list will be put into this back vector.
///
/// @author Dave Baur
/// @date 7/16/2013
inline void OutputMgr::pushActiveOutputters()
{
  activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::popActiveOutputters
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// @brief Pop active outputter stack
///
/// This method removes the "back" element of the outputter stack,
/// restoring the previous set to "active"
///
/// @author Dave Baur
/// @date 7/16/2013
inline void OutputMgr::popActiveOutputters()
{
  activeOutputterStack_.pop_back();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addActiveOutputter
// 
// Purpose       : ERK: sets previously allocated outputters (which are in the 
//                 outputterMap object) as "active".  In this context "Active"
//                 means that this outputter (or set of outputters) is
//                 currently being used.  
//
//                 This is useful in the following scenarios:
//                  - if netlist has multiple analysis types specified, and
//                    multiple corresponding .PRINT statements.  If, for 
//                    example, the netlist has both .AC and .TRAN, then 
//                    the .PRINT TRAN outputter will only be active during 
//                    the transient analysis and will be inactive during 
//                    the AC analysis.
//
// Comments from Dave Baur:
//
// The mapping from .PRINT <print-type> to outputters and then activiting the
// outputters at the appropriate time in the appropriate context is the
// trick here.
//
// The user communicates the output variables wanted via the .PRINT
// <print-type>.  The .PRINT <print-type> creates a mapping from
// print-type to output-type where each output-type represents a set of
// parameters to the written when the print-type is active.  The
// preparePrintParameters() function creates a mapping from .PRINT
// <print-type> to outputters, using the output-types to identify the
// parameters.  The application informs the output manager
// what needs to be written when called via the ActiveOutput sentry
// class and its functions.
//
// The ActiveOutput sentry class manages the stack of active
// outputters.  The construction pushes an empty vector of outputters
// onto the stack and the destructor pops the vector off, restoring the
// previous vector.
//
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// @brief Add all outputters of the selected type to the current active set
///
/// Sets previously allocated outputters (which are in the outputterMap
/// object) as "active".  In this context "Active" means that this
/// outputter (or set of outputters) is currently being used.
///
/// This is useful in the following scenarios:
/// - if netlist has multiple analysis types specified, and
///   multiple corresponding .PRINT statements.  If, for 
///   example, the netlist has both .AC and .TRAN, then 
///   the .PRINT TRAN outputter will only be active during 
///   the transient analysis and will be inactive during 
///   the AC analysis.
///
/// @param[in] print_type  .print line type to select
/// @param[in] analysis_mode  what analysis mode we are running
///
/// @author Dave Baur
/// @date 7/16/2013
inline void OutputMgr::addActiveOutputter(
    PrintType::PrintType print_type, Analysis::Mode analysis_mode)
{
  OutputterMap::iterator find_it = outputterMap_.find(print_type);

  if (find_it == outputterMap_.end())
  {
    // ERK:  Seems like there should be an error thrown here?  Or possibly 
    // it simply needs to be allocated?  As it is, it is a no-op.
  }

  if (!activeOutputterStack_.empty() && find_it != outputterMap_.end()) 
  {
    OutputterMap::mapped_type::iterator it  = (*find_it).second.begin(); 
    OutputterMap::mapped_type::iterator end = (*find_it).second.end();
    for ( ; it != end; ++it)
    {
      (*it)->setAnalysisMode(analysis_mode);
    }

    activeOutputterStack_.back().insert(activeOutputterStack_.back().end(), 
        (*find_it).second.begin(), (*find_it).second.end());
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addOutputter
// Purpose       : adds an allocated outputter to the outputterMap_ object.
//
// Special Notes : The outputterMap_ object contains all the allocated 
//                 outputters.  Some of these will be active, some may not,
//                 depending on the current analysis type.
//
//                 The map is indexed by the print type.  eg, .PRINT TRAN,
//                 where "TRAN" is the print type.  Each print type will
//                 have a vector of outputter pointers, meaning that there can
//                 be a bunch of separate print statements for a given type.
//                 (ie, you can have .PRINT TRAN format=tecplot and .PRINT TRAN
//                 format=probe, and they'll both be members of a vector that
//                 is indexed by PrintType::TRAN).
//
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
/// Add an outputter to the outputterMap_ object
///
/// @note This function is called by all of the "enable<type>Output" functions
/// in the Xyce::IO::Outputter namespace, each of which allocates
/// outputter objects for each print parameters found that matches
/// the correct type.  I (TVR) have not found any evidence that these
/// objects are ever deleted, nor is the outputterMap_ ever cleared.
/// That means that we are implicitly assuming that only one analysis
/// will ever be performed, and that each enable function will only be
/// called once.  If we ever modify Xyce so that multiple analyses can be
/// run in a single netlist run, then this will lead to leaks where
/// outputters that were allocated for a previous analysis won't be discarded.
/// And heaven forbid that we ever allow one to do more than one analysis of
/// a given type in the same run --- this would cause us to re-allocate
/// outputters for the same set of print parameters and add them a second
/// time to the outputter map.  Here there be tygers?
///
inline void OutputMgr::addOutputter(PrintType::PrintType print_type, Outputter::Interface *outputter)
{
  outputterMap_[print_type].push_back(outputter);
}

//-----------------------------------------------------------------------------
// Function      : operator<
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool operator<(const OutputMgr::OutputterKey &lhs, const OutputMgr::OutputterKey &rhs)
{
  return (lhs.analysisMode_ < rhs.analysisMode_)
    || (!(rhs.analysisMode_ < lhs.analysisMode_) && lhs.outputType_ < rhs.outputType_)
    || (!(rhs.analysisMode_ < lhs.analysisMode_ && rhs.outputType_ < lhs.outputType_) && lhs.format_ < rhs.format_);
}

bool registerPkgOptionsMgr(OutputMgr & output_manager, PkgOptionsMgr &options_manager);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputMgr_h
