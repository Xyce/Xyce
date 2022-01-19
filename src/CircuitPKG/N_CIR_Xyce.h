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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_CIRCUIT_h
#define Xyce_N_CIR_CIRCUIT_h

#include <string>
#include <map>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_ERH_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_NLS_fwd.h>
#include <N_TIA_fwd.h>

#include <N_UTL_ReportHandler.h>
#include <N_UTL_Stats.h>
#include <N_UTL_JSON.h>
#include <N_DEV_ADC.h>
#include <N_DEV_DAC.h>

#include <N_IO_CmdParse.h>

namespace Xyce {

namespace Circuit {

//-----------------------------------------------------------------------------
// Class         : Xyce
// Purpose       : This is the main "top level" class for Xyce.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class Simulator
{
 public:
  enum RunState {START, PARALLEL_INIT, PARSE_COMMAND_LINE, CHECK_NETLIST, OPEN_LOGSTREAM, ALLOCATE_SUBSYSTEMS, PARSE_NETLIST, SETUP_TOPOLOGY,
                 INSTANTIATE_DEVICES, SETUP_MATRIX_STRUCTURE, INITIALIZE_SYSTEM};
  enum RunStatus {ERROR, SUCCESS, DONE};

  Simulator(Parallel::Machine comm = MPI_COMM_NULL);

  virtual ~Simulator();

  virtual Analysis::AnalysisManager *newAnalysisManager(
    const IO::CmdParse &                command_line,
    IO::RestartMgr &                    restart_manager,
    Analysis::OutputMgrAdapter &        output_manager_adapter,
    Stats::Stat                         analysis_stat);

  Device::DeviceMgr &getDeviceManager() {
    return *deviceManager_;
  }

  IO::OutputMgr &getOutputManager() {
    return *outputManager_;
  }

  Analysis::AnalysisManager &getAnalysisManager() {
    return *analysisManager_;
  }

  Nonlinear::Manager &getNonlinearManager() {
    return *nonlinearManager_;
  }

  Linear::System &getLinearSystem() {
    return *linearSystem_;
  }

  Loader::CktLoader &getCircuitLoader();
  
  // These are all the API calls that we are suppose to be making available
  // for external programs and/or other objects

  //---------------------------------------------------------------------------
  // Function      : setNetlistParameters
  // Purpose       : This passes a vector of pairs "key", "value" that will
  //                 be substituted during the processing of the netlist.  This
  //                 more easily allows Dakota to change any netlist parameter
  //                 during netlist setup.
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical and MEMS Modeling
  // Creation Date : 10/9/2008
  //---------------------------------------------------------------------------
  void setNetlistParameters( const std::vector< std::pair< std::string, std::string > > & externalParams );


  //---------------------------------------------------------------------------
  // Function      : setNetlistParameters
  // Purpose       : Call through to the output manager to set the suffix to
  //                 be used on the output file, as in circuit + suffix + prn
  //                 This is useful in Dakota controlled runs to keep each
  //                 simulation from overwritting the last one.
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical and MEMS Modeling
  // Creation Date : 10/9/2008
  //---------------------------------------------------------------------------
  void setOutputFileSuffix( const std::string newSuffix );

  //---------------------------------------------------------------------------
  // Function      : run
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 02/19/01
  //---------------------------------------------------------------------------
  RunStatus run(int argc, char **argv);

  //---------------------------------------------------------------------------
  // Function      : initialize
  // Purpose       : To initialize Xyce to be driven by the SAX
  //                 simulation backplane. This includes the following:
  //                    Set up and register the Parallel Manager.
  //                    Parse the command-line arguments.
  //                    Redirect the output stream of processor 1,
  //                    if requested.
  //                    Read in the Netlist.
  //                    Allocate and register the external packages.
  //                    Set up the representation of the circuit topology.
  //                    Set up the matrix structures.
  //                    Initialize the solvers.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  RunStatus initialize(int argc, char **argv);

  /// First initialization call, intializes up to point where
  /// topology needs to query devices for number variables and jacstamp.
  /// This include reading netlist and instantiating devices.
  RunStatus initializeEarly(int argc, char **argv);

  /// Second initialization call, intialization steps starting from
  /// where topology needs to query devices for number variables and jacstamp
  RunStatus initializeLate();

  //---------------------------------------------------------------------------
  // Function      : runSimulation
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 5/27/00
  //---------------------------------------------------------------------------
  RunStatus runSimulation();

  bool getDeviceNames(const std::string & modelGroupName, std::vector<std::string> & deviceNames);
  bool getAllDeviceNames(std::vector<std::string> & deviceNames);
  bool getDeviceParamVal(const std::string full_param_name, double& val) const;
  bool checkDeviceParamName(const std::string full_param_name) const;

  bool getAdjGIDsForDevice(const std::string deviceName, std::vector<int> & adj_GIDs) const;
  bool getNumAdjNodesForDevice(const std::string deviceName, int & numAdjNodes) const;

  //---------------------------------------------------------------------------
  // Function      : getDACDeviceNames
  // Purpose       : Gets the (stripped) names of the DAC devices
  //                 in the circuit.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/13/03
  //---------------------------------------------------------------------------
  bool getDACDeviceNames(std::vector< std::string >& dacNames);

  //---------------------------------------------------------------------------
  // Function      : getADCMap
  // Purpose       : Gets the (stripped) names of the ADC devices
  //                 in the circuit(as key of map) and map of parameters
  //                 (keyed by parameter name) for each device
  // Special Notes :
  // Scope         : public
  // Creator       : Tom Russo, SNL, Component Information and Models
  // Creation Date : 05/06/2004
  //---------------------------------------------------------------------------
  bool getADCMap(std::map<std::string,std::map<std::string,double> >& ADCMap);

  //---------------------------------------------------------------------------
  // Function      : updateTimeVoltagePairs
  // Purpose       : Update the DAC devices in a circuit by adding the set
  //                 of time and voltage pairs built up on the "digital side"
  //                 since the last update and by removing the time-voltage
  //                 pairs for times that pre-date the given simulation time.
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 06/10/03
  //---------------------------------------------------------------------------
  bool updateTimeVoltagePairs(
        std::map< std::string, std::vector< std::pair<double,double> >* > const&
        timeVoltageUpdateMap);

  //---------------------------------------------------------------------------
  // Function      : getTimeVoltagePairs
  // Purpose       : query the DAC devices in a circuit for the set
  //                 of time and voltage pairs
  // Special Notes : Calling this function clears the ADC devices time-voltage
  //                 pair vector.  
  // Scope         : public
  // Creator       : Tom Russo, SNL ComponentInformation and Models
  // Creation Date : 05/10/2004
  //---------------------------------------------------------------------------
  bool getTimeVoltagePairs(
        std::map< std::string, std::vector< std::pair<double,double> > > &
        timeVoltageUpdateMap);
  
  //---------------------------------------------------------------------------
  // Function      : getTimeVoltagePairsSz
  // Purpose       : Returns the largest size of the TV Pairs data in the ADCs
  // Special Notes :
  // Scope         : public
  // Creator       : 
  // Creation Date : 11/11/2021
  //---------------------------------------------------------------------------
  bool getTimeVoltagePairsSz(int &maximumSize);

  //---------------------------------------------------------------------------
  // Function      : getTimeStatePairs
  // Purpose       : query the DAC devices in a circuit for the set
  //                 of time and state pairs
  // Special Notes :
  // Scope         : public
  // Creator       : Pete Sholander, SNL
  // Creation Date : 11/13/2018
  //---------------------------------------------------------------------------
  bool getTimeStatePairs(
        std::map< std::string, std::vector< std::pair<double,int> > > &
        timeStateUpdateMap);

  //----------------------------------------------------------------------------
  // Function       : setADCWidths
  // Purpose        : Update the ADC devices in a circuit by informing them
  //                  of the width of their bitvector output on the
  //                  "digital side"
  // Special Notes  :
  // Scope          :
  // Creator        : Tom Russo
  // Creation Date  : 05/07/2004
  //----------------------------------------------------------------------------
  bool setADCWidths(std::map< std::string, int > const& ADCWidthMap);

  //----------------------------------------------------------------------------
  // Function       : getADCWidth
  // Purpose        : get the width of the specified ADC devices bitvector output 
  //                  on the "digital side"
  // Special Notes  :
  // Scope          :
  // Creator        : Pete Sholander
  // Creation Date  : 11/16/2018
  //----------------------------------------------------------------------------
  bool getADCWidths(std::map< std::string, int > & ADCWidthMap);

  //---------------------------------------------------------------------------
  // Function      : simulateUntil
  // Purpose       : To continue the existing analog circuit simulation
  //                 until either the given <requestedUntilTime> is reached
  //                 or the simulation termination criterion is met.
  //                 Return a Boolean indicating whether the simulation
  //                 run was successful. (Note that the run is successful
  //                 even when the given <requestedUntilTime> is not reached,
  //                 so long as the run completed normally.)
  // Special Notes : The time variables are in units of seconds.
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/28/03
  //---------------------------------------------------------------------------
  bool simulateUntil(double requestedUntilTime, double& completedUntilTime);

  //---------------------------------------------------------------------------
  // Function      : finalize
  // Purpose       : To clean up after driving Xyce with the SIMBUS
  //                 simulation backplane. This includes the following:
  //                    Free any dynamically allocated memory...
  // Special Notes :
  // Scope         : public
  // Creator       : Lon Waters, Lisa Renee Maynes
  // Creation Date : 05/29/03
  //---------------------------------------------------------------------------
  RunStatus finalize();

  void reportTotalElapsedTime ();

  bool checkResponseVar(const std::string &variable_name) const;
  bool obtainResponse(const std::string& variable_name, double &result) const;

  // // report on whether simulation is finished or not
  bool simulationComplete();

 private:
    bool doAllocations_();
    bool doInitializations_();
    bool doRegistrations_();
    RunStatus setupTopology( unordered_map< std::string, std::string >& aliasMap );
    bool setUpMatrixStructure_();
    bool runSolvers_();

    Device::ADC::Instance *getADCInstance_(const std::string &deviceName);
    Device::DAC::Instance *getDACInstance_(const std::string &deviceName);
    void processParamOrDoc_(std::string & optionName, std::string & deviceName,
                            int modelLevel, bool printModel, bool printInstance);
    void finalizeLeadCurrentSetup_();

 private:
  RunState                              runState_;
  Parallel::Machine                     comm_;
  IO::ParsingMgr *                      parsingManager_;                ///< Parsing Manager
  Device::DeviceMgr *                   deviceManager_;                 ///< Device Manager

  Analysis::AnalysisCreatorRegistry *   analysisRegistry_;
  Analysis::ProcessorCreatorRegistry *  processorRegistry_;
  Topo::Topology *                      topology_;
  Linear::System *                      linearSystem_;                  ///< Linear algebra system
  Linear::Builder *                     builder_;                       ///< Linear algebra system component builder
  Analysis::AnalysisManager *           analysisManager_;               ///< Analysis manager
  Loader::CktLoader *                   circuitLoader_;                 ///< Circuit loader
  Analysis::OutputMgrAdapter *          outputManagerAdapter_;          ///< Output manager adapter
  Nonlinear::Manager *                  nonlinearManager_;              ///< Nonlinear solver manager
  Parallel::Manager *                   parallelManager_;               ///< Parallel distribution manager
  Util::Op::BuilderManager *            opBuilderManager_;              ///< Op builder manager
  IO::OutputMgr *                       outputManager_;                 ///< Output manager
  IO::Measure::Manager *                measureManager_;                ///< Measure manager
  IO::FourierMgr *                      fourierManager_;                ///< Fourier manager
  IO::FFTMgr *                          fftManager_;                    ///< FFT manager
  IO::InitialConditionsManager *        initialConditionsManager_;      ///< Initial conditions manager
  IO::OutputResponse *                  outputResponse_;                ///< Response output
  IO::RestartMgr *                      restartManager_;                ///< Restart manager
  IO::LoadManager *                     loadManager_;                   ///< .LOAD processing manager
  IO::PkgOptionsMgr *                   optionsManager_;                ///< Package options manager
  Stats::Stat                           rootStat_;                      ///< Stats collection root
  Stats::Stat                           analysisStat_;                  ///< Analysis stats
  Util::JSON                            auditJSON_;                     ///< Audit JSON structure
  Util::Timer *                         XyceTimerPtr_;                  ///< Xyce solver timing utility
  Util::Timer *                         ElapsedTimerPtr_;               ///< Elapsed time from beginning of run
  unordered_set<std::string> device_names_;

  protected:
  IO::CmdParse                          commandLine_;

  private:
  IO::HangingResistor                   hangingResistor_;

  // if the user is providing an external file with parameters in it
  std::vector< std::pair< std::string, std::string> > externalNetlistParams_;
  std::map<std::string, Device::DAC::Instance *>        dacDeviceMap_;
  std::map<std::string, Device::ADC::Instance *>        adcDeviceMap_;

  REH previousReportHandler_;
};

//-----------------------------------------------------------------------------
// Function      : Xyce_exit
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/26/04
//-----------------------------------------------------------------------------
 void Xyce_exit( int code , bool asymmetric);


} // namespace Circuit
} // namespace Xyce

typedef Xyce::Circuit::Simulator N_CIR_Xyce;

#endif // Xyce_N_CIR_CIRCUIT_h
