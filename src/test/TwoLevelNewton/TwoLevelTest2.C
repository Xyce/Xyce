// 
// Simple test of the two level API but to transmission lines in the circuit
//


#include <Xyce_config.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Math.h>
#include <N_CIR_SecondLevelSimulator.h>
#include <N_DEV_ExternalSimulationData.h>
#include <N_TIA_TwoLevelError.h>


//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  Xyce::lout() << "TwoLevelTest2: Starting Xyce as library test" << std::endl;

  bool xyceSuccess = false;
  // Set divide by zero, and invalid operation handling on linux
  // Set out of memory detection on all systems
  //std::set_new_handler (&_new_handler);

  //std::string netlist = "resInnerTran.cir";
  //Xyce::lout() << "---------------\nRunning the 2-level calculation:\n---------------" << netlist << std::endl;
  //topLevelNewton circuitCalculation(netlist);
  //bool bsuccess = circuitCalculation.runTran(iargs,cargs);

  auto xyceSimulator = new Xyce::Circuit::SecondLevelSimulator(); // don't need a comm object (optional argument)
  if( xyceSimulator == NULL)
  {
    // return a fail signal to the OS so the test fails
    Xyce::dout() << "Failed in call to allocate xyce second level simulator object" << std::endl;
    return -1;
  }
  const std::vector<std::string> argsToXyce = {"compInner2.cir"};
  
  xyceSuccess = xyceSimulator->initialize(argsToXyce);
  if( !xyceSuccess)
  {
    Xyce::dout() << "Failed in call to xyce initialize" << std::endl;
    return -1;
  }
  
  xyceSuccess = xyceSimulator->startupSolvers();
  if( !xyceSuccess)
  {
    Xyce::dout() << "Failed in call to xyce setup solvers" << std::endl;
    return -1;
  }
  
  // extSimData is a struct to hold lots of time step info for Xyce
  Xyce::Device::ExternalSimulationData  extSimData;
  // get time step data from inner problem that Xyce has loaded
  xyceSimulator->endTimeStep(extSimData);
  
  // transient source sine source vars
  double V0 = 0.0;
  double VA = 1000000.0;
  double FREQ = 2.5e+10;
  
  // voltageInputMap is passed to Xyce to set the transient inputs on the devices
  // vconnectXXXX
  std::map<std::string,double> voltageInputMap;
  voltageInputMap["vconnect0000"] = V0 + VA * std::sin(2.0*M_PI*(FREQ*0.0));;
  voltageInputMap["vconnect0001"] = 0.0;
  
  // vector and jacobian for communicating results from inner to outer problem
  // the output is the current the inner-problem needs from the vconnectXXXX devices
  // to support the boundary condition set by vconnectXXXX
  std::vector<double> outputVector;
  outputVector.resize(2,0.0);
  
  std::vector< std::vector<double> > innerJacobianStamp;
  // inner jacobian stamp is just for the two voltage source vconnectXXXX -- not the entire inner problem.
  std::vector<double> row(2,0.0);
  innerJacobianStamp.push_back(row);
  innerJacobianStamp.push_back(row);
  
  // Do the DC operating point calculation if needed:
  Xyce::lout() << "Begin DCOP step (inner solve)"  << std::endl;
  extSimData.is_transient = true;
  extSimData.current_time = 0.0;
  extSimData.final_time = extSimData.finalTime;
  extSimData.current_time_step_size = extSimData.nextTimeStep;
  extSimData.previous_time_step_size = extSimData.currTimeStep;
  extSimData.time_step_number = 0;
  extSimData.forceOrder=true;
  extSimData.imposedTimeIntegrationOrder=1;
  extSimData.forceBeginningIntegration=true;
  extSimData.imposedBeginningIntegration=true;
  
  
  Xyce::TimeIntg::TwoLevelError tlError;  // data structure related to global error control
  
  // signals to Xyce the state of the outer to inner junction.  Only important at DC op.
  bool initJctFlag = false;
  
  // note for the DC op we are not passing in voltageInputMap as that is set by the 
  // initial conditions on the inner circuit problem (i.e. either DC conditions on 
  // the voltage source set to zero by using NOOP)
  xyceSuccess =  xyceSimulator->startTimeStep(extSimData);
  xyceSuccess = xyceSimulator->simulateStep(initJctFlag, voltageInputMap, outputVector, innerJacobianStamp, tlError);
  if( !xyceSuccess )
  {
    Xyce::dout() << "Failed in call to xyce startTimeStep" << std::endl;
    return -1;    
  }
  
  if (xyceSuccess)
  {
    // process "successful" time step
    xyceSimulator->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP);
    xyceSimulator->endTimeStep(extSimData);
  }
  else
  {
    Xyce::dout() << "Initial Xyce newtonSolve step failed to converge. Exiting" << std::endl;
    return -1;
  }

  // now try transient stepping
  // set up a max time step
  double max_dt = extSimData.finalTime/10.0;

  //---------------------------------------------------------------------------
  // transient calculation. (post DCOP)
  extSimData.beginIntegrationFlag=true;
  int iTimeStep=0;
  bool finished = false;

  // time loop
  while(!finished)
  {
    // before resetting the sim data, check against max_dt.
    if ( extSimData.nextTimeStep > max_dt)
    {
      extSimData.nextTime -= extSimData.nextTimeStep;
      extSimData.nextTime += max_dt;
      extSimData.nextTimeStep = max_dt;
    }

    //resetSimData(extSimData);
    extSimData.is_transient = true;
    extSimData.current_time = extSimData.nextTime;
    extSimData.final_time = extSimData.finalTime;

    extSimData.current_time_step_size = extSimData.nextTimeStep;
    extSimData.previous_time_step_size = extSimData.currTimeStep;
    extSimData.time_step_number = extSimData.timeStepNumber;

    extSimData.forceOrder=true;
    extSimData.imposedTimeIntegrationOrder= extSimData.currentOrder;

    extSimData.forceBeginningIntegration=true;
    extSimData.imposedBeginningIntegration = extSimData.beginIntegrationFlag;

    Xyce::dout() 
      << "Begin Time step (outer solve) # " 
      << iTimeStep+1 
      << " current,step = " 
      << extSimData.nextTime
      << "," 
      << extSimData.nextTimeStep
      << std::endl;

    // update the sinewave source:
    double time = (extSimData.nextTime);
    voltageInputMap["vconnect0000"] = V0 + VA * std::sin(2.0*M_PI*(FREQ*time));

    xyceSuccess = xyceSimulator->startTimeStep(extSimData);
    xyceSuccess = xyceSimulator->simulateStep(initJctFlag, voltageInputMap, outputVector, innerJacobianStamp, tlError);
    
    if (xyceSuccess)
    {
      // process "successful" time step
      xyceSimulator->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT); 
      Xyce::dout() << "Current draw from inner problem:"  << outputVector[0] << ", " << outputVector[1] << std::endl;
    }
    else
    {
      Xyce::dout() << "Time step failed to converge. Exiting" << std::endl;
      return -1;
    }

    xyceSimulator->endTimeStep(extSimData);

    if (extSimData.current_time >= extSimData.finalTime) 
    { 
      finished = true; 
    }

    //Xyce::dout() << extSimData;

    extSimData.beginIntegrationFlag=false;
    ++iTimeStep;
  }

  
  xyceSuccess = xyceSimulator->finishSolvers();
  if( !xyceSuccess)
  {
    Xyce::dout() << "Failed in call to xyce finish solvers" << std::endl;
    return -1;
  }
  
  delete xyceSimulator;
  return 0;
}
