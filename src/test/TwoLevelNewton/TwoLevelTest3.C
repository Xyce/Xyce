// 
// Simple test of the two level API but to transmission lines in the circuit
// 
// In the outer problem use a  non-ideal power source (a discharging capacitor)
// inner problem uses an RLC transmission line model
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
  Xyce::lout() << "TwoLevelTest3: Starting Xyce as library test" << std::endl;

  bool xyceSuccess = false;

  auto xyceSimulator = new Xyce::Circuit::SecondLevelSimulator(); // don't need a comm object (optional argument)
  if( xyceSimulator == NULL)
  {
    // return a fail signal to the OS so the test fails
    Xyce::dout() << "Failed in call to allocate xyce second level simulator object" << std::endl;
    return -1;
  }
  const std::vector<std::string> argsToXyce = {"compInner3.cir"};
  
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
  
  // voltageInputMap is passed to Xyce to set the transient inputs on the devices
  // vconnectXXXX
  // I = C dV/dt
  // set Initial voltage to 1e5
  const double outerCap = 1e-12;
  double V0 = 1.0e5;
  
  std::map<std::string,double> voltageInputMap;
  voltageInputMap["vconnect0000"] = 0.0;
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
    // I = C dV/dt
    // I = C (V(t1) - V(t0)) / dt 
    // V(t1) = I*dt/C + V(t0)
    // simple euler update scheme.
    voltageInputMap["vconnect0000"] = outputVector[0]*extSimData.nextTimeStep/outerCap + V0;

    xyceSuccess = xyceSimulator->startTimeStep(extSimData);
    xyceSuccess = xyceSimulator->simulateStep(initJctFlag, voltageInputMap, outputVector, innerJacobianStamp, tlError);
    
    if (xyceSuccess)
    {
      // process "successful" time step
      xyceSimulator->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT); 
      Xyce::dout() << "Current draw from inner problem:"  << outputVector[0] << ", " << outputVector[1] << std::endl;
      V0 = voltageInputMap["vconnect0000"];
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
