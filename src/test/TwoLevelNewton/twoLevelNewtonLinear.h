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
// Filename       : twoLevelNewtonLinear.h
//
// Purpose        : 2-level Newton code used by several tests.  
//                  This hardwires a specific circuit for the "top level" 
//                  problem.  That circuit is described in detail in the
//                  comments.
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 12/5/2018
//
//-----------------------------------------------------------------------------

#ifndef Xyce_twoLevelNewtonLinear_h
#define Xyce_twoLevelNewtonLinear_h

#include <Xyce_config.h>

#include <baseNewton.h>

#include <N_CIR_SecondLevelSimulator.h>
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_TIA_TwoLevelError.h>
#include <N_ANP_fwd.h>
#include <N_DEV_ExternalSimulationData.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <vector>
#include <algorithm>

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class topLevelNewton : public baseNewton
{
  public:
    // variables: 
     Xyce::Circuit::SecondLevelSimulator * simulator_;

  public:

     // functions:
     // constructor
  private:
     topLevelNewton();

  public:
     topLevelNewton(std::string & netlist);

     // destructor
    ~topLevelNewton()
    {
      delete simulator_;
    }

    // this function runs the two-level problem, for a transient case.
    bool runTran( int iargs, char *cargs[]);

    // this function runs the two-level problem, for a DC sweep
    bool runDC( int iargs, char *cargs[]);

    // this function runs the whole inner problem.  No 2-level.
    bool runXyce()
    {
      if (!simulator_) return false;

      int argc = 2;
      char *argv[3];
      argv[0] = strdup("Xyce");
      argv[1] = strdup(netlistFilename_.c_str());
      argv[2] = 0;

      simulator_->run(argc, argv);

      delete[] argv[0];
      delete[] argv[1];

      return true;
    }

    void interfaceFunctions (); // documentation function. not meant to be called.
    void allocateLinearSystem(int size);
    void setupNonInnerMatrix();
    void updateInnerProblemInputs();
    void updateLinearSystem();

    bool checkConvergence();

    // variables:
    double R;
    double G;

    double v1;
    double v2;
    double v3;

    double Ir;
    double Iv0;
    double Iv1;

    double vconnect000;
    double vconnect001;

    std::map<std::string,double> voltageInputMap;
    std::vector<double> outputVector;

    Teuchos::SerialDenseMatrix< int, double > JacobianMatrix_exceptForInner;
    Teuchos::SerialDenseMatrix< int, double > RHS_exceptForInner;

    std::vector< std::vector<double> > innerJacobianStamp;

    Xyce::TimeIntg::TwoLevelError tlError;  // data structure related to global error control
    bool initJctFlag;
};

//-----------------------------------------------------------------------------
// constructor
//-----------------------------------------------------------------------------
inline topLevelNewton::topLevelNewton(std::string & netlist): 
    baseNewton (netlist),
    initJctFlag(true), R(1.0e-2), G(1/R),
    v1(0.0), v2(0.0), v3(0.0),
    Ir(0.0), Iv0(0.0), Iv1(0.0),
    vconnect000(1.0),
    vconnect001(2.0)
{
  Xyce::lout() << "netlistFilename_ = " << netlistFilename_ <<std::endl;

  int argc = 2;
  char* argv[3];
  argv[0] = strdup("Xyce");
  argv[1] = strdup(netlist.c_str());
  argv[2] = 0;

  simulator_ = new Xyce::Circuit::SecondLevelSimulator(); // don't need a comm object (optional argument)
  for(int i=0;i<argc;++i)
  {
    std::string arg(argv[i]);
    Xyce::lout() << "argv["<<i<<"] = " << arg << std::endl;
  }


  if (  ! simulator_->initialize(argc, argv) )
  {
    Xyce::dout() << "Failed to initialize Xyce for netlist " << netlistFilename_ << std::endl;
    exit(0); // fix this later
  }
  else
  {
    simulator_->startupSolvers();
  }

  free(argv[0]);
  free(argv[1]);

  // two level details
  voltageInputMap["vconnect0000"] = vconnect000;
  voltageInputMap["vconnect0001"] = vconnect001;

  outputVector.resize(2,0.0);

  std::vector<double> row(2,0.0);
  innerJacobianStamp.push_back(row);
  innerJacobianStamp.push_back(row);
}

//-----------------------------------------------------------------------------
// this function simply documents some of the API functions.  It isn't meant
// to be called.
inline void topLevelNewton::interfaceFunctions ()
{
}

//-----------------------------------------------------------------------------
inline void topLevelNewton::allocateLinearSystem(int size)
{
  // the "outer loop" linear system objects are owned and allocated/initialized
  // in the base object.
  baseNewton::allocateLinearSystem(size);

  // the rest are handled here:

  // "shape" is equivalent to resize for these objects
  JacobianMatrix_exceptForInner.shape(size,size);
  RHS_exceptForInner.shape(size,1);

  // initialize everything to zero
  JacobianMatrix_exceptForInner.putScalar(0.0);
  RHS_exceptForInner.putScalar(0.0);
}

//-----------------------------------------------------------------------------
//
// This is the problem we are solving:
//
//
//           Inner 
//     1     problem    2       R        3
//     *----/\/\/\/\----*----/\/\/\/\----*
//     |                                 |
//     O vconnect000                     O  vconnect001
//     |                                 |
//    _|_                               _|_
//    \/                                \/
//
//
//  Two independent voltage sources (Vsrcs) connected to ground, an inner "resistor" 
//  represented by the Xyce object, and an outer resistor that is handled as part
//  of the outer problem.  Everything is in series.
//
//  The variables to be solved are:
//
//  v1  voltage at node 1
//  v2  voltage at node 2
//  v3  voltage at node 3
//  i0  current thru vconnect0000
//  i1  current thru vconnect0001
//
//  The "conductance" provided by the inner problem will be a 2x2 sub-matrix:
//
//       v1   v2          rhs  (-f)
//  v1  +Gi  -Gi          -Ii
//  v2  -Gi  +Gi          +Ii
//
//  The "conductance" provided by the R resistor will likewise be a 2x2 sub-matrix:
//
//        v2   v3         rhs  (-f)
//  v2   +Gr  -Gr         -Ir
//  v3   -Gr  +Gr         +Ir
//
//  The matrix stamp (non-conductance) provided by a Vsrc will have the classic "1's pair" 
//  pattern, minus the gounded node entries.  So, for vconnect000:
//
//     v1    i0            rhs (-f)
// v1       +1.0           -Iv0
// i0 +1.0                 V0    <--- voltage value of vconnect000 source
//
// Same for vconnect0001:
//
//     v3    i1            rhs (-f)
// v3       +1.0           -Iv1
// i1 +1.0                 V1    <--- voltage value of vconnect001 source
//
//  So the total Jacobian matrix looks like this:
//
//        v1     v2     v3    i0    i1
//    v1 +Gi    -Gi           1.0
//    v2 -Gi   +Gi+Gr  -Gr
//    v3        -Gr    +Gr          1.0
//    i0  1.0
//    i1                1.0
//
//  And the non-inner parts of the Jacobian matrix look like this:
//
//        v1     v2     v3    i0    i1
//    v1                      1.0
//    v2        +Gr    -Gr
//    v3        -Gr    +Gr          1.0
//    i0  1.0
//    i1                1.0
//
// The matrix is constant, so might as well create a constant matrix for all the 
// parts that don't come from the Xyce inner problem.
// The Xyce inner problem could be linear (and thus have constant matrix 
// entries) or not.  From here we don't know.
inline void topLevelNewton::setupNonInnerMatrix()
{
  // sum these in, device by device:
  // inner problem "contributions", which for now are zero
  JacobianMatrix_exceptForInner(0,0) += 0.0;
  JacobianMatrix_exceptForInner(1,0) += 0.0;
  JacobianMatrix_exceptForInner(0,1) += 0.0;
  JacobianMatrix_exceptForInner(1,1) += 0.0;

  // R resistor contributions:

  JacobianMatrix_exceptForInner(1,1) += +G;
  JacobianMatrix_exceptForInner(2,1) += -G;
  JacobianMatrix_exceptForInner(1,2) += -G;
  JacobianMatrix_exceptForInner(2,2) += +G;

  // vconnect000
  JacobianMatrix_exceptForInner(0,3) += 1.0;
  JacobianMatrix_exceptForInner(3,0) += 1.0;

  // vconnect001
  JacobianMatrix_exceptForInner(2,4) += 1.0;
  JacobianMatrix_exceptForInner(4,2) += 1.0;
}


//-----------------------------------------------------------------------------
inline void topLevelNewton::updateInnerProblemInputs()
{
  voltageInputMap["vconnect0000"] = X(0,0);
  voltageInputMap["vconnect0001"] = X(1,0);
}

//-----------------------------------------------------------------------------
inline void topLevelNewton::updateLinearSystem()
{

  // evaluate the inner Xyce solve to obtain its reduced currents and conductances.
  updateInnerProblemInputs();
  bool bsuccess = simulator_->simulateStep(initJctFlag, voltageInputMap, outputVector, innerJacobianStamp, tlError);
  if (!bsuccess)
  {
    Xyce::dout() << "Inner solve failed" <<std::endl;
    bsuccess = false;
    exit(0);
  }

  // if we get here, that means the "inner" solver has succeeded.

  // global matrix:
  JacobianMatrix.putScalar(0.0);
  JacobianMatrix = JacobianMatrix_exceptForInner; // the rest of the matrix was set up 1x, as it is linear.

  JacobianMatrix(0,0) += innerJacobianStamp[0][0];
  JacobianMatrix(1,0) += innerJacobianStamp[1][0];
  JacobianMatrix(0,1) += innerJacobianStamp[0][1];
  JacobianMatrix(1,1) += innerJacobianStamp[1][1];

  //Xyce::lout() << "Inner problem conductance is: " << innerJacobianStamp[0][0] <<std::endl;

  // rhs:
  // compute the relevant currents from the solution vector.
  v1 = X(0,0);
  v2 = X(1,0);
  v3 = X(2,0);
  Iv0 = X(3,0);
  Iv1 = X(4,0);

  Ir = G*(v2-v3);

  // The RHS needs to contain -f.
  // The output of the inner solve (outputVector) is already -f (rather than f).  So sum those in "as-is"
  // Everything else (from outer problem) multiply by -1.
  RHS.putScalar(0.0);

  // KCL for node 1
  RHS(0,0) += -Iv0+outputVector[0];

  // KCL for node 2
  RHS(1,0) += -Ir+outputVector[1];

  // KCL for node 3
  RHS(2,0) += +Ir-Iv1;

  // voltage drop for vconnect000   = value of source minus the solution variable, since attached to gnd
  RHS(3,0) += (vconnect000-v1);

  // voltage drop for vconnect001   = value of source minus the solution variable, since attached to gnd
  RHS(4,0) += (vconnect001-v3);

}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
inline bool topLevelNewton::checkConvergence()
{
  double dvTol = 1.0e-6;
  std::vector<double> dVvec;
  dVvec.push_back(fabs(dX(0,0)));
  dVvec.push_back(fabs(dX(1,0)));
  dVvec.push_back(fabs(dX(2,0)));

  double max = *std::max_element(dVvec.begin(),dVvec.begin()+3);
  bool converged = (max < dvTol);

  if (converged) 
    Xyce::lout() << "Outer Newton loop converged.  maxNorm = " << max <<std::endl;
  else 
    Xyce::lout() << "Outer Newton loop failed.  maxNorm = " << max <<std::endl;

  return converged;
}

//-----------------------------------------------------------------------------
inline void setupSimDataTRANOP(Xyce::Device::ExternalSimulationData  & extSimData)
{
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
}

//-----------------------------------------------------------------------------
inline void resetSimData(Xyce::Device::ExternalSimulationData  & extSimData)
{
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
}

//-----------------------------------------------------------------------------
inline void setupSimDataDCOP(Xyce::Device::ExternalSimulationData  & extSimData)
{
  // most of the "time" related data isn't used for DC calculations.  So, just
  // set these to nominal values and move on.
  extSimData.is_transient = false;
  extSimData.current_time = 0.0;
  extSimData.final_time = 1.0;

  extSimData.current_time_step_size = 1.0e-5;
  extSimData.previous_time_step_size = 1.0e-5;

  extSimData.time_step_number = 0;

  extSimData.forceOrder=true;
  extSimData.imposedTimeIntegrationOrder=1;

  extSimData.forceBeginningIntegration=true;
  extSimData.imposedBeginningIntegration=true;
}

//-----------------------------------------------------------------------------
inline void resetSimDataDC(Xyce::Device::ExternalSimulationData  & extSimData)
{
  // most of the "time" related data isn't used for DC calculations, so just 
  // set these to equal the output values.
  extSimData.is_transient = false;
  extSimData.current_time = extSimData.nextTime;
  extSimData.final_time = extSimData.finalTime;

  extSimData.current_time_step_size = extSimData.nextTimeStep;
  extSimData.previous_time_step_size = extSimData.currTimeStep;
  extSimData.time_step_number = extSimData.timeStepNumber;

  extSimData.forceOrder=true;
  extSimData.imposedTimeIntegrationOrder= extSimData.currentOrder;

  extSimData.forceBeginningIntegration=true;
  extSimData.imposedBeginningIntegration = extSimData.beginIntegrationFlag;
}

#endif

