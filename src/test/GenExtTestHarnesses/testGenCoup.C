//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electronic Simulator.
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
// Purpose        : Test harness for General External Device coupling
//                  feature
//
// Special Notes  : Used by the GenCoupAPI test cases
//
// Creator        : Tom Russo
//
// Creation Date  : 13 March 2017
//
//-------------------------------------------------------------------------
// See comments at top of "main" function for details.

#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_DEV_VectorComputeInterface.h>
#include <N_IO_ExtOutInterface.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <Sacado_No_Kokkos.hpp>

/// Implement a simple resistor
///
/// This class implements a vector loader for a simple resistor device.
/// It does not specify a jacobian stamp, and so the device will use
/// a dense stamp by default
///
/// One can specify the resistance in the netlist using
///
///DPARAM={NAME=R VALUE=value}
///
///on the instance line.
///
class resistorVectorCompute : public Xyce::Device::VectorComputeInterface
{
public:
  resistorVectorCompute(double R) : R_(R) {};

  virtual bool computeXyceVectors(std::vector<double> & sV, double t,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdX,
                          std::vector<std::vector<double> > & dQdX)
  {
    int numVars=sV.size();

    F.resize(numVars);
    Q.clear();
    B.clear();
    dQdX.clear();
    dFdX.resize(numVars);
    for (int i=0; i<numVars; i++)
      dFdX[i].resize(numVars);

    std::vector<Sacado::Fad::DFad<double> > indepVars,Fcontribs;
    indepVars.resize(numVars);
    Fcontribs.resize(numVars);
    for (int i=0;i<numVars; i++)
    {
      indepVars[i]=sV[i];
      indepVars[i].diff(i,numVars);
    }

    Fcontribs[0] = (indepVars[0]-indepVars[1])/R_;
    Fcontribs[1] = -Fcontribs[0];

    for (int i=2; i<Fcontribs.size();i++)
      Fcontribs[i]=0;

    for (int i=0; i<numVars; i++)
    {
      F[i] = Fcontribs[i].val();
      for (int j=0; j<numVars; j++)
        dFdX[i][j] = Fcontribs[i].dx(j);
    }
    return true;
  };

  void processDoubleParams(std::vector<std::string> &pNames,
                           std::vector<double> &pValues)
  {
    for (int i = 0; i<pNames.size(); i++)
    {
      if (pNames[i] == "R")
      {
        R_ = pValues[i];
      }
      else
      {
        std::cerr << "Invalid parameter name " << pNames[i]
                  << " passed to resistorVectorCompute" << std::endl;
      }
    }
  };

private:
  double R_;

};

/// Implement a simple three-terminal device
///
/// This class exists so that we have one example of a 3-terminal device
/// which can have non-trivial lead currents at each port.
/// We don't bother with a jacobian stamp and let Xyce use a fully dense
/// stamp, even though we could use a more sparse stamp.
///
/// One can specify the resistances in the netlist using
///
/// DPARAM={NAME=R1,R2 VALUE=<R1value>,<R2value>}
///
/// on the instance line
class threeTerminalVectorCompute : public Xyce::Device::VectorComputeInterface
{
public:
  threeTerminalVectorCompute(double R1, double R2) : R1_(R1), R2_(R2) {};

  virtual bool computeXyceVectors(std::vector<double> & sV, double t,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdX,
                          std::vector<std::vector<double> > & dQdX)
  {
    int numVars=sV.size();

    F.resize(numVars);
    Q.clear();
    B.clear();
    dQdX.clear();
    dFdX.resize(numVars);
    for (int i=0; i<numVars; i++)
      dFdX[i].resize(numVars);

    std::vector<Sacado::Fad::DFad<double> > indepVars,Fcontribs;
    indepVars.resize(numVars);
    Fcontribs.resize(numVars);
    for (int i=0;i<numVars; i++)
    {
      indepVars[i]=sV[i];
      indepVars[i].diff(i,numVars);
    }

    Fcontribs[0] = (indepVars[0]-indepVars[1])/R1_
      + (indepVars[0]-indepVars[2])/R2_;
    Fcontribs[1] = -(indepVars[0]-indepVars[1])/R1_;
    Fcontribs[2] = -(indepVars[0]-indepVars[2])/R2_;

    for (int i=3; i<Fcontribs.size();i++)
      Fcontribs[i]=0;

    for (int i=0; i<numVars; i++)
    {
      F[i] = Fcontribs[i].val();
      for (int j=0; j<numVars; j++)
        dFdX[i][j] = Fcontribs[i].dx(j);
    }
    return true;
  };

  void processDoubleParams(std::vector<std::string> &pNames,
                           std::vector<double> &pValues)
  {
    for (int i = 0; i<pNames.size(); i++)
    {
      if (pNames[i] == "R1")
      {
        R1_ = pValues[i];
      }
      else if (pNames[i] == "R2")
      {
        R2_ = pValues[i];
      }
      else
      {
        std::cerr << "Invalid parameter name " << pNames[i]
                  << " passed to threeTerminalVectorCompute" << std::endl;
      }
    }
  };

private:
  double R1_, R2_;

};

/// implement a series RLC device with a sparse jacobian stamp
///
/// This vector loader demonstrates how to implement a device that
/// has a sparse jacobian stamp.
///
/// This device also has an inductor branch current, and an internal
/// node.  Thus it demonstrates how to use a device that has extra
/// internal variables (over and above the external connections that
/// are obvious from the netlist).
///
/// One can specify the resistance in the netlist using
///
/// DPARAM={NAME=R,L,C VALUE=<Rvalue>,<Lvalue>,<CValue>}
///
/// on the instance line

class rlcVectorCompute : public Xyce::Device::VectorComputeInterface
{
public:
  rlcVectorCompute(double R, double L, double C)
    : R_(R),
      L_(L),
      C_(C)
  {
    // set up our jacobian stamp
    // The jacobian stam establishes the pattern of nonzero entries
    // that this device contributes to the jacobian.
    // in this example we're just keeping jacStamp around as a public
    // variable that can be pulled out of the object.  It could have been done
    // in other ways.
    // There are four variables, the two external and two internal:
    jacStamp.resize(4);

    // The first variable is the positive external node.  The only
    // contribution to this is from the branch current, so we only
    // have one entry on this row
    jacStamp[node1].resize(1);
    jacStamp[node1][0] = branch;

    // The second variable is the negative external node.  The equation
    // for this variable depends on the node 2 value and also the internal
    // node value.  We set the sparsity pattern accordingly.
    jacStamp[node2].resize(2);
    jacStamp[node2][0]=node2;
    jacStamp[node2][1]=nodeInt;

    // The third variable is the internal node.  The equation for this
    // variable depends on the node 2 value and the internal node
    // value, as well as the branch current.  We set the sparsity
    // pattern accordingly.
    jacStamp[nodeInt].resize(3);
    jacStamp[nodeInt][0]=node2;
    jacStamp[nodeInt][1]=nodeInt;
    jacStamp[nodeInt][2]=branch;

    // The fourth variable is the branch current.  The equation for this
    // variable depends on the node 1 value and the internal node
    // value, as well as the branch current.  We set the sparsity
    // pattern accordingly.
    jacStamp[branch].resize(3);
    jacStamp[branch][0]=node1;
    jacStamp[branch][1]=nodeInt;
    jacStamp[branch][2]=branch;
  };

  virtual bool computeXyceVectors(std::vector<double> & sV, double t,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdX,
                          std::vector<std::vector<double> > & dQdX)
  {
    int numVars=sV.size();

    F.resize(numVars);
    Q.resize(numVars);
    B.clear();
    dFdX.resize(numVars);
    dQdX.resize(numVars);
    for (int i=0; i<numVars; i++)
    {
      dFdX[i].resize(numVars);
      dQdX[i].resize(numVars);
    }

    std::vector<Sacado::Fad::DFad<double> > indepVars,Fcontribs,Qcontribs;
    indepVars.resize(numVars);
    Fcontribs.resize(numVars);
    Qcontribs.resize(numVars);
    for (int i=0;i<numVars; i++)
    {
      indepVars[i]=sV[i];
      indepVars[i].diff(i,numVars);
      Fcontribs[i]=0;
      Qcontribs[i]=0;
    }

    // vars 0 and 1 are the external nodes
    // var 2 is the internal node
    // var 3 is the branch current
    // equation 0 and 1 are the equations for the external node
    // equation 2 is the equation for the internal node
    // equation 3 is the branch equation

    // Do the resistor and inductor as a single branch with just one
    // internal node (c.f rlc2.va in Xyce/utils/ADMS/examples/toys)
    // Branch equation is then
    //    (branch_current*R + L d(branch_current)/dt)-(v0-v2) = 0
    // Current between external node 0 and internal node 2 is just
    // the branch current
    Fcontribs[node1] = indepVars[branch];
    Fcontribs[nodeInt] = -indepVars[branch];

    Fcontribs[branch] = indepVars[branch]*R_
      - (indepVars[node1]-indepVars[nodeInt]);
    Qcontribs[branch] = L_*indepVars[branch];


    // Now the capacitor.  Current flowing from node 2 to node 1 is
    // just dQ/dt where Q=CV.
    Sacado::Fad::DFad<double> capCharge
      = C_*(indepVars[nodeInt]-indepVars[node2]);
    Qcontribs[nodeInt] = capCharge;
    Qcontribs[node2] = -capCharge;

    for (int i=0; i<numVars; i++)
    {
      F[i] = Fcontribs[i].val();
      Q[i] = Qcontribs[i].val();
      for (int j=0; j<numVars; j++)
      {
        dFdX[i][j] = Fcontribs[i].dx(j);
        dQdX[i][j] = Qcontribs[i].dx(j);
      }
    }
    return true;
  };

  void processDoubleParams(std::vector<std::string> &pNames,
                           std::vector<double> &pValues)
  {
    for (int i = 0; i<pNames.size(); i++)
    {
      if (pNames[i] == "R")
      {
        R_ = pValues[i];
      }
      else if (pNames[i] == "L")
      {
        L_ = pValues[i];
      }
      else if (pNames[i] == "C")
      {
        C_ = pValues[i];
      }
      else
      {
        std::cerr << "Invalid parameter name " << pNames[i]
                  << " passed to rlcVectorCompute" << std::endl;
      }
    }
  };

public:
  std::vector< std::vector<int> > jacStamp;
private:
  double R_;
  double L_;
  double C_;
  static const int node1=0;
  static const int node2=1;
  static const int nodeInt=2;
  static const int branch=3;
};

/// Implement a voltage-controlled current source
///
/// nodes 3 and 4 are the control nodes
/// current flows from positive (node 1) to negative (node2)
/// Illustrates how one could make a device that pulls in extra data
/// from the circuit without actually having a square jacobian stamp
///
/// One can set the VCCS's transconductance from the netlist using:
///
/// DPARAM={NAME=TRANSCONDUCTANCE VALUE=\<value\> }
///
/// on the instance line
class vccsVectorCompute : public Xyce::Device::VectorComputeInterface
{
public:
  vccsVectorCompute(double tC)
    : transConductance_(tC)
  {
    // set up our jacobian stamp
    // in this example we're just keeping jacStamp around as a public
    // variable that can be pulled out of the object.  It could have been done
    // in other ways.
    jacStamp.resize(4);

    jacStamp[nodePos].resize(2);
    jacStamp[nodePos][0] = controlNodePos;
    jacStamp[nodePos][1] = controlNodeNeg;
    jacStamp[nodeNeg].resize(2);
    jacStamp[nodeNeg][0] = controlNodePos;
    jacStamp[nodeNeg][1] = controlNodeNeg;

    // No loading happens to nodes 2 and 3, so those jacStamp rows are
    // zero length.
  };

  virtual bool computeXyceVectors(std::vector<double> & sV, double t,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdX,
                          std::vector<std::vector<double> > & dQdX)
  {
    int numVars=sV.size();

    std::cout << " In vccs computeXyceVectors, number vars = " << numVars
              << std::endl;

    F.resize(numVars);
    Q.clear();
    B.clear();
    dFdX.resize(numVars);
    dQdX.clear();
    for (int i=0; i<numVars; i++)
    {
      dFdX[i].resize(numVars);
    }

    std::vector<Sacado::Fad::DFad<double> > indepVars,Fcontribs;
    indepVars.resize(numVars);
    Fcontribs.resize(numVars);
    for (int i=0;i<numVars; i++)
    {
      indepVars[i]=sV[i];
      indepVars[i].diff(i,numVars);
      Fcontribs[i]=0;
    }


    Sacado::Fad::DFad<double> current
      = transConductance_*(indepVars[controlNodePos]-indepVars[controlNodeNeg]);
    Fcontribs[nodePos] = current;
    Fcontribs[nodeNeg] = -current;

    for (int i=0; i<numVars; i++)
    {
      F[i] = Fcontribs[i].val();
      for (int j=0; j<numVars; j++)
      {
        dFdX[i][j] = Fcontribs[i].dx(j);
      }
    }
    return true;
  };
  void processDoubleParams(std::vector<std::string> &pNames,
                           std::vector<double> &pValues)
  {
    for (int i = 0; i<pNames.size(); i++)
    {
      if (pNames[i] == "TRANSCONDUCTANCE")
      {
        transConductance_ = pValues[i];
      }
      else
      {
        std::cerr << "Invalid parameter name " << pNames[i]
                  << " passed to vccsVectorCompute" << std::endl;
      }
    }
  };

public:
  std::vector< std::vector<int> > jacStamp;
private:
  double transConductance_;

  static const int nodePos=0;
  static const int nodeNeg=1;
  static const int controlNodePos=2;
  static const int controlNodeNeg=3;
};

/// Implement a basic capacitor
///
/// This class serves to demonstrate how a device can be implemented
/// with both time and frequency domain loaders.
///
/// In time domain analyses (and small-signal AC) the
/// computeXyceVectors method will be called.  The
/// computeXyceFDVectors method will be called instead when doing
/// harmonic balance simulation.
///
/// One can set the capacitance using
///
/// DPARAM={NAME=C VALUE=\<value\> }
///
/// on the instance line
class capacitorVectorCompute : public Xyce::Device::VectorComputeInterface
{
public:
  capacitorVectorCompute(double C) : C_(C) {};

  virtual bool computeXyceVectors(std::vector<double> & sV, double t,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdX,
                          std::vector<std::vector<double> > & dQdX)
  {
    int numVars=sV.size();

    Q.resize(numVars);
    F.clear();
    B.clear();
    dFdX.clear();
    dQdX.resize(numVars);
    for (int i=0; i<numVars; i++)
      dQdX[i].resize(numVars);

    std::vector<Sacado::Fad::DFad<double> > indepVars,Qcontribs;
    indepVars.resize(numVars);
    Qcontribs.resize(numVars);
    for (int i=0;i<numVars; i++)
    {
      indepVars[i]=sV[i];
      indepVars[i].diff(i,numVars);
    }

    Qcontribs[0] = (indepVars[0]-indepVars[1])*C_;
    Qcontribs[1] = -Qcontribs[0];

    for (int i=2; i<Qcontribs.size();i++)
      Qcontribs[i]=0;

    for (int i=0; i<numVars; i++)
    {
      Q[i] = Qcontribs[i].val();
      for (int j=0; j<numVars; j++)
        dQdX[i][j] = Qcontribs[i].dx(j);
    }
    return true;
  };

  virtual bool computeXyceFDVectors(std::vector<std::complex<double> > &sV,
                                    double frequency,
                                    std::vector<std::complex<double> > &F,
                                    std::vector<std::complex<double> > &B,
                                    std::vector<std::vector<std::complex<double> > > &dFdX)
  {
    int numVars = sV.size();

    if (numVars != 2)
      return false;

    B.clear();
    F.resize(numVars);
    dFdX.resize(numVars);
    double pi = 3.14159265358979;

    for ( int i=0; i<numVars; i++)
    {
      dFdX[i].resize(numVars);
    }

    dFdX[0][0] = std::complex<double> ( 0, frequency*2*pi*C_);
    dFdX[0][1] = -dFdX[0][0];
    dFdX[1][0] = -dFdX[0][0];
    dFdX[1][1] = dFdX[0][0];

    F[0] = dFdX[0][0]*sV[0]+dFdX[0][1]*sV[1];
    F[1] = dFdX[1][0]*sV[0]+dFdX[1][1]*sV[1];

    return true;
  };

  virtual bool haveFDLoads() { return true; };

  void processDoubleParams(std::vector<std::string> &pNames,
                           std::vector<double> &pValues)
  {
    for (int i = 0; i<pNames.size(); i++)
    {
      if (pNames[i] == "C")
      {
        C_ = pValues[i];
      }
      else
      {
        std::cerr << "Invalid parameter name " << pNames[i]
                  << " passed to capacitorVectorCompute" << std::endl;
      }
    }
  };

private:
  double C_;

};

/// Implement the ExternalOutputInterface class
///
/// This class is used to test the ExternalOutputInterface mechanmism.
/// It is intended to mimic what Xyce would do with the data it receives,
/// creating a columnar output file like Xyce's "STD" format.
///
/// We implement the minimum subset of the ExternalOutputInterface
/// class, augmented a little with functions and member variables to suit
/// our current purposes.
class ioTestInterface : public Xyce::IO::ExternalOutputInterface
{
public:
  ioTestInterface(std::string filename,Xyce::IO::OutputType::OutputType theOutputType=Xyce::IO::OutputType::TRAN)
    :
    theOutputType_(theOutputType)
  {
    os_ =new std::ofstream(filename.c_str(),std::ios_base::out);
  };

  ~ioTestInterface()
  {
    delete os_;
  };

  Xyce::IO::OutputType::OutputType getOutputType()
  {
    return theOutputType_;
  }

  void addOutputString(std::string theString)
  {
    outputStrings_.push_back(theString);
  };

  void requestedOutputs(std::vector<std::string> &outputVars)
  {
    outputVars = outputStrings_;
  };

  void reportParseStatus(std::vector<bool> & statusVec)
  {
    for (int i=0;i<statusVec.size();i++)
    {
      std::cout << "String " << outputStrings_[i] << " is "
                <<  ( (statusVec[i])? " parseable" : " unparseable")
                << std::endl;
    }
  };

  void outputFieldNames(std::vector<std::string> & outputNames)
  {
    // The current "external output" feature doesn't do Xyce's "Expand
    // complex types" manipulation that some output formats do.  So
    // for frequency domain output, we *ALWAYS* get a single column
    // name, and always get complex data for each column.  Some of
    // these are nonsense, especially index and frequency, for which
    // we should only output the real part.  There are others like VR, VI, VDB
    // VM, VP and so on, but we are NOT handling those here.  This is a
    // demo code, not a fully functional standard print formatter.

    hideImag_.resize(outputNames.size(),false);
    intField_.resize(outputNames.size(),false);
    for (int i=0;i<outputNames.size();i++)
    {
      if (outputNames[i] == "INDEX")
      {
        *os_ << "Index" << "    ";
        hideImag_[i]=true;
        intField_[i]=true;
      }
      else if ( outputNames[i] == "FREQ")
      {
        *os_ << "FREQ" << "    ";
        hideImag_[i]=true;
      }
      else
      {
        // really need to test if any Frequency Domain output, but so far
        // AC is the only one implemented, so this will do for now.
        if (getOutputType() == Xyce::IO::OutputType::AC
            || getOutputType() == Xyce::IO::OutputType::HB_FD )
        {
          *os_ << "Re(" << outputNames[i] << ")    ";
          *os_ << "Im(" << outputNames[i] << ")    ";
        }
        else
        {
          *os_ << outputNames[i] << "    ";
        }
      }
    }
    *os_ << std::endl;
  };

  void outputReal(std::vector<double> & outputData)
  {
    for (int i=0;i<outputData.size();i++)
    {
      if (intField_[i])
      {
        *os_ <<  std::resetiosflags(std::ios_base::floatfield)
             << std::setiosflags(std::ios_base::fixed)
             << std::resetiosflags(std::ios_base::adjustfield)
             << std::setiosflags(std::ios_base::left)
             << std::setprecision(0) << std::setw(5);
      }
      else
      {
        *os_ <<  std::resetiosflags(std::ios_base::floatfield)
             << std::setiosflags(std::ios_base::scientific)
             << std::resetiosflags(std::ios_base::adjustfield)
             << std::setiosflags(std::ios_base::right)
             << std::setprecision(8) << std::setw(17);
      }
      *os_ << outputData[i] << "    ";
    }
    *os_ << std::endl;
  };

  void outputComplex(std::vector<std::complex<double> > & outputData)
  {
    for (int i=0;i<outputData.size();i++)
    {
      if (intField_[i])
      {
        *os_ <<  std::resetiosflags(std::ios_base::floatfield)
            << std::setiosflags(std::ios_base::fixed)
            << std::resetiosflags(std::ios_base::adjustfield)
            << std::setiosflags(std::ios_base::left)
            << std::setprecision(0) << std::setw(5);
      }
      else
      {
        *os_ <<  std::resetiosflags(std::ios_base::floatfield)
            << std::setiosflags(std::ios_base::scientific)
            << std::resetiosflags(std::ios_base::adjustfield)
            << std::setiosflags(std::ios_base::right)
            << std::setprecision(8) << std::setw(17);
      }
      *os_ << outputData[i].real() << "    ";
      if (!hideImag_[i])
        *os_ << outputData[i].imag() << "    " ;
    }
    *os_ << std::endl;
  };

  void finishOutput()
  {
    *os_ << "End of Xyce(TM) Simulation" << std::endl;
  };

private:
  std::vector<std::string> outputStrings_;
  std::ostream *os_;
  Xyce::IO::OutputType::OutputType theOutputType_;
  std::vector<bool>  hideImag_;
  std::vector<bool>  intField_;
};

/// Demo program for General External Coupling feature
///
/// This demo program, when linked with Xyce's shared libraries, will enable not
/// only running normal Xyce netlists
/// but also those that contain "General External" devices.
///
/// It only recognizes the General External devices if they're given specific
/// names, because we are not implenting a clever name scanner here.
///
/// The program will recognize these devices in the netlist based on their names.
///
/// YGENEXT R1 <node1> <node2>
///
/// YGENEXT C1 <node1> <node2>
///
/// YGENEXT RLC1 <node1> <node2>
///
/// YGENEXT VCCS1 <node+> <node-> <control +> <control -> <modelname>
///
/// Double-precision parameters may be added to the instance lines
/// in the netlist using the "DPARAMS" described in the individual class
/// documentation.
///
/// @note  A parser limitation of Xyce requires that a model be specified
/// any time any more than 2 nodes are given to the YGENEXT device.
/// This is true even though the YGenExt device has no model parameters.
/// A null model card like:
/// .model dummy genext ()
/// is sufficient to get past this limitation
///
/// YGENEXT threeterm1 <node1> <node2> <node3> <modelname>
///
/// You can only have at most one of each of these in your netlist.
///
/// Usage:
/// testGenCoup  <netlistname>
///
/// This version demonstrates a very simple resistor device, which depends on
/// the default behaviors of the GenExt device when neither numIntVars nor
/// jacobian stamp are set by the caller.
///
/// The RLC device has two internal variables and sets a sparse jacobian stamp.
int main(int argc, char **argv)
{
  bool debug=false;

  Xyce::Circuit::Simulator::RunStatus run_status=Xyce::Circuit::Simulator::SUCCESS;

  // Here is where we create the simulator object we will use to run Xyce
  Xyce::Circuit::GenCouplingSimulator xyce;

  bool iotest=false;
  ioTestInterface * theIOTestInterface=0;
  ioTestInterface * theIOTestInterface2=0;

  // Here we check for two special command line options (only one of which
  // will ever be respected) that will be used to turn on special output
  // capability testing.  If either is present as the first argument, we have
  // to set a flag, then readjust the args list to remove it.

  // iotest1 is used with the rlc_series_vccs test case in the GenCoupAPI test
  // suite.  iotest2 is used with the genextLeadCurr test.  iotest3 is used with the
  // genextAC test, and iotest4 with the genextDC test.  iotest5 will
  // be used with a harmonic balance test.
  //
  // Each is intended to show that it is possible to communicate the same data
  // through the external output interface as would be output by normal
  // .print lines.
  if ( argc > 1)
  {
    std::string arg1(argv[1]);
    if (arg1 == "-iotest1")
    {
      iotest=true;
      theIOTestInterface = new ioTestInterface("ioTest1.out");
      theIOTestInterface->addOutputString("Index");
      theIOTestInterface->addOutputString("TIME");
      theIOTestInterface->addOutputString("V(1)");
      theIOTestInterface->addOutputString("v(1a)");
      theIOTestInterface->addOutputString("i(v1)");
      theIOTestInterface->addOutputString("v(2)");
      theIOTestInterface2 = new ioTestInterface("ioTest1a.out");
      theIOTestInterface2->addOutputString("Index");
      theIOTestInterface2->addOutputString("TIME");
      theIOTestInterface2->addOutputString("V(*)");
      theIOTestInterface2->addOutputString("I(*)");
    }
    else if (arg1 == "-iotest2")
    {
      iotest=true;
      theIOTestInterface = new ioTestInterface("ioTest2.out");
      theIOTestInterface->addOutputString("Index");
      theIOTestInterface->addOutputString("TIME");
      theIOTestInterface->addOutputString("{I(vprobe1)-I1(YGenExt!R1)}");
      theIOTestInterface->addOutputString("{I(vprobe1)-I1(YGenExt!RLC1)}");
      theIOTestInterface->addOutputString("{I(vprobe1)-I1(YGenExt!ThreeTerm1)}");
      theIOTestInterface->addOutputString("{I(vprobe1)+I2(YGenExt!ThreeTerm1)+I3(YGenExt!ThreeTerm1)}");
      theIOTestInterface->addOutputString("{I1(YGenExt!ThreeTerm1)+I2(YGenExt!ThreeTerm1)+I3(YGenExt!ThreeTerm1)}");
    }
    else if (arg1 == "-iotest3")
    {
      iotest=true;
      theIOTestInterface = new ioTestInterface("ioTest3.out",
                                               Xyce::IO::OutputType::AC);
      theIOTestInterface->addOutputString("Index");
      theIOTestInterface->addOutputString("FREQ");
      theIOTestInterface->addOutputString("v(1)");
      theIOTestInterface->addOutputString("v(1a)");
      theIOTestInterface->addOutputString("I(v1)");
      theIOTestInterface->addOutputString("V(2)");
    }
    else if (arg1 == "-iotest4")
    {
      iotest = true;
      theIOTestInterface = new ioTestInterface("ioTest4.out",
                                               Xyce::IO::OutputType::DC);
      theIOTestInterface->addOutputString("Index");
      theIOTestInterface->addOutputString("v(1)");
      theIOTestInterface->addOutputString("I(v1)");
    }
    else if (arg1 == "-iotest5")
    {
      iotest = true;
      theIOTestInterface = new ioTestInterface("ioTest5.FD.out",
                                               Xyce::IO::OutputType::HB_FD);
      theIOTestInterface->addOutputString("Index");
      theIOTestInterface->addOutputString("FREQ");
      theIOTestInterface->addOutputString("v(1)");
      theIOTestInterface->addOutputString("v(2)");
      theIOTestInterface->addOutputString("I(v1)");
      theIOTestInterface2 = new ioTestInterface("ioTest5.TD.out",
                                               Xyce::IO::OutputType::HB_TD);
      theIOTestInterface2->addOutputString("Index");
      theIOTestInterface2->addOutputString("TIME");
      theIOTestInterface2->addOutputString("v(1)");
      theIOTestInterface2->addOutputString("v(2)");
      theIOTestInterface2->addOutputString("I(v1)");
    }

    if (iotest)
    {
      for (int argnum=2;argnum<argc;argnum++)
        argv[argnum-1]=argv[argnum];
      argc--;
    }
  }

  // Now we pass our command line arguments into Xyce's first initialization step
  run_status=xyce.initializeEarly(argc,argv);
  if (run_status==Xyce::Circuit::Simulator::ERROR) exit(1);
  if (run_status==Xyce::Circuit::Simulator::DONE) exit(0);

  // If we have external output to do, add the interfaces
  if (theIOTestInterface)
  {
    xyce.addOutputInterface(theIOTestInterface);
  }
  if (theIOTestInterface2)
  {
    xyce.addOutputInterface(theIOTestInterface2);
  }


  // We now must query Xyce for any general external devices present in the
  // netlist we just parsed.
  std::vector<std::string> deviceNames;
  bool bsuccess=xyce.getDeviceNames("YGENEXT",deviceNames);
  if (!bsuccess && debug)
  {
    std::cout << " No external devices found.. regular xyce netlist." << std::endl;
  }

  // Create instances of the vector loader classes we have defined
  resistorVectorCompute rvc(1);
  capacitorVectorCompute cvc(1);
  threeTerminalVectorCompute threeTermvc(1,2);
  vccsVectorCompute vccsvc(0);
  rlcVectorCompute rlcvc(1,1,1);


  // Now scan our list of general external devices, and find the ones we recognize.
  // For each recognized device, associate one of our vector loaders and set any
  // other special parameters
  for (int i =0; i<deviceNames.size(); i++)
  {
    std::cout << "GenExt Device " << i << " name is " << deviceNames[i] << std::endl;
    if (deviceNames[i] == "YGENEXT!R1")
    {
      // We found a resistor
      // Retrieve any double precision parameters that may be on the instance line
      std::vector<std::string> pNames;
      std::vector<double> pValues;
      bsuccess=xyce.getDParams(deviceNames[i],pNames,pValues);
      if (!bsuccess)
      {
        std::cerr << " getDParams failed on " << deviceNames[i] << std::endl;
        exit(1);
      }
      // pass those parameters to our object
      rvc.processDoubleParams(pNames,pValues);

      // Now tell Xyce to use the vector loader
      bsuccess=xyce.setVectorLoader(deviceNames[i],&rvc);
      if (!bsuccess)
      {
        std::cerr << " Failed to set vector loader for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }
    }
    if (deviceNames[i] == "YGENEXT!C1")
    {
      // We found a capacitor
      // Retrieve any double precision parameters that may be on the instance line
      std::vector<std::string> pNames;
      std::vector<double> pValues;
      bsuccess=xyce.getDParams(deviceNames[i],pNames,pValues);
      if (!bsuccess)
      {
        std::cerr << " getDParams failed on " << deviceNames[i] << std::endl;
        exit(1);
      }
      // pass those parameters to our object
      cvc.processDoubleParams(pNames,pValues);

      // Now tell Xyce to use the vector loader
      bsuccess=xyce.setVectorLoader(deviceNames[i],&cvc);
      if (!bsuccess)
      {
        std::cerr << " Failed to set vector loader for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }
    }
    if (deviceNames[i] == "YGENEXT!THREETERM1")
    {
      // We found a three-terminal device
      // Retrieve any double precision parameters that may be on the instance line
      std::vector<std::string> pNames;
      std::vector<double> pValues;
      bsuccess=xyce.getDParams(deviceNames[i],pNames,pValues);
      if (!bsuccess)
      {
        std::cerr << " getDParams failed on " << deviceNames[i] << std::endl;
        exit(1);
      }
      // pass those parameters to our object
      threeTermvc.processDoubleParams(pNames,pValues);


      // Now tell Xyce to use the vector loader
      bsuccess=xyce.setVectorLoader(deviceNames[i],&threeTermvc);
      if (!bsuccess)
      {
        std::cerr << " Failed to set vector loader for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }
    }

    if (deviceNames[i] == "YGENEXT!RLC1")
    {
      // We found a series RLC device
      // Retrieve any double precision parameters that may be on the instance line
      std::vector<std::string> pNames;
      std::vector<double> pValues;
      bsuccess=xyce.getDParams(deviceNames[i],pNames,pValues);
      if (!bsuccess)
      {
        std::cerr << " getDParams failed on " << deviceNames[i] << std::endl;
        exit(1);
      }
      // pass those parameters to our object
      rlcvc.processDoubleParams(pNames,pValues);

      // This device has two internal variables (an internal node and a branch current)
      bsuccess=xyce.setNumInternalVars(deviceNames[i],2);
      if (!bsuccess)
      {
        std::cerr << " Failed to set internal vars for" << deviceNames[i]
                  << std::endl;
        exit(1);
      }

      // We have implemented this device so it can use a sparse jacobian stamp.  Tell
      // Xyce what that stamp is.
      bsuccess=xyce.setJacStamp(deviceNames[i],rlcvc.jacStamp);
      if (!bsuccess)
      {
        std::cerr << " Failed to set Jacobian Stamp for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }

      // Now associate our vector loader with the device
      bsuccess=xyce.setVectorLoader(deviceNames[i],&rlcvc);
      if (!bsuccess)
      {
        std::cerr << " Failed to set vector loader for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }
    }

    if (deviceNames[i] == "YGENEXT!VCCS1")
    {
      std::vector<std::string> pNames;
      std::vector<double> pValues;
      bsuccess=xyce.getDParams(deviceNames[i],pNames,pValues);
      if (!bsuccess)
      {
        std::cerr << " getDParams failed on " << deviceNames[i] << std::endl;
        exit(1);
      }
      vccsvc.processDoubleParams(pNames,pValues);
      bsuccess=xyce.setJacStamp(deviceNames[i],vccsvc.jacStamp);
      if (!bsuccess)
      {
        std::cerr << " Failed to set Jacobian Stamp for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }

      bsuccess=xyce.setVectorLoader(deviceNames[i],&vccsvc);
      if (!bsuccess)
      {
        std::cerr << " Failed to set vector loader for " << deviceNames[i]
                  << std::endl;
        exit(1);
      }
    }

  }

  // All of our recognized general external devices are completely set up, it is now
  // safe to finish initialization
  run_status=xyce.initializeLate();
  if (run_status==Xyce::Circuit::Simulator::ERROR) exit(1);
  if (run_status==Xyce::Circuit::Simulator::DONE) exit(0);

  // Now we simply let Xyce run the simulation to completion.
  // Our objects will be called as needed throughout the run.
  // A more complex interaction could use the "simulateUntil"
  // method to make Xyce run only for a small time span and
  // return control to this calling program.
  try
  {
    run_status = xyce.runSimulation();
  }
  catch (std::exception &x)
  {
    run_status=Xyce::Circuit::Simulator::ERROR;
  }

  // the "finalize" function closes output files and
  // releases memory.
  if (run_status != Xyce::Circuit::Simulator::DONE)
  {
    try
    {
      xyce.finalize();
    }
    catch (std::exception &x)
    {
      run_status = Xyce::Circuit::Simulator::ERROR;
    }
  }
  exit ((run_status==Xyce::Circuit::Simulator::ERROR)?1:0);
}
