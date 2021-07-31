// Test harness for github issue 29/gitlab-ex issue 243
#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_IO_ExtOutInterface.h>
#include <iostream>
#include <fstream>
#include <iomanip>

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
    // this is less useful than it seems, and only reports ill-formatted
    // requests, not invalid requests for nonexistent quantities.
    // Let's leave it out and simplify our output.
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


int main(int argc, char **argv)
{
  bool debug=false;

  Xyce::Circuit::Simulator::RunStatus run_status=Xyce::Circuit::Simulator::SUCCESS;
  ioTestInterface * theIOTestInterface=0;
  Xyce::Circuit::GenCouplingSimulator xyce;
  run_status=xyce.initializeEarly(argc,argv);
  if (run_status==Xyce::Circuit::Simulator::ERROR) exit(1);
  if (run_status==Xyce::Circuit::Simulator::DONE) exit(0);

  theIOTestInterface = new ioTestInterface("ioTest.out");
  theIOTestInterface->addOutputString("Index");
  theIOTestInterface->addOutputString("TIME");
  theIOTestInterface->addOutputString("V(*)");
  theIOTestInterface->addOutputString("I(*)");

  xyce.addOutputInterface(theIOTestInterface);
  run_status=xyce.initializeLate();
  if (run_status==Xyce::Circuit::Simulator::ERROR) exit(1);
  if (run_status==Xyce::Circuit::Simulator::DONE) exit(0);
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

