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
// Purpose        : Base class for handling file output of simulation results.
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/05/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputFileBase_h
#define Xyce_N_IO_OutputFileBase_h

// ----------   Standard Includes   ----------
#include <string>
#include <list>
#include <vector>
#include <ostream>
#include <fstream>

// ----------   Xyce Includes   ----------
#include <N_UTL_Param.h>
#include <N_UTL_OptionBlock.h>
#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_ANP_SweepParam.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {

//-------------------------------------------------------------------------
// Class         : OutputFileBase
// Purpose       : Base class for output of simulation results to a file.
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class OutputFileBase
{
 
public:
  OutputFileBase();
  virtual ~OutputFileBase();
  
  // these functions have basic implementations that 
  // should work for most applications.  Thus they don't always 
  // need an override.
  virtual void openFile( std::string basename, std::string simulationSuffix="");
  void setFileSuffix( std::string newSuffix ) {fileSuffix_ = newSuffix;};
  virtual void closeOutput();

  void outputHOMOTOPYHeader( const std::vector<std::string> & paramNames);
  
  void outputHOMOTOPY( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals,
      Linear::Vector * solnVecPtr );
  
  // These functions will depend on the output format.  Thus, 
  // they emit errors if the base clase version is called.
  virtual void outputHeader() = 0;
  
  virtual void outputDC( 
      const int dcNumber, 
      const int maxDC, 
      const std::vector<Analysis::SweepParam> & dcParamVec1,
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr ) = 0;
      
  virtual void outputTran(
      const double & time, 
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr ) = 0;
    
  virtual void outputStep(
      const int stepNumber, 
      const int maxStep, 
      const std::vector<Analysis::SweepParam> & stepParamVec1,
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr ) = 0;
      
  virtual void outputAC( 
      const double & freq,
      Linear::Vector * freqDomainSolnVecReal,
      Linear::Vector * freqDomainSolnVecImaginary) = 0;
      
  virtual void outputMPDE(const double & time, Linear::Vector * solnVecPtr ) = 0;
  
  virtual void outputHB(
      const Linear::BlockVector & timeDomainSolnVec, 
      const Linear::BlockVector & freqDomainSolnVecReal,
      const Linear::BlockVector & freqDomainSolnVecImaginary) = 0;
  
  virtual void outputMOR() = 0;
  
  virtual void finishOutput() = 0;
  
  // these functions are intended to let Xyce re-read a simulation output 
  // file and then recalculate output metrics in .measure() statements 
  // without re-running the original simulation. 
  virtual bool openFileForRead( std::string filename );
  virtual bool getOutputVarNames( std::vector< std::string > & varNames );
  virtual int getOutputNextVarValuesParallel( Linear::Vector * varValues );
  virtual int getOutputNextVarValuesSerial( std::vector<double> * varValues );
  virtual bool closeFileForRead();
  
  void convertOutputNamesToSolVarNames( std::vector< std::string > & varNames,
                                        std::vector<int> & solSymIdx, 
                                        std::vector<int> & brVarSymIdx );
  
protected:  
  std::ostream * ostreamPtr_;           // the output stream to use.
  std::string outputFileBaseName_;      // the string name of the output file 
  std::string fileSuffix_;              // filename suffix if needed.
  std::string simulationSuffix_;        // any string that should be added as a suffix on the filename before the fileSuffix_ 
                                      
  std::string fullFileName_;            // should be  outputFileBaseName_ + simulationSuffixNumber_ + fileSuffix_
  
  std::string fileFormatName_;       // This is the name of the derived class file format. 
                                // Any virtual functions not overridden by the derived 
                                // class can make an error citing what wasn't supported.
                                // as in "Probe formated output does not currently support MPDE output."
  
  bool appendOutputFlag_;       // default to false, but if true, then append new data to output file.
  
  // for re-reading existing files
  std::ifstream * istreamPtr_;
  std::string inputFileBaseName_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputFileBase_h
