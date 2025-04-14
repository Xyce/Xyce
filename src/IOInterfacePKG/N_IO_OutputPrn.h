//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Class for re-reading Xyce file output, of simulation results,
//                  that is in FORMAT=<STD|CSV|NOINDEX>
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/06/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputPrn_h
#define Xyce_N_IO_OutputPrn_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include<N_IO_OutputFileBase.h>

// ---------- Forward Declarations ----------


namespace Xyce {
namespace IO {

class OutputPrn : public OutputFileBase
{
  public:

    OutputPrn();
    ~OutputPrn();

    // these functions are intended to let Xyce re-read a simulation output
    // file and then recalculate output metrics in .measure() statements
    // without re-running the original simulation.
    bool getOutputVarNames( std::vector< std::string > & varNames );
    int getOutputNextVarValuesParallel( Linear::Vector * varValues );
    int getOutputNextVarValuesSerial( std::vector<double> * varValues );

// These functions will depend on the output format.  Thus,
// they emit errors if the base clase version is called.
    virtual void outputHeader()
    {}

    virtual void outputDC(
      const int dcNumber,
      const int maxDC,
      const std::vector<Analysis::SweepParam> & dcParamVec1,
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr )
    {}

    virtual void outputTran(
      const double & time,
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr )
    {}

    virtual void outputStep(
      const int stepNumber,
      const int maxStep,
      const std::vector<Analysis::SweepParam> & stepParamVec1,
      Linear::Vector * solnVecPtr,
      Linear::Vector * stateVecPtr,
      Linear::Vector * storeVecPtr )
    {}

    virtual void outputAC(
      const double & freq,
      Linear::Vector * freqDomainSolnVecReal,
      Linear::Vector * freqDomainSolnVecImaginary)
    {}

    virtual void outputMPDE(const double & time, Linear::Vector * solnVecPtr )
    {}

    virtual void outputHB(
      const Linear::BlockVector & timeDomainSolnVec,
      const Linear::BlockVector & freqDomainSolnVecReal,
      const Linear::BlockVector & freqDomainSolnVecImaginary)
    {}

    virtual void outputMOR()
    {}

    virtual void finishOutput()
    {}
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputPrn_h
