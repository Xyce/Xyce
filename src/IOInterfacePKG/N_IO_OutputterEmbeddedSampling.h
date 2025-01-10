//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Outputter class for Embedded Sampling
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 7/26/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterEmbeddedSampling_h
#define Xyce_N_IO_OutputterEmbeddedSampling_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableEmbeddedSamplingOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode);


// Common functions used by the EmbeddedSamplingPrn, EmbeddedSamplingCSV and
// EmbeddedSamplingTecplot classes to output header info and data rows.
void makeEmbeddedSamplingColumnNames(
   const PrintParameters&       printParameters,
   std::vector<std::string>&    colNames,
   bool                         regressionPCEenable,
   bool                         projectionPCEenable,
   int                          numSamples,
   const std::vector<std::string> & regressionPCEcoeffs,
   const std::vector<std::string> & projectionPCEcoeffs,
   const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec);

void outputEmbeddedSamplingData(
   const PrintParameters&       printParameters,
   std::ostream *               os,
   const std::vector<complex>&  result_list,
   bool                         regressionPCEenable,
   bool                         projectionPCEenable,
   int                          numSamples,
   const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec);

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterEmbeddedSampling_h
