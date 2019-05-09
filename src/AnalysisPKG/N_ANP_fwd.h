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
// Purpose        : AC analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Ting Mei   
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_fwd_h
#define Xyce_N_ANP_fwd_h

#include <vector>

namespace Xyce {
namespace Analysis {

  // IMPORTANT:  If you add any modes to this enum, you MUST modify the
  // function "analysisModeName" so that the array of "mode_names" matches!
enum Mode
{
  ANP_MODE_INVALID,
  ANP_MODE_DC_OP,
  ANP_MODE_DC_SWEEP,
  ANP_MODE_DC_NLPOISSON,
  ANP_MODE_TRANSIENT,
  ANP_MODE_MPDE,
  ANP_MODE_HB,
  ANP_MODE_AC,
  ANP_MODE_NOISE,
  ANP_MODE_MOR
};

enum DCOPType {
  NL_POISSON,
  DRIFT_DIFFUSION,
  OFF
};

enum TwoLevelMode {TWO_LEVEL_MODE_TRANSIENT_DCOP = 0, TWO_LEVEL_MODE_TRANSIENT = 1, TWO_LEVEL_MODE_DC_SWEEP = 2};

class AC;
class AnalysisCreatorRegistry;
class AnalysisBase;
class AnalysisManager;
class DCSweep;
class Dakota;
class HB;
class MOR;
class MPDE;
class NOISE;
class NoiseData;
class OutputAdapter;
class OutputMgrAdapter;
class ProcessorBase;
class ProcessorCreatorRegistry;
class Step;
class Sampling;
class SweepParam;
class Transient;

typedef std::vector<SweepParam> SweepVector;

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_fwd_h
