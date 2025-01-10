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
// Purpose        : Dakota analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Dakota_h
#define Xyce_N_ANP_Dakota_h

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {
//-------------------------------------------------------------------------
// Class         : Dakota
// Purpose       : Dakota analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Dakota : public AnalysisBase
{
  public:
    Dakota(AnalysisManager &anaManagerPtr, AnalysisBase &anaType );

    virtual ~Dakota( );
  
    virtual bool getDCOPFlag () const {return true;}

protected:
    virtual void finalExpressionBasedSetup() {};
    virtual bool doRun();
    virtual bool doInit() { return true; }
    virtual bool doLoopProcess() { return true; }
    virtual bool doProcessSuccessfulStep() { return true; }
    virtual bool doProcessFailedStep() { return true; }
    virtual bool doFinish() { return true; }
    virtual bool doHandlePredictor() { return true; }

private:
    AnalysisBase &              mainAnalysis_;

  // these should be part of this class, but are still too entangled in
  // the AnalysisManager class to be moved
  //vector <N_TIA_SweepParam> stepParamVec_;
  //
  //bool stepLoopInitialized_;    // true if the step loop has been set up.
  //int stepLoopSize_;
  //int stepLoopIter_;

};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Dakota_h
