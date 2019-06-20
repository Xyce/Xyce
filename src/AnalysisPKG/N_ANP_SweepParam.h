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
// Purpose       : Handles the parameter data associated with sweeps.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 9/4/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_SweepParam_h
#define Xyce_N_ANP_SweepParam_h

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#if __cplusplus>=201103L
// note, this only works with C++11!
#include <random>
#else
#include <N_UTL_RandomNumbers.h>
#endif

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Param.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Class         : SweepParam
//
// Purpose       : This class contains basic parameter data for parameter
//                 sweeps, for a single parameter.  If there are multiple
//                 parameters in the sweep, each one gets a class like
//                 this.
//
// Special Notes : "Step" here refers to steps in a parameter sweep loop,
//                 not time steps or DC sweep steps.
//
//
// Creator       : Eric Keiter, SNL
// Creation Date : 10/31/03
//-----------------------------------------------------------------------------
class SweepParam
{
public:
  // Default constructor
  SweepParam () : 
   name(""),
   type("LIN"),
   startVal(0.0),
   stopVal(0.0),
   stepVal(0.0),
   stepMult(0.0),
   mean(0.0),
   stdDev(0.0),
   alpha(1.0),
   beta(1.0),

   upper_bound(0.0),
   lower_bound(0.0),
   upper_boundGiven(false),
   lower_boundGiven(false),

   currentVal(0.0),
   numSteps(0),
   count(-1),
   maxStep(0),
   interval(1),
   outerStepNumber(0),
   valList(0),
   dataSetName(""),
   sweepResetFlag_(false),
   lastLocalStepNumber_(-1)
   {}

  // Destructor
  ~SweepParam ()
  {}

  bool updateCurrentVal (int stepNumber);
  bool getSweepResetFlag() {return sweepResetFlag_;}

  std::string name;
  std::string type;

  double startVal;
  double stopVal;
  double stepVal;
  double stepMult;

  // for normal distribution
  double mean;
  double stdDev;

  // for gamma distribution
  double alpha;
  double beta;

  double upper_bound;
  double lower_bound;
  bool upper_boundGiven;
  bool lower_boundGiven;

  double currentVal;

  int numSteps;

  int count;
  int maxStep;
  int interval;

  int outerStepNumber;

  std::vector<double> valList;

  std::string dataSetName;

 private:
  bool sweepResetFlag_;
  int lastLocalStepNumber_;
};

typedef std::vector<SweepParam> SweepVector;

//-----------------------------------------------------------------------------
// Function      : SweepParam::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
/// @param os 
/// @param sp 
///
/// @return 
///
std::ostream & operator<<(std::ostream & os, const SweepParam & sp);


//-----------------------------------------------------------------------------
// Function      : parseSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Oct  1 09:19:28 2014
//-----------------------------------------------------------------------------
///
/// Populate the sweep params from the parameter list.
///
/// @invariant
///
/// @param sweep_param  sweep parameters objec t populate
/// @param first        begin iterator of the parameter list
/// @param last         end iterator of the parameter list
///
SweepParam parseSweepParams(Util::ParamList::const_iterator first, Util::ParamList::const_iterator last);

bool updateSweepParams(int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);
bool updateSweepParams(Loader::Loader &loader, int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end, bool overrideOriginal);
int setupSweepLoop(Parallel::Machine comm, Loader::Loader &loader, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);
int setSweepLoopVals(std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);

bool processDataStatements(
    const Util::OptionBlock & paramsBlock,
    std::map< std::string, std::vector<std::string> > & dataNamesMap,
    std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap);

bool convertData(
    SweepVector & stepSweepVector,
    const std::map< std::string, std::vector<std::string> > & dataNamesMap,
    const std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap
    );

bool isDataSpecified(const Util::OptionBlock & paramsBlock);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_SweepParam_h
