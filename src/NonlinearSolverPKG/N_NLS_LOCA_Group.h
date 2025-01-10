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

//-------------------------------------------------------------------------
//
// Purpose        : Interface to Xyce for LOCA groups.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, NLS, 9233
//
// Creation Date  : 02/17/03
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_LOCA_Group_h
#define Xyce_N_NLS_LOCA_Group_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>

// ----------   NOX Includes   ----------

#include "N_NLS_NOX_Group.h"       // base class
#include "LOCA_Abstract_Group.H"   // base class
#include "LOCA_Parameter_Vector.H" // data member
#include "LOCA_DerivUtils.H"       // data member
#include "N_LAS_Vector.h"          // data member
#include "Teuchos_RCP.hpp" // data member

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_LOCA {

//-----------------------------------------------------------------------------
// Class         : N_NLS::LOCA::Group
//
// Purpose       :
//
//      NOX Group Interface for Xyce
//
// Creator       : Roger Pawlowski, SNL, 9233
//
// Creation Date : 2/17/03
//-----------------------------------------------------------------------------

class Group : public N_NLS_NOX::Group, public LOCA::Abstract::Group {

public:

  //! Basic Constructor
  Group(Teuchos::RCP<LOCA::GlobalData> globalData,
        N_NLS_NOX::SharedSystem& s, Xyce::Loader::NonlinearEquationLoader& l, Xyce::IO::OutputMgr& o,
        Xyce::Analysis::AnalysisManager & t);

  //! Copy Constructor
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor
  ~Group();

  //! Assignment Operator
  NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  //! Assignment Operator
  N_NLS_NOX::Group& operator=(const N_NLS_NOX::Group& source);

  //! Assignment Operator
  LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

  //! Assignment Operator
  Group& operator=(const Group& source);

  //! Special LOCA assignment operator
  void 	copy (const NOX::Abstract::Group &source);

  //! Cloning function
  Teuchos::RCP<NOX::Abstract::Group>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Overloaded function evluation routine
  NOX::Abstract::Group::ReturnType computeF();

  //! Overloaded Jacobian evaluation routine
  NOX::Abstract::Group::ReturnType computeJacobian();

  //! Overloaded dfdp sensitivity calculation
  NOX::Abstract::Group::ReturnType computeDfDpMulti	(const std::vector< int > & paramIDs, NOX::Abstract::MultiVector & dfdp, bool isValidF);

  void setParams(const LOCA::ParameterVector& p);

  const LOCA::ParameterVector& getParams() const;

  void setParam(int paramID, double value);

  double getParam(int paramID) const;

  void setParam(std::string paramID, double value);

  double getParam(std::string paramID) const;

  void setScaleVec(const NOX::Abstract::Vector& s);

  const NOX::Abstract::Vector& getScaleVec() const;

  NOX::Abstract::Group::ReturnType
    augmentJacobianForHomotopy(double conParamValue);

  void printSolution (const double conParam) const;

  void printSolution (const NOX::Abstract::Vector &x,
		      const double conParam) const;

  void preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus);
  void postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus);

  // Pseudo Transient methods
  void setAugmentLinearSystem(bool enable,
		  const Teuchos::RCP<N_NLS_NOX::AugmentLinSys>& ls);

  // Continuation flag accessors
  void setNonContinuationFlag (bool value);
  bool getNonContinuationFlag ();

  void setSourceNamesAndScalings( const std::vector<std::string>& srcNames, const std::vector<double>& scales ) 
  { srcNames_ = srcNames; srcScales_ = scales; }

private:

  //! Keep a reference to the loader to set parameters.
  Xyce::Loader::NonlinearEquationLoader& loader;

  //! For output to a file we need xyce's output manager.
  Xyce::IO::OutputMgr& outputMgr;

  //! need xyce's time integration manager.
  Xyce::Analysis::AnalysisManager & anaInt;

  //! Parameter vector container
  LOCA::ParameterVector params;

  //! Utilities for computing derivatives
  LOCA::DerivUtils derivUtils;

  //! Temporary vector used for homotopy calculation
  Xyce::Linear::Vector* tmpVectorPtr;

  //! LOCA Scaling Vector
  const NOX::Abstract::Vector* scalingVecPtr;

  // Objects for Pseudo transient continuation
  bool useAugmentLinSys_;
  Teuchos::RCP<N_NLS_NOX::AugmentLinSys> augmentLSStrategy_;

  // Flag to indicate if this is a traditional newton solve or not.
  bool nonContinuationSolve_;

  // Source names and values being stepped during sequential source stepping
  std::vector<std::string> srcNames_;
  std::vector<double> srcScales_;
};

}}} // namespace N_NLS_LOCA

#endif // Xyce_N_NLS_LOCA_Group_h

