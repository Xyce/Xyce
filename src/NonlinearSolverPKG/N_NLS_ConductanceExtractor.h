//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/03/06
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ConductanceExtractor_h
#define Xyce_N_NLS_ConductanceExtractor_h

// ---------- Standard Includes ----------
#include <vector>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_LOA_NonlinearEquationLoader.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : ConductanceExtractor
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------

class ConductanceExtractor 
{
public:
  ConductanceExtractor (NonLinearSolver &nonlinear_solver, Topo::Topology & topology);

  ~ConductanceExtractor ();
  
  bool extract( 
    const std::map<std::string,double> &        inputMap,
    std::vector<double> &                       outputVector,
    std::vector< std::vector<double> > &        jacobian);
  
  bool extract ( 
      const std::string & isoName,
      std::vector< std::vector<double> > & jacobian );
  
  bool setOptions(const Util::OptionBlock& OB);
  
  void printJacobian (
      std::ostream &os,
     const std::map<std::string,double> & inputMap,
     std::vector< std::vector<double> > & jacobian);
  
  void print(
      std::ostream &os,
      const std::string & varName);

private:
  bool setupIDs_( const std::map<std::string, double> & inputMap);
  bool setup_dIdX_Vectors_();
  
  bool setupISO2_IDs_(const std::string & isoName);

private:
  int solutionSize_;

  // temporary stuff, for use with iso devices:
  std::map<std::string, double> varMap_;

  // GID variables
  bool gidsSetUpFlag_;
  std::vector<int> currentGIDs_;
  std::vector<int> currentLIDs_;
  std::vector<int> vsrcPosGIDs_;
  std::vector<int> vsrcPosLIDs_;

  // package references:
  NonLinearSolver & nls_;
  Topo::Topology & top_;
  
  // linear system data:
  Linear::System  * lasSysPtr_;
  Loader::NonlinearEquationLoader  * loaderPtr_;
  Linear::Vector  * rhsVectorPtr_;
  Linear::Vector  * dfdvVectorPtr_;
  Linear::Vector  * NewtonVectorPtr_;
  Linear::Vector  * dxdvVectorPtr_;
  Linear::Solver  * lasSolverPtr_;

  Linear::Vector  * matrixDiagonalPtr_;

  std::vector<Linear::Vector*> dIdxPtrVector_;

  Linear::Matrix  * jacobianMatrixPtr_;
  Linear::Vector  * savedRHSVectorPtr_;
  Linear::Vector  * savedNewtonVectorPtr_;
  Linear::Vector  * gradVectorPtr_;

  Linear::Vector  * columnVectorPtr_;
  const Parallel::ParMap  * columnMapPtr_;
};

} // namespace Nonlinear
} // namespace Xyce

#endif

