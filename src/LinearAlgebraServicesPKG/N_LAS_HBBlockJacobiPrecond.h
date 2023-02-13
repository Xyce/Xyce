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

//-----------------------------------------------------------------------------
//
// Purpose        : block jacobi preconditioner designed for harmonic balance
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 11/11/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_HBBlockJacobiPrecond_h
#define Xyce_N_LAS_HBBlockJacobiPrecond_h

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>

#include <N_LAS_Preconditioner.h>
#include <N_LAS_Problem.h>
#include <N_UTL_OptionBlock.h>

class Epetra_Operator;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_LinearProblem;
class Amesos_BaseSolver;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : HBBlockJacobiPrecond
// Purpose       : interface to block jacobi preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
class HBBlockJacobiPrecond : public Preconditioner
{

public:
  // Constructors
  HBBlockJacobiPrecond(Linear::Builder &builder);

  // Destructor
  virtual ~HBBlockJacobiPrecond() {}

  // Set the preconditioner options
  bool setOptions(const Util::OptionBlock & OB);
  bool setDefaultOptions();

  // Set individual preconditioner options
  bool setParam( const Util::Param & param );

  // Set the fast times being used in the HB analysis.
  void setFastTimes( const std::vector<double> & times )
    { times_ = times; }

  // Register the HB loader
  void registerHBLoader( const Teuchos::RCP<Loader::HBLoader> & hbLoaderPtr )
    { hbLoaderPtr_ = hbLoaderPtr; }

  void setHBFreqs( const std::vector<double> & freqs )
    { freqs_ = freqs; }

  void setHBOsc( const bool osc )
    { hbOsc_ = osc; }

  // Register the HB builder
  void registerHBBuilder( const Teuchos::RCP<HBBuilder> & hbBuilderPtr )
    { hbBuilderPtr_ = hbBuilderPtr; }

  // Register the linear system pointer
  void registerLinearSystem( const Teuchos::RCP< System >& lasSysPtr )
    { lasSysPtr_ = lasSysPtr; }

  // Set or reset the matrix pattern for the preconditioner using problem.
  bool initGraph( const Teuchos::RCP<Problem> & problem );

  // Set the matrix values for the preconditioner
  bool initValues( const Teuchos::RCP<Problem> & problem );

  // Compute the preconditioner using the current matrix values.
  bool compute();

  // Apply the preconditioner; y = M*x.
  int apply( MultiVector & x, MultiVector & y );

  // Return the preconditioner as an Epetra_Operator object.
  Teuchos::RCP<Epetra_Operator> epetraObj() { return epetraPrec_; }

private:

  bool isCorrected_;
  bool hbOsc_;

  // Fourier information.
  // N_ is the number of Fourier coefficients.
  // M_ is the number of positive Fourier coefficients, [0,1,...,M_,-M_,...,-1]
  int N_, beginN_, endN_, M_, maxRefNNZs_;

  // Fast times.
  std::vector<double> times_;

  std::vector<double> freqs_;

  // Harmonic Balance loader.
  Teuchos::RCP<Loader::HBLoader> hbLoaderPtr_;

  // Harmonic Balance builder.
  Teuchos::RCP<HBBuilder> hbBuilderPtr_;

  // Application builder.
  Linear::Builder &               builder_;

  // Linear system.
  Teuchos::RCP<System> lasSysPtr_;

  // Epetra Map for each linear problem's real equivalent form.
  Teuchos::RCP<Epetra_Map> epetraMap_;

  // Epetra CrsGraph for each linear problem's real equivalent form.
  Teuchos::RCP<Epetra_CrsGraph> epetraGraph_;

  // Amesos interface to Klu solver.
  std::vector<Teuchos::RCP<Amesos_BaseSolver> > amesosPtr_;

  // Epetra_CrsMatrix storage for each matrix.
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > epetraMatrix_;

  // Epetra_CrsMatrix storage for the correction matrices (if necessary).
  std::vector<Teuchos::RCP<FilteredMatrix> > diffCMatrix_, diffGMatrix_;

  // Epetra_MultiVector storage, can be reused for each linear system.
  Teuchos::RCP<Epetra_MultiVector> epetraRHS_, epetraSoln_;

  // Current problems being preconditioned.
  std::vector<Teuchos::RCP<Epetra_LinearProblem> > epetraProblem_;

  // Preconditioner as an Epetra_Operator object.
  Teuchos::RCP<Epetra_Operator> epetraPrec_;

  // Options block.
  Teuchos::RCP<const Util::OptionBlock> options_;

  // Serialized matrix objects.
  std::vector<Teuchos::RCP<Epetra_Map> > serialEpetraMap_; 
  std::vector<Teuchos::RCP<Epetra_Import> > serialImporter_;
  std::vector<Teuchos::RCP<Epetra_CrsGraph> > serialGraph_;
  Teuchos::RCP<Epetra_CrsMatrix> serialMatrix_;
  Teuchos::RCP<Epetra_MultiVector> serialVector_;

  // Matrix objects that are only on one processor, using split communicator.
  std::vector<int> gidList_;
  Teuchos::RCP<Epetra_Map> singleMap_;
  Teuchos::RCP<Epetra_CrsGraph> singleGraph_;
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > singleMatrix_;
  Teuchos::RCP<Epetra_MultiVector> singleRHS_, singleSoln_;
  
  // No copying
  HBBlockJacobiPrecond(const HBBlockJacobiPrecond & right);
  HBBlockJacobiPrecond & operator=(const HBBlockJacobiPrecond & right);

  // No comparison
  bool operator==(const HBBlockJacobiPrecond & right) const;
  bool operator!=(const HBBlockJacobiPrecond & right) const;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_HBBlockJacobiPrecond_h
