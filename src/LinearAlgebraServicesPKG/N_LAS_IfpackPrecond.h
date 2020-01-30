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
// Purpose        : interface to preconditioner
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_IfpackPrecond_h
#define Xyce_N_LAS_IfpackPrecond_h

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_LAS_Preconditioner.h>
#include <N_LAS_Problem.h>

// ----------  Fwd Declares     ----------

class Epetra_Operator;

class Ifpack_Preconditioner;
class Ifpack_CrsRiluk;
class Ifpack_IlukGraph;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : IfpackPrecond
// Purpose       : interface to ifpack preconditioner
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
class IfpackPrecond : public Preconditioner
{

public:
  // Constructors
  IfpackPrecond();

  // Destructor
  virtual ~IfpackPrecond() {}

  // Set the preconditioner options
  bool setOptions(const Util::OptionBlock & OB);
  bool setDefaultOptions();

  // Set individual preconditioner options
  bool setParam( const Util::Param & param );

  // Set or reset the matrix pattern for the preconditioner using problem.
  // \note The preconditioner will be recreated each time this is called.
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

  // Indicate whether the old or new Ifpack factory interface is being used.
  bool useFactory_;

  // If the factory is used, which preconditioner should be used.
  std::string ifpackType_;

  // Diagonal perturbation, if necessary, to compute better conditioned preconditioner.
  double diagPerturb_;

  // Overlap for Additive Schwarz preconditioner
  int overlap_;
  // Drop tolerance for LU or ILUT preconditioners
  double dropTol_;
  // Fill-factor for ILUT preconditioner
  double ilutFill_;
  // Diagonal shifting relative threshold
  double rThresh_;
  // Diagonal shifting absolute threshold
  double aThresh_;

  // Default preconditioner values
  static const bool useFactory_default_;
  static const std::string ifpackType_default_;
  static const double diagPerturb_default_;
  static const int overlap_default_;
  static const double dropTol_default_;
  static const double ilutFill_default_;
  static const double rThresh_default_;
  static const double aThresh_default_;

  // Old Ifpack preconditioning classes.
  Teuchos::RCP<Ifpack_IlukGraph> ilukGraph_;
  Teuchos::RCP<Ifpack_CrsRiluk>  rILUK_;

  // New Ifpack abstract interface to preconditioners.
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrecond_;

  // Current matrix being preconditioned.
  Teuchos::RCP<Problem> problem_;

  // Preconditioner as an Epetra_Operator object.
  Teuchos::RCP<Epetra_Operator> epetraPrec_;

  // Options block.
  Teuchos::RCP<Util::OptionBlock> options_;

  // No copying
  IfpackPrecond(const IfpackPrecond & right);
  IfpackPrecond & operator=(const IfpackPrecond & right);

  // No comparison
  bool operator==(const IfpackPrecond & right) const;
  bool operator!=(const IfpackPrecond & right) const;

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_IfpackPrecond_h
