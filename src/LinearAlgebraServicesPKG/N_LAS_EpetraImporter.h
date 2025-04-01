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
// Purpose        : Abstract interface to an importer between two different maps.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 12/19/20
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_EpetraImporter_h
#define Xyce_N_LAS_EpetraImporter_h

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>

#include <Epetra_Import.h>

#include <N_LAS_Importer.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraImporter
// Purpose       : Epetra interface to a importer object for two separate maps.
// Special Notes : 
// Creator       : Heidi Thornquist, SNL
// Creation Date : 12/19/20
//-----------------------------------------------------------------------------
class EpetraImporter : public Importer
{

public:

  // Basic constructor with from and to maps for the importer
  EpetraImporter( const Parallel::ParMap & target_map, const Parallel::ParMap & source_map );

  // Return Epetra importer object
  Epetra_Import& epetraObj() { return *importer_; }
 
  // Destructor
  ~EpetraImporter() { delete importer_; }

private: 

  Epetra_Import * importer_;

  // Copy constructor
  EpetraImporter( const EpetraImporter& import );

};

} // namespace Linear
} // namespace Xyce

#endif
