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
// Purpose        : Implementation file for importer interface 
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Computational Sciences
//
// Creation Date  : 12/18/20
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_EpetraImporter.h>
#include <N_PDS_EpetraParMap.h>
#include <Epetra_Map.h>

// ---------  Other Includes  -----------

namespace Xyce {
namespace Linear {

  // Basic constructor with from and to maps for the importer
  EpetraImporter::EpetraImporter( const Parallel::ParMap & target_map, const Parallel::ParMap & source_map )
  : Importer( target_map, source_map )
  {
    const Parallel::EpetraParMap& e_target_map = dynamic_cast<const Parallel::EpetraParMap &>( target_map );
    const Parallel::EpetraParMap& e_source_map = dynamic_cast<const Parallel::EpetraParMap &>( source_map );
  
    importer_ = new Epetra_Import( *e_target_map.petraMap(), *e_source_map.petraMap() );
  }

} // namespace Linear
} // namespace Xyce
