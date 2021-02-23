//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Structs to facilitate device evaluation and assembly
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 07/31/2017
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_AssemblyTypes_h
#define Xyce_N_UTL_AssemblyTypes_h

namespace Xyce{
  namespace Util{

///////////////////////////////////////////////////////////////////////////////
///
/// Holds informations for passing of frequency `F` and `B` entries
///
/// A device needs to insert the values of the `F` and `B` vectors into the
/// correct locations of the F and B matrices.  This struct makes up the
/// basic unit for storage of that information.  It will eventually be placed
/// in a vector for all the appropriate devices in the netlist.
///
/// \author Jason Verley
/// \date 7/13/17
///
///////////////////////////////////////////////////////////////////////////////
struct FreqVecEntry
{
  std::complex<double> val;   /// Value of the entry
  int lid;                /// Local ID of the positive node
};

///////////////////////////////////////////////////////////////////////////////
///
/// Holds informations for passing of frequency `dF/dx` entries
///
/// A device needs to insert the values of the `dF/dx` elements into the
/// correct locations of the Jacobian.  This struct makes up the basic unit for
/// storage of that information.  It will eventually be placed in a vector for
/// all the appropriate devices in the netlist.
///
/// \author Jason Verley
/// \date 7/13/17
///
///////////////////////////////////////////////////////////////////////////////
struct FreqMatEntry
{
  std::complex<double> val;   /// Value of the entry
  int row_lid;                /// Local row ID of the matrix element
  int col_lid;                /// Local row ID of the matrix element
};

}  // end namespace Util
}  // end namespace Xyce

#endif // Xyce_N_UTL_AssemblyTypes_h

