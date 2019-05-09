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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Thomas V. Russo, SNL, Component Information and Models
//
// Creation Date  : 08/19/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_TransportHelper_h
#define Xyce_N_DEV_TransportHelper_h

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : TransportHelper
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
class TransportHelper
{
public:
  TransportHelper () :
    flux_bc1(0.0),
    flux_bc2(0.0),
    bcScale1(1.0),
    bcScale2(1.0),
    D_specie  (3.6e-11),
    transportFlag(false)
  {};

  TransportHelper (double D, std::string & n) :
    name(n),
    flux_bc1(0.0),
    flux_bc2(0.0),
    bcScale1(1.0),
    bcScale2(1.0),
    D_specie  (D),
    transportFlag(false)
  {
    if (D!= 0.0) transportFlag = true;
  };

  std::string name;
  std::vector<int> regSubIndexVec;
  std::vector<double> fluxVec;

  double flux_bc1;
  double flux_bc2;
  double bcScale1;
  double bcScale2;

  std::vector<int> specie_id;
  double D_specie; // diffusion constant.

  bool transportFlag;
};

} // namespace Device
} // namespace Xyce

#endif

