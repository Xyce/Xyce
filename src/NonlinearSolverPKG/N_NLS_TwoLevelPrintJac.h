//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 02/21/2019
//
//
//-------------------------------------------------------------------------
#ifndef Xyce_N_NLS_TwoLevelPrintJac_h
#define Xyce_N_NLS_TwoLevelPrintJac_h

#include <vector>
#include <string>
#include <iomanip>
#include <iostream>

namespace Xyce {
namespace Nonlinear {

inline void printJacobian(
    std::ostream & os,
    const std::string idString,
    const std::vector<std::string> & names, 
    const std::vector< std::vector<double> > & jacobian)
{
  int colw=20;
  int idw=idString.size();
  os << idString << std::left << std::setw(25) << " ConArray:";
  os << std::right;

  int numElectrodes = names.size();

  for (int iE1 = 0; iE1 < numElectrodes; ++iE1)
  {
    os << std::setw(colw) << names[iE1];
  }
  os << std::endl;
  for (int iE1 = 0; iE1 < numElectrodes; ++iE1)
  {
    os << idString << " ConArray:"<< std::setw(15) << names[iE1];
    for (int iE2 = 0; iE2 < numElectrodes; ++iE2)
    {
      os << std::scientific << std::setw(colw) << std::setprecision(8) << jacobian[iE1][iE2];
    }
    os << std::endl;
  }
  os << std::endl;
}
} // Nonlinear
} //Xyce

#endif 


