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
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2022
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_OutputConductanceFile_h
#define Xyce_N_ANP_OutputConductanceFile_h

#include <sstream>
#include <iomanip>

#include <N_NLS_ConductanceExtractor.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : writeConductanceFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
inline void writeConductanceFile(
    const std::vector<std::string> &device_names, 
    Nonlinear::ConductanceExtractor &conductance_extractor, 
    const std::string &filename)
{
  std::map<std::string,double> inputMap;

  // load inputMap from tiaParam.device_names option
  for (std::vector<std::string>::const_iterator it = device_names.begin(), 
      end = device_names.end(); it != end; ++it) 
  {
    inputMap[*it] = 0.0;
  }

  int isize = inputMap.size();
  std::vector<double> outputVector(isize, 0.0);
  std::vector< std::vector<double> > jacobian(isize);
  for (int i = 0; i < isize; ++i)
  {
    jacobian[i].resize(isize, 0.0);
  }

  bool b1 = conductance_extractor.extract(inputMap, outputVector, jacobian);

  int iE1, iE2;
  int numElectrodes = isize;

  FILE *fp1;
  fp1 = fopen(filename.c_str(), "w");

  fprintf(fp1, "%s", "Conductance array: \n");
  fprintf(fp1,"%s", "              ");
  if (b1)
  {
    std::map<std::string,double>::const_iterator iterM = inputMap.begin();
    std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
    for (iE2 = 0; iE2 < numElectrodes; ++iE2, ++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1, "\t%14s", srcname.c_str());
    }
    fprintf(fp1, "%s", "\n");

    iterM = inputMap.begin();
    for (iE1 = 0; iE1 < numElectrodes; ++iE1, ++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1,"%14s",srcname.c_str());
      for (iE2 = 0; iE2 < numElectrodes; ++iE2)
      {
        fprintf(fp1,"\t%14.4e",jacobian[iE1][iE2]);
      }
      fprintf(fp1,"%s", "\n");
    }
    fprintf(fp1,"%s", "\n");
  }
  else
  {
    fprintf(fp1,"%s", "\nConductance calculation failed!\n");
  }

  fclose(fp1);
};

} // namespace Analysis
} // namespace Xyce

#endif
