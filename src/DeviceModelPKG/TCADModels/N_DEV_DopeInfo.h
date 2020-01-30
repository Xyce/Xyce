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
// Purpose        : This file contains classes related to the doping specificaiton.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DopeInfo_h
#define Xyce_N_DEV_DopeInfo_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>
#include <N_UTL_Interpolators.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_DopeInfo
// Purpose       : Doping region information.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/28/03
//-----------------------------------------------------------------------------
class DopeInfo : public CompositeParam
{
  friend class ParametricData<DopeInfo>;

public:
  static ParametricData<DopeInfo> &getParametricData();

  DopeInfo();
  bool processParam(Param & ndParam, std::string & param, DevicePDEInstance & di);
  void processParams();

  void setupInfo(
     std::vector<double> & CVec,
     std::vector<double> & CdonorVec,
     std::vector<double> & CacceptorVec,
     std::vector<double> & xVec);

  void setupInfo(
     std::vector<double> & CVec,
     std::vector<double> & CdonorVec,
     std::vector<double> & CacceptorVec,
     std::vector<double> & xVec,
     std::vector<bool>   & el2Vec);

  void setupInfo2d(
     std::vector<double> & CVec,
     std::vector<double> & CdonorVec,
     std::vector<double> & CacceptorVec,
     std::vector<double> & xVec, std::vector<double> &yVec);

  static double nsdep(double x, double W, double Dt);
  static double ngdep(double x, double y, double W, double ax, double ay);
  static double ngdep2(double x, double y, double ax, double ay);
  static double erf(double x);

  static void readDopingFile(std::string & filename, std::vector<double> & xloc,
                             std::vector<double> & nvec);

  static void readDopingFile(std::string & filename, std::vector<double> & xloc,
                             std::vector<double> & nvec, std::vector<double> & pvec);

public:
  std::string name;         // this is also the index into the map.
  std::string type;         // p-type or n-type
  std::string funcType;     // uniform, step, gaussian.
  std::string speciesName;  // boron, phosphorus, etc.
  std::string fileName;     // doping file.
  std::string exprString;   // if expression based, store here.

  // location params:
  double xmin;
  double xmax;
  bool xminGiven;
  bool xmaxGiven;

  double xloc;
  double xwidth;

  double ymin;
  double ymax;
  bool yminGiven;
  bool ymaxGiven;

  double yloc;
  double ywidth;

  double Nmax;     // maximum magnitude
  double Nmin;     // minimum magnitude

  int EL2Present;

  // sometimes having a really high doping hurts convergence.
  // This lets the user set truncate it.
  double Nmax_chop;
  bool Nmax_chopGiven;

  int flatX;
  int flatY;

  // arrays for spline fitting, if specified from file:
  std::vector<double> xlocVec;
  std::vector<double> dopeVec;
  std::vector<double> interpolatedDopeVec;
  Util::akima<double> dopeInterpolator;
};

// inline functions
//-----------------------------------------------------------------------------
// Function      : DopeInfo::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const DopeInfo & di)
{
  os << di.name << ":\n";
  os << "  type     = " << di.type << "\n";
  os << "  funcType = " << di.funcType << "\n";

  if (di.funcType == "expression")
  {
    os << "  exp. string = " << di.exprString << "\n";
  }

  os << "  xloc     = " << di.xloc     << "\n";
  os << "  xwidth   = " << di.xwidth   << "\n";
  os << "  yloc     = " << di.yloc     << "\n";
  os << "  ywidth   = " << di.ywidth   << "\n";

  os << "  Nmax     = " << di.Nmax     << "\n";
  os << "  Nmin     = " << di.Nmin     << "\n";

  os << "  flatX    = " << di.flatX << "\n";
  os << "  flatY    = " << di.flatY << "\n";

  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

#endif

