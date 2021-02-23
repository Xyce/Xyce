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
// Purpose        : This file contains the details of the dope info class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/04/08
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>
#include <fstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DevicePDEInstance.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Expression.h>

#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<DopeInfo>::ParametricData()
{
  addPar("NMAX", 1.0e+15, &DopeInfo::Nmax)
    .setDescription("Maximum value of impurity concentration")
    .setUnit(U_CMM3);
  addPar("NMIN", 0.0, &DopeInfo::Nmin)
    .setDescription("Minimum value of impurity concentration")
    .setUnit(U_CMM3);
  addPar("NMAXCHOP", 1.0e+20, &DopeInfo::Nmax_chop)
    .setGivenMember(&DopeInfo::Nmax_chopGiven)
    .setUnit(U_CMM3);
  addPar("XLOC", 0.0, &DopeInfo::xloc)
    .setDescription("Peak location of the doping in the x-direction")
    .setUnit(U_CM);
  addPar("XMIN", 0.0, &DopeInfo::xmin)
    .setGivenMember(&DopeInfo::xminGiven)
    .setUnit(U_CM);
  addPar("XMAX", 0.0, &DopeInfo::xmax)
    .setGivenMember(&DopeInfo::xmaxGiven)
    .setUnit(U_CM);
  addPar("XWIDTH", 1.0e-3, &DopeInfo::xwidth)
    .setDescription("Distance from nmax to nmin. This is only applicable for the function=gaussian case.")
    .setUnit(U_CM);
  addPar("YLOC", 0.0, &DopeInfo::yloc)
    .setDescription("2D ONLY: Peak location of the doping in the y-direction ()")
    .setUnit(U_CM);
  addPar("YMIN", 0.0, &DopeInfo::ymin)
    .setGivenMember(&DopeInfo::yminGiven)
    .setDescription("2D ONLY: ")
    .setUnit(U_CM);
  addPar("YMAX", 0.0, &DopeInfo::ymax)
    .setGivenMember(&DopeInfo::ymaxGiven)
    .setDescription("2D ONLY: ")
    .setUnit(U_CM);
  addPar("YWIDTH", 1.0e-3, &DopeInfo::ywidth)
    .setDescription("2D ONLY: Distance from nmax to nmin. This is only applicable for the function=gaussian case.")
    .setUnit(U_CM);

  // Set up map for non-double precision variables:
  addPar("NAME", std::string("none"), &DopeInfo::name);
  addPar("FUNCTION", std::string("uniform"), &DopeInfo::funcType)
    .setDescription("Functional form of doping region; options are uniform, gaussian, and step.");
  addPar("TYPE", std::string("ntype"), &DopeInfo::type)
    .setDescription("ntype or ptype");
  addPar("FLATX", 0, &DopeInfo::flatX)
    .setDescription("Determines the doping shape (half-gaussian or a full gaussian)");
  addPar("FLATY", 0, &DopeInfo::flatY)
    .setDescription("2D ONLY: Determines the doping shape (half-gaussian or a full gaussian)");
  addPar("SPECIES", std::string("none"), &DopeInfo::speciesName);
  addPar("FILE", std::string("none"), &DopeInfo::fileName);
  addPar("EXPRESSION", std::string("none"), &DopeInfo::exprString);
  addPar("EL2", 0, &DopeInfo::EL2Present);
}


ParametricData<DopeInfo> &DopeInfo::getParametricData() {
  static ParametricData<DopeInfo> parMap;

  return parMap;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::DopeInfo
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/07/05
// ----------------------------------------------------------------------------
DopeInfo::DopeInfo (const Xyce::Device::SolverState & ss)
  : CompositeParam (getParametricData()),
    name("reg0"),
    type("ntype"),
    funcType("uniform"),
    speciesName("none"),
    fileName("none"),

    xmin(0.0),
    xmax(0.0),
    xloc(0.0),
    xwidth(0.0),

    ymin(0.0),
    ymax(0.0),
    yloc(0.0),
    ywidth(0.0),

    Nmax(1.0e+15),
    Nmin(1.0e+11),

    EL2Present(0),

    Nmax_chop(1.0e+99),
    Nmax_chopGiven(false),
    flatX(0),
    flatY(0),
    solState_(ss)
{}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::processParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
// ----------------------------------------------------------------------------
bool DopeInfo::processParam(Param & ndParam, std::string & param, DevicePDEInstance & di)
{
  bool bsuccess = true;

  return bsuccess;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/31/03
// ----------------------------------------------------------------------------
void DopeInfo::processParams()
{
  {
    ParameterMap::const_iterator p_i = getParameterMap().find(std::string("FUNCTION"));
    const Descriptor &p = *(*p_i).second;

    ExtendedString tmp = getValue<std::string, DopeInfo>(*this, p);
    setValue<std::string, DopeInfo>(*this, p, static_cast<std::string>(tmp.toLower()));
  }

  {
    ParameterMap::const_iterator p_i = getParameterMap().find(std::string("TYPE"));
    const Descriptor &p = *(*p_i).second;

    ExtendedString tmp = getValue<std::string, DopeInfo>(*this, p);
    setValue<std::string, DopeInfo>(*this, p, static_cast<std::string>(tmp.toLower()));
  }
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::setupInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/20/08
// ----------------------------------------------------------------------------
void DopeInfo::setupInfo(
  std::vector<double> & CVec,
  std::vector<double> & CdonorVec,
  std::vector<double> & CacceptorVec,
  std::vector<double> & xVec)
{
  int i(0);
  int NX (CVec.size());
  interpolatedDopeVec.resize(NX,0.0);

  double sign = 0.0;
  if (type == "ptype" || type == "acceptor")
  {
    sign = -1.0;
  }
  else if (type == "ntype" || type == "donor")
  {
    sign = 1.0;
  }

  if (funcType == "uniform")
  {
    for (i=0;i<NX;++i)
    {
      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the range, skip
        {
          continue;
        }
      }

      CVec[i] += sign*Nmax;
      interpolatedDopeVec[i] = Nmax;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax;
      }
    }
  }
  else if (funcType == "gaussian")
  {
    double deltaX = fabs(xwidth);

    double ax = log(Nmax/Nmin)/(deltaX*deltaX);
    for (i=0;i<NX;++i)
    {
      double scalarX = 1.0;
      double abs_dx;

      if (given("XLOC") && given("XWIDTH") && deltaX != 0.0)
      {
        abs_dx = fabs(xVec[i]-xloc);

        // if true gaussian, x-section:
        if (flatX==0)
        {
          scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
        }
        else if (flatX>0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] > xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
        else if (flatX<0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] < xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
      }

      CVec[i] += sign*Nmax*scalarX;
      interpolatedDopeVec[i] = Nmax*scalarX;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax*scalarX;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax*scalarX;
      }
    }
  }
  else if (funcType == "step")
  {
    for (i=0;i<NX;++i)
    {
      double x  = xVec[i];
      bool regOn = true;

      if (given("XLOC"))
      {
        if (flatX ==  0) regOn = true;

        if (flatX == -1)
        {
          if (x > xloc) regOn = false;
          else          regOn = true;
        }

        if (flatX == +1)
        {
          if (x < xloc) regOn = false;
          else          regOn = true;
        }
      }

      CVec[i] += (regOn)?(sign*Nmax):(sign*Nmin);
      interpolatedDopeVec[i] = (regOn)?(Nmax):(Nmin);
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
    }
  }
  else if (funcType == "expression")
  {
    if (exprString == "none")
    {
      Report::UserFatal() <<  "Dope Region : "
                          <<  name
                          <<  " has specified the expression specification, but not provided an expression.";
    }
    else
    {
      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "DopeInfo::setupInfo: exprString = " << exprString << std::endl;
      }
      Util::Expression expr(solState_.expressionGroup_,exprString);
      //expr.set(exprString);

      for (i=0;i<NX;++i)
      {
        double dopeValue(0.0);

#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
        expr.set_var(std::string("#X"), xVec[i]);
#endif
        expr.evaluateFunction (dopeValue);
        CVec[i] += sign*dopeValue;
        interpolatedDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else if (funcType == "file")
  {
    if (fileName == "none")
    {
      Report::UserFatal() << "Dope Region : "
                          << name
                          << " has specified the file specification, but not specified a file name.";
    }
    else
    {
      readDopingFile (fileName, xlocVec, dopeVec);
      dopeInterpolator.clear(); 
      dopeInterpolator.init(xlocVec, dopeVec);

      // if the user has requested that this be truncated to a max value,
      // do it here:
      if (Nmax_chopGiven)
      {
        int dopeSize=dopeVec.size();
        for (int id=0;id<dopeSize;++id)
        {
          if (dopeVec[id] > Nmax_chop)
          {
            dopeVec[id] = Nmax_chop;
          }
        }
      }

      for (i=0;i<NX;++i)
      {
        double xtmp = xVec[i];
        double dopeValue(0.0);
        dopeInterpolator.eval(xlocVec, dopeVec, xtmp, dopeValue);
        CVec[i] += sign*dopeValue;
        interpolatedDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else
  {
    Report::UserFatal() << "Unrecognized Dope Region function type:  "
                        << funcType
                        << "  for region: "
                        << name;
  }
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::setupInfo2d
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/20/08
// ----------------------------------------------------------------------------
void DopeInfo::setupInfo2d(
  std::vector<double> & CVec,
  std::vector<double> & CdonorVec,
  std::vector<double> & CacceptorVec,
  std::vector<double> & xVec,
  std::vector<double> & yVec)
{
  int i(0);
  int numMeshPoints (CVec.size());

  double sign = 1.0;
  if (type == "ptype" || type == "acceptor") sign = -1.0;

  if (funcType == "uniform")
  {
    for (i=0;i<numMeshPoints;++i)
    {

      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the x-range, skip
        {
          continue;
        }
      }

      if (ymaxGiven && yminGiven)
      {
        if (yVec[i] > ymax || yVec[i] < ymin) // if outside the y-range, skip
        {
          continue;
        }
      }

      CVec[i] += sign*Nmax;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax;
      }
      else
      {
        CdonorVec[i] += Nmax;
      }
    }
  }
  else if (funcType == "gaussian")
  {
    double deltaX = fabs(xwidth);
    double deltaY = fabs(ywidth);

    double ax = 0.0;
    double ay = 0.0;

    if (deltaX!=0.0)
    {
      ax = log(Nmax/Nmin)/(deltaX*deltaX);
    }

    if (deltaY!=0.0)
    {
      ay = log(Nmax/Nmin)/(deltaY*deltaY);
    }

    for (i=0;i<numMeshPoints;++i)
    {
      double scalarX = 1.0;
      double scalarY = 1.0;
      double abs_dx, abs_dy;

      if (given("XLOC") && given("XWIDTH") && deltaX != 0.0)
      {
        abs_dx = fabs(xVec[i]-xloc);

        // if true gaussian, x-section:
        if (flatX==0)
        {
          scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
        }
        else if (flatX>0) // half-guassian, x-section:
        {
          bool flatReg = (xVec[i] > xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
        else if (flatX<0) // half-guassian, x-section:
        {
          bool flatReg = (xVec[i] < xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
      }

      if (given("YLOC") && given("YWIDTH") && deltaY != 0.0)
      {
        abs_dy = fabs(yVec[i]-yloc);
        // if true gaussian, y-section:
        if (flatY==0)
        {
          scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
        }
        else if (flatY>0) // half-guassian, y-section:
        {
          bool flatReg = (yVec[i] > yloc);

          if (!flatReg)
          {
            scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
          }
        }
        else if (flatY<0) // half-guassian, y-section:
        {
          bool flatReg = (yVec[i] < yloc);

          if (!flatReg)
          {
            scalarY *= ngdep2(0.0, abs_dy, 1.0, ay);
          }
        }
      }

      CVec[i] += sign*Nmax*scalarX*scalarY;

      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax*scalarX*scalarY;
      }
      else
      {
        CdonorVec[i] += Nmax*scalarX*scalarY;
      }
    }
  }
  else if (funcType == "step")
  {
    for (i=0;i<numMeshPoints;++i)
    {
      double x = xVec[i];
      double y = yVec[i];
      bool regOnX = true;
      bool regOnY = true;

      if (given("YLOC"))
      {
        if (flatY ==  0) regOnY = true;

        if (flatY == -1)
        {
          if(y > yloc) regOnY = false;
          else            regOnY = true;
        }

        if (flatY == +1)
        {
          if (y < yloc) regOnY = false;
          else             regOnY = true;
        }
      }

      if (given("XLOC"))
      {
        if (flatX ==  0) regOnX = true;

        if (flatX == -1)
        {
          if(x > xloc) regOnX = false;
          else            regOnX = true;
        }

        if (flatX == +1)
        {
          if (x < xloc) regOnX = false;
          else             regOnX = true;
        }
      }
      bool regOn = (regOnX && regOnY);

      CVec[i] += (regOn)?(sign*Nmax):(sign*Nmin);
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
      else
      {
        CdonorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
    }
  }
  else
  {
    Report::UserFatal() <<  "Unrecognized Dope Region function type:  "
                        << funcType
                        <<  "  for region: "
                        << name;
  }

}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::nsdep
// Purpose       : This function returns an approximate deposition profile
//                 of a step implant driven in an inert environment.
// Special Notes :
//                        1       W/2 + x          W/2 - x
//        nsdep(x,W,Dt) = - (erf(---------) + erf(----------))
//                        2      2*sqrt(Dt)       2*sqrt(Dt)
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::nsdep(double x, double W, double Dt)
{
  double D  = 2.0 * sqrt(Dt);
  double Wh = W / 2.0;
  return 0.5 * (erf((Wh + x)/D) + erf((Wh - x)/D));
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::ngdep
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
//                 Adapted from a similar function in SGF.
//
// Special Notes : This function is a little flakey, in that it assumes
//                 that we're using a cylindrical geomtry
//                 (x = radius, y = height.) It also assumes that the (0,0)
//                 origin is in the upper-left-hand corner of the mesh.
//
//                 So, y=0.0 is the upper surface of the device, and the
//                 implant is coming from the top (above y=0.0).  Hence,
//                 the (y>0) conditional.
//
//                 Also, the width parameter (W), is set up to be a
//                 diameter about x=0, which is the reason for the 0.5*W -
//                 half of this diameter will impact this radius.
//
//                 The implant is completely flat and constant in the
//                 x-direction, as long as fabs(x) is less than W/2.0.
//                 Beyond W/2.0, the gaussian profile kicks in.
//
//                 The parameters ax and ay are scaling parameters, and
//                 correspond to how much you want the doping to vary with
//                 space.  A typical value for either can be set up as:
//
//                 Ax = ln (Nhi / Nlo)/(Rx*Rx)
//
//                 where:
//
//                   Nhi = the max. level of doping
//                   Nlo = the min. level of doping
//                   Rx  = distance over which this doping should vary.
//
//                   Nhi/Nlo = 10^N, where N = # of orders of magnitude
//                   that should vary between x=0 and x=Rx.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::ngdep(double x, double y, double W, double ax, double ay)
{
  double xprime = fabs(x) - (0.5 * W);
  return ((xprime <= 0.0) ? 1.0 : exp(-ax*xprime*xprime))*
    ((y      >  0.0) ? 0.0 : exp(-ay*y*y));
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::ngdep2
//
// Purpose       : This function returns an approximate Gaussian deposition.
//
// Special Notes : This function is a modification of the original ngdep
//                 (see above), and I designed it to address some of the
//                 peculiarities of ngdep.
//
//                 (1) I've gotten rid of the width, W.  I'm just
//                 going to assume that whoever is calling this function
//                 can set that(the constant region) up on their own.
//
//                 (2) I've removed the conditionals cause things to
//                 be set to zero, or one, or whatever, if you are on one
//                 side or another of the suface.  I'm assuming that
//                 whoever calls this function can do that themselves, if
//                 they need to.
//
//                 (3) I've removed the stuff that sets the retVal to zero
//                 for y>0.  Again, this is the user's problem.
//
//                 ax and ay mean the same as they did for the original
//                 ngdep.  (see above).
//
//                 It is possible to use this for the 1D case, pretty
//                 easily.  Set the xflag to false, hold y fixed at zero,
//                 and have x correspond to the 1D mesh locations. (or, set
//                 xflag to true, hold x fixed at zero, and let y
//                 correspond to 1D mesh locations - either way).
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/28/03
// ----------------------------------------------------------------------------
double DopeInfo::ngdep2(double x, double y, double ax, double ay)
{
  double retVal = exp(-ax*x*x)* exp(-ay*y*y);
  return retVal;
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::erf
// Purpose       : This function returns the error function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/25/02
// ----------------------------------------------------------------------------
double DopeInfo::erf(double x)
{
  double t1 = 1.0 / (1.0 + 0.3275911 * fabs(x));
  double t2 = t1 * t1;
  double t3 = t2 * t1;
  double t4 = t3 * t1;
  double t5 = t4 * t1;
  double result = 1.0 - (0.254829592*t1 - 0.284496736*t2 + 1.421413741*t3 -
                         1.453152027*t4 + 1.061405429*t5) * exp(-x*x);
  return (x < 0.0) ? -result : result;
}

//-----------------------------------------------------------------------------
// Function      : DopeInfo::readDopingFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/07
//-----------------------------------------------------------------------------
void DopeInfo::readDopingFile(
  std::string &         filename,
  std::vector<double> & xloc,
  std::vector<double> & nvec)
{
  std::ifstream input;
  double x_loc(0.0);
  double value(0.0);
  xloc.clear();
  nvec.clear();

  // Error out if the user-specified DOPING_FILE does not exist, cannot
  // be opened, or is a directory name rather than a file name.  See SON 
  // Bug 785 and SRN Bug 2100 for more details.
  if ( !(Util::checkIfValidFile(filename)) )
  {
    Report::UserFatal() << "Error: Cannot find doping file: " << filename;
  }

  input.open( filename.c_str(), std::ios::in );
  if ( input.good() )
  {
    bool endOfFile = input.eof();
    while (!endOfFile)
    {
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> x_loc;
      }
      else
      {
        break;
      }

      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value;
      }
      else
      {
        break;
      }
      xloc.push_back(x_loc);
      nvec.push_back(value);
    }
    input.close();
  }
  else
  {
    Report::UserFatal() << "Error: Cannot open doping file: " << filename;
  }
}

//-----------------------------------------------------------------------------
// Function      : DopeInfo::readDopingFile
// Purpose       :
// Special Notes : This version assumes 2 dopants are in the file, P and N.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/07
//-----------------------------------------------------------------------------
void DopeInfo::readDopingFile(
  std::string &         filename,
  std::vector<double> & xloc,
  std::vector<double> & nvec,
  std::vector<double> & pvec)
{
  std::ifstream input;
  double x_loc(0.0);
  double value1(0.0);
  double value2(0.0);
  xloc.clear();
  nvec.clear();
  pvec.clear();

  // Error out if the user-specified DOPING_FILE does not exist, cannot
  // be opened, or is a directory name rather than a file name.  See SON 
  // Bug 785 and SRN Bug 2100 for more details.
  if ( !(Util::checkIfValidFile(filename)) )
  {
    Report::UserFatal() << "Error: Cannot find doping file: " << filename;
  }

  input.open( filename.c_str(), std::ios::in );
  if ( input.good() )
  {
    bool endOfFile = input.eof();
    while (!endOfFile)
    {
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> x_loc;
      }
      else
      {
        break;
      }

      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value1;
      }
      else
      {
        break;
      }
      endOfFile = input.eof();
      if (!endOfFile)
      {
        input >> value2;
      }
      else
      {
        break;
      }

      xloc.push_back(x_loc);
      nvec.push_back(value1);
      pvec.push_back(value2);
    }
    input.close();
  }
  else
  {
    Report::UserFatal() << "Error: Cannot open doping file: " << filename;
  }
}

// ----------------------------------------------------------------------------
// Function      : DopeInfo::setupInfo
// Purpose       :
// Special Notes : Hack for support of EL2
// Scope         : public
// Creator       : Jason Verley, SNL
// Creation Date : 11/20/08
// ----------------------------------------------------------------------------
void DopeInfo::setupInfo(
  std::vector<double> & CVec,
  std::vector<double> & CdonorVec,
  std::vector<double> & CacceptorVec,
  std::vector<double> & xVec,
  std::vector<bool>   & el2Vec)
{
  int i(0);
  int NX (CVec.size());
  interpolatedDopeVec.resize(NX,0.0);

  double sign = 0.0;
  if (type == "ptype" || type == "acceptor")
  {
    sign = -1.0;
  }
  else if (type == "ntype" || type == "donor")
  {
    sign = 1.0;
  }

  if (EL2Present > 0)
  {
    for (i=0;i<NX;++i)
    {
      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the range, skip
        {
          continue;
        }
      }
      el2Vec[i] = true;
    }
  }

  if (funcType == "uniform")
  {
    for (i=0;i<NX;++i)
    {
      if (xmaxGiven && xminGiven)
      {
        if (xVec[i] > xmax || xVec[i] < xmin) // if outside the range, skip
        {
          continue;
        }
      }

      CVec[i] += sign*Nmax;
      interpolatedDopeVec[i] = Nmax;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax;
      }
    }
  }
  else if (funcType == "gaussian")
  {
    double deltaX = fabs(xwidth);

    double ax = log(Nmax/Nmin)/(deltaX*deltaX);
    for (i=0;i<NX;++i)
    {
      double scalarX = 1.0;
      double abs_dx;

      if (given("XLOC") && given("XWIDTH") && deltaX != 0.0)
      {
        abs_dx = fabs(xVec[i]-xloc);

        // if true gaussian, x-section:
        if (flatX==0)
        {
          scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
        }
        else if (flatX>0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] > xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
        else if (flatX<0) // half-gaussian, x-section:
        {
          bool flatReg = (xVec[i] < xloc);

          if (!flatReg)
          {
            scalarX *= ngdep2(abs_dx, 0.0, ax, 1.0);
          }
        }
      }

      CVec[i] += sign*Nmax*scalarX;
      interpolatedDopeVec[i] = Nmax*scalarX;
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += Nmax*scalarX;
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += Nmax*scalarX;
      }
    }
  }
  else if (funcType == "step")
  {
    for (i=0;i<NX;++i)
    {
      double x  = xVec[i];
      bool regOn = true;

      if (given("XLOC"))
      {
        if (flatX ==  0) regOn = true;

        if (flatX == -1)
        {
          if (x > xloc) regOn = false;
          else          regOn = true;
        }

        if (flatX == +1)
        {
          if (x < xloc) regOn = false;
          else          regOn = true;
        }
      }

      CVec[i] += (regOn)?(sign*Nmax):(sign*Nmin);
      interpolatedDopeVec[i] = (regOn)?(Nmax):(Nmin);
      if (type == "ptype" || type == "acceptor")
      {
        CacceptorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
      else if (type == "ntype" || type == "donor")
      {
        CdonorVec[i] += (regOn)?(Nmax):(sign*Nmin);
      }
    }
  }
  else if (funcType == "expression")
  {
    if (exprString == "none")
    {
      Report::UserFatal() << "Dope Region : "
                          << name
                          << " has specified the expression specification, but not provided an expression.";
    }
    else
    {
      if (DEBUG_DEVICE)
      {
        Xyce::dout() << "DopeInfo::setupInfo: exprString = " << exprString << std::endl;
      }
      Util::Expression expr(solState_.expressionGroup_,exprString);
      //expr.set(exprString);

      for (i=0;i<NX;++i)
      {
        double dopeValue(0.0);

#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
        expr.set_var(std::string("#X"), xVec[i]);
#endif
        expr.evaluateFunction (dopeValue);
        CVec[i] += sign*dopeValue;
        interpolatedDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else if (funcType == "file")
  {
    if (fileName == "none")
    {
      Report::UserFatal() <<  "Dope Region : "
                          << name
                          << " has specified the file specification, but not specified a file name.";
    }
    else
    {
      readDopingFile (fileName, xlocVec, dopeVec);
      dopeInterpolator.clear(); 
      dopeInterpolator.init(xlocVec, dopeVec);

      // if the user has requested that this be truncated to a max value,
      // do it here:
      if (Nmax_chopGiven)
      {
        int dopeSize=dopeVec.size();
        for (int id=0;id<dopeSize;++id)
        {
          if (dopeVec[id] > Nmax_chop)
          {
            dopeVec[id] = Nmax_chop;
          }
        }
      }

      for (i=0;i<NX;++i)
      {
        double xtmp = xVec[i];
        double dopeValue(0.0);

        dopeInterpolator.eval(xlocVec, dopeVec, xtmp, dopeValue);

        CVec[i] += sign*dopeValue;
        interpolatedDopeVec[i] = dopeValue;
        if (type == "ptype" || type == "acceptor")
        {
          CacceptorVec[i] += dopeValue;
        }
        else if (type == "ntype" || type == "donor")
        {
          CdonorVec[i] += dopeValue;
        }
      }
    }
  }
  else
  {
    Report::UserFatal() << "Unrecognized Dope Region function type:  "
                        << funcType
                        << "  for region: "
                        << name;
  }
}

} // namespace Device
} // namespace Xyce
