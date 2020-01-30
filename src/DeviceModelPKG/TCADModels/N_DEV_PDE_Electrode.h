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

//-------------------------------------------------------------------------
//
// Purpose        : This is the class for mesh processing/ownership.
//                  of two dimensional meshes.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_PDE_Electrode__h
#define Xyce_N_DEV_PDE_Electrode__h

#include <N_DEV_CompositeParam.h>


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : PDE_Electrode
//
// Purpose       : This class contains the base class for user
//                 specifications of electrodes.
//
// Special Notes : In general, this class will ONLY contain info that came
//                 from netlist user-specified vector-composites.  It
//                 does not contain much else.
//
//                 For both the 1D and 2D device, there are other classes
//                 (such as bcData for 1D and deviceInterfaceNode for 2D),
//                 which contain a lot of other information, such as the
//                 names of the circuit nodes, and all the indexing
//                 information.
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_Electrode : public CompositeParam
{
public:
  PDE_Electrode (ParametricData<void> &parametric_data)
    : CompositeParam(parametric_data),
      name   ("ANODE"),// got to call it something...
      nodeName("node1"),
      bcName("bc1"),
      material ("neutral"),
      materialGiven(false),
      oxideBndryFlag(false),
      oxthick(0.0),
      oxcharge(0.0)
  {};

  virtual ~PDE_Electrode () {}

private:
  PDE_Electrode(const PDE_Electrode &);

public:
  virtual void processParams () {}

public:
  std::string name;     // name of the electrode.
  std::string nodeName; // name of the ckt node.
  std::string bcName;   // name of the bc.
  std::string material;
  bool   materialGiven;
  bool   oxideBndryFlag;
  double oxthick;
  double oxcharge;

};

//-----------------------------------------------------------------------------
// Class         : PDE_1DElectrode
// Purpose       : This class contains user specification of a 1D electrode.
//
// Special Notes :
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_1DElectrode : public PDE_Electrode
{
  friend class ParametricData<PDE_1DElectrode>;

public:
  static ParametricData<PDE_1DElectrode> &getParametricData();

  PDE_1DElectrode ();

  virtual ~PDE_1DElectrode () {}
  virtual void processParams ();

private:
  PDE_1DElectrode(const PDE_1DElectrode &);

public:
  double area;
  bool areaGiven;
  double location;
  bool sideGiven;
  std::string side;    // this class implicitly assumes that the device is
  // one dimensional.  The options for side in this
  // class are: left (x=0), middle (0<x<xmax), right (x=xmax)

};

//-----------------------------------------------------------------------------
// Class         : PDE_2DElectrode
// Purpose       : This class contains user specification of a 2D electrode.
//
// Special Notes :
//
// Creator       : Eric Keiter
// Creation Date : 04/17/03
//-----------------------------------------------------------------------------
class PDE_2DElectrode : public PDE_Electrode
{
  friend class ParametricData<PDE_2DElectrode>;

public:
  static ParametricData<PDE_2DElectrode> &getParametricData();

  PDE_2DElectrode ();

  virtual ~PDE_2DElectrode () {}
  virtual void processParams ();

private:
  PDE_2DElectrode(const PDE_2DElectrode &);

public:
  double start;     // beginning location.
  double end;       // ending location.

  bool startGiven;  // beginning location.
  bool endGiven;    // ending location.

  bool sideGiven;
  std::string side;    // this class implicitly assumes that
  // the device is a 4-sided parallelogram.  Any
  // electrode, therefore, is on the top, bottom, left,
  // or right side.

  int iA, iB;
  int uLabel;   // label index
};

} // namespace Device
} // namespace Xyce

#endif
