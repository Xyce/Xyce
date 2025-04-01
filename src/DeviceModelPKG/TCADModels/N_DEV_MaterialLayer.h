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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MaterialLayer_h
#define Xyce_N_DEV_MaterialLayer_h

#include <string>

#include <N_DEV_CompositeParam.h>
#include <N_DEV_MaterialSupport.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MaterialLayer
// Purpose       :
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 7/14/11
//-----------------------------------------------------------------------------
class MaterialLayer : public CompositeParam
{
  friend class ParametricData<MaterialLayer>;

public:
  static ParametricData<MaterialLayer> &getParametricData();

  MaterialLayer (
                 std::string materialName = std::string("gaas"),
                 double donor = 0.0, 
                 double acceptor = 5.00E+19
      );
  virtual ~MaterialLayer ()
  {}

  friend std::ostream & operator<<(std::ostream & os, const MaterialLayer & ml);

public:
  std::string name;
  bool nameGiven;
  std::string material;
  bool materialGiven;
  int NX;
  bool NXGiven;
  int LX;
  int begin;  // beginning mesh point
  int end;    // end mesh point +1

  //////////////////////////////////////
  double diel;  // dielectric constant
  bool dielGiven;  // dielectric constant given flag

  double Ec; // conduction band edge
  bool EcGiven; // conduction band edge given flag
  double Ev; // valance band edge
  bool EvGiven ; // valance band edge given flag
  double EcEff; // conduction band edge, including BGN
  double EvEff; // valance band edge, including BGN

  double bg;  // bandgap
  double bgEff; // effective bandgap (including band-gap narrowing, bgn)

  double Cdonor; // n doping concentration
  bool CdonorGiven; // n doping concentration given flag
  double Cacceptor; // p doping concentration
  bool CacceptorGiven; // p doping concentration given flag

  double narco; // band gap narrowing of conduction band
  bool narcoGiven; // band gap narrowing of conduction band given flag
  double narva; // band gap narrowing of valence band
  bool narvaGiven; // band gap narrowing of valence band given flag

  double dnco; // conduction band density of states multiplier = (md*/mo)^3/2
  double dnva;  // valence band density of states multiplier = (md*/mo)^3/2

  double Nc; // conduction band DOS
  bool NcGiven;
  double Nv; // valance band DOS
  bool NvGiven;

  double emass; // electron DOS effective mass
  bool emassGiven; // electron DOS effective mass given flag
  double hmass; // hole DOS effective mass
  bool hmassGiven; // hole DOS effective mass given flag

  double elmob0; // zero field mobility for electrons (cm2/Vs)
  bool elmob0Given; // zero field mobility for electrons (cm2/Vs) given flag

  double elvsat; // saturation veocity for electrons (cm/s)
  bool elvsatGiven; // saturation veocity for electrons (cm/s) given flag
  double eleo;   // Eo(V/cm) in mobility field dependence

  double homob0; // zero field mobility for holes (cm2/Vs)
  bool homob0Given; // zero field mobility for holes (cm2/Vs) given flag

  double hovsat; // saturation veocity for holes (cm/s)
  bool hovsatGiven; // saturation veocity for holes (cm/s) given flag

  double dir; // direct recombination rate coefficient (cm3/s)
  bool dirGiven; // direct recombination rate coefficient (cm3/s) given flag

  double augnpp; // Auger recombination rate coefficient for npp (cm3/s)
  double augpnn; // Auger recombination rate coefficient for pnn (cm3/s)
  bool augnppGiven; // Auger recombination rate coefficient for npp (cm3/s) given flag
  bool augpnnGiven; // Auger recombination rate coefficient for pnn (cm3/s) given flag

  double srh; // SRH rate coeff (inverse lifetime)
  double srhdet; // energy shift from midgap for SRH

  double Ni; // intrinsic concentration
  bool NiGiven; // intrinsic concentration given flag
  double NiEff; // effective intrinsic concentration, including BGN
  double width;
  bool widthGiven;

  double gradedLayerWidth;
  bool gradedLayerWidthGiven;

  double temperature;

  double electronThermalV;
  bool electronThermalVGiven;
  double holeThermalV;
  bool holeThermalVGiven;
  double latticeConstant;
  bool latticeConstantGiven;
  double defectReactionRadius;
  bool defectReactionRadiusGiven;

  void processParams ();
};

} // namespace Device
} // namespace Xyce

#endif

