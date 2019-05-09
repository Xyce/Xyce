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


//----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 07/01/11
//
//
//
//
//----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_MaterialLayer.h>
#include <N_DEV_Const.h>


namespace Xyce {
namespace Device {

template<>
ParametricData<MaterialLayer>::ParametricData()
{
    // Set up map for double precision variables:
  addPar ("DIEL", 13.1, &MaterialLayer::diel)
    .setGivenMember(&MaterialLayer::dielGiven);
    //
  addPar ("CON", 1.422482, &MaterialLayer::Ec)
    .setGivenMember(&MaterialLayer::EcGiven);

  addPar ("VAL", 0.0, &MaterialLayer::Ev)
    .setGivenMember(&MaterialLayer::EvGiven);

  addPar ("NDOPE", 0.0, &MaterialLayer::Cdonor)
    .setGivenMember(&MaterialLayer::CdonorGiven);

  addPar ("PDOPE", 5.0e+19, &MaterialLayer::Cacceptor)
    .setGivenMember(&MaterialLayer::CacceptorGiven);

    ////////
  addPar ("NARCO", 0.047, &MaterialLayer::narco)
    .setGivenMember(&MaterialLayer::narcoGiven);

  addPar ("NARVA", 0.047, &MaterialLayer::narva)
    .setGivenMember(&MaterialLayer::narvaGiven);

    //////
  addPar ("EMASS", 0.067, &MaterialLayer::emass)
    .setGivenMember(&MaterialLayer::emassGiven);

  addPar ("HMASS", 0.5, &MaterialLayer::hmass)
    .setGivenMember(&MaterialLayer::hmassGiven);

    /////
  addPar ("ELMOB0", 2240.0, &MaterialLayer::elmob0)
    .setGivenMember(&MaterialLayer::elmob0Given);

  addPar ("HOMOB0", 30.0, &MaterialLayer::homob0)
    .setGivenMember(&MaterialLayer::homob0Given);

    /////
  addPar ("ELVSAT", 7.7e+6, &MaterialLayer::elvsat)
    .setGivenMember(&MaterialLayer::elvsatGiven);

  addPar ("HOVSAT", 7.7e+6, &MaterialLayer::hovsat)
    .setGivenMember(&MaterialLayer::hovsatGiven);

  addPar ("NI", 1.79e+6, &MaterialLayer::Ni)
    .setGivenMember(&MaterialLayer::NiGiven);

  addPar ("WIDTH", 1.0e-6, &MaterialLayer::width)
    .setGivenMember(&MaterialLayer::widthGiven);

  addPar ("GRADEDWIDTH", 0.0, &MaterialLayer::gradedLayerWidth)
    .setGivenMember(&MaterialLayer::gradedLayerWidthGiven),

    // non-doubles
    addPar ("MATERIAL", "gaas", &MaterialLayer::material)
    .setGivenMember(&MaterialLayer::materialGiven);

  addPar ("NAME", "EMITTER", &MaterialLayer::name)
    .setGivenMember(&MaterialLayer::nameGiven);

  addPar ("NX", 25, &MaterialLayer::NX)
    .setGivenMember(&MaterialLayer::NXGiven);

  addPar ("ValenceBandDOS", 2.66e19, &MaterialLayer::Nv)
    .setGivenMember(&MaterialLayer::NvGiven);

  addPar ("ConductionBandDOS", 2.89e19, &MaterialLayer::Nc)
    .setGivenMember(&MaterialLayer::NcGiven);

}

ParametricData<MaterialLayer> &MaterialLayer::getParametricData() {
  static ParametricData<MaterialLayer> parMap;

  return parMap;
}

template <typename ScalarT> ScalarT Xycemax ( ScalarT f1, ScalarT f2) { return f1 > f2 ? f1 : f2; }
template <typename ScalarT> ScalarT Xycemin ( ScalarT f1, ScalarT f2) { return f1 < f2 ? f1 : f2; }

// Class MaterialLayer
//-----------------------------------------------------------------------------
// Function      : MaterialLayer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 07/01/11
//-----------------------------------------------------------------------------
  MaterialLayer::MaterialLayer(std::string materialName,
                               double donor, double acceptor)
  : CompositeParam(getParametricData()),
    name("EMITTER"),
    nameGiven(false),
    material(materialName),
    materialGiven(false),
    NX(25),
    NXGiven(false),
    LX(24),
    begin(0),  // beginning mesh point
    end(25),    // end mesh point (+1)
    diel(13.1),  // dielectric constant
    dielGiven(false),  // dielectric constant given flag

    Ec(1.422482), // conduction band edge
    EcGiven(false), // conduction band edge given flag
    Ev(0), // valance band edge
    EvGiven(false), // valance band edge given flag
    EcEff(1.375482), // conduction band edge, including BGN
    EvEff(0.047), // valance band edge, including BGN
    bg(1.422482),  // bandgap
    bgEff(1.328482), // effective bandgap (including band-gap narrowing, bgn)
    Cdonor(donor), // n doping concentration
    CdonorGiven(true), // n doping concentration given flag
    Cacceptor(acceptor), // p doping concentration
    CacceptorGiven(true), // p doping concentration given flag


    // these defaults assume band-gap narrowing is turned off:
    narco(0.0), // band gap narrowing of conduction band
    narcoGiven(false), // band gap narrowing of conduction band given flag
    narva(0.0), // band gap narrowing of valence band
    narvaGiven(false), // band gap narrowing of valence band given flag


    dnco(0.01734), // conduction band density of states multiplier ((md*/mo)^3/2
    dnva(0.4318),  // valence band density of states multiplier ((md*/mo)^3/2
    Nc(1.0), // conduction band DOS
    NcGiven(false),
    Nv(1.0), // valance band DOS
    NvGiven(false),
    emass(0.067), // electron DOS effective mass 
    emassGiven(false), // electron DOS effective mass given flag
    hmass(0.5), // hole DOS effective mass
    hmassGiven(false), // hole DOS effective mass given flag
    elmob0(2240), // zero field mobility for electrons (cm2/Vs)
    elmob0Given(false), // zero field mobility for electrons (cm2/Vs) given flag

    elvsat (7.70E+06), // saturation velocity for electrons (cm/s)
    elvsatGiven (false), // saturation velocity for electrons (cm/s) given flag
    eleo(4000), // Eo(V/cm) in mobility field dependence
    homob0 (30), // zero field mobility for holes (cm2/Vs)
    homob0Given(false), // zero field mobility for holes (cm2/Vs) given flag
    hovsat (7.70E+06), // saturation velocity for holes (cm/s)
    hovsatGiven (false), // saturation velocity for holes (cm/s) given flag

    dir(1.00E-10), // direct recombination rate coefficient (cm3/s)
    dirGiven(false), // direct recombination rate coefficient (cm3/s) given flag
    augnpp(1.00E-30), // Auger recombination rate coefficient for npp (cm3/s)
    augpnn(1.00E-30), // Auger recombination rate coefficient for pnn (cm3/s)
    augnppGiven(false), // Auger recombination rate coefficient for npp (cm3/s) given flag
    augpnnGiven(false), // Auger recombination rate coefficient for pnn (cm3/s) given flag
    srh(0), // SRH rate coeff (inverse lifetime)
    srhdet(0), // energy shift from midgap for SRH
    Ni(1.79e+6),
    NiGiven(false),
    NiEff(1.79e+6), // Ni including BGN
    width(1.0e-5),
    widthGiven(false),
    gradedLayerWidth(0.0),
    gradedLayerWidthGiven(false),
    temperature(CONSTREFTEMP),
    electronThermalV(2.3e7),
    electronThermalVGiven(false),
    holeThermalV(1.9e7),
    holeThermalVGiven(false),
    latticeConstant(4.25e-7),
    latticeConstantGiven(false),
    defectReactionRadius(4.25e-7),
    defectReactionRadiusGiven(false)
{
  ExtendedString mN = material;
  mN.toLower();
  material = mN;

  if (NXGiven)
  {
    end = begin + NX;
  }

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : MaterialLayer::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/28/12
//-----------------------------------------------------------------------------
void MaterialLayer::processParams ()
{
  // set some material parameters to appropriate values
  ExtendedString tmpMat = material;
  tmpMat.toUpper ();

  if (!emassGiven) { emass = MaterialSupport::get_DOS_EffectiveMassN(tmpMat); }
  if (!hmassGiven) { hmass = MaterialSupport::get_DOS_EffectiveMassP(tmpMat); }

  dnco = pow(emass,1.5);
  dnva = pow(hmass,1.5);

  if(!NcGiven)
    Nc = MaterialSupport::getNc(tmpMat, temperature);
  if(!NvGiven)
    Nv = MaterialSupport::getNv(tmpMat, temperature);

  if (!narcoGiven) 
  {
    narco =  0.047;
  }
  if (!narvaGiven) 
  {
    narva =  0.047;
  }

  if (tmpMat == "GAAS")
    {
      if (!dielGiven) { diel = 13.1;}
      if (!elmob0Given) { elmob0 = 2240;}
      if (!homob0Given) { homob0 = 30;}
      if (!NiGiven) { Ni = 1.79e+6; }
      if (!EcGiven) {Ec =  1.422482;}
      if (!EvGiven) {Ev =  0;}
      if (!elvsatGiven) {elvsat =  7.7e+6;}
      if (!hovsatGiven) {hovsat =  7.7e+6;}
      if (!electronThermalVGiven) {electronThermalV = 4.5e7;}
      if (!holeThermalVGiven) {holeThermalV = 4.5e7;}
      if (!defectReactionRadiusGiven) {defectReactionRadius = 4.25e-7;}
      if (!latticeConstantGiven) {latticeConstant = 5.63e-8;}
      if (!narcoGiven) { narco =  0.047;}
      if (!narvaGiven) { narva =  0.047;}
      if (!CdonorGiven) { Cdonor =  0.0;}
      if (!CacceptorGiven) { Cacceptor =  5.0e+19;}
      if (!augpnnGiven) { augpnn = 1.0e-30;}
      if (!augnppGiven) { augnpp = 1.0e-30;}
      if (!dirGiven) { dir = 1.0e-10;}
    }
  else if(tmpMat == "INGAP")
  {
    if (!dielGiven) { diel=11.8;}
    if (!elmob0Given) { elmob0 = 850;}
    if (!homob0Given) { homob0 = 30;}
    if (!NiGiven) { Ni = 1.2e+3;}
    if (!EcGiven) {Ec =  1.45248;}
    if (!EvGiven) {Ev =  -0.4085;}
    if (!elvsatGiven) {elvsat =  7.7e+6;}
    if (!hovsatGiven) {hovsat =  7.7e+6;}
    if (!narcoGiven) { narco =  0.013;}
    if (!narvaGiven) { narva =  0.013;}
    if (!CdonorGiven) { Cdonor =  5.0e+17;}
    if (!CacceptorGiven) { Cacceptor =  0.0;}
  }
  else if(tmpMat == "ALGAAS")
  {
    if (!dielGiven) { diel=12.95;}
    if (!elmob0Given) { elmob0 = 600;}
    if (!homob0Given) { homob0 = 80;}
    if (!NiGiven) { Ni = 2.1e+3;}
    if (!EcGiven) {Ec =  1.45248;}
    if (!EvGiven) {Ev =  -0.4085;}
    if (!elvsatGiven) {elvsat =  7.7e+6;}
    if (!hovsatGiven) {hovsat =  7.7e+6;}
    if (!narcoGiven) { narco =  0.010;}
    if (!narvaGiven) { narva =  0.010;}
    if (!CdonorGiven) { Cdonor =  0.0;}
    if (!CacceptorGiven) { Cacceptor =  5.0e+17;}
  }
  else if (tmpMat == "INALAS")
  {
    if (!dielGiven) { diel = 12.5;}
    if (!elmob0Given) { elmob0 = 1000;}
    if (!homob0Given) { homob0 = 30;}
    if (!NiGiven) { Ni = 1.79e+6; }
    if (!EcGiven) {Ec =  1.46;}
    if (!EvGiven) {Ev =  0;}
    if (!elvsatGiven) {elvsat =  4.7e+6;}
    if (!hovsatGiven) {hovsat =  3.0e+6;}
    if (!narcoGiven) { narco =  0.0;}
    if (!narvaGiven) { narva =  0.0;}
    if (!CdonorGiven) { Cdonor =  0.0;}
    if (!CacceptorGiven) { Cacceptor =  5.0e+19;}
  }
  else if (tmpMat == "INGAAS")
  {
    if (!dielGiven) { diel = 14;}
    if (!elmob0Given) { elmob0 = 3000;}
    if (!homob0Given) { homob0 = 76;}
    if (!NiGiven) { Ni = 1.79e+6; }
    if (!EcGiven) {Ec =  0.96;}
    if (!EvGiven) {Ev =  0.21;}
    if (!elvsatGiven) {elvsat =  8.4e+6;}
    if (!hovsatGiven) {hovsat =  4.8e+6;}
    if (!narcoGiven) { narco =  0.036;}
    if (!narvaGiven) { narva =  0.036;}
    if (!CdonorGiven) { Cdonor =  0.0;}
    if (!CacceptorGiven) { Cacceptor =  5.0e+19;}
  }
  else if (tmpMat == "INGAAP")
  {
    if (!dielGiven) { diel = 14;}
    if (!elmob0Given) { elmob0 = 3500;}
    if (!homob0Given) { homob0 = 50;}
    if (!NiGiven) { Ni = 1.79e+6; }
    if (!EcGiven) {Ec =  0.75;}
    if (!EvGiven) {Ev =  0;}
    if (!elvsatGiven) {elvsat =  1.3e+7;}
    if (!hovsatGiven) {hovsat =  6.6e+6;}
    if (!narcoGiven) { narco =  0.047;}
    if (!narvaGiven) { narva =  0.047;}
    if (!CdonorGiven) { Cdonor =  0.0;}
    if (!CacceptorGiven) { Cacceptor =  5.0e+19;}
  }
  else if (tmpMat == "INT")
  {
    if (!dielGiven) { diel = 12.6;}
    if (!elmob0Given) { elmob0 = 4100;}
    if (!homob0Given) { homob0 = 150;}
    if (!NiGiven) { Ni = 1.79e+6; }
    if (!EcGiven) {Ec =  1.21;}
    if (!EvGiven) {Ev =  -0.14;}
    if (!elvsatGiven) {elvsat =  1.3e+7;}
    if (!hovsatGiven) {hovsat =  6.6e+6;}
    if (!narcoGiven) { narco =  0.047;}
    if (!narvaGiven) { narva =  0.047;}
    if (!CdonorGiven) { Cdonor =  0.0;}
    if (!CacceptorGiven) { Cacceptor =  2.0e+16;}
  }
  else if(tmpMat == "GAN")
  {
    if (!dielGiven) { diel=8.9;}
    if (!elmob0Given) { elmob0 = 1200;}
    if (!homob0Given) { homob0 = 11;}
    if (!NiGiven) { Ni = 1.0e-10;}
    //if (!EcGiven) {Ec =  1.45248;}
    //if (!EvGiven) {Ev =  -0.4085;}
    if (!elvsatGiven) {elvsat =  7.7e+6;}
    if (!hovsatGiven) {hovsat =  7.7e+6;}
    if (!CdonorGiven) { Cdonor =  5.0e+17;}
    if (!CacceptorGiven) { Cacceptor =  0.0;}
  }
  else if(tmpMat == "SI")
  {
    double Eg = MaterialSupport::bandgap(tmpMat,300.0);
    double chi = MaterialSupport::affin(tmpMat);

    if (!dielGiven) { diel=11.8;}
    if (!elmob0Given) { elmob0 = 850;} // check this! copied from ingap
    if (!homob0Given) { homob0 = 30;} // check this! copied from ingap
    if (!NiGiven) { Ni = 1.25e+10;}

    // Eg = 1.12 according to http://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html
    if (!EcGiven) {Ec =  1.12;} 
    if (!EvGiven) {Ev =  0.0;} 

    if (!elvsatGiven) {elvsat =  7.7e+6;} // check this!
    if (!hovsatGiven) {hovsat =  7.7e+6;} // check this!

    // band-gap narrowing: this is doping dependent, so should be computing.  
    // By default is zero for typical Si devices.
    //if (!narcoGiven) { narco =  0.013;} //check this!
    //if (!narvaGiven) { narva =  0.013;} //check this!
    if (!narcoGiven) { narco =  0.0;} //check this!
    if (!narvaGiven) { narva =  0.0;} //check this!

    if (!CdonorGiven) { Cdonor =  5.0e+17;}
    if (!CacceptorGiven) { Cacceptor =  0.0;}
    if (!augpnnGiven) { augpnn = 2.8e-31;}
    if (!augnppGiven) { augnpp = 9.9e-32;}
    if (!dirGiven) { dir = 0.0;}
    if (!electronThermalVGiven) {electronThermalV = 2.3e7;}
    if (!holeThermalVGiven) {holeThermalV = 1.9e7;}
    if (!defectReactionRadiusGiven) {defectReactionRadius = 1.4e-4;}
    if (!latticeConstantGiven) {latticeConstant = 5.00e-8;}  //This should be 5.43e-8

    // from Larry to match XPD
    if (!NcGiven){Nc = 2.86e19;}
    if (!NvGiven){Nv = 2.66e19;}
  }

  EcEff = Ec-narco;
  EvEff = Ev-narva;
}

//-----------------------------------------------------------------------------
// Function      : MaterialLayer::operator<<
// Purpose       : "<<" operator
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/28/12
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const MaterialLayer & ml)
{
  os << " Material Layer Data: name = " << ml.name ;
  os << " material = " << ml.material;

  os << " NX = " << ml.NX<<std::endl;
  os << " LX = " << ml.LX<<std::endl;
  os << " begin = " << ml.begin<<std::endl;
  os << " end = " << ml.end<<std::endl;
  os << " diel = " << ml.diel<<std::endl;
  os << " Ec = " << ml.Ec<<std::endl;
  os << " Ev = " << ml.Ev<<std::endl;
  os << " Cdonor = " << ml.Cdonor<<std::endl;
  os << " Cacceptor = " << ml.Cacceptor<<std::endl;
  os << " narco = " << ml.narco<<std::endl;
  os << " narva = " << ml.narva<<std::endl;
  os << " dnco = " << ml.dnco<<std::endl;
  os << " dnva = " << ml.dnva<<std::endl;
  os << " Nc = " << ml.Nc<<std::endl;
  os << " Nv = " << ml.Nv<<std::endl;
  os << " emass = " << ml.emass<<std::endl;
  os << " hmass = " << ml.hmass<<std::endl;
  os << " elmob0 = " << ml.elmob0<<std::endl;
  os << " elvsat = " << ml.elvsat<<std::endl;
  os << " eleo = " << ml.eleo<<std::endl;
  os << " homob0 = " << ml.homob0<<std::endl;
  os << " hovsat = " << ml.hovsat<<std::endl;
  os << " dir = " << ml.dir<<std::endl;
  os << " augnpp = " << ml.augnpp<<std::endl;
  os << " srh = " << ml.srh<<std::endl;
  os << " Ni = " << ml.Ni<<std::endl;
  os << " width = " << ml.width<<std::endl;
  os << " gradedLayerWidth = " << ml.gradedLayerWidth<<std::endl;
  os << " temperature = " << ml.temperature<<std::endl;

  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce
