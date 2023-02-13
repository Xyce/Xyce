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
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 07/27/2006
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_DEV_Specie_H
#define N_DEV_Specie_H

#include <string>
#include <vector>

#include <N_DEV_Const.h>
#include <iostream>

#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {

class Specie
{
public:
  Specie(std::string name, double diff_prefac, double act_energy, 
         int charge_state, int index=0)
    : Name(name),
      DiffusionPrefactor(diff_prefac),
      ActivationEnergy(act_energy),
      ChargeState(charge_state),
      carrierIndex(-1),
      sigma(0.0),
      hopLength(0.0),
      thermalVelocity(0.0),
      enhancedDiffusion(false),
      myIndex(index-3),
      BCCarrierCharge(0)
  {
    //Note that myIndex is being assigned with 3 subtracted from the count.  This is because the carriers, the
    //first two species in the list, are treated as "Constant" species and are not in the concentrations arrays
    //that hold defect concentrations that have to be used in the reactions.  Also, the count starts at 1 so an
    //there's an extra -- to get to the zero index
  } ;

  inline const std::string & getName() const {return (Name); };
  inline void setName(std::string &name) {Name = name;};
  inline int getChargeState() {return (ChargeState);};
  inline void setChargeState(int chargestate) {ChargeState=chargestate;};
  inline double getDiffPrefactor() {return(DiffusionPrefactor);} ; 
  inline void setDiffPrefactor(double p) {DiffusionPrefactor=p;};
  inline double getActEnergy() {return(ActivationEnergy);};    
  inline void setActEnergy(double Energy) {ActivationEnergy=Energy;};
  inline bool getEnhancedDiffusion() {return enhancedDiffusion;}
  void setBCEnhancedDiffusion(int cI, double sigma, int BCCC, double hopLength);
  inline int getBCCarrierIndex(){return carrierIndex;}
  inline double getBCHopLength(){return hopLength;}
  inline double getBCSigma(){return sigma;}
  inline int getMyIndex(){ return myIndex;}
  inline int getBCCarrierCharge(){ return BCCarrierCharge;}
  inline void setBCThermalVelocity(double TV){ thermalVelocity = TV;}
  template <class ScalarT>
    ScalarT getDiffusionCoefficient(ScalarT Temperature);
  template <class ScalarT>
    ScalarT getDiffusionCoefficient(ScalarT Temperature, 
                                    std::vector<ScalarT> &concs,std::vector<ScalarT> &constant_vec);
private:
  std::string Name;
  double DiffusionPrefactor;
  double ActivationEnergy;
  int ChargeState;
  //Some species have diffusion enhanced by carrier capture
  int carrierIndex;
  double sigma;
  double hopLength;
  double thermalVelocity;
  bool enhancedDiffusion;
  int myIndex;
  int BCCarrierCharge;
  
};

//-----------------------------------------------------------------------------
// Function      : Specie::getDiffusionCoefficient
// Purpose       : Accessor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
template <class ScalarT>
inline ScalarT Specie::getDiffusionCoefficient(ScalarT Temperature)
{

  return DiffusionPrefactor*exp(-ActivationEnergy/(CONSTboltz*Temperature/CONSTQ));
}


//-----------------------------------------------------------------------------
// Function      : Specie::getDiffusionCoefficient
// Purpose       : compute the diffusion coefficient included bourgoin corbett enhancement
// Special Notes : 
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/16/2014
//-----------------------------------------------------------------------------
 template <class ScalarT>
   ScalarT Specie::getDiffusionCoefficient(ScalarT Temperature, 
                                           std::vector<ScalarT> &concs,std::vector<ScalarT> &constant_vec)
   {
    
     ScalarT DF = DiffusionPrefactor*exp(-ActivationEnergy/(CONSTboltz*Temperature/CONSTQ));

     if(enhancedDiffusion)
       DF += sigma*thermalVelocity*constant_vec[carrierIndex]*hopLength*hopLength/6.0;

     return DF;

   }



} // namespace Device
} // namespace Xyce

#endif
