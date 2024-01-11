//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : This file implements the BSIM4 MOSFET model.  It
//                  is intended to be compatible with the Berkeley SPICE
//                  (3f5) version, BSIM4 version 4.8.2 and implements only
//                  those functions that differ from version to version.
//
//
// Special Notes  :
//
// Creator        : Tom Russo
//
// Creation Date  : 14 Sep 2022
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MOSFET1.h>
#include <N_DEV_MOSFET_B4.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_Math.h>
#include <N_ANP_NoiseData.h>

// ---------- BSIM4 constants ---------------
// Many constants are obtained from N_DEV_Const.h
// A few overrides (2) are below, to make this device as close
// to the original spice3f5 code as possible.

#define CONSTEPS0 8.85418e-12

#define Charge_q     (1.60219e-19) // electron charge, used in the
                                   // updateTemperature function instead of
                                   // of CONSTQ.  Having 2 constants for
                                   // the same quantity (with different
                                   // precision) doesn't make sense, but
                                   // it is what is in the spice3f5 bsim4.

#define CONSTKboQ 8.617087e-5  // another updateTemperature constant, which is
                               // not really necessary but I am keeping for
                               // compatibility with the spice3f5 bsim4.

#define CONSTvt0     (CONSTboltz * (27.0 +CONSTCtoK)/CONSTQ)

#define CONSTMM  3  // smooth coeff

#define DEXP(A,B,C) {                                                         \
        if (A > CONSTEXP_THRESHOLD) {                                         \
            B = CONSTMAX_EXP*(1.0+(A)-CONSTEXP_THRESHOLD);                    \
            C = CONSTMAX_EXP;                                                 \
        } else if (A < -CONSTEXP_THRESHOLD)  {                                \
            B = CONSTMIN_EXP;                                                 \
            C = 0;                                                            \
        } else   {                                                            \
            B = exp(A);                                                       \
            C = B;                                                            \
        }                                                                     \
    }


#define DELTA  1.0E-9
#define DEXP2(A,B) {                                                       \
        if (A > CONSTEXP_THRESHOLD) {                                      \
            B = CONSTMAX_EXP*(1.0+(A)-CONSTEXP_THRESHOLD);                 \
        } else if (A < -CONSTEXP_THRESHOLD)  {                             \
            B = CONSTMIN_EXP;                                              \
        } else   {                                                         \
            B = exp(A);                                                    \
        }                                                                  \
    }


namespace Xyce {
namespace Device {

namespace MOSFET_B4 {

//-----------------------------------------------------------------------------
// Function      : Instance::processParams4p82
// Purpose       : Version-specific process params for 4.8.2
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
//-----------------------------------------------------------------------------
bool Instance::processParams4p82_ ()
{
  double Rtot;

  // Set any non-constant parameter defaults:
  if (!RBDBgiven)
    rbdb = model_.rbdb;
  if (!RBSBgiven)
    rbsb = model_.rbsb;
  if (!RBPBgiven)
    rbpb = model_.rbpb;
  if (!RBPSgiven)
    rbps = model_.rbps;
  if (!RBPDgiven)
    rbpd = model_.rbpd;
  if (!XGWgiven)
    xgw = model_.xgw;
  if (!NGCONgiven)
    ngcon = model_.ngcon;
  if (!SDgiven)
    sd = 2.0 * model_.dmcg;
  if (!TEMPgiven)
    temp = getDeviceOptions().temp.getImmutableValue<double>();
  if (!drainAreaGiven)
    drainArea = getDeviceOptions().defad;
  if (!sourceAreaGiven)
    sourceArea = getDeviceOptions().defas;

  // Process instance (*M_iter) selectors, some
  // may override their global counterparts
  //
  if (!RBODYMODgiven)
  {
    rbodyMod = model_.rbodyMod;
  }
  else if ((rbodyMod != 0) && (rbodyMod != 1) && (rbodyMod != 2))
  {
    rbodyMod = model_.rbodyMod;
    UserWarning(*this) << "rbodyMod has been set to its global value: " << model_.rbodyMod;
  }

  if (!RGATEMODgiven)
  {
    rgateMod = model_.rgateMod;
  }
  else if ((rgateMod != 0) && (rgateMod != 1) && (rgateMod != 2) && (rgateMod != 3))
  {
    rgateMod = model_.rgateMod;
    UserWarning(*this) << "rgateMod has been set to its global value: " << model_.rgateMod;
  }

  if (!GEOMODgiven)
  {
    geoMod = model_.geoMod;
  }

  if (!RGEOMODgiven)
  {
    rgeoMod = model_.rgeoMod;
  }
  else if ((rgeoMod != 0) && (rgeoMod != 1))
  {
    rgeoMod = model_.rgeoMod;
    UserWarning(*this) << "rgeoMod has been set to its global value: " << model_.rgeoMod;
  }

  if (!TRNQSMODgiven)
  {
    trnqsMod = model_.trnqsMod;
  }
  else if ((trnqsMod != 0) && (trnqsMod != 1))
  {
    trnqsMod = model_.trnqsMod;
    UserWarning(*this) << "trnqsMod has been set to its global value: ";
  }

  if (!ACNQSMODgiven)
  {
    acnqsMod = model_.acnqsMod;
  }
  else if ((acnqsMod != 0) && (acnqsMod != 1))
  {
    acnqsMod = model_.acnqsMod;
    UserWarning(*this) << "acnqsMod has been set to its global value: ";
  }

  bool noiseAnalGiven=getSolverState().earlyNoiseFlag_;

  // process drain series resistance
  int createNode = 0;
  if ( (model_.rdsMod != 0) || (model_.tnoiMod == 1 && noiseAnalGiven))
  {
    createNode = 1;
  }
  else if (model_.sheetResistance > 0)
  {
    if (drainSquaresGiven && drainSquares > 0)
    {
      createNode = 1;
    }
    else if (!drainSquaresGiven && (rgeoMod != 0))
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
              w, model_.sheetResistance,
              DMCGeff, DMCIeff, DMDGeff, 0, Rtot);

      if(Rtot > 0)
      {
        createNode = 1;
      }
    }
  }

  if ( createNode != 0 )
  {
    drainMOSFET_B4Exists = true;
  }
  else
  {
    drainMOSFET_B4Exists = false;
  }

  // process source series resistance
  createNode = 0;
  if ( (model_.rdsMod != 0) || (model_.tnoiMod == 1 && noiseAnalGiven))
  {
    createNode = 1;
  }
  else if (model_.sheetResistance > 0)
  {
    if (sourceSquaresGiven && sourceSquares > 0)
    {
      createNode = 1;
    }
    else if (!sourceSquaresGiven && (rgeoMod != 0))
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
              w, model_.sheetResistance,
              DMCGeff, DMCIeff, DMDGeff, 1, Rtot);

      if(Rtot > 0)
      {
        createNode = 1;
      }
    }
  }

  if ( createNode != 0 )
  {
    sourceMOSFET_B4Exists = true;
  }
  else
  {
    sourceMOSFET_B4Exists = false;
  }

  // now set the temperature related stuff.
  updateTemperature(temp);


  // set up numIntVars:
  numIntVars = 0;

  if (drainMOSFET_B4Exists) ++ numIntVars;
  if (sourceMOSFET_B4Exists) ++ numIntVars;

  if (rgateMod == 1 || rgateMod == 2) ++numIntVars;
  else if (rgateMod == 3) numIntVars+=2;

  if ( trnqsMod ) ++numIntVars;
  if ( rbodyMod ) numIntVars+=3;

  if (icVBSGiven) ++numIntVars;
  if (icVDSGiven) ++numIntVars;
  if (icVGSGiven) ++numIntVars;

  // set up numStateVars
  numStateVars = 3;
  setNumStoreVars(22);

  if (rgateMod == 3)
  {
    numStateVars += 1;
  }

  // parasitic capacitors:
  if (rbodyMod)
  {
    numStateVars += 2;
  }

  if (trnqsMod)
  {
    numStateVars += 2;
  }

  // If there are any time dependent parameters, set their values at for
  // the current time.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature4p82_
// Purpose       : This updates all the instance-owned paramters which
//                 are temperature dependent.
//
// Special Notes : Annoyingly, some model-owned parameters need to be
//                 tweaked here because of how the SPICE code is set up.
//
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
bool Instance::updateTemperature4p82_ (const double & temp_tmp)
{
  std::string msg="";

  double tmp(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), Eg(0.0), Eg0(0.0), ni,epssub;
  double T0(0.0), T1(0.0);
  double T2(0.0), T3(0.0), T4(0.0), T5(0.0), T6(0.0), T7(0.0), T8(0.0), T9(0.0), Lnew(0.0), Wnew(0.0);
  double delTemp(0.0), TRatio(0.0), Inv_L(0.0), Inv_W(0.0), Inv_LW(0.0), Vtm0, Tnom(0.0);
  double dumPs(0.0), dumPd(0.0), dumAs(0.0), dumAd(0.0), PowWeffWr(0.0);
  double Nvtms(0.0), Nvtmd(0.0), SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double T10(0.0), T11(0.0);
  double Inv_saref(0.0), Inv_sbref(0.0), Inv_sa(0.0), Inv_sb(0.0), rho(0.0), Ldrn(0.0), dvth0_lod(0.0);
  double W_tmp(0.0), Inv_ODeff(0.0), OD_offset(0.0), dk2_lod(0.0), deta0_lod(0.0);
  double lnl(0.0), lnw(0.0), lnnf(0.0), rbpbx(0.0), rbpby(0.0), rbsbx(0.0), rbsby(0.0), rbdbx(0.0), rbdby(0.0),bodymode(0.0);
  double kvsat(0.0), wlod(0.0), sceff(0.0), Wdrn(0.0);
  double V0, lt1, ltw, Theta0, Delt_vth, TempRatio, Vth_NarrowW, Lpe_Vb;
  //double Vth; // converted to instance variable
  double n, n0, Vgsteff, Vgs_eff, toxpf, toxpi, Tcen, toxe, epsrox, vddeot;
  double vtfbphi2eot, phieot, TempRatioeot, Vtm0eot, Vtmeot,vbieot;

  int niter;

  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << "Instance::updateTemperature\n";
    Xyce::dout() << "name = " << getName() << std::endl;
  }

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;

  Tnom = model_.tnom;
  TRatio = temp/Tnom;

  ///////////////////////////////////////////////////////////////////////////////
  // Model-specific stuff:
  // This is kludgey - the model-specific stuff should be handled in a model function.
  // Some of this used to be in the model class's processParams, but was
  // moved back here because it makes updating from new spice BSIM4 code
  // less painful, even though it makes no sense here because an instance
  // function is modifying model data structures.
  //
  // This stuff all comes from the SPICE code's "BSIM4temp" function which
  // is a mélange of instance and model setting.  While all the model setting
  // stuff technically should belong in Xyce's Model::processParams,
  // keeping it here, where other code from BSIM4temp is taken, is why having
  // it here makes things easier when updating or adding new versions.

  if (model_.mtrlMod == 0)
  {
    if ((model_.toxeGiven) && (model_.toxpGiven) &&
        (model_.dtoxGiven) &&
        (model_.toxe != (model_.toxp +model_.dtox)))
    {
      UserWarning(*this) << "toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored";
    }
    else if ((model_.toxeGiven) && (!model_.toxpGiven))
    {
      model_.toxp = model_.toxe - model_.dtox;
    }
    else if ((!model_.toxeGiven) && (model_.toxpGiven))
    {
      model_.toxe = model_.toxp + model_.dtox;
      if (!model_.toxmGiven)
      {
        model_.toxm = model_.toxe;
      }
    }
    if (!model_.cfGiven)
      model_.cf = 2.0 * model_.epsrox * CONSTEPS0/M_PI
        * log(1.0 + 0.4E-6 / model_.toxe);
  }
  else if (model_.mtrlCompatMod != 0)
  {
    T0 = model_.epsrox / 3.9;
    if ((model_.eotGiven) && (model_.toxpGiven) && (model_.dtoxGiven)
        && (fabs(model_.eot * T0 - (model_.toxp + model_.dtox)) > 1.0e-20))
    {
      UserWarning(*this) << "eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored";
    }
    else if ((model_.eotGiven) && (!model_.toxpGiven))
    {
      model_.toxp = T0 * model_.eot - model_.dtox;
    }
    else if ((!model_.eotGiven) && (model_.toxpGiven))
    {
      model_.eot = (model_.toxp + model_.dtox) / T0;
      if (!model_.toxmGiven)
      {
        model_.toxm = model_.eot;
      }
    }
  }

  if (model_.mtrlMod)
  {
    epsrox = 3.9;
    toxe = model_.eot;
    epssub = CONSTEPS0 * model_.epsrsub;
  }
  else
  {
    epsrox = model_.epsrox;
    toxe = model_.toxe;
    epssub = CONSTEPSSI;
  }


  model_.coxe = epsrox * CONSTEPS0 / toxe;
  if (model_.mtrlMod == 0 || model_.mtrlCompatMod != 0)
    model_.coxp = model_.epsrox * CONSTEPS0 / model_.toxp;

  if (!model_.cgdoGiven)
  {
    if (model_.dlcGiven && (model_.dlc > 0.0))
    {
      model_.cgdo = model_.dlc * model_.coxe - model_.cgdl ;
    }
    else
    {
      model_.cgdo = 0.6 * model_.xj * model_.coxe;
    }
  }

  if (!model_.cgsoGiven)
  {
    if (model_.dlcGiven && (model_.dlc > 0.0))
    {
      model_.cgso = model_.dlc * model_.coxe - model_.cgsl ;
    }
    else
    {
      model_.cgso = 0.6 * model_.xj * model_.coxe;
    }
  }
  if (!model_.cgboGiven)
  {
    model_.cgbo = 2.0 * model_.dwc * model_.coxe;
  }

  model_.vcrit = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
  model_.factor1 = sqrt(epssub / (epsrox * CONSTEPS0) * toxe);

  Vtm0 = model_.vtm0 = CONSTKboQ * Tnom;

  if (model_.mtrlMod == 0)
  {
    Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
    ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15)
      * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
  }
  else
  {
    Eg0 = model_.bg0sub - model_.tbgasub * Tnom * Tnom
      / (Tnom + model_.tbgbsub);
    T0 =  model_.bg0sub - model_.tbgasub * 90090.0225
      / (300.15 + model_.tbgbsub);
    ni = model_.ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15)
      * exp((T0 - Eg0) / (2.0 * Vtm0));
  }

  model_.Eg0 = Eg0;
  model_.vtm = CONSTKboQ * temp;
  if (model_.mtrlMod == 0)
  {
    Eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
  }
  else
  {
    Eg = model_.bg0sub - model_.tbgasub * temp * temp
      / (temp + model_.tbgbsub);
  }

  if (temp != Tnom)
  {
    T0 = Eg0 / Vtm0 - Eg / model_.vtm;
    T1 = log(temp / Tnom);
    T2 = T0 + model_.SjctTempExponent * T1;
    T3 = exp(T2 / model_.SjctEmissionCoeff);
    model_.SjctTempSatCurDensity = model_.SjctSatCurDensity * T3;

    model_.SjctSidewallTempSatCurDensity
         = model_.SjctSidewallSatCurDensity * T3;

    model_.SjctGateSidewallTempSatCurDensity
                     = model_.SjctGateSidewallSatCurDensity * T3;

    T2 = T0 + model_.DjctTempExponent * T1;
    T3 = exp(T2 / model_.DjctEmissionCoeff);

    model_.DjctTempSatCurDensity = model_.DjctSatCurDensity * T3;

    model_.DjctSidewallTempSatCurDensity
           = model_.DjctSidewallSatCurDensity * T3;

    model_.DjctGateSidewallTempSatCurDensity
           = model_.DjctGateSidewallSatCurDensity * T3;

  }
  else
  {
   model_.SjctTempSatCurDensity = model_.SjctSatCurDensity;
   model_.SjctSidewallTempSatCurDensity = model_.SjctSidewallSatCurDensity;

   model_.SjctGateSidewallTempSatCurDensity
                      = model_.SjctGateSidewallSatCurDensity;

   model_.DjctTempSatCurDensity = model_.DjctSatCurDensity;

   model_.DjctSidewallTempSatCurDensity
                      = model_.DjctSidewallSatCurDensity;

   model_.DjctGateSidewallTempSatCurDensity
                      = model_.DjctGateSidewallSatCurDensity;
  }

  if (model_.SjctTempSatCurDensity < 0.0)
     model_.SjctTempSatCurDensity = 0.0;

  if (model_.SjctSidewallTempSatCurDensity < 0.0)
     model_.SjctSidewallTempSatCurDensity = 0.0;

  if (model_.SjctGateSidewallTempSatCurDensity < 0.0)
     model_.SjctGateSidewallTempSatCurDensity = 0.0;

  if (model_.DjctTempSatCurDensity < 0.0)
     model_.DjctTempSatCurDensity = 0.0;

  if (model_.DjctSidewallTempSatCurDensity < 0.0)
     model_.DjctSidewallTempSatCurDensity = 0.0;

  if (model_.DjctGateSidewallTempSatCurDensity < 0.0)
     model_.DjctGateSidewallTempSatCurDensity = 0.0;

  // Temperature dependence of D/B and S/B diode capacitance begins
  delTemp = temp - Tnom;
  T0 = model_.tcj * delTemp;
  if (T0 >= -1.0)
  {   model_.SunitAreaTempJctCap = model_.SunitAreaJctCap *(1.0 + T0); //bug_fix -JX
     model_.DunitAreaTempJctCap = model_.DunitAreaJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitAreaJctCap > 0.0)
   {
     model_.SunitAreaTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjs to be negative. Cjs is clamped to zero";

   }
   if (model_.DunitAreaJctCap > 0.0)
   {
     model_.DunitAreaTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjd to be negative. Cjd is clamped to zero";
   }
  }

  T0 = model_.tcjsw * delTemp;
  if (model_.SunitLengthSidewallJctCap < 0.0)
  {
   model_.SunitLengthSidewallJctCap = 0.0;
   UserWarning(*this) << "CJSWS is negative. Cjsws is clamped to zero";
  }
  if (model_.DunitLengthSidewallJctCap < 0.0)
  {
   model_.DunitLengthSidewallJctCap = 0.0;
   UserWarning(*this) << "CJSWD is negative. Cjswd is clamped to zero";
  }
  if (T0 >= -1.0)
  {
   model_.SunitLengthSidewallTempJctCap = model_.SunitLengthSidewallJctCap *(1.0 + T0);
   model_.DunitLengthSidewallTempJctCap = model_.DunitLengthSidewallJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitLengthSidewallJctCap > 0.0)
   {
     model_.SunitLengthSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero";
   }
   if (model_.DunitLengthSidewallJctCap > 0.0)
   {
     model_.DunitLengthSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero";
   }
  }

  T0 = model_.tcjswg * delTemp;
  if (T0 >= -1.0)
  {
   model_.SunitLengthGateSidewallTempJctCap = model_.SunitLengthGateSidewallJctCap *(1.0 + T0);
   model_.DunitLengthGateSidewallTempJctCap = model_.DunitLengthGateSidewallJctCap *(1.0 + T0);
  }
  else
  {
   if (model_.SunitLengthGateSidewallJctCap > 0.0)
   {
     model_.SunitLengthGateSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero";
   }
   if (model_.DunitLengthGateSidewallJctCap > 0.0)
   {
     model_.DunitLengthGateSidewallTempJctCap = 0.0;
     UserWarning(*this) << "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero";
   }
  }

  model_.PhiBS = model_.SbulkJctPotential - model_.tpb * delTemp;

  if (model_.PhiBS < 0.01)
  {
   model_.PhiBS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01";
  }

  model_.PhiBD = model_.DbulkJctPotential - model_.tpb * delTemp;
  if (model_.PhiBD < 0.01)
  {
   model_.PhiBD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01";
  }

  model_.PhiBSWS = model_.SsidewallJctPotential - model_.tpbsw * delTemp;
  if (model_.PhiBSWS <= 0.01)
  {
   model_.PhiBSWS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01";
  }

  model_.PhiBSWD = model_.DsidewallJctPotential - model_.tpbsw * delTemp;
  if (model_.PhiBSWD <= 0.01)
  {
   model_.PhiBSWD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01";
  }

  model_.PhiBSWGS = model_.SGatesidewallJctPotential - model_.tpbswg * delTemp;
  if (model_.PhiBSWGS <= 0.01)
  {
   model_.PhiBSWGS = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01";
  }

  model_.PhiBSWGD = model_.DGatesidewallJctPotential - model_.tpbswg * delTemp;
  if (model_.PhiBSWGD <= 0.01)
  {
   model_.PhiBSWGD = 0.01;
   UserWarning(*this) << "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01";
  } // End of junction capacitance


  if (model_.ijthdfwd <= 0.0)
  {
   model_.ijthdfwd = 0.0;
   UserWarning(*this) << "Ijthdfwd reset to " << model_.ijthdfwd;
  }
  if (model_.ijthsfwd <= 0.0)
  {
   model_.ijthsfwd = 0.0;
   UserWarning(*this) << "Ijthsfwd reset to " << model_.ijthsfwd;
  }
  if (model_.ijthdrev <= 0.0)
  {
   model_.ijthdrev = 0.0;
   UserWarning(*this) << "Ijthdrev reset to " << model_.ijthdrev;
  }
  if (model_.ijthsrev <= 0.0)
  {
   model_.ijthsrev = 0.0;
   UserWarning(*this) << "Ijthsrev reset to " << model_.ijthsrev;
  }

  if ((model_.xjbvd <= 0.0) && (model_.dioMod == 2))
  {
   model_.xjbvd = 0.0;
   UserWarning(*this) << "Xjbvd reset to " << model_.xjbvd;
  }
  else if ((model_.xjbvd < 0.0) && (model_.dioMod == 0))
  {
   model_.xjbvd = 0.0;
   UserWarning(*this) << "Xjbvd reset to " << model_.xjbvd;
  }

  if (model_.bvd <= 0.0)
  {
   model_.bvd = 0.0;
   UserWarning(*this) << "BVD reset to " << model_.bvd;
  }

  if ((model_.xjbvs <= 0.0) && (model_.dioMod == 2))
  {
   model_.xjbvs = 0.0;
   UserWarning(*this) << "Xjbvs reset to " << model_.xjbvs;
  }
  else if ((model_.xjbvs < 0.0) && (model_.dioMod == 0))
  {
   model_.xjbvs = 0.0;
   UserWarning(*this) << "Xjbvs reset to " << model_.xjbvs;
  }

  if (model_.bvs <= 0.0)
  {
   model_.bvs = 0.0;
   UserWarning(*this) << "BVS reset to " << model_.bvs;
  }


  ///////////////////////////////////////////////////////////////////////////////
  // Instance stuff:
  // (loop through all the instances of the model)

  // stress effect
  Ldrn = l;
  Wdrn = w / nf;

  // This next block determines whether  or not to use a previously allocated
  // set of size dependent parameters. These are stored in a list that is
  // owned by the model.  If the values for length and width match those of
  // a previously allocated set, then use the old set.  If not, allocate a new set.

  std::list<SizeDependParam*>::iterator it_dpL =
    model_.sizeDependParamList.begin();
  std::list<SizeDependParam*>::iterator end_dpL =
    model_.sizeDependParamList.end();

  paramPtr = NULL;

  for( ; it_dpL != end_dpL; ++it_dpL )
  {
    if( ((*it_dpL)->Length  == l)
     && ((*it_dpL)->Width   == w)
     && ((*it_dpL)->NFinger == nf)
     && ((*it_dpL)->referenceTemperature == temp_tmp))
    {
      paramPtr = (*it_dpL);
    }
  }

  // This was inside of the "Size_Not_Found" if-statement, but that
  // won't work here - it winds up being uninitialized whenever the
  // size pointer is found
  Lnew = l  + model_.xl ;
  Wnew = w / nf + model_.xw;

  if ( paramPtr != NULL )
  {
  }
  else
  {
    paramPtr = new SizeDependParam ();

    model_.sizeDependParamList.push_back( paramPtr );
    paramPtr->referenceTemperature = temp_tmp;

    //paramPtr->pNext = NULL;

    paramPtr->Length = l;
    paramPtr->Width = w;
    paramPtr->NFinger = nf;
    //Lnew = l  + model_.xl ;
    //Wnew = w / nf + model_.xw;

    T0 = pow(Lnew, model_.Lln);
    T1 = pow(Wnew, model_.Lwn);
    tmp1 = model_.Ll / T0 + model_.Lw / T1
         + model_.Lwl / (T0 * T1);
    paramPtr->dl = model_.Lint + tmp1;
    tmp2 = model_.Llc / T0 + model_.Lwc / T1
         + model_.Lwlc / (T0 * T1);
    paramPtr->dlc = model_.dlc + tmp2;

    T2 = pow(Lnew, model_.Wln);
    T3 = pow(Wnew, model_.Wwn);
    tmp1 = model_.Wl / T2 + model_.Ww / T3
         + model_.Wwl / (T2 * T3);
    paramPtr->dw = model_.Wint + tmp1;
    tmp2 = model_.Wlc / T2 + model_.Wwc / T3
         + model_.Wwlc / (T2 * T3);
    paramPtr->dwc = model_.dwc + tmp2;
    paramPtr->dwj = model_.dwj + tmp2;

    paramPtr->leff = Lnew - 2.0 * paramPtr->dl;

    if (paramPtr->leff <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel length <= 0";
    }

    paramPtr->weff = Wnew - 2.0 * paramPtr->dw;
    if (paramPtr->weff <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width <= 0";
    }

    paramPtr->leffCV = Lnew - 2.0 * paramPtr->dlc;
    if (paramPtr->leffCV <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel length for C-V <= 0";
    }

    paramPtr->weffCV = Wnew - 2.0 * paramPtr->dwc;
    if (paramPtr->weffCV <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width for C-V <= 0";
    }

    paramPtr->weffCJ = Wnew - 2.0 * paramPtr->dwj;
    if (paramPtr->weffCJ <= 0.0)
    {
      UserError(*this) << "mosfet " << getName() << " model " << model_.getName()
                        << "  Effective channel width for S/D junctions <= 0";
    }


    if (model_.binUnit == 1)
    {
      Inv_L = 1.0e-6 / paramPtr->leff;
      Inv_W = 1.0e-6 / paramPtr->weff;
      Inv_LW = 1.0e-12 / (paramPtr->leff * paramPtr->weff);
    }
    else
    {
      Inv_L = 1.0 / paramPtr->leff;
      Inv_W = 1.0 / paramPtr->weff;
      Inv_LW = 1.0 / (paramPtr->leff * paramPtr->weff);
    }
    paramPtr->cdsc = model_.cdsc
                  + model_.lcdsc * Inv_L
                  + model_.wcdsc * Inv_W
                  + model_.pcdsc * Inv_LW;
    paramPtr->cdscb = model_.cdscb
                  + model_.lcdscb * Inv_L
                  + model_.wcdscb * Inv_W
                  + model_.pcdscb * Inv_LW;

    paramPtr->cdscd = model_.cdscd
                  + model_.lcdscd * Inv_L
                  + model_.wcdscd * Inv_W
                  + model_.pcdscd * Inv_LW;

    paramPtr->cit = model_.cit
                  + model_.lcit * Inv_L
                  + model_.wcit * Inv_W
                  + model_.pcit * Inv_LW;
    paramPtr->nfactor = model_.nfactor
                  + model_.lnfactor * Inv_L
                  + model_.wnfactor * Inv_W
                  + model_.pnfactor * Inv_LW;
    paramPtr->tnfactor = model_.tnfactor
                  + model_.ltnfactor * Inv_L
                  + model_.wtnfactor * Inv_W
                  + model_.ptnfactor * Inv_LW;
    paramPtr->xj = model_.xj
                  + model_.lxj * Inv_L
                  + model_.wxj * Inv_W
                  + model_.pxj * Inv_LW;
    paramPtr->vsat = model_.vsat
                  + model_.lvsat * Inv_L
                  + model_.wvsat * Inv_W
                  + model_.pvsat * Inv_LW;
    paramPtr->at = model_.at
                  + model_.lat * Inv_L
                  + model_.wat * Inv_W
                  + model_.pat * Inv_LW;
    paramPtr->a0 = model_.a0
                  + model_.la0 * Inv_L
                  + model_.wa0 * Inv_W
                  + model_.pa0 * Inv_LW;

    paramPtr->ags = model_.ags
                  + model_.lags * Inv_L
                  + model_.wags * Inv_W
                  + model_.pags * Inv_LW;

    paramPtr->a1 = model_.a1
                  + model_.la1 * Inv_L
                  + model_.wa1 * Inv_W
                  + model_.pa1 * Inv_LW;
    paramPtr->a2 = model_.a2
                  + model_.la2 * Inv_L
                  + model_.wa2 * Inv_W
                  + model_.pa2 * Inv_LW;
    paramPtr->keta = model_.keta
                  + model_.lketa * Inv_L
                  + model_.wketa * Inv_W
                  + model_.pketa * Inv_LW;
    paramPtr->nsub = model_.nsub
                  + model_.lnsub * Inv_L
                  + model_.wnsub * Inv_W
                  + model_.pnsub * Inv_LW;
    paramPtr->ndep = model_.ndep
                  + model_.lndep * Inv_L
                  + model_.wndep * Inv_W
                  + model_.pndep * Inv_LW;
    paramPtr->nsd = model_.nsd
                     + model_.lnsd * Inv_L
                     + model_.wnsd * Inv_W
                     + model_.pnsd * Inv_LW;
    paramPtr->phin = model_.phin
                      + model_.lphin * Inv_L
                      + model_.wphin * Inv_W
                      + model_.pphin * Inv_LW;
    paramPtr->ngate = model_.ngate
                      + model_.lngate * Inv_L
                      + model_.wngate * Inv_W
                      + model_.pngate * Inv_LW;
    paramPtr->gamma1 = model_.gamma1
                      + model_.lgamma1 * Inv_L
                      + model_.wgamma1 * Inv_W
                      + model_.pgamma1 * Inv_LW;
    paramPtr->gamma2 = model_.gamma2
                      + model_.lgamma2 * Inv_L
                      + model_.wgamma2 * Inv_W
                      + model_.pgamma2 * Inv_LW;
    paramPtr->vbx = model_.vbx
                      + model_.lvbx * Inv_L
                      + model_.wvbx * Inv_W
                      + model_.pvbx * Inv_LW;
    paramPtr->vbm = model_.vbm
                      + model_.lvbm * Inv_L
                      + model_.wvbm * Inv_W
                      + model_.pvbm * Inv_LW;
    paramPtr->xt = model_.xt
                      + model_.lxt * Inv_L
                      + model_.wxt * Inv_W
                      + model_.pxt * Inv_LW;
    paramPtr->vfb = model_.vfb
                     + model_.lvfb * Inv_L
                     + model_.wvfb * Inv_W
                     + model_.pvfb * Inv_LW;
    paramPtr->k1 = model_.k1
                      + model_.lk1 * Inv_L
                      + model_.wk1 * Inv_W
                      + model_.pk1 * Inv_LW;
    paramPtr->kt1 = model_.kt1
                      + model_.lkt1 * Inv_L
                      + model_.wkt1 * Inv_W
                      + model_.pkt1 * Inv_LW;
    paramPtr->kt1l = model_.kt1l
                      + model_.lkt1l * Inv_L
                      + model_.wkt1l * Inv_W
                      + model_.pkt1l * Inv_LW;
    paramPtr->k2 = model_.k2
                      + model_.lk2 * Inv_L
                      + model_.wk2 * Inv_W
                      + model_.pk2 * Inv_LW;
    paramPtr->kt2 = model_.kt2
                      + model_.lkt2 * Inv_L
                      + model_.wkt2 * Inv_W
                      + model_.pkt2 * Inv_LW;
    paramPtr->k3 = model_.k3
                      + model_.lk3 * Inv_L
                      + model_.wk3 * Inv_W
                      + model_.pk3 * Inv_LW;
    paramPtr->k3b = model_.k3b
                      + model_.lk3b * Inv_L
                      + model_.wk3b * Inv_W
                      + model_.pk3b * Inv_LW;
    paramPtr->w0 = model_.w0
                      + model_.lw0 * Inv_L
                      + model_.ww0 * Inv_W
                      + model_.pw0 * Inv_LW;
    paramPtr->lpe0 = model_.lpe0
                      + model_.llpe0 * Inv_L
                      + model_.wlpe0 * Inv_W
                      + model_.plpe0 * Inv_LW;
    paramPtr->lpeb = model_.lpeb
                      + model_.llpeb * Inv_L
                      + model_.wlpeb * Inv_W
                      + model_.plpeb * Inv_LW;
    paramPtr->dvtp0 = model_.dvtp0
                       + model_.ldvtp0 * Inv_L
                       + model_.wdvtp0 * Inv_W
                       + model_.pdvtp0 * Inv_LW;
    paramPtr->dvtp1 = model_.dvtp1
                       + model_.ldvtp1 * Inv_L
                       + model_.wdvtp1 * Inv_W
                       + model_.pdvtp1 * Inv_LW;
    paramPtr->dvtp2 = model_.dvtp2
                       + model_.ldvtp2 * Inv_L
                       + model_.wdvtp2 * Inv_W
                       + model_.pdvtp2 * Inv_LW;
    paramPtr->dvtp3 = model_.dvtp3
                       + model_.ldvtp3 * Inv_L
                       + model_.wdvtp3 * Inv_W
                       + model_.pdvtp3 * Inv_LW;
    paramPtr->dvtp4 = model_.dvtp4
                       + model_.ldvtp4 * Inv_L
                       + model_.wdvtp4 * Inv_W
                       + model_.pdvtp4 * Inv_LW;
    paramPtr->dvtp5 = model_.dvtp5
                       + model_.ldvtp5 * Inv_L
                       + model_.wdvtp5 * Inv_W
                       + model_.pdvtp5 * Inv_LW;
    paramPtr->dvt0 = model_.dvt0
                      + model_.ldvt0 * Inv_L
                      + model_.wdvt0 * Inv_W
                      + model_.pdvt0 * Inv_LW;
    paramPtr->dvt1 = model_.dvt1
                      + model_.ldvt1 * Inv_L
                      + model_.wdvt1 * Inv_W
                      + model_.pdvt1 * Inv_LW;
    paramPtr->dvt2 = model_.dvt2
                      + model_.ldvt2 * Inv_L
                      + model_.wdvt2 * Inv_W
                      + model_.pdvt2 * Inv_LW;
    paramPtr->dvt0w = model_.dvt0w
                      + model_.ldvt0w * Inv_L
                      + model_.wdvt0w * Inv_W
                      + model_.pdvt0w * Inv_LW;
    paramPtr->dvt1w = model_.dvt1w
                      + model_.ldvt1w * Inv_L
                      + model_.wdvt1w * Inv_W
                      + model_.pdvt1w * Inv_LW;
    paramPtr->dvt2w = model_.dvt2w
                      + model_.ldvt2w * Inv_L
                      + model_.wdvt2w * Inv_W
                      + model_.pdvt2w * Inv_LW;
    paramPtr->drout = model_.drout
                      + model_.ldrout * Inv_L
                      + model_.wdrout * Inv_W
                      + model_.pdrout * Inv_LW;
    paramPtr->dsub = model_.dsub
                      + model_.ldsub * Inv_L
                      + model_.wdsub * Inv_W
                      + model_.pdsub * Inv_LW;
    paramPtr->vth0 = model_.vth0
                      + model_.lvth0 * Inv_L
                      + model_.wvth0 * Inv_W
                      + model_.pvth0 * Inv_LW;
    paramPtr->ua = model_.ua
                      + model_.lua * Inv_L
                      + model_.wua * Inv_W
                      + model_.pua * Inv_LW;
    paramPtr->ua1 = model_.ua1
                      + model_.lua1 * Inv_L
                      + model_.wua1 * Inv_W
                      + model_.pua1 * Inv_LW;
    paramPtr->ub = model_.ub
                      + model_.lub * Inv_L
                      + model_.wub * Inv_W
                      + model_.pub * Inv_LW;
    paramPtr->ub1 = model_.ub1
                      + model_.lub1 * Inv_L
                      + model_.wub1 * Inv_W
                      + model_.pub1 * Inv_LW;
    paramPtr->uc = model_.uc
                      + model_.luc * Inv_L
                      + model_.wuc * Inv_W
                      + model_.puc * Inv_LW;
    paramPtr->uc1 = model_.uc1
                      + model_.luc1 * Inv_L
                      + model_.wuc1 * Inv_W
                      + model_.puc1 * Inv_LW;
    paramPtr->ud = model_.ud
                      + model_.lud * Inv_L
                      + model_.wud * Inv_W
                      + model_.pud * Inv_LW;
    paramPtr->ud1 = model_.ud1
                      + model_.lud1 * Inv_L
                      + model_.wud1 * Inv_W
                      + model_.pud1 * Inv_LW;
    paramPtr->up = model_.up
                      + model_.lup * Inv_L
                      + model_.wup * Inv_W
                      + model_.pup * Inv_LW;
    paramPtr->lp = model_.lp
                      + model_.llp * Inv_L
                      + model_.wlp * Inv_W
                      + model_.plp * Inv_LW;
    paramPtr->eu = model_.eu
                    + model_.leu * Inv_L
                    + model_.weu * Inv_W
                    + model_.peu * Inv_LW;
    paramPtr->u0 = model_.u0
                    + model_.lu0 * Inv_L
                    + model_.wu0 * Inv_W
                    + model_.pu0 * Inv_LW;
    paramPtr->ute = model_.ute
                    + model_.lute * Inv_L
                    + model_.wute * Inv_W
                    + model_.pute * Inv_LW;
    paramPtr->ucs = model_.ucs
                    + model_.lucs * Inv_L
                    + model_.wucs * Inv_W
                    + model_.pucs * Inv_LW;
    paramPtr->ucste = model_.ucste
                    + model_.lucste * Inv_L
                    + model_.wucste * Inv_W
                    + model_.pucste * Inv_LW;
    paramPtr->voff = model_.voff
                    + model_.lvoff * Inv_L
                    + model_.wvoff * Inv_W
                    + model_.pvoff * Inv_LW;
    paramPtr->tvoff = model_.tvoff
                    + model_.ltvoff * Inv_L
                    + model_.wtvoff * Inv_W
                    + model_.ptvoff * Inv_LW;
    paramPtr->minv = model_.minv
                      + model_.lminv * Inv_L
                      + model_.wminv * Inv_W
                      + model_.pminv * Inv_LW;
    paramPtr->minvcv = model_.minvcv
                      + model_.lminvcv * Inv_L
                      + model_.wminvcv * Inv_W
                      + model_.pminvcv * Inv_LW;
    paramPtr->fprout = model_.fprout
                       + model_.lfprout * Inv_L
                       + model_.wfprout * Inv_W
                       + model_.pfprout * Inv_LW;
    paramPtr->pdits = model_.pdits
                       + model_.lpdits * Inv_L
                       + model_.wpdits * Inv_W
                       + model_.ppdits * Inv_LW;
    paramPtr->pditsd = model_.pditsd
                        + model_.lpditsd * Inv_L
                        + model_.wpditsd * Inv_W
                        + model_.ppditsd * Inv_LW;
    paramPtr->delta = model_.delta
                        + model_.ldelta * Inv_L
                        + model_.wdelta * Inv_W
                        + model_.pdelta * Inv_LW;
    paramPtr->rdsw = model_.rdsw
                        + model_.lrdsw * Inv_L
                        + model_.wrdsw * Inv_W
                        + model_.prdsw * Inv_LW;
    paramPtr->rdw = model_.rdw
                      + model_.lrdw * Inv_L
                      + model_.wrdw * Inv_W
                      + model_.prdw * Inv_LW;
    paramPtr->rsw = model_.rsw
                      + model_.lrsw * Inv_L
                      + model_.wrsw * Inv_W
                      + model_.prsw * Inv_LW;
    paramPtr->prwg = model_.prwg
                      + model_.lprwg * Inv_L
                      + model_.wprwg * Inv_W
                      + model_.pprwg * Inv_LW;
    paramPtr->prwb = model_.prwb
                      + model_.lprwb * Inv_L
                      + model_.wprwb * Inv_W
                      + model_.pprwb * Inv_LW;
    paramPtr->prt = model_.prt
                      + model_.lprt * Inv_L
                      + model_.wprt * Inv_W
                      + model_.pprt * Inv_LW;
    paramPtr->eta0 = model_.eta0
                      + model_.leta0 * Inv_L
                      + model_.weta0 * Inv_W
                      + model_.peta0 * Inv_LW;
    paramPtr->teta0 = model_.teta0
                      + model_.lteta0 * Inv_L
                      + model_.wteta0 * Inv_W
                      + model_.pteta0 * Inv_LW;
    paramPtr->tvoffcv = model_.tvoffcv    /* v4.8.0  */
                      + model_.ltvoffcv * Inv_L
                      + model_.wtvoffcv * Inv_W
                      + model_.ptvoffcv * Inv_LW;
    paramPtr->etab = model_.etab
                      + model_.letab * Inv_L
                      + model_.wetab * Inv_W
                      + model_.petab * Inv_LW;
    paramPtr->pclm = model_.pclm
                      + model_.lpclm * Inv_L
                      + model_.wpclm * Inv_W
                      + model_.ppclm * Inv_LW;
    paramPtr->pdibl1 = model_.pdibl1
                      + model_.lpdibl1 * Inv_L
                      + model_.wpdibl1 * Inv_W
                      + model_.ppdibl1 * Inv_LW;
    paramPtr->pdibl2 = model_.pdibl2
                      + model_.lpdibl2 * Inv_L
                      + model_.wpdibl2 * Inv_W
                      + model_.ppdibl2 * Inv_LW;
    paramPtr->pdiblb = model_.pdiblb
                      + model_.lpdiblb * Inv_L
                      + model_.wpdiblb * Inv_W
                      + model_.ppdiblb * Inv_LW;
    paramPtr->pscbe1 = model_.pscbe1
                      + model_.lpscbe1 * Inv_L
                      + model_.wpscbe1 * Inv_W
                      + model_.ppscbe1 * Inv_LW;
    paramPtr->pscbe2 = model_.pscbe2
                      + model_.lpscbe2 * Inv_L
                      + model_.wpscbe2 * Inv_W
                      + model_.ppscbe2 * Inv_LW;
    paramPtr->pvag = model_.pvag
                      + model_.lpvag * Inv_L
                      + model_.wpvag * Inv_W
                      + model_.ppvag * Inv_LW;
    paramPtr->wr = model_.wr
                      + model_.lwr * Inv_L
                      + model_.wwr * Inv_W
                      + model_.pwr * Inv_LW;
    paramPtr->dwg = model_.dwg
                      + model_.ldwg * Inv_L
                      + model_.wdwg * Inv_W
                      + model_.pdwg * Inv_LW;
    paramPtr->dwb = model_.dwb
                      + model_.ldwb * Inv_L
                      + model_.wdwb * Inv_W
                      + model_.pdwb * Inv_LW;
    paramPtr->b0 = model_.b0
                      + model_.lb0 * Inv_L
                      + model_.wb0 * Inv_W
                      + model_.pb0 * Inv_LW;
    paramPtr->b1 = model_.b1
                      + model_.lb1 * Inv_L
                      + model_.wb1 * Inv_W
                      + model_.pb1 * Inv_LW;
    paramPtr->alpha0 = model_.alpha0
                      + model_.lalpha0 * Inv_L
                      + model_.walpha0 * Inv_W
                      + model_.palpha0 * Inv_LW;
    paramPtr->alpha1 = model_.alpha1
                        + model_.lalpha1 * Inv_L
                        + model_.walpha1 * Inv_W
                        + model_.palpha1 * Inv_LW;
    paramPtr->beta0 = model_.beta0
                        + model_.lbeta0 * Inv_L
                        + model_.wbeta0 * Inv_W
                        + model_.pbeta0 * Inv_LW;
    paramPtr->agidl = model_.agidl
                       + model_.lagidl * Inv_L
                       + model_.wagidl * Inv_W
                       + model_.pagidl * Inv_LW;
    paramPtr->bgidl = model_.bgidl
                       + model_.lbgidl * Inv_L
                       + model_.wbgidl * Inv_W
                       + model_.pbgidl * Inv_LW;
    paramPtr->cgidl = model_.cgidl
                       + model_.lcgidl * Inv_L
                       + model_.wcgidl * Inv_W
                       + model_.pcgidl * Inv_LW;
    paramPtr->egidl = model_.egidl
                       + model_.legidl * Inv_L
                       + model_.wegidl * Inv_W
                       + model_.pegidl * Inv_LW;
    paramPtr->rgidl = model_.rgidl
                       + model_.lrgidl * Inv_L
                       + model_.wrgidl * Inv_W
                       + model_.prgidl * Inv_LW;
    paramPtr->kgidl = model_.kgidl
                       + model_.lkgidl * Inv_L
                       + model_.wkgidl * Inv_W
                       + model_.pkgidl * Inv_LW;
    paramPtr->fgidl = model_.fgidl
                       + model_.lfgidl * Inv_L
                       + model_.wfgidl * Inv_W
                       + model_.pfgidl * Inv_LW;
    paramPtr->agisl = model_.agisl
                         + model_.lagisl * Inv_L
                         + model_.wagisl * Inv_W
                         + model_.pagisl * Inv_LW;
    paramPtr->bgisl = model_.bgisl
                         + model_.lbgisl * Inv_L
                         + model_.wbgisl * Inv_W
                         + model_.pbgisl * Inv_LW;
    paramPtr->cgisl = model_.cgisl
                         + model_.lcgisl * Inv_L
                         + model_.wcgisl * Inv_W
                         + model_.pcgisl * Inv_LW;
    paramPtr->egisl = model_.egisl
                         + model_.legisl * Inv_L
                         + model_.wegisl * Inv_W
                         + model_.pegisl * Inv_LW;
    paramPtr->rgisl = model_.rgisl
                         + model_.lrgisl * Inv_L
                         + model_.wrgisl * Inv_W
                         + model_.prgisl * Inv_LW;
    paramPtr->kgisl = model_.kgisl
                         + model_.lkgisl * Inv_L
                         + model_.wkgisl * Inv_W
                         + model_.pkgisl * Inv_LW;
    paramPtr->fgisl = model_.fgisl
                         + model_.lfgisl * Inv_L
                         + model_.wfgisl * Inv_W
                         + model_.pfgisl * Inv_LW;
    paramPtr->aigc = model_.aigc
                       + model_.laigc * Inv_L
                       + model_.waigc * Inv_W
                       + model_.paigc * Inv_LW;
    paramPtr->bigc = model_.bigc
                       + model_.lbigc * Inv_L
                       + model_.wbigc * Inv_W
                       + model_.pbigc * Inv_LW;
    paramPtr->cigc = model_.cigc
                       + model_.lcigc * Inv_L
                       + model_.wcigc * Inv_W
                       + model_.pcigc * Inv_LW;
    paramPtr->aigsd = model_.aigsd
                       + model_.laigsd * Inv_L
                       + model_.waigsd * Inv_W
                       + model_.paigsd * Inv_LW;
    paramPtr->bigsd = model_.bigsd
                       + model_.lbigsd * Inv_L
                       + model_.wbigsd * Inv_W
                       + model_.pbigsd * Inv_LW;
    paramPtr->cigsd = model_.cigsd
                       + model_.lcigsd * Inv_L
                       + model_.wcigsd * Inv_W
                       + model_.pcigsd * Inv_LW;
    paramPtr->aigs = model_.aigs
                       + model_.laigs * Inv_L
                       + model_.waigs * Inv_W
                       + model_.paigs * Inv_LW;
    paramPtr->bigs = model_.bigs
                       + model_.lbigs * Inv_L
                       + model_.wbigs * Inv_W
                       + model_.pbigs * Inv_LW;
    paramPtr->cigs = model_.cigs
                       + model_.lcigs * Inv_L
                       + model_.wcigs * Inv_W
                       + model_.pcigs * Inv_LW;
    paramPtr->aigd = model_.aigd
                       + model_.laigd * Inv_L
                       + model_.waigd * Inv_W
                       + model_.paigd * Inv_LW;
    paramPtr->bigd = model_.bigd
                       + model_.lbigd * Inv_L
                       + model_.wbigd * Inv_W
                       + model_.pbigd * Inv_LW;
    paramPtr->cigd = model_.cigd
                       + model_.lcigd * Inv_L
                       + model_.wcigd * Inv_W
                       + model_.pcigd * Inv_LW;
    paramPtr->aigbacc = model_.aigbacc
                         + model_.laigbacc * Inv_L
                         + model_.waigbacc * Inv_W
                         + model_.paigbacc * Inv_LW;
    paramPtr->bigbacc = model_.bigbacc
                         + model_.lbigbacc * Inv_L
                         + model_.wbigbacc * Inv_W
                         + model_.pbigbacc * Inv_LW;
    paramPtr->cigbacc = model_.cigbacc
                         + model_.lcigbacc * Inv_L
                         + model_.wcigbacc * Inv_W
                         + model_.pcigbacc * Inv_LW;
    paramPtr->aigbinv = model_.aigbinv
                         + model_.laigbinv * Inv_L
                         + model_.waigbinv * Inv_W
                         + model_.paigbinv * Inv_LW;
    paramPtr->bigbinv = model_.bigbinv
                         + model_.lbigbinv * Inv_L
                         + model_.wbigbinv * Inv_W
                         + model_.pbigbinv * Inv_LW;
    paramPtr->cigbinv = model_.cigbinv
                         + model_.lcigbinv * Inv_L
                         + model_.wcigbinv * Inv_W
                         + model_.pcigbinv * Inv_LW;
    paramPtr->nigc = model_.nigc
                         + model_.lnigc * Inv_L
                         + model_.wnigc * Inv_W
                         + model_.pnigc * Inv_LW;
    paramPtr->nigbacc = model_.nigbacc
                         + model_.lnigbacc * Inv_L
                         + model_.wnigbacc * Inv_W
                         + model_.pnigbacc * Inv_LW;
    paramPtr->nigbinv = model_.nigbinv
                         + model_.lnigbinv * Inv_L
                         + model_.wnigbinv * Inv_W
                         + model_.pnigbinv * Inv_LW;
    paramPtr->ntox = model_.ntox
                      + model_.lntox * Inv_L
                      + model_.wntox * Inv_W
                      + model_.pntox * Inv_LW;
    paramPtr->eigbinv = model_.eigbinv
                         + model_.leigbinv * Inv_L
                         + model_.weigbinv * Inv_W
                         + model_.peigbinv * Inv_LW;
    paramPtr->pigcd = model_.pigcd
                       + model_.lpigcd * Inv_L
                       + model_.wpigcd * Inv_W
                       + model_.ppigcd * Inv_LW;
    paramPtr->poxedge = model_.poxedge
                         + model_.lpoxedge * Inv_L
                         + model_.wpoxedge * Inv_W
                         + model_.ppoxedge * Inv_LW;
    paramPtr->xrcrg1 = model_.xrcrg1
                        + model_.lxrcrg1 * Inv_L
                        + model_.wxrcrg1 * Inv_W
                        + model_.pxrcrg1 * Inv_LW;
    paramPtr->xrcrg2 = model_.xrcrg2
                        + model_.lxrcrg2 * Inv_L
                        + model_.wxrcrg2 * Inv_W
                        + model_.pxrcrg2 * Inv_LW;
    paramPtr->lambda = model_.lambda
                        + model_.llambda * Inv_L
                        + model_.wlambda * Inv_W
                        + model_.plambda * Inv_LW;
    paramPtr->vtl = model_.vtl
                        + model_.lvtl * Inv_L
                        + model_.wvtl * Inv_W
                        + model_.pvtl * Inv_LW;
    paramPtr->xn = model_.xn
                        + model_.lxn * Inv_L
                        + model_.wxn * Inv_W
                        + model_.pxn * Inv_LW;
    paramPtr->vfbsdoff = model_.vfbsdoff
                        + model_.lvfbsdoff * Inv_L
                        + model_.wvfbsdoff * Inv_W
                        + model_.pvfbsdoff * Inv_LW;
    paramPtr->tvfbsdoff = model_.tvfbsdoff
                        + model_.ltvfbsdoff * Inv_L
                        + model_.wtvfbsdoff * Inv_W
                        + model_.ptvfbsdoff * Inv_LW;

    paramPtr->cgsl = model_.cgsl
                        + model_.lcgsl * Inv_L
                        + model_.wcgsl * Inv_W
                        + model_.pcgsl * Inv_LW;
    paramPtr->cgdl = model_.cgdl
                        + model_.lcgdl * Inv_L
                        + model_.wcgdl * Inv_W
                        + model_.pcgdl * Inv_LW;
    paramPtr->ckappas = model_.ckappas
                        + model_.lckappas * Inv_L
                        + model_.wckappas * Inv_W
                        + model_.pckappas * Inv_LW;
    paramPtr->ckappad = model_.ckappad
                         + model_.lckappad * Inv_L
                         + model_.wckappad * Inv_W
                         + model_.pckappad * Inv_LW;
    paramPtr->cf = model_.cf
                        + model_.lcf * Inv_L
                        + model_.wcf * Inv_W
                        + model_.pcf * Inv_LW;
    paramPtr->clc = model_.clc
                        + model_.lclc * Inv_L
                        + model_.wclc * Inv_W
                        + model_.pclc * Inv_LW;
    paramPtr->cle = model_.cle
                        + model_.lcle * Inv_L
                        + model_.wcle * Inv_W
                        + model_.pcle * Inv_LW;
    paramPtr->vfbcv = model_.vfbcv
                        + model_.lvfbcv * Inv_L
                        + model_.wvfbcv * Inv_W
                        + model_.pvfbcv * Inv_LW;
    paramPtr->acde = model_.acde
                      + model_.lacde * Inv_L
                      + model_.wacde * Inv_W
                      + model_.pacde * Inv_LW;
    paramPtr->moin = model_.moin
                      + model_.lmoin * Inv_L
                      + model_.wmoin * Inv_W
                      + model_.pmoin * Inv_LW;
    paramPtr->noff = model_.noff
                      + model_.lnoff * Inv_L
                      + model_.wnoff * Inv_W
                      + model_.pnoff * Inv_LW;
    paramPtr->voffcv = model_.voffcv
                        + model_.lvoffcv * Inv_L
                        + model_.wvoffcv * Inv_W
                        + model_.pvoffcv * Inv_LW;
    paramPtr->kvth0we = model_.kvth0we
                        + model_.lkvth0we * Inv_L
                        + model_.wkvth0we * Inv_W
                        + model_.pkvth0we * Inv_LW;
    paramPtr->k2we = model_.k2we
                        + model_.lk2we * Inv_L
                        + model_.wk2we * Inv_W
                        + model_.pk2we * Inv_LW;
    paramPtr->ku0we = model_.ku0we
                        + model_.lku0we * Inv_L
                        + model_.wku0we * Inv_W
                        + model_.pku0we * Inv_LW;

    paramPtr->abulkCVfactor = 1.0 + pow((paramPtr->clc / paramPtr->leffCV), paramPtr->cle);

    T0 = (TRatio - 1.0);

    PowWeffWr = pow(paramPtr->weffCJ * 1.0e6, paramPtr->wr) * nf;

    T1 = T2 = T3 = T4 = 0.0;
    paramPtr->ucs = paramPtr->ucs * pow(TRatio, paramPtr->ucste);
    if(model_.tempMod == 0)
    {
      paramPtr->ua = paramPtr->ua + paramPtr->ua1 * T0;
      paramPtr->ub = paramPtr->ub + paramPtr->ub1 * T0;
      paramPtr->uc = paramPtr->uc + paramPtr->uc1 * T0;
      paramPtr->ud = paramPtr->ud + paramPtr->ud1 * T0;
        paramPtr->vsattemp = paramPtr->vsat - paramPtr->at * T0;
      T10 = paramPtr->prt * T0;
      if(model_.rdsMod)
      {
          // External Rd(V)
          T1 = paramPtr->rdw + T10;
            T2 = model_.rdwmin + T10;
          // External Rs(V)
          T3 = paramPtr->rsw + T10;
            T4 = model_.rswmin + T10;
      }
      // Internal Rds(V) in IV
      paramPtr->rds0 = (paramPtr->rdsw + T10) * nf / PowWeffWr;
      paramPtr->rdswmin = (model_.rdswmin + T10) * nf / PowWeffWr;
    }
    else
    {
      if (model_.tempMod == 3)
      {
        paramPtr->ua = paramPtr->ua * pow(TRatio, paramPtr->ua1) ;
        paramPtr->ub = paramPtr->ub * pow(TRatio, paramPtr->ub1);
        paramPtr->uc = paramPtr->uc * pow(TRatio, paramPtr->uc1);
        paramPtr->ud = paramPtr->ud * pow(TRatio, paramPtr->ud1);
      }
      else
      {
        // tempMod = 1, 2
        paramPtr->ua = paramPtr->ua * (1.0 + paramPtr->ua1 * delTemp) ;
        paramPtr->ub = paramPtr->ub * (1.0 + paramPtr->ub1 * delTemp);
        paramPtr->uc = paramPtr->uc * (1.0 + paramPtr->uc1 * delTemp);
        paramPtr->ud = paramPtr->ud * (1.0 + paramPtr->ud1 * delTemp);
      }
      paramPtr->vsattemp = paramPtr->vsat * (1.0 - paramPtr->at * delTemp);
      T10 = 1.0 + paramPtr->prt * delTemp;
      if(model_.rdsMod)
      {
        // External Rd(V)
        T1 = paramPtr->rdw * T10;
        T2 = model_.rdwmin * T10;
        // External Rs(V)
        T3 = paramPtr->rsw * T10;
        T4 = model_.rswmin * T10;
      }
      // Internal Rds(V) in IV
      paramPtr->rds0 = paramPtr->rdsw * T10 * nf / PowWeffWr;
      paramPtr->rdswmin = model_.rdswmin * T10 * nf / PowWeffWr;
    }
    if (T1 < 0.0)
    {
      T1 = 0.0;
      UserWarning(*this) << "Rdw at current temperature is negative; set to 0";
    }
    if (T2 < 0.0)
    {
      T2 = 0.0;
      UserWarning(*this) << "Rdwmin at current temperature is negative; set to 0";
    }
    paramPtr->rd0 = T1 / PowWeffWr;
    paramPtr->rdwmin = T2 / PowWeffWr;
    if (T3 < 0.0)
    {
      T3 = 0.0;
      UserWarning(*this) << "Rsw at current temperature is negative; set to 0";
    }
    if (T4 < 0.0)
    {
      T4 = 0.0;
      UserWarning(*this) << "Rswmin at current temperature is negative; set to 0";
    }
    paramPtr->rs0 = T3 / PowWeffWr;
    paramPtr->rswmin = T4 / PowWeffWr;

    if (paramPtr->u0 > 1.0)
        paramPtr->u0 = paramPtr->u0 / 1.0e4;

    // mobility channel length dependence
    T5 = 1.0 - paramPtr->up * exp( - paramPtr->leff / paramPtr->lp);
    paramPtr->u0temp = paramPtr->u0 * T5
    * pow(TRatio, paramPtr->ute);
    if (paramPtr->eu < 0.0)
    {
      paramPtr->eu = 0.0;
      UserWarning(*this) << "eu has been negative; reset to 0.0";
    }
    if (paramPtr->ucs < 0.0)
    {
      paramPtr->ucs = 0.0;
      UserWarning(*this) << "ucs has been negative; reset to 0.0";
    }

    paramPtr->vfbsdoff = paramPtr->vfbsdoff * (1.0 + paramPtr->tvfbsdoff * delTemp);
    paramPtr->voff = paramPtr->voff * (1.0 + paramPtr->tvoff * delTemp);

    paramPtr->nfactor = paramPtr->nfactor + paramPtr->tnfactor * delTemp / Tnom;
    paramPtr->voffcv = paramPtr->voffcv * (1.0 + paramPtr->tvoffcv * delTemp);
    paramPtr->eta0 = paramPtr->eta0 + paramPtr->teta0 * delTemp / Tnom;

    // Source End Velocity Limit
    if((model_.vtlGiven) && (model_.vtl > 0.0) )
    {
       if(model_.lc < 0.0) paramPtr->lc = 0.0;
       else   paramPtr->lc = model_.lc ;
       T0 = paramPtr->leff / (paramPtr->xn * paramPtr->leff + paramPtr->lc);
       paramPtr->tfactor = (1.0 - T0) / (1.0 + T0 );
    }

    paramPtr->cgdo = (model_.cgdo + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgso = (model_.cgso + paramPtr->cf) * paramPtr->weffCV;
    paramPtr->cgbo = model_.cgbo * paramPtr->leffCV * nf;

    if (!model_.ndepGiven && model_.gamma1Given)
    {
      T0 = paramPtr->gamma1 * model_.coxe;
      paramPtr->ndep = 3.01248e22 * T0 * T0;
    }

    paramPtr->phi = Vtm0 * log(paramPtr->ndep / ni) + paramPtr->phin + 0.4;

    paramPtr->sqrtPhi = sqrt(paramPtr->phi);
    paramPtr->phis3 = paramPtr->sqrtPhi * paramPtr->phi;

    paramPtr->Xdep0 = sqrt(2.0 * epssub / (Charge_q
                        * paramPtr->ndep * 1.0e6))
                        * paramPtr->sqrtPhi;

    paramPtr->sqrtXdep0 = sqrt(paramPtr->Xdep0);

    if (model_.mtrlMod == 0)
    {
      paramPtr->litl = sqrt(3.0 * 3.9 / epsrox * paramPtr->xj * toxe);
    }
    else
    {
      paramPtr->litl = sqrt(model_.epsrsub/epsrox * paramPtr->xj * toxe);
    }

    paramPtr->vbi = Vtm0 * log(paramPtr->nsd * paramPtr->ndep / (ni * ni));

    if (model_.mtrlMod == 0)
    {
      if (paramPtr->ngate > 0.0)
      {   paramPtr->vfbsd = Vtm0 * log(paramPtr->ngate
                                          / paramPtr->nsd);
      }
      else
        paramPtr->vfbsd = 0.0;
    }
    else
    {
      T0 = Vtm0 * log(paramPtr->nsd/ni);
      T1 = 0.5 * Eg0;
      if(T0 > T1)
        T0 = T1;
      T2 = model_.easub + T1 - model_.dtype * T0;
      paramPtr->vfbsd = model_.phig - T2;
    }

    paramPtr->cdep0 = sqrt(Charge_q * epssub
                        * paramPtr->ndep * 1.0e6 / 2.0 / paramPtr->phi);

    paramPtr->ToxRatio = exp(paramPtr->ntox
                            * log(model_.toxref / toxe))
                            / toxe / toxe;

    paramPtr->ToxRatioEdge = exp(paramPtr->ntox
                              * log(model_.toxref
                              / (toxe * paramPtr->poxedge)))
                              / toxe / toxe
                              / paramPtr->poxedge / paramPtr->poxedge;

    paramPtr->Aechvb = (model_.dtype == CONSTNMOS) ? 4.97232e-7 : 3.42537e-7;
    paramPtr->Bechvb = (model_.dtype == CONSTNMOS) ? 7.45669e11 : 1.16645e12;
    if (model_.versionDouble <= 4.80)
    {
      paramPtr->AechvbEdgeS = paramPtr->Aechvb * paramPtr->weff
                              * model_.dlcig * paramPtr->ToxRatioEdge;
      paramPtr->AechvbEdgeD = paramPtr->Aechvb * paramPtr->weff
                              * model_.dlcigd * paramPtr->ToxRatioEdge;
    }
    else
    {
      if (model_.dlcig < 0.0)
      {
        UserWarning(*this) << "Warning: dlcig = " << model_.dlcig <<
          " is negative. Set to zero.";
        model_.dlcig = 0.0;
      }
      paramPtr->AechvbEdgeS = paramPtr->Aechvb * paramPtr->weff
                             * model_.dlcig * paramPtr->ToxRatioEdge;
      if (model_.dlcigd < 0.0)
      {
        UserWarning(*this) << "Warning: dlcigd = " << model_.dlcigd <<
          " is negative. Set to zero.";
        model_.dlcigd = 0.0;
      }
      paramPtr->AechvbEdgeD = paramPtr->Aechvb * paramPtr->weff
                             * model_.dlcigd * paramPtr->ToxRatioEdge;
    }
    paramPtr->BechvbEdge = -paramPtr->Bechvb
                            * model_.toxe * paramPtr->poxedge;
    paramPtr->Aechvb *= paramPtr->weff * paramPtr->leff
                          * paramPtr->ToxRatio;
    paramPtr->Bechvb *= -toxe;


    paramPtr->mstar = 0.5 + atan(paramPtr->minv) / M_PI;
    paramPtr->mstarcv = 0.5 + atan(paramPtr->minvcv) / M_PI;
    paramPtr->voffcbn =  paramPtr->voff + model_.voffl / paramPtr->leff;
    paramPtr->voffcbncv =  paramPtr->voffcv + model_.voffcvl / paramPtr->leff;

    paramPtr->ldeb = sqrt(epssub * Vtm0 / (Charge_q
                                               * paramPtr->ndep * 1.0e6)) / 3.0;
    paramPtr->acde *= pow((paramPtr->ndep / 2.0e16), -0.25);


    if (model_.k1Given || model_.k2Given)
    {
      if (!model_.k1Given)
      {
        UserWarning(*this) << "k1 should be specified with k2";
        paramPtr->k1 = 0.53;
      }
      if (!model_.k2Given)
      {
        UserWarning(*this) << "k2 should be specified with k1";
        paramPtr->k2 = -0.0186;
      }
      if (model_.nsubGiven)
      {
        UserWarning(*this) << "nsub is ignored because k1 or k2 is given";
      }
      if (model_.xtGiven)
      {
        UserWarning(*this) << "xt is ignored because k1 or k2 is given";
      }
      if (model_.vbxGiven)
      {
        UserWarning(*this) << "vbx is ignored because k1 or k2 is given";
      }
      if (model_.gamma1Given)
      {
        UserWarning(*this) << "gamma1 is ignored because k1 or k2 is given";
      }
      if (model_.gamma2Given)
      {
        UserWarning(*this) << "gamma2 is ignored because k1 or k2 is given";
      }
    }
    else
    {
      if (!model_.vbxGiven)
            paramPtr->vbx = paramPtr->phi - 7.7348e-4
                             * paramPtr->ndep * paramPtr->xt * paramPtr->xt;

      if (paramPtr->vbx > 0.0)
            paramPtr->vbx = -paramPtr->vbx;
      if (paramPtr->vbm > 0.0)
            paramPtr->vbm = -paramPtr->vbm;

      if (!model_.gamma1Given)
            paramPtr->gamma1 = 5.753e-12
                                * sqrt(paramPtr->ndep) / model_.coxe;

      if (!model_.gamma2Given)
            paramPtr->gamma2 = 5.753e-12
                              * sqrt(paramPtr->nsub) / model_.coxe;

      T0 = paramPtr->gamma1 - paramPtr->gamma2;

      T1 = sqrt(paramPtr->phi - paramPtr->vbx) - paramPtr->sqrtPhi;

      T2 = sqrt(paramPtr->phi * (paramPtr->phi - paramPtr->vbm)) - paramPtr->phi;

      paramPtr->k2 = T0 * T1 / (2.0 * T2 + paramPtr->vbm);

      paramPtr->k1 = paramPtr->gamma2 - 2.0
                    * paramPtr->k2 * sqrt(paramPtr->phi - paramPtr->vbm);
    }

    if (!model_.vfbGiven)
    {
      if (model_.vth0Given)
      {   paramPtr->vfb = model_.dtype * paramPtr->vth0
                             - paramPtr->phi - paramPtr->k1
                             * paramPtr->sqrtPhi;
      }
      else
      {
        if ((model_.mtrlMod) && (model_.phigGiven) &&
            (model_.nsubGiven))
        {
          T0 = Vtm0 * log(paramPtr->nsub/ni);
          T1 = 0.5 * Eg0;
          if(T0 > T1)
            T0 = T1;
          T2 = model_.easub + T1 + model_.dtype * T0;
          paramPtr->vfb = model_.phig - T2;
        }
        else
        {
          paramPtr->vfb = -1.0;
        }
      }
    }
    if (!model_.vth0Given)
    {
      paramPtr->vth0 = model_.dtype * (paramPtr->vfb
                          + paramPtr->phi + paramPtr->k1
                          * paramPtr->sqrtPhi);
    }

    paramPtr->k1ox = paramPtr->k1 * toxe / model_.toxm;

    tmp = sqrt(epssub / (epsrox * CONSTEPS0)
        * toxe * paramPtr->Xdep0);

    T0 = paramPtr->dsub * paramPtr->leff / tmp;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      paramPtr->theta0vb0 = T1 / T4;
    }
    else
    {
      paramPtr->theta0vb0 = 1.0 / (CONSTMAX_EXP - 2.0);
    }

    T0 = paramPtr->drout * paramPtr->leff / tmp;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T5 = T1 / T4;
    }
    else
        T5 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted

    paramPtr->thetaRout = paramPtr->pdibl1 * T5 + paramPtr->pdibl2;

    tmp = sqrt(paramPtr->Xdep0);
    tmp1 = paramPtr->vbi - paramPtr->phi;
    tmp2 = model_.factor1 * tmp;

    T0 = paramPtr->dvt1w * paramPtr->weff * paramPtr->leff / tmp2;

    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T8 = T1 / T4;
    }
    else
    {
      T8 = 1.0 / (CONSTMAX_EXP - 2.0);
    }

    T0 = paramPtr->dvt0w * T8;
    T8 = T0 * tmp1;

    T0 = paramPtr->dvt1 * paramPtr->leff / tmp2;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T9 = T1 / T4;
    }
    else
      T9 = 1.0 / (CONSTMAX_EXP - 2.0);
    T9 = paramPtr->dvt0 * T9 * tmp1;

    T4 = toxe * paramPtr->phi
       / (paramPtr->weff + paramPtr->w0);

    T0 = sqrt(1.0 + paramPtr->lpe0 / paramPtr->leff);
    if((model_.tempMod == 1) || (model_.tempMod == 0))
      T3 = (paramPtr->kt1 + paramPtr->kt1l / paramPtr->leff)
          * (TRatio - 1.0);
    if((model_.tempMod == 2) || (model_.tempMod == 3))
          T3 = - paramPtr->kt1 * (TRatio - 1.0);

    T5 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi + T3;

    paramPtr->vfbzbfactor
      = - T8 - T9 + paramPtr->k3 * T4 + T5 - paramPtr->phi - paramPtr->k1 * paramPtr->sqrtPhi;

    // stress effect

    wlod = model_.wlod;
    if (model_.wlod < 0.0)
    {
      UserWarning(*this) << "WLOD =is less than 0. 0.0 is used";
      wlod = 0.0;
    }
    T0 = pow(Lnew, model_.llodku0);
    W_tmp = Wnew + wlod;
    T1 = pow(W_tmp, model_.wlodku0);
    tmp1 = model_.lku0 / T0 + model_.wku0 / T1
           + model_.pku0 / (T0 * T1);
    paramPtr->ku0 = 1.0 + tmp1;

    T0 = pow(Lnew, model_.llodvth);
    T1 = pow(W_tmp, model_.wlodvth);
    tmp1 = model_.lkvth0 / T0 + model_.wkvth0 / T1
         + model_.pkvth0 / (T0 * T1);
    paramPtr->kvth0 = 1.0 + tmp1;
    paramPtr->kvth0 = sqrt(paramPtr->kvth0*paramPtr->kvth0 + DELTA);

    T0 = (TRatio - 1.0);
    paramPtr->ku0temp = paramPtr->ku0 * (1.0 + model_.tku0 *T0) + DELTA;

    Inv_saref = 1.0/(model_.saref + 0.5*Ldrn);
    Inv_sbref = 1.0/(model_.sbref + 0.5*Ldrn);
    paramPtr->inv_od_ref = Inv_saref + Inv_sbref;
    paramPtr->rho_ref = model_.ku0 / paramPtr->ku0temp * paramPtr->inv_od_ref;

    /*high k*/
    /*Calculate VgsteffVth for mobMod=3*/
    if(model_.mobMod==3)
    { /*Calculate n @ Vbs=Vds=0*/
      lt1 = model_.factor1* paramPtr->sqrtXdep0;
      T0 = paramPtr->dvt1 * paramPtr->leff / lt1;
      if (T0 < CONSTEXP_THRESHOLD)
      {
        T1 = exp(T0);
        T2 = T1 - 1.0;
        T3 = T2 * T2;
        T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
        Theta0 = T1 / T4;
      }
      else
        Theta0 = 1.0 / (CONSTMAX_EXP - 2.0);

      tmp1 = epssub / paramPtr->Xdep0;
      tmp2 = paramPtr->nfactor * tmp1;
      tmp3 = (tmp2 + paramPtr->cdsc * Theta0 + paramPtr->cit) / model_.coxe;
      if (tmp3 >= -0.5)
        n0 = 1.0 + tmp3;
      else
      {
        T0 = 1.0 / (3.0 + 8.0 * tmp3);
        n0 = (1.0 + 3.0 * tmp3) * T0;
      }

      T0 = n0 * model_.vtm;
      T1 = paramPtr->voffcbn;
      T2 = T1/T0;
      if (T2 < -CONSTEXP_THRESHOLD)
      {   T3 = model_.coxe * CONSTMIN_EXP / paramPtr->cdep0;
        T4 = paramPtr->mstar + T3 * n0;
      }
      else if (T2 > CONSTEXP_THRESHOLD)
      {   T3 = model_.coxe * CONSTMAX_EXP / paramPtr->cdep0;
        T4 = paramPtr->mstar + T3 * n0;
      }
      else
      {  T3 = exp(T2)* model_.coxe / paramPtr->cdep0;
        T4 = paramPtr->mstar + T3 * n0;
      }
      paramPtr->VgsteffVth = T0 * log(2.0)/T4;
    }

    /* New DITS term added in 4.7 */
    T0 = -paramPtr->dvtp3 * log(paramPtr->leff);
    DEXP2(T0, T1);
    paramPtr->dvtp2factor = paramPtr->dvtp5 + paramPtr->dvtp2 * T1;
  } // End of size if-statement

  //  stress effect
  if( (sa > 0.0) && (sb > 0.0) &&
      ((nf == 1.0) || ((nf > 1.0) && (sd > 0.0))) )
  {
    Inv_sa = 0;
    Inv_sb = 0;

    kvsat = model_.kvsat;
    if (model_.kvsat < -1.0 )
    {
      UserWarning(*this) << "KVSAT is too small; -1.0 is used";
      kvsat = -1.0;
    }
    if (model_.kvsat > 1.0)
    {
      UserWarning(*this) << "KVSAT is too big; 1.0 is used";
      kvsat = 1.0;
    }

    int i=0;
    for(i = 0; i < nf; i++)
    {
      T0 = 1.0 / nf / (sa + 0.5*Ldrn + i * (sd +Ldrn));
      T1 = 1.0 / nf / (sb + 0.5*Ldrn + i * (sd +Ldrn));
      Inv_sa += T0;
      Inv_sb += T1;
    }
    Inv_ODeff = Inv_sa + Inv_sb;
    rho = model_.ku0 / paramPtr->ku0temp * Inv_ODeff;
    T0 = (1.0 + rho)/(1.0 + paramPtr->rho_ref);
    u0temp = paramPtr->u0temp * T0;

    T1 = (1.0 + kvsat * rho)/(1.0 + kvsat * paramPtr->rho_ref);
    vsattemp = paramPtr->vsattemp * T1;

    OD_offset = Inv_ODeff - paramPtr->inv_od_ref;
    dvth0_lod = model_.kvth0 / paramPtr->kvth0 * OD_offset;
    dk2_lod = model_.stk2 / pow(paramPtr->kvth0, model_.lodk2) *
                           OD_offset;
    deta0_lod = model_.steta0 / pow(paramPtr->kvth0, model_.lodeta0) *
                             OD_offset;
    vth0 = paramPtr->vth0 + dvth0_lod;

    eta0 = paramPtr->eta0 + deta0_lod;
    k2 = paramPtr->k2 + dk2_lod;
  }
  else
  {
    u0temp = paramPtr->u0temp;
    vth0 = paramPtr->vth0;
    vsattemp = paramPtr->vsattemp;
    eta0 = paramPtr->eta0;
    k2 = paramPtr->k2;
  }

  //  Well Proximity Effect
  if (model_.wpemod)
  {
    if( (!scaGiven) && (!scbGiven) && (!sccGiven) )
    {
      if((scGiven) && (sc > 0.0) )
      {
        T1 = sc + Wdrn;
        T2 = 1.0 / model_.scref;
        sca = model_.scref * model_.scref / (sc * T1);
        scb = ( (0.1 * sc + 0.01 * model_.scref)
            * exp(-10.0 * sc * T2)
            - (0.1 * T1 + 0.01 * model_.scref)
            * exp(-10.0 * T1 * T2) ) / Wdrn;
        scc = ( (0.05 * sc + 0.0025 * model_.scref)
                                  * exp(-20.0 * sc * T2)
                                  - (0.05 * T1 + 0.0025 * model_.scref)
                                  * exp(-20.0 * T1 * T2) ) / Wdrn;
      }
      else
      {
        UserWarning(*this) << "No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive";
      }
    }
    if (sca < 0.0)
    {
      UserWarning(*this) << "SCA = " << sca << " is negative. Set to 0.0.";
      sca = 0.0;
    }
    if (scb < 0.0)
    {
      UserWarning(*this) << "SCB = " << scb << " is negative. Set to 0.0.";
      scb = 0.0;
    }
    if (scc < 0.0)
    {
      UserWarning(*this) << "SCC = " << scc << " is negative. Set to 0.0.";
      scc = 0.0;
    }
    if (sc < 0.0)
    {
      UserWarning(*this) << "SC = " << sc << " is negative. Set to 0.0.";
      sc = 0.0;
    }
    sceff = sca + model_.web * scb
                + model_.wec * scc;
    vth0 += paramPtr->kvth0we * sceff;
    k2 +=  paramPtr->k2we * sceff;
    T3 =  1.0 + paramPtr->ku0we * sceff;
    if (T3 <= 0.0)
    {
      T3 = 0.0;
      UserWarning(*this) << "ku0we = %g is negatively too high. Negative mobility!";
    }
    u0temp *= T3;
  }

  // adding delvto
  vth0 += delvto;
  vfb = paramPtr->vfb + model_.dtype * delvto;

  // Instance variables calculation
  T3 = model_.dtype * vth0
         - vfb - paramPtr->phi;
  T4 = T3 + T3;
  T5 = 2.5 * T3;
  vtfbphi1 = (model_.dtype == CONSTNMOS) ? T4 : T5;
  if (vtfbphi1 < 0.0)
      vtfbphi1 = 0.0;

  vtfbphi2 = 4.0 * T3;
  if (vtfbphi2 < 0.0)
      vtfbphi2 = 0.0;

  if (k2 < 0.0)
  {
    T0 = 0.5 * paramPtr->k1 / k2;
    vbsc = 0.9 * (paramPtr->phi - T0 * T0);
    if (vbsc > -3.0)
              vbsc = -3.0;
    else if (vbsc < -30.0)
              vbsc = -30.0;
  }
  else
    vbsc = -30.0;
  if (vbsc > paramPtr->vbm)
      vbsc = paramPtr->vbm;
  k2ox = k2 * toxe
                        / model_.toxm;

  vfbzb = paramPtr->vfbzbfactor
                    +  model_.dtype * vth0 ;

  cgso = paramPtr->cgso;
  cgdo = paramPtr->cgdo;

  lnl = log(paramPtr->leff * 1.0e6);
  lnw = log(paramPtr->weff * 1.0e6);
  lnnf = log(nf);

  bodymode = 5;
  if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    bodymode = 1;
  else
    if( (!model_.rbsbx0Given && !model_.rbsby0Given) ||
        (!model_.rbdbx0Given && !model_.rbdby0Given) )
    bodymode = 3;

  if(rbodyMod == 2)
  {
    if (bodymode == 5)
    {
      rbsbx =  model_.rbsbx0 * exp(model_.rbsdbxl*lnl +
                                   model_.rbsdbxw*lnw + model_.rbsdbxnf*lnnf);
      rbsby =  model_.rbsby0 * exp(model_.rbsdbyl*lnl +
                                   model_.rbsdbyw*lnw + model_.rbsdbynf*lnnf);
      rbsb = rbsbx * rbsby / (rbsbx + rbsby);

      rbdbx =  model_.rbdbx0 * exp(model_.rbsdbxl*lnl +
                                   model_.rbsdbxw*lnw + model_.rbsdbxnf*lnnf);
      rbdby =  model_.rbdby0 * exp(model_.rbsdbyl*lnl +
                                   model_.rbsdbyw*lnw + model_.rbsdbynf*lnnf);

      rbdb = rbdbx * rbdby / (rbdbx + rbdby);
    }

    if ((bodymode == 3)|| (bodymode == 5))
    {
      rbps = model_.rbps0 * exp( model_.rbpsl * lnl +
                                 model_.rbpsw * lnw + model_.rbpsnf * lnnf );
      rbpd = model_.rbpd0 * exp( model_.rbpdl * lnl +
                                 model_.rbpdw * lnw + model_.rbpdnf * lnnf );
    }

    rbpbx =  model_.rbpbx0 * exp(model_.rbpbxl*lnl +
                                 model_.rbpbxw*lnw + model_.rbpbxnf*lnnf );
    rbpby =  model_.rbpby0 * exp(model_.rbpbyl*lnl +
                                 model_.rbpbyw*lnw + model_.rbpbynf*lnnf );
    rbpb = rbpbx*rbpby/(rbpbx + rbpby);
  }


  if ((rbodyMod == 1 ) || ((rbodyMod == 2 ) && (bodymode == 5)) )
  {
    if (rbdb < 1.0e-3) grbdb = 1.0e3; // in mho
    else               grbdb = model_.gbmin + 1.0 / rbdb;

    if (rbpb < 1.0e-3) grbpb = 1.0e3;
    else               grbpb = model_.gbmin + 1.0 / rbpb;

    if (rbps < 1.0e-3) grbps = 1.0e3;
    else               grbps = model_.gbmin + 1.0 / rbps;

    if (rbsb < 1.0e-3) grbsb = 1.0e3;
    else               grbsb = model_.gbmin + 1.0 / rbsb;

    if (rbpd < 1.0e-3) grbpd = 1.0e3;
    else               grbpd = model_.gbmin + 1.0 / rbpd;

  }

  if((rbodyMod == 2) && (bodymode == 3))
  {
    grbdb = grbsb = model_.gbmin;

    if (rbpb < 1.0e-3) grbpb = 1.0e3;
    else               grbpb = model_.gbmin + 1.0 / rbpb;

    if (rbps < 1.0e-3) grbps = 1.0e3;
    else               grbps = model_.gbmin + 1.0 / rbps;

    if (rbpd < 1.0e-3) grbpd = 1.0e3;
    else               grbpd = model_.gbmin + 1.0 / rbpd;
  }

  if((rbodyMod == 2) && (bodymode == 1))
  {
    grbdb = grbsb = model_.gbmin;
    grbps = grbpd = 1.0e3;
    if (rbpb < 1.0e-3)
    {
      grbpb = 1.0e3;
    }
    else
    {
      grbpb = model_.gbmin + 1.0 / rbpb;
    }
  }

  // Process geomertry dependent parasitics
  grgeltd = model_.rshg * (xgw
        + paramPtr->weffCJ / 3.0 / ngcon) /
        (ngcon * nf * (Lnew - model_.xgl));

  if (grgeltd > 0.0)
  {
     grgeltd = 1.0 / grgeltd;
  }
  else
  {
   grgeltd = 1.0e3; // mho
   if (rgateMod != 0)
   {
     UserWarning(*this) << "The gate conductance reset to 1.0e3 mho";
   }
  }

  DMCGeff = model_.dmcg - model_.dmcgt;
  DMCIeff = model_.dmci;
  DMDGeff = model_.dmdg - model_.dmcgt;

  // New Diode Model v4.7
  if (sourcePerimeterGiven)
  {
    // given
    if (sourcePerimeter == 0.0)
      Pseff = 0.0;
    else if (sourcePerimeter < 0.0)
    {
      UserWarning(*this) << "Source Perimeter is specified as negative, it is set to zero.";
      Pseff = 0.0;
    }
    else
    {
      if (model_.perMod == 0)
        Pseff = sourcePerimeter;
      else
        Pseff = sourcePerimeter - paramPtr->weffCJ * nf;
    }
  }
  else
  {
    // not given
    PAeffGeo(nf, geoMod, min,
             paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
             (Pseff), dumPd, dumAs, dumAd);
  }
  if (Pseff < 0.0)
  {
    // v4.7 final check
    Pseff = 0.0;
    UserWarning(*this) << "Pseff is negative, it is set to zero.";
  }

  if (drainPerimeterGiven)
  {
    // given
    if (drainPerimeter == 0.0)
      Pdeff = 0.0;
    else if (drainPerimeter < 0.0)
    {
      UserWarning(*this) << "Drain Perimeter is specified as negative, it is set to zero.";
      Pdeff = 0.0;
    }
    else
    {
      if (model_.perMod == 0)
        Pdeff = drainPerimeter;
      else
        Pdeff = drainPerimeter - paramPtr->weffCJ * nf;
    }
  }
  else
  {
    // not given
    PAeffGeo(nf, geoMod, min,
             paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
             dumPs, (Pdeff), dumAs, dumAd);
  }
  if (Pdeff < 0.0)
  {
    Pdeff = 0.0;
    UserWarning(*this) << "Pdeff is negative, it is set to zero.";
  }

  if (sourceAreaGiven)
   Aseff = sourceArea;
  else
   PAeffGeo(nf, geoMod, min,
              paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
              dumPs, dumPd, (Aseff), dumAd);
  if (Aseff < 0.0)
  {
    Aseff = 0.0;
    UserWarning(*this) << "Aseff is negative, it is set to zero.";
  }

  if (drainAreaGiven)
    Adeff = drainArea;
  else
    PAeffGeo(nf, geoMod, min,
              paramPtr->weffCJ, DMCGeff, DMCIeff, DMDGeff,
              dumPs, dumPd, dumAs, (Adeff));
  if (Adeff < 0.0)
  {
    Adeff = 0.0;
    UserWarning(*this) << "Adeff is negative, it is set to zero.";
  }

  // Processing S/D resistance and conductance below
  if (sourceMOSFET_B4Exists)
  {
    sourceConductance = 0.0;
    if(sourceSquaresGiven)
    {
      sourceConductance = model_.sheetResistance * sourceSquares;
    }
    else if (rgeoMod > 0)
    {
      RdseffGeo(nf, geoMod,
        rgeoMod, min,
        paramPtr->weffCJ, model_.sheetResistance,
        DMCGeff, DMCIeff, DMDGeff, 1, (sourceConductance));
    }
    else
    {
      sourceConductance = 0.0;
    }

    if (sourceConductance > 0.0)
    {
      sourceConductance = 1.0 / sourceConductance;
    }
    else
    {
      sourceConductance = 1.0e3; // mho
      UserWarning(*this) << "Source conductance reset to 1.0e3 mho";
    }
  }
  else
  {
    sourceConductance = 0.0;
  }

  if (drainMOSFET_B4Exists)
  {
    drainConductance = 0.0;
    if(drainSquaresGiven)
    {
      drainConductance = model_.sheetResistance * drainSquares;
    }
    else if (rgeoMod > 0)
    {
      RdseffGeo(nf, geoMod, rgeoMod, min,
        paramPtr->weffCJ, model_.sheetResistance,
        DMCGeff, DMCIeff, DMDGeff, 0, (drainConductance));
    }
    else
    {
      drainConductance = 0.0;
    }

    if (drainConductance > 0.0)
    {
      drainConductance = 1.0 / drainConductance;
    }
    else
    {
      drainConductance = 1.0e3; // mho
      UserWarning(*this) << "Drain conductance reset to 1.0e3 mho";
    }
  }
  else
  {
    drainConductance = 0.0;
  }

  // End of Rsd processing

  Nvtms = model_.vtm * model_.SjctEmissionCoeff;
  if ((Aseff <= 0.0) && (Pseff <= 0.0))
  {
    SourceSatCurrent = 0.0;
  }
  else
  {
    SourceSatCurrent = Aseff * model_.SjctTempSatCurDensity
                       + Pseff * model_.SjctSidewallTempSatCurDensity
                       + paramPtr->weffCJ * nf
                       * model_.SjctGateSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent > 0.0)
  {
    switch(model_.dioMod)
    {
      case 0:
        if ((model_.bvs / Nvtms) > CONSTEXP_THRESHOLD)
          XExpBVS = model_.xjbvs * CONSTMIN_EXP;
        else
          XExpBVS = model_.xjbvs * exp(-model_.bvs / Nvtms);
        break;
      case 1:
        DioIjthVjmEval(Nvtms, model_.ijthsfwd, SourceSatCurrent,
                    0.0, (vjsmFwd));
        IVjsmFwd = SourceSatCurrent * exp(vjsmFwd / Nvtms);
        break;
      case 2:
        if ((model_.bvs / Nvtms) > CONSTEXP_THRESHOLD)
        {
          XExpBVS = model_.xjbvs * CONSTMIN_EXP;
          tmp = CONSTMIN_EXP;
        }
        else
        {
          XExpBVS = exp(-model_.bvs / Nvtms);
          tmp = XExpBVS;
          XExpBVS *= model_.xjbvs;
        }

        DioIjthVjmEval(Nvtms, model_.ijthsfwd, SourceSatCurrent,
                                XExpBVS, (vjsmFwd));
        T0 = exp(vjsmFwd / Nvtms);
        IVjsmFwd = SourceSatCurrent * (T0 - XExpBVS / T0
                            + XExpBVS - 1.0);
        SslpFwd = SourceSatCurrent
                             * (T0 + XExpBVS / T0) / Nvtms;

        T2 = model_.ijthsrev / SourceSatCurrent;
        if (T2 < 1.0)
        {
          T2 = 10.0;
          UserWarning(*this) << "ijthsrev too small and set to 10 times IsbSat";
        }
        vjsmRev = -model_.bvs
                             - Nvtms * log((T2 - 1.0) / model_.xjbvs);
        T1 = model_.xjbvs * exp(-(model_.bvs
            + vjsmRev) / Nvtms);
        IVjsmRev = SourceSatCurrent * (1.0 + T1);
        SslpRev = -SourceSatCurrent * T1 / Nvtms;
        break;
      default:
        UserError(*this) << "Specified dioMod not matched.  dioMod = " << model_.dioMod;
    }
  }

  Nvtmd = model_.vtm * model_.DjctEmissionCoeff;
  if ((Adeff <= 0.0) && (Pdeff <= 0.0))
  {
    DrainSatCurrent = 0.0;
  }
  else
  {
    DrainSatCurrent = Adeff * model_.DjctTempSatCurDensity
                    + Pdeff * model_.DjctSidewallTempSatCurDensity
                    + paramPtr->weffCJ * nf
                    * model_.DjctGateSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent > 0.0)
  {
    switch(model_.dioMod)
    {
      case 0:
        if ((model_.bvd / Nvtmd) > CONSTEXP_THRESHOLD)
          XExpBVD = model_.xjbvd * CONSTMIN_EXP;
        else
          XExpBVD = model_.xjbvd * exp(-model_.bvd / Nvtmd);
        break;
      case 1:
        DioIjthVjmEval(Nvtmd, model_.ijthdfwd, DrainSatCurrent,
                            0.0, (vjdmFwd));
        IVjdmFwd = DrainSatCurrent * exp(vjdmFwd / Nvtmd);
        break;
      case 2:
        if ((model_.bvd / Nvtmd) > CONSTEXP_THRESHOLD)
        {   XExpBVD = model_.xjbvd * CONSTMIN_EXP;
            tmp = CONSTMIN_EXP;
        }
        else
        {   XExpBVD = exp(-model_.bvd / Nvtmd);
            tmp = XExpBVD;
            XExpBVD *= model_.xjbvd;
        }

        DioIjthVjmEval(Nvtmd, model_.ijthdfwd, DrainSatCurrent,
                            XExpBVD, (vjdmFwd));
        T0 = exp(vjdmFwd / Nvtmd);
        IVjdmFwd = DrainSatCurrent * (T0 - XExpBVD / T0 + XExpBVD - 1.0);
        DslpFwd = DrainSatCurrent * (T0 + XExpBVD / T0) / Nvtmd;

        T2 = model_.ijthdrev / DrainSatCurrent;
        if (T2 < 1.0)
        {
          T2 = 10.0;
          UserWarning(*this) << "ijthdrev too small and set to 10 times IdbSat";
        }
        vjdmRev = -model_.bvd
                           - Nvtmd * log((T2 - 1.0) / model_.xjbvd); // bugfix
        T1 = model_.xjbvd * exp(-(model_.bvd
           + vjdmRev) / Nvtmd);
        IVjdmRev = DrainSatCurrent * (1.0 + T1);
        DslpRev = -DrainSatCurrent * T1 / Nvtmd;
        break;
      default:
        UserError(*this) << "Specified dioMod not matched.  dioMod = " << model_.dioMod;
    }
  }

  // GEDL current reverse bias
  T0 = (TRatio - 1.0);
  model_.njtsstemp = model_.njts * (1.0 + model_.tnjts * T0);
  model_.njtsswstemp = model_.njtssw * (1.0 + model_.tnjtssw * T0);
  model_.njtsswgstemp = model_.njtsswg * (1.0 + model_.tnjtsswg * T0);
  model_.njtsdtemp = model_.njtsd * (1.0 + model_.tnjtsd * T0);
  model_.njtsswdtemp = model_.njtsswd * (1.0 + model_.tnjtsswd * T0);
  model_.njtsswgdtemp = model_.njtsswgd * (1.0 + model_.tnjtsswgd * T0);

  T7 = Eg0 / model_.vtm * T0;

  T9 = model_.xtss * T7;
  DEXP2(T9, T1);
  T9 = model_.xtsd * T7;
  DEXP2(T9, T2);
  T9 = model_.xtssws * T7;
  DEXP2(T9, T3);
  T9 = model_.xtsswd * T7;
  DEXP2(T9, T4);
  T9 = model_.xtsswgs * T7;
  DEXP2(T9, T5);
  T9 = model_.xtsswgd * T7;
  DEXP2(T9, T6);

  // IBM TAT
  if (model_.jtweff < 0.0)
  {
    model_.jtweff = 0.0;
    UserWarning(*this) << "TAT width dependence effect is negative. Jtweff is clamped to zero.";
  }

  T11 = sqrt(model_.jtweff / paramPtr->weffCJ) + 1.0;

  T10 = paramPtr->weffCJ * nf;
  SjctTempRevSatCur = T1 * Aseff * model_.jtss;
  DjctTempRevSatCur = T2 * Adeff * model_.jtsd;
  SswTempRevSatCur = T3 * Pseff * model_.jtssws;
  DswTempRevSatCur = T4 * Pdeff * model_.jtsswd;
  SswgTempRevSatCur = T5 * T10 * T11 * model_.jtsswgs;
  DswgTempRevSatCur = T6 * T10 * T11 * model_.jtsswgd;

  if ((model_.mtrlMod != 0) && (model_.mtrlCompatMod == 0))
  {
    /* Calculate TOXP from EOT */

    /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
    Vtm0eot = CONSTKboQ * model_.tempeot;
    Vtmeot = Vtm0eot;
    vbieot = Vtm0eot * log(paramPtr->nsd * paramPtr->ndep / (ni * ni));
    phieot = Vtm0eot * log(paramPtr->ndep / ni) + paramPtr->phin + 0.4;
    tmp2 = vfb + phieot;
    vddeot = model_.dtype * model_.vddeot;
    T0 = model_.epsrgate * CONSTEPS0;
    if ((paramPtr->ngate > 1.0e18) && (paramPtr->ngate < 1.0e25)
        && (vddeot > tmp2) && (T0!=0))
    {
      T1 = 1.0e6 * CONSTQ * T0 * paramPtr->ngate /
        (model_.coxe * model_.coxe);
      T8 = vddeot - tmp2;
      T4 = sqrt(1.0 + 2.0 * T8 / T1);
      T2 = 2.0 * T8 / (T4 + 1.0);
      T3 = 0.5 * T2 * T2 / T1;
      T7 = 1.12 - T3 - 0.05;
      T6 = sqrt(T7 * T7 + 0.224);
      T5 = 1.12 - 0.5 * (T7 + T6);
      Vgs_eff = vddeot - T5;
    }
    else
      Vgs_eff = vddeot;

    /* Calculate Vth @ Vds=Vbs=0 */
    V0 = vbieot - phieot;
    lt1 = model_.factor1* paramPtr->sqrtXdep0;
    ltw = lt1;
    T0 = paramPtr->dvt1 * model_.leffeot / lt1;
    if (T0 < CONSTEXP_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      Theta0 = T1 / T4;
    }
    else
      Theta0 = 1.0 / (CONSTMAX_EXP - 2.0);
    Delt_vth = paramPtr->dvt0 * Theta0 * V0;
    T0 = paramPtr->dvt1w * model_.weffeot * model_.leffeot / ltw;
    if (T0 < CONSTEXP_THRESHOLD)
    {   T1 = exp(T0);
      T2 = T1 - 1.0;
      T3 = T2 * T2;
      T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
      T5 = T1 / T4;
    }
    else
      T5 = 1.0 / (CONSTMAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
    T2 = paramPtr->dvt0w * T5 * V0;
    TempRatioeot =  model_.tempeot / model_.tnom - 1.0;
    T0 = sqrt(1.0 + paramPtr->lpe0 / model_.leffeot);
    T1 = paramPtr->k1ox * (T0 - 1.0) * sqrt(phieot)
      + (paramPtr->kt1 + paramPtr->kt1l / model_.leffeot) * TempRatioeot;
    Vth_NarrowW = toxe * phieot
      / (model_.weffeot + paramPtr->w0);
    Lpe_Vb = sqrt(1.0 + paramPtr->lpeb / model_.leffeot);
    Vth = model_.dtype * vth0 +
      (paramPtr->k1ox - paramPtr->k1) * sqrt(phieot) * Lpe_Vb
      - Delt_vth - T2 + paramPtr->k3 * Vth_NarrowW + T1;

    /* Calculate n */
    tmp1 = epssub / paramPtr->Xdep0;
    tmp2 = paramPtr->nfactor * tmp1;
    tmp3 = (tmp2 + paramPtr->cdsc * Theta0 + paramPtr->cit) / model_.coxe;
    if (tmp3 >= -0.5)
      n = 1.0 + tmp3;
    else
    {
      T0 = 1.0 / (3.0 + 8.0 * tmp3);
			n = (1.0 + 3.0 * tmp3) * T0;
    }

    /* Vth correction for Pocket implant */
    if (paramPtr->dvtp0 > 0.0)
    {
      T3 = model_.leffeot + paramPtr->dvtp0 * 2.0;
      if (model_.tempMod < 2)
        T4 = Vtmeot * log(model_.leffeot / T3);
      else
        T4 = Vtm0eot * log(model_.leffeot / T3);
      Vth -= n * T4;
    }
    Vgsteff = Vgs_eff-Vth;
    /* calculating Toxp */
    T3 = model_.dtype * vth0 - vfb - phieot;
    T4 = T3 + T3;
    T5 = 2.5 * T3;

    vtfbphi2eot = 4.0 * T3;
    if (vtfbphi2eot < 0.0)
      vtfbphi2eot = 0.0;

    niter = 0;
    toxpf = toxe;
    do
    {
      toxpi = toxpf;
      tmp2 = 2.0e8 * toxpf;
      T0 = (Vgsteff + vtfbphi2eot) / tmp2;
      T1 = 1.0 + exp(model_.bdos * 0.7 * log(T0));
      Tcen = model_.ados * 1.9e-9 / T1;
      toxpf = toxe - epsrox/model_.epsrsub * Tcen;
      niter++;
    } while ((niter<=4)&&(fabs(toxpf-toxpi)>1e-12));
    toxp = toxpf;
    coxp = epsrox * CONSTEPS0 / model_.toxp;
  }
  else
  {
    toxp=model_.toxp;
    coxp=model_.coxp;
  }

  ///////////////////////////////////////////////////////////////////////////////

  updateTemperatureCalled_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars4p82_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars4p82_ ()
{
  bool bsuccess = true;

  // begin the b4ld.c parameters:
  double dgstot_dvd(0.0), dgstot_dvg(0.0), dgstot_dvs(0.0), dgstot_dvb(0.0);
  double dgdtot_dvd(0.0), dgdtot_dvg(0.0), dgdtot_dvs(0.0), dgdtot_dvb(0.0);
  double Rs(0.0), Rd(0.0);
  double dRs_dvg(0.0), dRd_dvg(0.0),  dRs_dvb(0.0), dRd_dvb(0.0);
  double dT0_dvg(0.0), dT1_dvb(0.0), dT3_dvg(0.0), dT3_dvb(0.0);

  double vgd_old(0.0), vsbd_old(0.0), vsbd(0.0);
  double SourceSatCurrent(0.0), DrainSatCurrent(0.0);
  double VgstNVt(0.0), ExpVgst(0.0);
  double czbd(0.0), czbdsw(0.0), czbdswg(0.0),
         czbs(0.0), czbssw(0.0), czbsswg(0.0);
  double evbd(0.0), evbs(0.0),  arg(0.0), sarg(0.0);
  double Vfbeff(0.0), dVfbeff_dVg(0.0), dVfbeff_dVb(0.0), V3(0.0), V4(0.0);
  double MJD(0.0), MJSWD(0.0), MJSWGD(0.0);
  double MJS(0.0), MJSWS(0.0), MJSWGS(0.0);
  double Ggidld(0.0), Ggidlg(0.0), Ggidlb(0.0);
  double Voxacc(0.0), dVoxacc_dVg(0.0), dVoxacc_dVb(0.0);
  double Voxdepinv(0.0);
  double dVoxdepinv_dVg(0.0), dVoxdepinv_dVd(0.0), dVoxdepinv_dVb(0.0);
  double VxNVt(0.0), ExpVxNVt(0.0);
  double Vaux(0.0), dVaux_dVg(0.0), dVaux_dVd(0.0), dVaux_dVb(0.0);
  double Igc(0.0), dIgc_dVg(0.0), dIgc_dVd(0.0), dIgc_dVb(0.0);
  double dIgcs_dVg(0.0), dIgcs_dVd(0.0), dIgcs_dVb(0.0);
  double  dIgcd_dVg(0.0), dIgcd_dVd(0.0), dIgcd_dVb(0.0);
  double dIgs_dVg(0.0), dIgs_dVs(0.0), dIgd_dVg(0.0), dIgd_dVd(0.0);
  double Igbacc(0.0), dIgbacc_dVg(0.0);
  double dIgbacc_dVb(0.0);
  double Igbinv(0.0), dIgbinv_dVg(0.0), dIgbinv_dVd(0.0), dIgbinv_dVb(0.0);
  double Pigcd(0.0), dPigcd_dVg(0.0), dPigcd_dVd(0.0), dPigcd_dVb(0.0);

  double Vgs_eff(0.0), Vfb(0.0);
  double Vth_NarrowW(0.0);
  double Phis(0.0), dPhis_dVb(0.0), sqrtPhis(0.0), dsqrtPhis_dVb(0.0);
  //double Vth(0.0);  // converted to instance variable
  double dVth_dVb(0.0), dVth_dVd(0.0);
  double Vgst(0.0), dVgst_dVg(0.0), dVgst_dVb(0.0), dVgs_eff_dVg(0.0);
  double Nvtms(0.0), Nvtmd(0.0);
  double Vtm(0.0), Vtm0(0.0);
  double n(0.0), dn_dVb(0.0), dn_dVd(0.0), voffcv (0.0);
  double noff(0.0), dnoff_dVd(0.0), dnoff_dVb(0.0);
  double CoxWLcen(0.0), QovCox(0.0), LINK(0.0), V0(0.0);
  double DeltaPhi(0.0), dDeltaPhi_dVg(0.0), VgDP(0.0), dVgDP_dVg(0.0);
  double Cox(0.0), Tox(0.0);
  double Tcen(0.0), dTcen_dVg(0.0), dTcen_dVd(0.0), dTcen_dVb(0.0);
  double Ccen(0.0);
  //double Coxeff(0.0); // made into an instance variable
  double dCoxeff_dVd(0.0), dCoxeff_dVg(0.0), dCoxeff_dVb(0.0);
  double Denomi(0.0), dDenomi_dVg(0.0), dDenomi_dVd(0.0), dDenomi_dVb(0.0);
  double dueff_dVg(0.0), dueff_dVd(0.0), dueff_dVb(0.0);
  double Esat(0.0);
  //double Vdsat(0.0); // made into an instance variable
  double dEsatL_dVg(0.0), dEsatL_dVd(0.0), dEsatL_dVb(0.0);
  double dVdsat_dVg(0.0), dVdsat_dVb(0.0);
  double dVdsat_dVd(0.0), Vasat(0.0), dAlphaz_dVg(0.0), dAlphaz_dVb(0.0);
  double dVasat_dVg(0.0), dVasat_dVb(0.0);
  double dVasat_dVd(0.0), Va(0.0), dVa_dVd(0.0), dVa_dVg(0.0), dVa_dVb(0.0);
  double Vbseff(0.0), dVbseff_dVb(0.0), VbseffCV(0.0), dVbseffCV_dVb(0.0);
  double VgsteffVth(0.0), dT11_dVg(0.0);
  double Arg1(0.0), One_Third_CoxWL(0.0), Two_Third_CoxWL(0.0), Alphaz(0.0);

  double T0,dT0_dVg(0.0), dT0_dVd(0.0), dT0_dVb(0.0);
  double T1,dT1_dVg(0.0), dT1_dVd(0.0), dT1_dVb(0.0);

  double T2(0.0), dT2_dVg(0.0), dT2_dVd(0.0), dT2_dVb(0.0);
  double T3(0.0), dT3_dVg(0.0), dT3_dVd(0.0), dT3_dVb(0.0);
  double T4(0.0),               dT4_dVd(0.0), dT4_dVb(0.0);
  double T5(0.0), dT5_dVg(0.0), dT5_dVd(0.0), dT5_dVb(0.0);
  double T6(0.0), dT6_dVg(0.0), dT6_dVd(0.0), dT6_dVb(0.0);
  double T7(0.0), dT7_dVg(0.0), dT7_dVd(0.0), dT7_dVb(0.0);
  double T8(0.0), dT8_dVg(0.0), dT8_dVd(0.0), dT8_dVb(0.0);
  double T9(0.0), dT9_dVg(0.0), dT9_dVd(0.0), dT9_dVb(0.0);
  double T10(0.0), dT10_dVg(0.0), dT10_dVb(0.0), dT10_dVd(0.0);
  double T11(0.0), T12(0.0), T13(0.0), T14(0.0);
  double tmp(0.0);
  double dAbulk_dVb(0.0), Abulk0(0.0), dAbulk0_dVb(0.0);
  double Cclm(0.0), dCclm_dVg(0.0), dCclm_dVd(0.0), dCclm_dVb(0.0);
  double FP(0.0), dFP_dVg(0.0);
  double PvagTerm(0.0), dPvagTerm_dVg(0.0);
  double dPvagTerm_dVd(0.0), dPvagTerm_dVb(0.0);
  double VADITS(0.0), dVADITS_dVg(0.0), dVADITS_dVd(0.0);
  double Lpe_Vb(0.0);
  double dDITS_Sft_dVb(0.0), dDITS_Sft_dVd(0.0);
  double DITS_Sft2, dDITS_Sft2_dVd;
  double VACLM(0.0), dVACLM_dVg(0.0), dVACLM_dVd(0.0), dVACLM_dVb(0.0);
  double VADIBL(0.0), dVADIBL_dVg(0.0), dVADIBL_dVd(0.0), dVADIBL_dVb(0.0);
  double Xdep(0.0), dXdep_dVb(0.0);
  double lt1(0.0), dlt1_dVb(0.0), ltw(0.0), dltw_dVb(0.0);
  double Delt_vth(0.0), dDelt_vth_dVb(0.0);
  double Theta0(0.0), dTheta0_dVb(0.0);

  double TempRatio(0.0), tmp1(0.0), tmp2(0.0), tmp3(0.0), tmp4(0.0);
  double DIBL_Sft(0.0), dDIBL_Sft_dVd(0.0);
  double Lambda(0.0), dLambda_dVg(0.0);
  double a1(0.0);

  double dVgsteff_dVg(0.0), dVgsteff_dVd(0.0), dVgsteff_dVb(0.0);
  double dVdseff_dVg(0.0), dVdseff_dVd(0.0), dVdseff_dVb(0.0);
  double VdseffCV(0.0),
         dVdseffCV_dVg(0.0), dVdseffCV_dVd(0.0), dVdseffCV_dVb(0.0);
  double diffVds(0.0);
  double dAbulk_dVg(0.0);
  double beta(0.0), dbeta_dVg(0.0), dbeta_dVd(0.0), dbeta_dVb(0.0);
  double gche(0.0), dgche_dVg(0.0), dgche_dVd(0.0), dgche_dVb(0.0);
  double fgche1(0.0), dfgche1_dVg(0.0), dfgche1_dVd(0.0), dfgche1_dVb(0.0);
  double fgche2(0.0), dfgche2_dVg(0.0), dfgche2_dVd(0.0), dfgche2_dVb(0.0);
  double Idl(0.0), dIdl_dVg(0.0), dIdl_dVd(0.0), dIdl_dVb(0.0);
  double Idsa(0.0), dIdsa_dVg(0.0), dIdsa_dVd(0.0), dIdsa_dVb(0.0);
  double Ids(0.0), Gmb(0.0);
  double devbs_dvb(0.0), devbd_dvb(0.0);
  double Isub(0.0), Gbd(0.0), Gbg(0.0), Gbb(0.0), Gds(0.0);
  double VASCBE(0.0), dVASCBE_dVg(0.0), dVASCBE_dVd(0.0), dVASCBE_dVb(0.0);
  double CoxeffWovL(0.0);
  double Rds(0.0), dRds_dVg(0.0), dRds_dVb(0.0), WVCox(0.0), WVCoxRds(0.0);
  double Vgst2Vtm(0.0), VdsatCV(0.0);
  double Leff(0.0), Weff(0.0), dWeff_dVg(0.0), dWeff_dVb(0.0);
  double AbulkCV(0.0), dAbulkCV_dVb(0.0);

  double Cgg1(0.0), Cgb1(0.0), Cgd1(0.0), Cbg1(0.0), Cbb1(0.0), Cbd1(0.0);
  double Qac0(0.0), Qsub0(0.0);
  double dQac0_dVg(0.0), dQac0_dVb(0.0);
  double dQsub0_dVg(0.0), dQsub0_dVd(0.0), dQsub0_dVb(0.0);
  double Ggislg(0.0), Ggislb(0.0), Ggisls(0.0);
  double Nvtmrss(0.0), Nvtmrssws(0.0), Nvtmrsswgs(0.0);
  double Nvtmrsd(0.0), Nvtmrsswd(0.0), Nvtmrsswgd(0.0);

  double vs(0.0), Fsevl(0.0);
  double dvs_dVg(0.0), dvs_dVd(0.0), dvs_dVb(0.0), dFsevl_dVg(0.0);
  double dFsevl_dVd(0.0), dFsevl_dVb(0.0);
  double vgdx(0.0), vgsx(0.0),epssub(0.0),toxe(0.0),epsrox(0.0);
  double von_local(0.0);

  // end b4ld.c parameters

  ScalingFactor = 1.0e-9;

  // Don't do charge computations in DC sweeps.
  if (getSolverState().tranopFlag || getSolverState().acopFlag || getSolverState().transientFlag)
  {
    ChargeComputationNeeded = true;
  }
  else
  {
    ChargeComputationNeeded = false;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag && isActive(Diag::TIME_STEP) )
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::updateIntermediateVars\n";
    Xyce::dout() << "  name = " << getName();
    Xyce::dout() << "  model name = " << model_.getName();
    Xyce::dout() <<"   dtype is " << model_.dtype << std::endl;
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "  " << std::endl;
  }

  int Check = 0;
  int Check1 = 0;
  int Check2 = 0;

  limitedFlag=false;

  // The first block of code in b4ld.c basically sets up, locally,
  // what the load function should use as values for the various solution
  // variables.  There is a series of IF statements which are dependent
  // upon the mode.  (transient, initializing transient, operating point,
  // small signal, etc.).  Xyce treats the operating point and transient
  // calculation in the same way, from the device's point of view, and
  // we don't support any of the other modes.  Therefore most of these
  // mode options are not here - only the transient mode stuff.

  // First get some of the needed solution variables:
  Vd     = 0.0;
  Vs     = 0.0;
  Vb     = 0.0;
  Vsp    = 0.0;
  Vdp    = 0.0;
  Vgp    = 0.0;
  Vbp    = 0.0;
  Vge    = 0.0;
  Vgm    = 0.0;
  Vdb    = 0.0;
  Vsb    = 0.0;
  Qtotal = 0.0;

  Vd  = (extData.nextSolVectorRawPtr)[li_Drain] ;
  Vs  = (extData.nextSolVectorRawPtr)[li_Source] ;
  Vb  = (extData.nextSolVectorRawPtr)[li_Body] ;
  Vsp = (extData.nextSolVectorRawPtr)[li_SourcePrime] ;
  Vdp = (extData.nextSolVectorRawPtr)[li_DrainPrime] ;

  Vgp = (extData.nextSolVectorRawPtr)[li_GatePrime] ;
  Vbp = (extData.nextSolVectorRawPtr)[li_BodyPrime] ;
  Vge = (extData.nextSolVectorRawPtr)[li_GateExt] ;

  if (li_GateMid >= 0) // only true for rgateMod==3
  {
    Vgm = (extData.nextSolVectorRawPtr)[li_GateMid] ;
  }

  Vdb = (extData.nextSolVectorRawPtr)[li_DrainBody] ;
  Vsb = (extData.nextSolVectorRawPtr)[li_SourceBody] ;

  if (trnqsMod)
  {
    Qtotal = (extData.nextSolVectorRawPtr)[li_Charge];
  }
  else
  {
    Qtotal = 0.0;
  }

  Vddp  = Vd   - Vdp;
  Vssp  = Vs   - Vsp;
  //Vbsp  = Vb   - Vsp;
  //Vbdp  = Vb   - Vdp;
  //Vgsp  = Vg   - Vsp;
  //Vgdp  = Vg   - Vdp;
  //Vgb   = Vg   - Vb;

  //Vdpsp = Vdp  - Vsp;

  // substrate network:
  Vdbb  = Vdb - Vb;
  Vdbbp = Vdb - Vbp;
  Vsbb  = Vsb - Vb;
  Vsbbp = Vsb - Vbp;
  Vbpb  = Vbp - Vb;

  // modified from b4ld:
  vds  = model_.dtype * (Vdp - Vsp);
  vgs  = model_.dtype * (Vgp - Vsp);
  vbs  = model_.dtype * (Vbp - Vsp);
  vges = model_.dtype * (Vge - Vsp);
  vgms = model_.dtype * (Vgm - Vsp);
  vdbs = model_.dtype * (Vdb - Vsp);
  vsbs = model_.dtype * (Vsb - Vsp);
  vses = model_.dtype * (Vs - Vsp);
  vdes = model_.dtype * (Vd - Vsp);
  qdef = model_.dtype * (Qtotal);

  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;
  vged = vges - vds;
  vgmd = vgms - vds;
  vgmb = vgms - vbs;
  vdbd = vdbs - vds;

  vbs_jct = (!rbodyMod) ? vbs : vsbs;
  vbd_jct = (!rbodyMod) ? vbd : vdbd;

  // Set up the linear resistors.  We need the exact drops, not munged by
  // type.
  Vgegp = Vge - Vgp;
  Vgegm = Vge - Vgm;
  Vgmgp = Vgm - Vgp;

  origFlag = 1;

  vbd_orig = vbd;
  vbs_orig = vbs;
  vgs_orig = vgs;
  vds_orig = vds;
  vgd_orig = vgd;
  vges_orig = vges;
  vgms_orig = vgms;
  vdes_orig = vdes;
  vses_orig = vses;
  vdbs_orig = vdbs;
  vsbs_orig = vsbs;
  vdbd_orig = vdbd;
  vged_orig = vged;
  vgmd_orig = vgmd;
  vbs_jct_orig = vbs_jct;
  vbd_jct_orig = vbd_jct;
  vgmb_orig = vgmb;
  vgb_orig = vgb;

  // What follows is a block of code designed to impose some  limits,
  //  or initial conditions on the junction voltages.  Initial conditions
  //  should only be imposed on the first Newton step of an operating point.
  //
  // The first possible limit on the  junction voltages has to do with
  // limiting the percent change of junction voltages between  Newton
  // iterations.  The second has to do with avoiding extra floating point
  // operations in the event that the device has in some sense converged
  // (aka BYPASS).  Although the primary point of BYPASS is to reduce
  // neccessary work, it also seems to reduce the number of Newton iterations.
  //
  // NOTE:  We do not support BYPASS.
  //
  // The "old" variables should be the values for the previous
  // Newton iteration, if indeed there was a previous Newton
  // iteration.  If not, just set the  old values equal to
  // the current ones.
  //

  // set an initial condition if appropriate:
  if (getSolverState().initJctFlag_ && !OFF && getDeviceOptions().voltageLimiterFlag)
  {
    if (getSolverState().inputOPFlag)
    {
      Linear::Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
      if ((*flagSolVectorPtr)[li_Drain]       == 0 ||
          (*flagSolVectorPtr)[li_GateExt]     == 0 ||
          (*flagSolVectorPtr)[li_Source]      == 0 ||
          (*flagSolVectorPtr)[li_Body]        == 0 ||
          (*flagSolVectorPtr)[li_DrainPrime]  == 0 ||
          (*flagSolVectorPtr)[li_GatePrime]   == 0 ||
          (*flagSolVectorPtr)[li_GateMid]     == 0 ||
          (*flagSolVectorPtr)[li_SourcePrime] == 0 ||
          (*flagSolVectorPtr)[li_BodyPrime]   == 0 ||
          (*flagSolVectorPtr)[li_DrainBody]   == 0 ||
          (*flagSolVectorPtr)[li_SourceBody]  == 0 ||
          (*flagSolVectorPtr)[li_Charge]      == 0 )
      {
        vds = 0.1;
        vdes = 0.11;
        vses = -0.01;
        vgs = vges = vgms = model_.dtype * vth0 + 0.1;
        origFlag = 0;
      }
    }
    else
    {
      vds = 0.1;
      vdes = 0.11;
      vses = -0.01;
      vgs = vges = vgms = model_.dtype * vth0 + 0.1;
      origFlag = 0;
    }
    vbs = vdbs = vsbs = 0.0;
    vbd = vbs - vds;
    vdbd = vdbs - vds;
    vgd = vgs - vds;
    vged = vges-vds;
    vgmd = vgms-vds;
    //origFlag = 0;
  }
  else if ((getSolverState().initFixFlag || getSolverState().initJctFlag_) && OFF)
  {
    vds = vgs = vbs = vges = vgms = 0.0;
    vds = vsbs = vdes = vses = qdef = 0.0;
  }


  if (getSolverState().newtonIter == 0)
  {

    if (!getSolverState().dcopFlag || (getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    // ie, first newton step of a transient time step or DCOP continuation step.
    {
      vbd_old = (extData.currStoVectorRawPtr)[li_store_vbd];
      vbs_old = (extData.currStoVectorRawPtr)[li_store_vbs];
      vgs_old = (extData.currStoVectorRawPtr)[li_store_vgs];
      vds_old = (extData.currStoVectorRawPtr)[li_store_vds];
      vges_old = (extData.currStoVectorRawPtr)[li_store_vges];
      vgms_old = (extData.currStoVectorRawPtr)[li_store_vgms];
      vdes_old = (extData.currStoVectorRawPtr)[li_store_vdes];
      vses_old = (extData.currStoVectorRawPtr)[li_store_vses];
      vdbs_old = (extData.currStoVectorRawPtr)[li_store_vdbs];
      vsbs_old = (extData.currStoVectorRawPtr)[li_store_vsbs];
      vdbd_old = (extData.currStoVectorRawPtr)[li_store_vdbd];
      vged_old = (extData.currStoVectorRawPtr)[li_store_vged];
      vgmd_old = (extData.currStoVectorRawPtr)[li_store_vgmd];
    }
    else
    {  // no history
      vbd_old = vbd;
      vbs_old = vbs;
      vgs_old = vgs;
      vds_old = vds;
      vges_old = vges;
      vgms_old = vgms;
      vdes_old = vdes;
      vses_old = vses;
      vdbs_old = vdbs;
      vsbs_old = vsbs;
      vdbd_old = vdbd;
      vged_old = vged;
      vgmd_old = vgmd;
    }
  }
  else
  {
    vbd_old = (extData.nextStoVectorRawPtr)[li_store_vbd];
    vbs_old = (extData.nextStoVectorRawPtr)[li_store_vbs];
    vgs_old = (extData.nextStoVectorRawPtr)[li_store_vgs];
    vds_old = (extData.nextStoVectorRawPtr)[li_store_vds];
    vges_old = (extData.nextStoVectorRawPtr)[li_store_vges];
    vgms_old = (extData.nextStoVectorRawPtr)[li_store_vgms];
    vdes_old = (extData.nextStoVectorRawPtr)[li_store_vdes];
    vses_old = (extData.nextStoVectorRawPtr)[li_store_vses];
    vdbs_old = (extData.nextStoVectorRawPtr)[li_store_vdbs];
    vsbs_old = (extData.nextStoVectorRawPtr)[li_store_vsbs];
    vdbd_old = (extData.nextStoVectorRawPtr)[li_store_vdbd];
    vged_old = (extData.nextStoVectorRawPtr)[li_store_vged];
    vgmd_old = (extData.nextStoVectorRawPtr)[li_store_vgmd];
  }

  vgd_old = vgs_old - vds_old;

  // This next block performs checks on the junction voltages and
  // imposes limits on them if they are too big.
  // Note:  In the level=1 von is multiplied by dtype.  Here it is not.  They
  // are both right.

  if (getDeviceOptions().voltageLimiterFlag && !(getSolverState().initFixFlag && OFF))
  {
    // only do this if we are beyond the first Newton iteration.  On the
    // first newton iteration, the "old" values are from a previous time
    // step.
    von_local=von;
    if (getSolverState().newtonIter >= 0 && !(getSolverState().initJctFlag_))
    {
      if (vds_old >= 0.0)
      {
        vgs = devSupport.fetlim(vgs, vgs_old, von_local);
        vds = vgs - vgd;
        vds = devSupport.limvds(vds, vds_old);
        vgd = vgs - vds;
        if (rgateMod == 3)
        {
         vges = devSupport.fetlim(vges, vges_old, von_local);
         vgms = devSupport.fetlim(vgms, vgms_old, von_local);
         vged = vges - vds;
         vgmd = vgms - vds;
        }
        else if ((rgateMod == 1) || (rgateMod == 2))
        {
          vges = devSupport.fetlim(vges, vges_old, von_local);
          vged = vges - vds;
        }

        if (model_.rdsMod)
        {
          vdes = devSupport.limvds(vdes, vdes_old);
          vses = -devSupport.limvds(-vses, -vses_old);
        }
      }
      else
      {
        vgd = devSupport.fetlim(vgd, vgd_old, von_local);
        vds = vgs - vgd;
        vds = -devSupport.limvds(-vds, -vds_old);
        vgs = vgd + vds;

        if (rgateMod == 3)
        {
          vged = devSupport.fetlim(vged, vged_old, von_local);
          vges = vged + vds;
          vgmd = devSupport.fetlim(vgmd, vgmd_old, von_local);
          vgms = vgmd + vds;
        }
        if ((rgateMod == 1) || (rgateMod == 2))
        {
          vged = devSupport.fetlim(vged, vged_old, von_local);
          vges = vged + vds;
        }

        if (model_.rdsMod)
        {
          vdes = -devSupport.limvds(-vdes, -vdes_old);
          vses = devSupport.limvds(vses, vses_old);
        }
      }

      if (vds >= 0.0)
      {
        vbs = devSupport.pnjlim(vbs, vbs_old,
                 CONSTvt0, model_.vcrit, &Check);
        vbd = vbs - vds;
        if (rbodyMod)
        {
          vdbs = devSupport.pnjlim(vdbs, vdbs_old,
                      CONSTvt0, model_.vcrit, &Check1);
          vdbd = vdbs - vds;
          vsbs = devSupport.pnjlim(vsbs, vsbs_old,
                      CONSTvt0, model_.vcrit, &Check2);
          if ((Check1 != 0) || (Check2 != 0))
          {
            Check = 1;
          }
        }
      }
      else
      {
        vbd = devSupport.pnjlim(vbd, vbd_old,
                   CONSTvt0, model_.vcrit, &Check);
        vbs = vbd + vds;
        if (rbodyMod)
        {
          vdbd = devSupport.pnjlim(vdbd, vdbd_old,
                        CONSTvt0, model_.vcrit, &Check1);
          vdbs = vdbd + vds;
          vsbd_old = vsbs_old - vds_old;
          vsbd = vsbs - vds;
          vsbd = devSupport.pnjlim(vsbd, vsbd_old, CONSTvt0, model_.vcrit, &Check2);
          vsbs = vsbd + vds;
          if ((Check1 != 0) || (Check2 != 0))
          {
            Check = 1;
          }
        }
      }
    }

    // for convergence testing:
    if (Check == 1) limitedFlag=true;

  } // getDeviceOptions().voltageLimiterFlag

  // Calculate DC currents and their derivatives
  vbd = vbs - vds;
  vgd = vgs - vds;
  vgb = vgs - vbs;
  vged = vges - vds;
  vgmd = vgms - vds;
  vgmb = vgms - vbs;
  vdbd = vdbs - vds;

  vbs_jct = (!rbodyMod) ? vbs : vsbs;
  vbd_jct = (!rbodyMod) ? vbd : vdbd;

  if (getDeviceOptions().voltageLimiterFlag)
  {
    double machprec= Util::MachineDependentParams::MachinePrecision();

    if (
       fabs( vbs_orig - vbs) > machprec ||
       fabs( vgs_orig - vgs) > machprec ||
       fabs( vds_orig - vds) > machprec ||
       fabs( vges_orig - vges) > machprec ||
       fabs( vgms_orig - vgms) > machprec ||
       fabs( vdes_orig - vdes) > machprec ||
       fabs( vses_orig - vses) > machprec ||
       fabs( vdbs_orig - vdbs) > machprec ||
       fabs( vsbs_orig - vsbs) > machprec 
       )
    {
      origFlag = 0;
    }
  }


  if (DEBUG_DEVICE && getDeviceOptions().voltageLimiterFlag)
  {
    if (isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      if (!origFlag)
      {
        double machprec = Util::MachineDependentParams::MachinePrecision();
#define SET_UPDATE(VAR) double VAR ## _update = VAR ## _orig -  VAR;
#define ORIG_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR ## _orig>=0)?"+":"")  << VAR ## _orig;
#define LIMIT_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR >=0)?"+":"")  << VAR;
#define DIFF_OUTPUT(VAR) if ( fabs(VAR ## _update) > machprec ) Xyce::dout()<<" "#VAR ":"	  << ((VAR ## _update >=0)?"+":"")  << VAR ## _update;

        SET_UPDATE(vbd) SET_UPDATE(vbs) SET_UPDATE(vgs) SET_UPDATE(vds) SET_UPDATE(vgd) SET_UPDATE(vges)
        SET_UPDATE(vgms) SET_UPDATE(vdes) SET_UPDATE(vses) SET_UPDATE(vdbs) SET_UPDATE(vsbs) SET_UPDATE(vdbd)
        SET_UPDATE(vbs_jct) SET_UPDATE(vbd_jct) SET_UPDATE(vgmb) SET_UPDATE(vged) SET_UPDATE(vgmd)
        SET_UPDATE(vgb)

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Blim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        LIMIT_OUTPUT(vbd) LIMIT_OUTPUT(vbs) LIMIT_OUTPUT(vgs) LIMIT_OUTPUT(vds) LIMIT_OUTPUT(vgd) LIMIT_OUTPUT(vges)
        LIMIT_OUTPUT(vgms) LIMIT_OUTPUT(vdes) LIMIT_OUTPUT(vses) LIMIT_OUTPUT(vdbs) LIMIT_OUTPUT(vsbs) LIMIT_OUTPUT(vdbd)
        LIMIT_OUTPUT(vbs_jct) LIMIT_OUTPUT(vbd_jct) LIMIT_OUTPUT(vgmb) LIMIT_OUTPUT(vged) LIMIT_OUTPUT(vgmd)
        LIMIT_OUTPUT(vgb)
        Xyce::dout() << std::endl;

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Alim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        LIMIT_OUTPUT(vbd) LIMIT_OUTPUT(vbs) LIMIT_OUTPUT(vgs) LIMIT_OUTPUT(vds) LIMIT_OUTPUT(vgd) LIMIT_OUTPUT(vges)
        LIMIT_OUTPUT(vgms) LIMIT_OUTPUT(vdes) LIMIT_OUTPUT(vses) LIMIT_OUTPUT(vdbs) LIMIT_OUTPUT(vsbs) LIMIT_OUTPUT(vdbd)
        LIMIT_OUTPUT(vbs_jct) LIMIT_OUTPUT(vbd_jct) LIMIT_OUTPUT(vgmb) LIMIT_OUTPUT(vged) LIMIT_OUTPUT(vgmd)
        LIMIT_OUTPUT(vgb)
        Xyce::dout()  <<std::endl;

        Xyce::dout().width(3);
        Xyce::dout() << getSolverState().newtonIter;
        Xyce::dout().width(5);Xyce::dout() << getName();
        Xyce::dout() << " Dlim:";
        Xyce::dout().width(12); Xyce::dout().precision(4); Xyce::dout().setf(std::ios::scientific);
        DIFF_OUTPUT(vbd) DIFF_OUTPUT(vbs) DIFF_OUTPUT(vgs) DIFF_OUTPUT(vds) DIFF_OUTPUT(vgd) DIFF_OUTPUT(vges)
        DIFF_OUTPUT(vgms) DIFF_OUTPUT(vdes) DIFF_OUTPUT(vses) DIFF_OUTPUT(vdbs) DIFF_OUTPUT(vsbs) DIFF_OUTPUT(vdbd)
        DIFF_OUTPUT(vbs_jct) DIFF_OUTPUT(vbd_jct) DIFF_OUTPUT(vgmb) DIFF_OUTPUT(vged) DIFF_OUTPUT(vgmd)
        DIFF_OUTPUT(vgb)
        Xyce::dout()  <<std::endl;
      }
    }
  }

  // Source/drain junction diode DC model begins
  Nvtms = model_.vtm * model_.SjctEmissionCoeff;
  if ((Aseff <= 0.0) && (Pseff <= 0.0))
  {
    SourceSatCurrent = 0.0;
  }
  else
  {
    SourceSatCurrent = Aseff * model_.SjctTempSatCurDensity
                     + Pseff * model_.SjctSidewallTempSatCurDensity
                     + paramPtr->weffCJ * nf
                     * model_.SjctGateSidewallTempSatCurDensity;
  }

  if (SourceSatCurrent <= 0.0)
  {
    gbs = getDeviceOptions().gmin;
    cbs = gbs * vbs_jct;
  }
  else
  {
    switch(model_.dioMod)
    {
      case 0:
          evbs = exp(vbs_jct / Nvtms);
          T1 = model_.xjbvs * exp(-(model_.bvs + vbs_jct) / Nvtms);
          // WDLiu: Magic T1 in this form; different from  beta.
          gbs = SourceSatCurrent * (evbs + T1) / Nvtms + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs + XExpBVS
                                    - T1 - 1.0) + getDeviceOptions().gmin * vbs_jct;
          break;
      case 1:
        T2 = vbs_jct / Nvtms;
        if (T2 < -CONSTEXP_THRESHOLD)
        {
          gbs = getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (CONSTMIN_EXP - 1.0)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else if (vbs_jct <= vjsmFwd)
        {
          evbs = exp(T2);
          gbs = SourceSatCurrent * evbs / Nvtms + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs - 1.0)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else
        {
          T0 = IVjsmFwd / Nvtms;
          gbs = T0 + getDeviceOptions().gmin;
          cbs = IVjsmFwd - SourceSatCurrent + T0
                * (vbs_jct - vjsmFwd) + getDeviceOptions().gmin * vbs_jct;
        }
        break;
      case 2:
        if (vbs_jct < vjsmRev)
        {
          T0 = vbs_jct / Nvtms;
          if (T0 < -CONSTEXP_THRESHOLD)
          {
            evbs = CONSTMIN_EXP;
            devbs_dvb = 0.0;
          }
          else
          {
            evbs = exp(T0);
            devbs_dvb = evbs / Nvtms;
          }

          T1 = evbs - 1.0;
          T2 = IVjsmRev + SslpRev * (vbs_jct - vjsmRev);
          gbs = devbs_dvb * T2 + T1 * SslpRev + getDeviceOptions().gmin;
          cbs = T1 * T2 + getDeviceOptions().gmin * vbs_jct;
        }
        else if (vbs_jct <= vjsmFwd)
        {
          T0 = vbs_jct / Nvtms;
          if (T0 < -CONSTEXP_THRESHOLD)
          {
            evbs = CONSTMIN_EXP;
            devbs_dvb = 0.0;
          }
          else
          {
            evbs = exp(T0);
            devbs_dvb = evbs / Nvtms;
          }

          T1 = (model_.bvs + vbs_jct) / Nvtms;
          if (T1 > CONSTEXP_THRESHOLD)
          {
            T2 = CONSTMIN_EXP;
            T3 = 0.0;
          }
          else
          {
            T2 = exp(-T1);
            T3 = -T2 /Nvtms;
          }
          gbs = SourceSatCurrent * (devbs_dvb - model_.xjbvs * T3)
                + getDeviceOptions().gmin;
          cbs = SourceSatCurrent * (evbs + XExpBVS - 1.0
                                    - model_.xjbvs * T2)
                + getDeviceOptions().gmin * vbs_jct;
        }
        else
        {
          gbs = SslpFwd + getDeviceOptions().gmin;
          cbs = IVjsmFwd + SslpFwd * (vbs_jct
                                      - vjsmFwd) + getDeviceOptions().gmin * vbs_jct;
        }
        break;
    default: break;
    }
  }

  Nvtmd = model_.vtm * model_.DjctEmissionCoeff;
  if ((Adeff <= 0.0) && (Pdeff <= 0.0))
  {
    DrainSatCurrent = 0.0;
  }
  else
  {
    DrainSatCurrent = Adeff * model_.DjctTempSatCurDensity
                    + Pdeff * model_.DjctSidewallTempSatCurDensity
                    + paramPtr->weffCJ * nf
                    * model_.DjctGateSidewallTempSatCurDensity;
  }

  if (DrainSatCurrent <= 0.0)
  {
    gbd = getDeviceOptions().gmin;
    cbd = gbd * vbd_jct;
  }
  else
  {
    switch(model_.dioMod)
    {
      case 0:
          evbd = exp(vbd_jct / Nvtmd);
          T1 = model_.xjbvd * exp(-(model_.bvd + vbd_jct) / Nvtmd);
          // WDLiu: Magic T1 in this form; different from  beta.
          gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + getDeviceOptions().gmin;
          cbd = DrainSatCurrent * (evbd + XExpBVD
                                   - T1 - 1.0) + getDeviceOptions().gmin * vbd_jct;
          break;
      case 1:
          T2 = vbd_jct / Nvtmd;
          if (T2 < -CONSTEXP_THRESHOLD)
          {
            gbd = getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (CONSTMIN_EXP - 1.0)
                  + getDeviceOptions().gmin * vbd_jct;
          }
          else if (vbd_jct <= vjdmFwd)
          {
            evbd = exp(T2);
            gbd = DrainSatCurrent * evbd / Nvtmd + getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (evbd - 1.0)
                  + getDeviceOptions().gmin * vbd_jct;
          }
          else
          {
            T0 = IVjdmFwd / Nvtmd;
            gbd = T0 + getDeviceOptions().gmin;
            cbd = IVjdmFwd - DrainSatCurrent + T0
                  * (vbd_jct - vjdmFwd) + getDeviceOptions().gmin * vbd_jct;
          }
          break;
      case 2:
          if (vbd_jct < vjdmRev)
          {
            T0 = vbd_jct / Nvtmd;
            if (T0 < -CONSTEXP_THRESHOLD)
            {
              evbd = CONSTMIN_EXP;
              devbd_dvb = 0.0;
            }
            else
            {
              evbd = exp(T0);
              devbd_dvb = evbd / Nvtmd;
            }

            T1 = evbd - 1.0;
            T2 = IVjdmRev + DslpRev * (vbd_jct - vjdmRev);
            gbd = devbd_dvb * T2 + T1 * DslpRev + getDeviceOptions().gmin;
            cbd = T1 * T2 + getDeviceOptions().gmin * vbd_jct;
          }
          else if (vbd_jct <= vjdmFwd)
          {
            T0 = vbd_jct / Nvtmd;
            if (T0 < -CONSTEXP_THRESHOLD)
            {
              evbd = CONSTMIN_EXP;
              devbd_dvb = 0.0;
            }
            else
            {
              evbd = exp(T0);
              devbd_dvb = evbd / Nvtmd;
            }

            T1 = (model_.bvd + vbd_jct) / Nvtmd;
            if (T1 > CONSTEXP_THRESHOLD)
            {
              T2 = CONSTMIN_EXP;
              T3 = 0.0;
            }
            else
            {
              T2 = exp(-T1);
              T3 = -T2 /Nvtmd;
            }
            gbd = DrainSatCurrent * (devbd_dvb - model_.xjbvd * T3)
                  + getDeviceOptions().gmin;
            cbd = DrainSatCurrent * (evbd + XExpBVD - 1.0
                                     - model_.xjbvd * T2) + getDeviceOptions().gmin * vbd_jct;
          }
          else
          {
            gbd = DslpFwd + getDeviceOptions().gmin;
            cbd = IVjdmFwd + DslpFwd * (vbd_jct
                                        - vjdmFwd) + getDeviceOptions().gmin * vbd_jct;
          }
          break;
      default: break;
    }
  }

  /* trap-assisted tunneling and recombination current for reverse bias  */
  Nvtmrssws = model_.vtm0 * model_.njtsswstemp;
  Nvtmrsswgs = model_.vtm0 * model_.njtsswgstemp;
  Nvtmrss = model_.vtm0 * model_.njtsstemp;
  Nvtmrsswd = model_.vtm0 * model_.njtsswdtemp;
  Nvtmrsswgd = model_.vtm0 * model_.njtsswgdtemp;
  Nvtmrsd = model_.vtm0 * model_.njtsdtemp;

  if ((model_.vtss - vbs_jct) < (model_.vtss * 1e-3))
  {
    T9 = 1.0e3;
    T0 = - vbs_jct / Nvtmrss * T9;
    DEXP(T0, T1, T10);
    dT1_dVb = T10 / Nvtmrss * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtss - vbs_jct);
    T0 = -vbs_jct / Nvtmrss * model_.vtss * T9;
    dT0_dVb = model_.vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T1, T10);
    dT1_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsd - vbd_jct) < (model_.vtsd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsd * T9;
    DEXP(T0, T2, T10);
    dT2_dVb = T10 / Nvtmrsd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsd * model_.vtsd * T9;
    dT0_dVb = model_.vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T2, T10);
    dT2_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtssws - vbs_jct) < (model_.vtssws * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbs_jct / Nvtmrssws * T9;
    DEXP(T0, T3, T10);
    dT3_dVb = T10 / Nvtmrssws * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtssws - vbs_jct);
    T0 = -vbs_jct / Nvtmrssws * model_.vtssws * T9;
    dT0_dVb = model_.vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T3, T10);
    dT3_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswd - vbd_jct) < (model_.vtsswd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsswd * T9;
    DEXP(T0, T4, T10);
    dT4_dVb = T10 / Nvtmrsswd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsswd * model_.vtsswd * T9;
    dT0_dVb = model_.vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T4, T10);
    dT4_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswgs - vbs_jct) < (model_.vtsswgs * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbs_jct / Nvtmrsswgs * T9;
    DEXP(T0, T5, T10);
    dT5_dVb = T10 / Nvtmrsswgs * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswgs - vbs_jct);
    T0 = -vbs_jct / Nvtmrsswgs * model_.vtsswgs * T9;
    dT0_dVb = model_.vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9) ;
    DEXP(T0, T5, T10);
    dT5_dVb = T10 * dT0_dVb;
  }

  if ((model_.vtsswgd - vbd_jct) < (model_.vtsswgd * 1e-3) )
  {
    T9 = 1.0e3;
    T0 = -vbd_jct / Nvtmrsswgd * T9;
    DEXP(T0, T6, T10);
    dT6_dVb = T10 / Nvtmrsswgd * T9;
  }
  else
  {
    T9 = 1.0 / (model_.vtsswgd - vbd_jct);
    T0 = -vbd_jct / Nvtmrsswgd * model_.vtsswgd * T9;
    dT0_dVb = model_.vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9) ;
    DEXP(T0, T6, T10);
    dT6_dVb = T10 * dT0_dVb;
  }

  gbs += SjctTempRevSatCur * dT1_dVb
        + SswTempRevSatCur * dT3_dVb
        + SswgTempRevSatCur * dT5_dVb;
  cbs -= SjctTempRevSatCur * (T1 - 1.0)
        + SswTempRevSatCur * (T3 - 1.0)
        + SswgTempRevSatCur * (T5 - 1.0);
  gbd += DjctTempRevSatCur * dT2_dVb
        + DswTempRevSatCur * dT4_dVb
        + DswgTempRevSatCur * dT6_dVb;
  cbd -= DjctTempRevSatCur * (T2 - 1.0)
        + DswTempRevSatCur * (T4 - 1.0)
        + DswgTempRevSatCur * (T6 - 1.0);

  // End of diode DC model

  if (vds >= 0.0)
  {
    mode = 1;
    Vds = vds;
    Vgs = vgs;
    Vbs = vbs;
    Vdb = vds - vbs;  // WDLiu: for GIDL
  }
  else
  {
    mode = -1;
    Vds = -vds;
    Vgs = vgd;
    Vbs = vbd;
    Vdb = -vbs;
  }

  // dunga
  if(model_.mtrlMod)
  {
    epsrox = 3.9;
    toxe = model_.eot;
    epssub = CONSTEPS0 * model_.epsrsub;
  }
  else
  {
    epsrox = model_.epsrox;
    toxe = model_.toxe;
    epssub = CONSTEPSSI;
  }

  /////////////////////////////////////////////////////////////////////////////
  // mosfet continuation.
  // This idea is based, loosely, on a paper by Jaijeet
  // Rosychowdhury.  If the artificial parameter flag has been enabled,
  // modify Vds and Vgs.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    if ( getSolverState().dcopFlag || getSolverState().tranopFlag )
    {
      Xyce::dout() << "HOMOTOPY INFO: gainscale   = " << getSolverState().gainScale_ << std::endl
                   << "HOMOTOPY INFO: before vds  = " << Vds << std::endl
                   << "HOMOTOPY INFO: before vgst = " << Vgs << std::endl;
    }
  }
  if (getSolverState().artParameterFlag_)
  {

    double alpha = getSolverState().gainScale_;
    double vgstConst = getDeviceOptions().vgstConst;

    Vds = devSupport.contVds (Vds,getSolverState().nltermScale_, getDeviceOptions().vdsScaleMin);
    Vgs = devSupport.contVgst(Vgs, alpha, vgstConst);
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag && isActive(Diag::TIME_STEP))
  {
    if ( getSolverState().dcopFlag || getSolverState().tranopFlag )
    {
      Xyce::dout() << "HOMOTOPY INFO: after vds   = " << Vds << std::endl;
      Xyce::dout() << "HOMOTOPY INFO: after vgst  = " << Vgs << std::endl;
    }
  }
  // end of mosfet continuation block.
  /////////////////////////////////////////////////////////////////////////////

  T0 = Vbs - vbsc - 0.001;
  T1 = sqrt(T0 * T0 - 0.004 * vbsc);
  if (T0 >= 0.0)
  {
    Vbseff = vbsc + 0.5 * (T0 + T1);
    dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
  }
  else
  {
    T2 = -0.002 / (T1 - T0);
    Vbseff = vbsc * (1.0 + T2);
    dVbseff_dVb = T2 * vbsc / T1;
  }

  // JX: Correction to forward body bias
  T9 = 0.95 * paramPtr->phi;
  T0 = T9 - Vbseff - 0.001;
  T1 = sqrt(T0 * T0 + 0.004 * T9);
  Vbseff = T9 - 0.5 * (T0 + T1);
  dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);

  Phis = paramPtr->phi - Vbseff;
  dPhis_dVb = -1.0;
  sqrtPhis = sqrt(Phis);
  dsqrtPhis_dVb = -0.5 / sqrtPhis;

  Xdep = paramPtr->Xdep0 * sqrtPhis / paramPtr->sqrtPhi;
  dXdep_dVb = (paramPtr->Xdep0 / paramPtr->sqrtPhi) * dsqrtPhis_dVb;

  Leff = paramPtr->leff;
  Vtm = model_.vtm;
  Vtm0 = model_.vtm0;

  // Vth Calculation
  T3 = sqrt(Xdep);
  V0 = paramPtr->vbi - paramPtr->phi;

  T0 = paramPtr->dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2;
  }
  else
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2 * T4 * T4;
  }
  lt1 = model_.factor1 * T3 * T1;
  dlt1_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = paramPtr->dvt2w;
  }
  else
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = paramPtr->dvt2w * T4 * T4;
  }
  ltw = model_.factor1 * T3 * T1;
  dltw_dVb = model_.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

  T0 = paramPtr->dvt1 * Leff / lt1;
  if (T0 < CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 - 1.0;
    T3 = T2 * T2;
    T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
    Theta0 = T1 / T4;
    dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
    dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + CONSTMIN_EXP)) / T4 / T4;
  }
  else
  {
    Theta0 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted
    dTheta0_dVb = 0.0;
  }
  thetavth = paramPtr->dvt0 * Theta0;
  Delt_vth = thetavth * V0;
  dDelt_vth_dVb = paramPtr->dvt0 * dTheta0_dVb * V0;

  T0 = paramPtr->dvt1w * paramPtr->weff * Leff / ltw;
  if (T0 < CONSTEXP_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 - 1.0;
    T3 = T2 * T2;
    T4 = T3 + 2.0 * T1 * CONSTMIN_EXP;
    T5 = T1 / T4;
    dT1_dVb = -T0 * T1 * dltw_dVb / ltw;
    dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + CONSTMIN_EXP)) / T4 / T4;
  }
  else
  {
    T5 = 1.0 / (CONSTMAX_EXP - 2.0); // 3.0 * CONSTMIN_EXP omitted
    dT5_dVb = 0.0;
  }

  T0 = paramPtr->dvt0w * T5;
  T2 = T0 * V0;
  dT2_dVb = paramPtr->dvt0w * dT5_dVb * V0;

  TempRatio =  temp / model_.tnom - 1.0;
  T0 = sqrt(1.0 + paramPtr->lpe0 / Leff);
  T1 = paramPtr->k1ox * (T0 - 1.0) * paramPtr->sqrtPhi
     + (paramPtr->kt1 + paramPtr->kt1l / Leff
     + paramPtr->kt2 * Vbseff) * TempRatio;
  Vth_NarrowW = toxe * paramPtr->phi
              / (paramPtr->weff + paramPtr->w0);

  T3 = eta0 + paramPtr->etab * Vbseff;
  if (T3 < 1.0e-4)
  {
    T9 = 1.0 / (3.0 - 2.0e4 * T3);
    T3 = (2.0e-4 - T3) * T9;
    T4 = T9 * T9;
  }
  else
  {
    T4 = 1.0;
  }
  dDIBL_Sft_dVd = T3 * paramPtr->theta0vb0;
  DIBL_Sft = dDIBL_Sft_dVd * Vds;

  Lpe_Vb = sqrt(1.0 + paramPtr->lpeb / Leff);

  Vth = model_.dtype * vth0 + (paramPtr->k1ox * sqrtPhis
      - paramPtr->k1 * paramPtr->sqrtPhi) * Lpe_Vb
      - k2ox * Vbseff - Delt_vth - T2 + (paramPtr->k3
      + paramPtr->k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

  dVth_dVb = Lpe_Vb * paramPtr->k1ox * dsqrtPhis_dVb - k2ox
           - dDelt_vth_dVb - dT2_dVb + paramPtr->k3b * Vth_NarrowW
           - paramPtr->etab * Vds * paramPtr->theta0vb0 * T4
           + paramPtr->kt2 * TempRatio;
  dVth_dVd = -dDIBL_Sft_dVd;


  // Calculate n
  tmp1 = epssub / Xdep;
  nstar = model_.vtm / CONSTQ * (model_.coxe + tmp1 + paramPtr->cit);
  tmp2 = paramPtr->nfactor * tmp1;
  tmp3 = paramPtr->cdsc + paramPtr->cdscb * Vbseff
       + paramPtr->cdscd * Vds;
  tmp4 = (tmp2 + tmp3 * Theta0 + paramPtr->cit) / model_.coxe;
  if (tmp4 >= -0.5)
  {
    n = 1.0 + tmp4;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
           + paramPtr->cdscb * Theta0) / model_.coxe;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.coxe;
  }
  else
  {
    T0 = 1.0 / (3.0 + 8.0 * tmp4);
    n = (1.0 + 3.0 * tmp4) * T0;
    T0 *= T0;
    dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
           + paramPtr->cdscb * Theta0) / model_.coxe * T0;
    dn_dVd = paramPtr->cdscd * Theta0 / model_.coxe * T0;
  }


  // Vth correction for Pocket implant
  if (paramPtr->dvtp0 > 0.0)
  {
    T0 = -paramPtr->dvtp1 * Vds;
    if (T0 < -CONSTEXP_THRESHOLD)
    {
      T2 = CONSTMIN_EXP;
      dT2_dVd = 0.0;
    }
    else
    {
      T2 = exp(T0);
      dT2_dVd = -paramPtr->dvtp1 * T2;
    }

    T3 = Leff + paramPtr->dvtp0 * (1.0 + T2);
    dT3_dVd = paramPtr->dvtp0 * dT2_dVd;
    if (model_.tempMod < 2)
    {
      T4 = Vtm * log(Leff / T3);
      dT4_dVd = -Vtm * dT3_dVd / T3;
    }
    else
    {
      T4 = model_.vtm0 * log(Leff / T3);
      dT4_dVd = -model_.vtm0 * dT3_dVd / T3;
    }
    dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
    dDITS_Sft_dVb = T4 * dn_dVb;

    Vth -= n * T4;
    dVth_dVd -= dDITS_Sft_dVd;
    dVth_dVb -= dDITS_Sft_dVb;
  }

  // v4.7 DITS_SFT2
  if ((paramPtr->dvtp4 == 0.0) || (paramPtr->dvtp2factor == 0.0))
  {
    T0 = 0.0;
    DITS_Sft2 = 0.0;
  }
  else
  {
    T1 = 2.0 * paramPtr->dvtp4 * Vds;
    DEXP(T1, T0, T10);
    DITS_Sft2 = paramPtr->dvtp2factor * (T0 - 1) / (T0 + 1);
    dDITS_Sft2_dVd = paramPtr->dvtp2factor * paramPtr->dvtp4 * 4.0 * T10 / ((T0 + 1) * (T0 + 1));
    Vth -= DITS_Sft2;
    dVth_dVd -= dDITS_Sft2_dVd;
  }

  von = Vth;

  // Poly Gate Si Depletion Effect
  T0 = vfb + paramPtr->phi;
  if(model_.mtrlMod == 0)
    T1 = CONSTEPSSI;
  else
    T1 = model_.epsrgate * CONSTEPS0;

  polyDepletion(T0, paramPtr->ngate, T1, model_.coxe,
                vgs, vgs_eff, dvgs_eff_dvg);

  polyDepletion(T0, paramPtr->ngate, T1, model_.coxe,
                vgd, vgd_eff, dvgd_eff_dvg);

  if(mode>0)
  {
    Vgs_eff = vgs_eff;
    dVgs_eff_dVg = dvgs_eff_dvg;
  }
  else
  {
    Vgs_eff = vgd_eff;
    dVgs_eff_dVg = dvgd_eff_dvg;
  }

  Vgst = Vgs_eff - Vth;

  // Calculate Vgsteff
  T0 = n * Vtm;
  T1 = paramPtr->mstar * Vgst;
  T2 = T1 / T0;
  if (T2 > CONSTEXP_THRESHOLD)
  {
    T10 = T1;
    dT10_dVg = paramPtr->mstar * dVgs_eff_dVg;
    dT10_dVd = -dVth_dVd * paramPtr->mstar;
    dT10_dVb = -dVth_dVb * paramPtr->mstar;
  }
  else if (T2 < -CONSTEXP_THRESHOLD)
  {
    T10 = Vtm * log(1.0 + CONSTMIN_EXP);
    dT10_dVg = 0.0;
    dT10_dVd = T10 * dn_dVd;
    dT10_dVb = T10 * dn_dVb;
    T10 *= n;
  }
  else
  {
    ExpVgst = exp(T2);
    T3 = Vtm * log(1.0 + ExpVgst);
    T10 = n * T3;
    dT10_dVg = paramPtr->mstar * ExpVgst / (1.0 + ExpVgst);
    dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
    dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
    dT10_dVg *= dVgs_eff_dVg;
  }

  T1 = paramPtr->voffcbn - (1.0 - paramPtr->mstar) * Vgst;
  T2 = T1 / T0;

  if (T2 < -CONSTEXP_THRESHOLD)
  {
    T3 = model_.coxe * CONSTMIN_EXP / paramPtr->cdep0;
    T9 = paramPtr->mstar + T3 * n;
    dT9_dVg = 0.0;
    dT9_dVd = dn_dVd * T3;
    dT9_dVb = dn_dVb * T3;
  }
  else if (T2 > CONSTEXP_THRESHOLD)
  {
    T3 = model_.coxe * CONSTMAX_EXP / paramPtr->cdep0;
    T9 = paramPtr->mstar + T3 * n;
    dT9_dVg = 0.0;
    dT9_dVd = dn_dVd * T3;
    dT9_dVb = dn_dVb * T3;
  }
  else
  {
    ExpVgst = exp(T2);
    T3 = model_.coxe / paramPtr->cdep0;
    T4 = T3 * ExpVgst;
    T5 = T1 * T4 / T0;
    T9 = paramPtr->mstar + n * T4;
    dT9_dVg = T3 * (paramPtr->mstar - 1.0) * ExpVgst / Vtm;
    dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
    dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
    dT9_dVg *= dVgs_eff_dVg;
  }

  Vgsteff = T10 / T9;
  Vgsteff_forNoise = Vgsteff;
  T11 = T9 * T9;
  dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
  dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
  dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

  // Calculate Effective Channel Geometry
  T9 = sqrtPhis - paramPtr->sqrtPhi;
  Weff = paramPtr->weff - 2.0 * (paramPtr->dwg * Vgsteff
      + paramPtr->dwb * T9);
  dWeff_dVg = -2.0 * paramPtr->dwg;
  dWeff_dVb = -2.0 * paramPtr->dwb * dsqrtPhis_dVb;

  if (Weff < 2.0e-8) // to avoid the discontinuity problem due to Weff
  {
    T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
    Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
    T0 *= T0 * 4.0e-16;
    dWeff_dVg *= T0;
    dWeff_dVb *= T0;
  }

  if (model_.rdsMod == 1)
  {
    Rds = dRds_dVg = dRds_dVb = 0.0;
  }
  else
  {
    T0 = 1.0 + paramPtr->prwg * Vgsteff;
    dT0_dVg = -paramPtr->prwg / T0 / T0;
    T1 = paramPtr->prwb * T9;
    dT1_dVb = paramPtr->prwb * dsqrtPhis_dVb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01); // 0.01 = 4.0 * 0.05 * 0.05
    dT3_dVg = 1.0 + T2 / (T3 - T2);
    dT3_dVb = dT3_dVg * dT1_dVb;
    dT3_dVg *= dT0_dVg;

    T4 = paramPtr->rds0 * 0.5;
    Rds = paramPtr->rdswmin + T3 * T4;
    dRds_dVg = T4 * dT3_dVg;
    dRds_dVb = T4 * dT3_dVb;

    if (Rds > 0.0)
    {
     grdsw = 1.0 / Rds * nf;
    }
    else
    {
     grdsw = 0.0;
    }
  }

  // Calculate Abulk
  T9 = 0.5 * paramPtr->k1ox * Lpe_Vb / sqrtPhis;
  T1 = T9 + k2ox - paramPtr->k3b * Vth_NarrowW;
  dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

  T9 = sqrt(paramPtr->xj * Xdep);
  tmp1 = Leff + 2.0 * T9;
  T5 = Leff / tmp1;
  tmp2 = paramPtr->a0 * T5;
  tmp3 = paramPtr->weff + paramPtr->b1;
  tmp4 = paramPtr->b0 / tmp3;
  T2 = tmp2 + tmp4;
  dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
  T6 = T5 * T5;
  T7 = T5 * T6;

  Abulk0 = 1.0 + T1 * T2;
  dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

  T8 = paramPtr->ags * paramPtr->a0 * T7;
  dAbulk_dVg = -T1 * T8;
  Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
  dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb + 3.0 * T1 * dT2_dVb);

  if (Abulk0 < 0.1) // added to avoid the problems caused by Abulk0
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk0);
    Abulk0 = (0.2 - Abulk0) * T9;
    dAbulk0_dVb *= T9 * T9;
  }

  if (Abulk < 0.1)
  {
    T9 = 1.0 / (3.0 - 20.0 * Abulk);
    Abulk = (0.2 - Abulk) * T9;
    T10 = T9 * T9;
    dAbulk_dVb *= T10;
    dAbulk_dVg *= T10;
  }

  Abulk_forNoise = Abulk;  // ERK. this is a bit screwy.  But in
                           // spice3/ngspice, Abulk is saved as an
                           // instance variable and used later in
                           // noise calculations. But, within the load
                           // function the "local" copy of Abulk gets
                           // modified further afterwards.  So it is
                           // then "wrong" for noise at that point.
                           // As we've made Abulk a class variable,
                           // there needs to be an extra copy for
                           // noise, saved at the right time.

  T2 = paramPtr->keta * Vbseff;
  if (T2 >= -0.9)
  {
    T0 = 1.0 / (1.0 + T2);
    dT0_dVb = -paramPtr->keta * T0 * T0;
  }
  else
  {
    T1 = 1.0 / (0.8 + T2);
    T0 = (17.0 + 20.0 * T2) * T1;
    dT0_dVb = -paramPtr->keta * T1 * T1;
  }

  dAbulk_dVg *= T0;
  dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
  dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
  Abulk *= T0;
  Abulk0 *= T0;

  // Mobility calculation
  if (model_.mtrlMod && (model_.mtrlCompatMod == 0))
    T14 = 2.0 * model_.dtype *(model_.phig - model_.easub - 0.5*model_.Eg0 + 0.45);
  else
    T14 = 0.0;

  if (model_.mobMod == 0)
  { T0 = Vgsteff + Vth + Vth - T14;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T3 * (T2 + paramPtr->ub * T3) + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = (T2 + 2.0 * paramPtr->ub * T3) / toxe;
    T13 = 2.0 * (dDenomi_dVg + T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + paramPtr->uc * T3;
    dDenomi_dVg+= T7;
  }
  else if (model_.mobMod == 1)
  {   T0 = Vgsteff + Vth + Vth - T14;
    T2 = 1.0 + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T4 = T3 * (paramPtr->ua + paramPtr->ub * T3);
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T4 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = (paramPtr->ua + 2.0 * paramPtr->ub * T3) * T2 / toxe;
    T13 = 2.0 * (dDenomi_dVg + T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + paramPtr->uc * T4;
    dDenomi_dVg+= T7;
  }
  else if (model_.mobMod == 2)
  {   T0 = (Vgsteff + vtfbphi1) / toxe;
    T1 = exp(paramPtr->eu * log(T0));
    dT1_dVg = T1 * paramPtr->eu / T0 / toxe;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T12 = sqrt(Vth * Vth + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * Vth;
    T6 = T8 * Vth;
    T5 = T1 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    T11 = T7 * Vth/T12;
    dDenomi_dVg = T2 * dT1_dVg + T7;
    T13 = 2.0 * (T11 + T8);
    dDenomi_dVd = T13 * dVth_dVd;
    dDenomi_dVb = T13 * dVth_dVb + T1 * paramPtr->uc;
  }
  else if (model_.mobMod == 4) // Synopsys 08/30/2013 add
  {
    T0 = Vgsteff + vtfbphi1 - T14;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T12 = sqrt(vtfbphi1*vtfbphi1 + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * vtfbphi1;
    T6 = T8 * vtfbphi1;
    T5 = T3 * (T2 + paramPtr->ub * T3) + T6;
    T7 = - 2.0 * T6 * T9;
    dDenomi_dVg = (T2 + 2.0 * paramPtr->ub * T3) / toxe;
    dDenomi_dVd = 0.0;
    dDenomi_dVb = paramPtr->uc * T3;
    dDenomi_dVg+= T7;
  }
  else if (model_.mobMod == 5) // Synopsys 08/30/2013 add
  {
    T0 = Vgsteff + vtfbphi1 - T14;
    T2 = 1.0 + paramPtr->uc * Vbseff;
    T3 = T0 / toxe;
    T4 = T3 * (paramPtr->ua + paramPtr->ub * T3);
    T12 = sqrt(vtfbphi1 * vtfbphi1 + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * vtfbphi1;
    T6 = T8 * vtfbphi1;
    T5 = T4 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    dDenomi_dVg = (paramPtr->ua + 2.0 * paramPtr->ub * T3) * T2
      / toxe;
    dDenomi_dVd = 0.0;
    dDenomi_dVb = paramPtr->uc * T4;
    dDenomi_dVg+= T7;
  }
  else if (model_.mobMod == 6) // Synopsys 08/30/2013 modify
  {
    T0 = (Vgsteff + vtfbphi1) / toxe;
    T1 = exp(paramPtr->eu * log(T0));
    dT1_dVg = T1 * paramPtr->eu / T0 / toxe;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;

    T12 = sqrt(vtfbphi1 * vtfbphi1 + 0.0001);
    T9 = 1.0/(Vgsteff + 2*T12);
    T10 = T9*toxe;
    T8 = paramPtr->ud * T10 * T10 * vtfbphi1;
    T6 = T8 * vtfbphi1;
    T5 = T1 * T2 + T6;
    T7 = - 2.0 * T6 * T9;
    dDenomi_dVg = T2 * dT1_dVg + T7;
    dDenomi_dVd = 0;
    dDenomi_dVb = T1 * paramPtr->uc;
  }
  else
  {
    // high K mobility

    // universal mobility
    T0 = (Vgsteff + vtfbphi1) * 1.0e-8 / toxe / 6.0;
    T1 = exp(paramPtr->eu * log(T0));
    dT1_dVg = T1 * paramPtr->eu * 1.0e-8 / T0 / toxe / 6.0;
    T2 = paramPtr->ua + paramPtr->uc * Vbseff;

    // Coulombic
    VgsteffVth = paramPtr->VgsteffVth;

    T10 = exp(paramPtr->ucs * log(0.5 + 0.5 * Vgsteff / VgsteffVth));
    T11 = paramPtr->ud / T10;
    dT11_dVg = -0.5 * paramPtr->ucs * T11 / (0.5 + 0.5 * Vgsteff / VgsteffVth) / VgsteffVth;

    dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
    dDenomi_dVd = 0.0;
    dDenomi_dVb = T1 * paramPtr->uc;

    T5 = T1 * T2 + T11;
  }

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
    dDenomi_dVg *= T9;
    dDenomi_dVd *= T9;
    dDenomi_dVb *= T9;
  }

  ueff = ueff = u0temp / Denomi;
  T9 = -ueff / Denomi;
  dueff_dVg = T9 * dDenomi_dVg;
  dueff_dVd = T9 * dDenomi_dVd;
  dueff_dVb = T9 * dDenomi_dVb;

  // Saturation Drain Voltage  Vdsat
  WVCox = Weff * vsattemp * model_.coxe;
  WVCoxRds = WVCox * Rds;

  Esat = 2.0 * vsattemp / ueff;
  EsatL = EsatL = Esat * Leff;
  T0 = -EsatL /ueff;
  dEsatL_dVg = T0 * dueff_dVg;
  dEsatL_dVd = T0 * dueff_dVd;
  dEsatL_dVb = T0 * dueff_dVb;

  // Sqrt()
  a1 = paramPtr->a1;
  if (a1 == 0.0)
  {
    Lambda = paramPtr->a2;
    dLambda_dVg = 0.0;
  }
  else if (a1 > 0.0)
  {
    T0 = 1.0 - paramPtr->a2;
    T1 = T0 - paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * T0);
    Lambda = paramPtr->a2 + T0 - 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }
  else
  {
    T1 = paramPtr->a2 + paramPtr->a1 * Vgsteff - 0.0001;
    T2 = sqrt(T1 * T1 + 0.0004 * paramPtr->a2);
    Lambda = 0.5 * (T1 + T2);
    dLambda_dVg = 0.5 * paramPtr->a1 * (1.0 + T1 / T2);
  }

  Vgst2Vtm = Vgsteff + 2.0 * Vtm;
  if (Rds > 0)
  {
    tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
    tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
  }
  else
  {
    tmp2 = dWeff_dVg / Weff;
    tmp3 = dWeff_dVb / Weff;
  }
  if ((Rds == 0.0) && (Lambda == 1.0))
  {
    T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
    tmp1 = 0.0;
    T1 = T0 * T0;
    T2 = Vgst2Vtm * T0;
    T3 = EsatL * Vgst2Vtm;
    Vdsat = T3 * T0;

    dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
    dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
    dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

    dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
    dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
    dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
  }
  else
  {
    tmp1 = dLambda_dVg / (Lambda * Lambda);
    T9 = Abulk * WVCoxRds;
    T8 = Abulk * T9;
    T7 = Vgst2Vtm * T9;
    T6 = Vgst2Vtm * WVCoxRds;
    T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
    dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
     + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

    dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
     + (1.0 / Lambda - 1.0) * dAbulk_dVb);
    dT0_dVd = 0.0;
    T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

    dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
     + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
     + T7 * tmp2 + T6 * dAbulk_dVg);
    dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
           + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
    dT1_dVd = Abulk * dEsatL_dVd;

    T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
    dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
       + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
    dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
    dT2_dVd = Vgst2Vtm * dEsatL_dVd;

    T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
    Vdsat = (T1 - T3) / T0;

    dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
     / T3;
    dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
     / T3;
    dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
     / T3;

    dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
      - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
    dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
      - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
    dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
  }
  vdsat = Vdsat;

  // Calculate Vdseff
  T1 = Vdsat - Vds - paramPtr->delta;
  dT1_dVg = dVdsat_dVg;
  dT1_dVd = dVdsat_dVd - 1.0;
  dT1_dVb = dVdsat_dVb;

  T2 = sqrt(T1 * T1 + 4.0 * paramPtr->delta * Vdsat);
  T0 = T1 / T2;
  T9 = 2.0 * paramPtr->delta;
  T3 = T9 / T2;
  dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
  dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
  dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

  if (T1 >= 0.0)
  {
    Vdseff = Vdsat - 0.5 * (T1 + T2);
    dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
    dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
    dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
  }
  else
  {
    T4 = T9 / (T2 - T1);
    T5 = 1.0 - T4;
    T6 = Vdsat * T4 / (T2 - T1);
    Vdseff = Vdsat * T5;
    dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
    dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
    dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
  }

  if (Vds == 0.0)
  {
    Vdseff = 0.0;
    dVdseff_dVg = 0.0;
    dVdseff_dVb = 0.0;
  }

  if (Vdseff > Vds)
  {
    Vdseff = Vds;
  }

  diffVds = Vds - Vdseff;
  Vdseff_forNoise = Vdseff;

  // Velocity Overshoot
  if((model_.lambdaGiven) && (model_.lambda > 0.0) )
  {
    T1 =  Leff * ueff;
    T2 = paramPtr->lambda / T1;
    T3 = -T2 / T1 * Leff;
    dT2_dVd = T3 * dueff_dVd;
    dT2_dVg = T3 * dueff_dVg;
    dT2_dVb = T3 * dueff_dVb;
    T5 = 1.0 / (Esat * paramPtr->litl);
    T4 = -T5 / EsatL;
    dT5_dVg = dEsatL_dVg * T4;
    dT5_dVd = dEsatL_dVd * T4;
    dT5_dVb = dEsatL_dVb * T4;
    T6 = 1.0 + diffVds  * T5;
    dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
    dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
    dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
    T7 = 2.0 / (T6 * T6 + 1.0);
    T8 = 1.0 - T7;
    T9 = T6 * T7 * T7;
    dT8_dVg = T9 * dT6_dVg;
    dT8_dVd = T9 * dT6_dVd;
    dT8_dVb = T9 * dT6_dVb;
    T10 = 1.0 + T2 * T8;
    dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
    dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
    dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
    if(T10 == 1.0)
    {
      dT10_dVg = dT10_dVd = dT10_dVb = 0.0;
    }

    dEsatL_dVg *= T10;
    dEsatL_dVg += EsatL * dT10_dVg;
    dEsatL_dVd *= T10;
    dEsatL_dVd += EsatL * dT10_dVd;
    dEsatL_dVb *= T10;
    dEsatL_dVb += EsatL * dT10_dVb;
    EsatL *= T10;
    Esat = EsatL / Leff;
  }

  // Calculate Vasat
  tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff;
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

  T7 = 2.0 * WVCoxRds * tmp4;
  dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
          - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
          + Vdsat * dAbulk_dVg);

  dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
          - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
  dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

  T9 = WVCoxRds * Abulk;
  T1 = 2.0 / Lambda - 1.0 + T9;
  dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
  dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

  Vasat = T0 / T1;
  dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
  dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
  dVasat_dVd = dT0_dVd / T1;

  // Calculate Idl first
  tmp1 = vtfbphi2;
  tmp2 = 2.0e8 * toxp;
  dT0_dVg = 1.0 / tmp2;
  T0 = (Vgsteff + tmp1) * dT0_dVg;

  tmp3 = exp(model_.bdos * 0.7 * log(T0));
  T1 = 1.0 + tmp3;
  T2 = model_.bdos * 0.7 * tmp3 / T0;
  Tcen = model_.ados * 1.9e-9 / T1;
  dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

  Coxeff = epssub * coxp / (epssub + coxp * Tcen);
  dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

  CoxeffWovL = Coxeff * Weff / Leff;
  beta = ueff * CoxeffWovL;
  T3 = ueff / Leff;
  dbeta_dVg = CoxeffWovL * dueff_dVg + T3
          * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
  dbeta_dVd = CoxeffWovL * dueff_dVd;
  dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

  AbovVgst2Vtm = Abulk / Vgst2Vtm;
  T0 = 1.0 - 0.5 * Vdseff * AbovVgst2Vtm;
  dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
          - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
  dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
  dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff) / Vgst2Vtm;

  fgche1 = Vgsteff * T0;
  dfgche1_dVg = Vgsteff * dT0_dVg + T0;
  dfgche1_dVd = Vgsteff * dT0_dVd;
  dfgche1_dVb = Vgsteff * dT0_dVb;

  T9 = Vdseff / EsatL;
  fgche2 = 1.0 + T9;
  dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
  dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
  dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

  gche = beta * fgche1 / fgche2;
  dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
          - gche * dfgche2_dVg) / fgche2;
  dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
          - gche * dfgche2_dVd) / fgche2;
  dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
          - gche * dfgche2_dVb) / fgche2;

  T0 = 1.0 + gche * Rds;
  Idl = gche / T0;
  T1 = (1.0 - Idl * Rds) / T0;
  T2 = Idl * Idl;
  dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
  dIdl_dVd = T1 * dgche_dVd;
  dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

  // Calculate degradation factor due to pocket implant
  if (paramPtr->fprout <= 0.0)
  {
    FP = 1.0;
    dFP_dVg = 0.0;
  }
  else
  {
    T9 = paramPtr->fprout * sqrt(Leff) / Vgst2Vtm;
    FP = 1.0 / (1.0 + T9);
    dFP_dVg = FP * FP * T9 / Vgst2Vtm;
  }

  // Calculate VACLM
  T8 = paramPtr->pvag / EsatL;
  T9 = T8 * Vgsteff;
  if (T9 > -0.9)
  {
    PvagTerm = 1.0 + T9;
    dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
    dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
    dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
  }
  else
  {
    T4 = 1.0 / (17.0 + 20.0 * T9);
    PvagTerm = (0.8 + T9) * T4;
    T4 *= T4;
    dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
    T9 *= T4 / EsatL;
    dPvagTerm_dVb = -T9 * dEsatL_dVb;
    dPvagTerm_dVd = -T9 * dEsatL_dVd;
  }

  if ((paramPtr->pclm > CONSTMIN_EXP) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 + Rds * Idl;
    dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
    dT0_dVd = Rds * dIdl_dVd;
    dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

    T2 = Vdsat / Esat;
    T1 = Leff + T2;
    dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
    dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
    dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

    Cclm = FP * PvagTerm * T0 * T1 / (paramPtr->pclm * paramPtr->litl);
    dCclm_dVg = Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm
              + dT0_dVg / T0 + dT1_dVg / T1);
    dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0
              + dT1_dVb / T1);
    dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0
              + dT1_dVd / T1);
    VACLM = Cclm * diffVds;

    dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
    dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
    dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
  }
  else
  {
    VACLM = Cclm = CONSTMAX_EXP;
    dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
    dCclm_dVd = dCclm_dVg = dCclm_dVb = 0.0;
  }

  // Calculate VADIBL
  if (paramPtr->thetaRout > CONSTMIN_EXP)
  {
    T8 = Abulk * Vdsat;
    T0 = Vgst2Vtm * T8;
    dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
      + Vgst2Vtm * Vdsat * dAbulk_dVg;
    dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
    dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

    T1 = Vgst2Vtm + T8;
    dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
    dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
    dT1_dVd = Abulk * dVdsat_dVd;

    T9 = T1 * T1;
    T2 = paramPtr->thetaRout;
    VADIBL = (Vgst2Vtm - T0 / T1) / T2;
    dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
    dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
    dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

    T7 = paramPtr->pdiblb * Vbseff;
    if (T7 >= -0.9)
    {
      T3 = 1.0 / (1.0 + T7);
      VADIBL *= T3;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = (dVADIBL_dVb - VADIBL * paramPtr->pdiblb) * T3;
      dVADIBL_dVd *= T3;
    }
    else
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      dVADIBL_dVg *= T3;
      dVADIBL_dVb = dVADIBL_dVb * T3
        - VADIBL * paramPtr->pdiblb * T4 * T4;
      dVADIBL_dVd *= T3;
      VADIBL *= T3;
    }

    dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
    dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
    dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
    VADIBL *= PvagTerm;
  }
  else
  {
    VADIBL = CONSTMAX_EXP;
    dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
  }

  // Calculate Va
  Va = Vasat + VACLM;
  dVa_dVg = dVasat_dVg + dVACLM_dVg;
  dVa_dVb = dVasat_dVb + dVACLM_dVb;
  dVa_dVd = dVasat_dVd + dVACLM_dVd;

  // Calculate VADITS
  T0 = paramPtr->pditsd * Vds;
  if (T0 > CONSTEXP_THRESHOLD)
  {
    T1 = CONSTMAX_EXP;
    dT1_dVd = 0;
  }
  else
  {
    T1 = exp(T0);
    dT1_dVd = T1 * paramPtr->pditsd;
  }

  if (paramPtr->pdits > CONSTMIN_EXP)
  {
    T2 = 1.0 + model_.pditsl * Leff;
    VADITS = (1.0 + T2 * T1) / paramPtr->pdits;
    dVADITS_dVg = VADITS * dFP_dVg;
    dVADITS_dVd = FP * T2 * dT1_dVd / paramPtr->pdits;
    VADITS *= FP;
  }
  else
  {
    VADITS = CONSTMAX_EXP;
    dVADITS_dVg = dVADITS_dVd = 0;
  }

  // Calculate VASCBE
  if ((paramPtr->pscbe2 > 0.0) && (paramPtr->pscbe1 >= 0.0))
  {
    if (diffVds > paramPtr->pscbe1 * paramPtr->litl / CONSTEXP_THRESHOLD)
    {
      T0 =  paramPtr->pscbe1 * paramPtr->litl / diffVds;
      VASCBE = Leff * exp(T0) / paramPtr->pscbe2;
      T1 = T0 * VASCBE / diffVds;
      dVASCBE_dVg = T1 * dVdseff_dVg;
      dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
      dVASCBE_dVb = T1 * dVdseff_dVb;
    }
    else
    {
      VASCBE = CONSTMAX_EXP * Leff/paramPtr->pscbe2;
      dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
    }
  }
  else
  {
    VASCBE = CONSTMAX_EXP;
    dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
  }

  // Add DIBL to Ids
  T9 = diffVds / VADIBL;
  T0 = 1.0 + T9;
  Idsa = Idl * T0;
  dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
  dIdsa_dVd = T0 * dIdl_dVd + Idl
          * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
  dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

  // Add DITS to Ids
  T9 = diffVds / VADITS;
  T0 = 1.0 + T9;
  dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
  dIdsa_dVd = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
  dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
  Idsa *= T0;

  // Add CLM to Ids
  T0 = log(Va / Vasat);
  dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
  dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
  dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
  T1 = T0 / Cclm;
  T9 = 1.0 + T1;
  dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
  dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
  dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

  dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
  dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
  dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
  Idsa *= T9;

  // Substrate current begins
  tmp = paramPtr->alpha0 + paramPtr->alpha1 * Leff;
  if ((tmp <= 0.0) || (paramPtr->beta0 <= 0.0))
  {
    Isub = Gbd = Gbb = Gbg = 0.0;
  }
  else
  {
    T2 = tmp / Leff;
    if (diffVds > paramPtr->beta0 / CONSTEXP_THRESHOLD)
    {
      T0 = -paramPtr->beta0 / diffVds;
      T1 = T2 * diffVds * exp(T0);
      T3 = T1 / diffVds * (T0 - 1.0);
      dT1_dVg = T3 * dVdseff_dVg;
      dT1_dVd = T3 * (dVdseff_dVd - 1.0);
      dT1_dVb = T3 * dVdseff_dVb;
    }
    else
    {
      T3 = T2 * CONSTMIN_EXP;
      T1 = T3 * diffVds;
      dT1_dVg = -T3 * dVdseff_dVg;
      dT1_dVd = T3 * (1.0 - dVdseff_dVd);
      dT1_dVb = -T3 * dVdseff_dVb;
    }
    T4 = Idsa * Vdseff;
    Isub = T1 * T4;
    Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg)
        + T4 * dT1_dVg;
    Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd)
        + T4 * dT1_dVd;
    Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb)
        + T4 * dT1_dVb;

    Gbd += Gbg * dVgsteff_dVd;
    Gbb += Gbg * dVgsteff_dVb;
    Gbg *= dVgsteff_dVg;
    Gbb *= dVbseff_dVb;
  }
  csub = Isub;
  gbbs = Gbb;
  gbgs = Gbg;
  gbds = Gbd;

  // Add SCBE to Ids
  T9 = diffVds / VASCBE;
  T0 = 1.0 + T9;
  Ids = Idsa * T0;

  Gm = T0 * dIdsa_dVg - Idsa
   * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
  Gds = T0 * dIdsa_dVd + Idsa
    * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
  Gmb = T0 * dIdsa_dVb - Idsa
    * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;


  tmp1 = Gds + Gm * dVgsteff_dVd;
  tmp2 = Gmb + Gm * dVgsteff_dVb;
  tmp3 = Gm;

  Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
  Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd)
    + Vdseff * tmp1;
  Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb)
    + Vdseff * tmp2) * dVbseff_dVb;

  cdrain = Ids * Vdseff;

  // Source End Velocity Limit
  if((model_.vtlGiven) && (model_.vtl > 0.0) )
  {
    T12 = 1.0 / Leff / CoxeffWovL;
    T11 = T12 / Vgsteff;
    T10 = -T11 / Vgsteff;
    vs = cdrain * T11; // vs
    dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
    dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
    dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
    T0 = 2 * CONSTMM;
    T1 = vs / (paramPtr->vtl * paramPtr->tfactor);

    if(T1 > 0.0)
    {
      T2 = 1.0 + exp(T0 * log(T1));
      T3 = (T2 - 1.0) * T0 / vs;
      Fsevl = 1.0 / exp(log(T2)/ T0);
      dT2_dVg = T3 * dvs_dVg;
      dT2_dVd = T3 * dvs_dVd;
      dT2_dVb = T3 * dvs_dVb;
      T4 = -1.0 / T0 * Fsevl / T2;
      dFsevl_dVg = T4 * dT2_dVg;
      dFsevl_dVd = T4 * dT2_dVd;
      dFsevl_dVb = T4 * dT2_dVb;
    }
    else
    {
      Fsevl = 1.0;
      dFsevl_dVg = 0.0;
      dFsevl_dVd = 0.0;
      dFsevl_dVb = 0.0;
    }
    Gm *=Fsevl;
    Gm += cdrain * dFsevl_dVg;
    Gmb *=Fsevl;
    Gmb += cdrain * dFsevl_dVb;
    Gds *=Fsevl;
    Gds += cdrain * dFsevl_dVd;

    cdrain *= Fsevl;
  }

  gds = Gds;
  gm = Gm;
  gmbs = Gmb;
  IdovVds = Ids;
  if( IdovVds <= model_.idovvdsc) IdovVds = model_.idovvdsc;

  // Calculate Rg
  if ((rgateMod > 1) || (trnqsMod != 0) || (acnqsMod != 0))
  {
    T9 = paramPtr->xrcrg2 * model_.vtm;
    T0 = T9 * beta;
    dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
    dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
    dT0_dVg = dbeta_dVg * T9;

    gcrg = paramPtr->xrcrg1 * ( T0 + Ids);
    gcrgd = paramPtr->xrcrg1 * (dT0_dVd + tmp1);
    gcrgb = paramPtr->xrcrg1 * (dT0_dVb + tmp2)
                     * dVbseff_dVb;
    gcrgg = paramPtr->xrcrg1 * (dT0_dVg + tmp3)
                     * dVgsteff_dVg;

    if (nf != 1.0)
    {
      gcrg *= nf;
      gcrgg *= nf;
      gcrgd *= nf;
      gcrgb *= nf;
    }

    if (rgateMod == 2)
    {
      T10 = grgeltd * grgeltd;
      T11 = grgeltd + gcrg;
      gcrg = grgeltd * gcrg / T11;
      T12 = T10 / T11 / T11;
      gcrgg *= T12;
      gcrgd *= T12;
      gcrgb *= T12;
    }
    gcrgs = -(gcrgg + gcrgd + gcrgb);
  }

  // Calculate bias-dependent external S/D resistance
  if (model_.rdsMod)
  {   // Rs(V)
    T0 = vgs - paramPtr->vfbsd;
    T1 = sqrt(T0 * T0 + 1.0e-4);
    vgs_eff = 0.5 * (T0 + T1);
    dvgs_eff_dvg = vgs_eff / T1;

    T0 = 1.0 + paramPtr->prwg * vgs_eff;
    dT0_dvg = -paramPtr->prwg / T0 / T0 * dvgs_eff_dvg;
    T1 = -paramPtr->prwb * vbs;
    dT1_dvb = -paramPtr->prwb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01);
    dT3_dvg = T3 / (T3 - T2);
    dT3_dvb = dT3_dvg * dT1_dvb;
    dT3_dvg *= dT0_dvg;

    T4 = paramPtr->rs0 * 0.5;
    Rs = paramPtr->rswmin + T3 * T4;
    dRs_dvg = T4 * dT3_dvg;
    dRs_dvb = T4 * dT3_dvb;

    T0 = 1.0 + sourceConductance * Rs;
    gstot = sourceConductance / T0;
    T0 = -gstot * gstot;
    dgstot_dvd = 0.0; // place holder
    dgstot_dvg = T0 * dRs_dvg;
    dgstot_dvb = T0 * dRs_dvb;
    dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

    // Rd(V)
    T0 = vgd - paramPtr->vfbsd;
    T1 = sqrt(T0 * T0 + 1.0e-4);
    vgd_eff = 0.5 * (T0 + T1);
    dvgd_eff_dvg = vgd_eff / T1;

    T0 = 1.0 + paramPtr->prwg * vgd_eff;
    dT0_dvg = -paramPtr->prwg / T0 / T0 * dvgd_eff_dvg;
    T1 = -paramPtr->prwb * vbd;
    dT1_dvb = -paramPtr->prwb;

    T2 = 1.0 / T0 + T1;
    T3 = T2 + sqrt(T2 * T2 + 0.01);
    dT3_dvg = T3 / (T3 - T2);
    dT3_dvb = dT3_dvg * dT1_dvb;
    dT3_dvg *= dT0_dvg;

    T4 = paramPtr->rd0 * 0.5;
    Rd = paramPtr->rdwmin + T3 * T4;
    dRd_dvg = T4 * dT3_dvg;
    dRd_dvb = T4 * dT3_dvb;

    T0 = 1.0 + drainConductance * Rd;
    gdtot = drainConductance / T0;
    T0 = -gdtot * gdtot;
    dgdtot_dvs = 0.0;
    dgdtot_dvg = T0 * dRd_dvg;
    dgdtot_dvb = T0 * dRd_dvb;
    dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

    gstotd = vses * dgstot_dvd;
    gstotg = vses * dgstot_dvg;
    gstots = vses * dgstot_dvs;
    gstotb = vses * dgstot_dvb;

    T2 = vdes - vds;
    gdtotd = T2 * dgdtot_dvd;
    gdtotg = T2 * dgdtot_dvg;
    gdtots = T2 * dgdtot_dvs;
    gdtotb = T2 * dgdtot_dvb;
  }
  else // WDLiu: for bypass
  {
    gstot = gstotd = gstotg = 0.0;
    gstots = gstotb = 0.0;
    gdtot = gdtotd = gdtotg = 0.0;
    gdtots = gdtotb = 0.0;
  }

  // GIDL/GISL Models
  if (model_.mtrlMod == 0)
  {
    T0 = 3.0 * toxe;
  }
  else
  {
    T0 = model_.epsrsub * toxe / epsrox;
  }

  // Calculate GIDL current
  if (model_.gidlMod == 0)
  {
    if(model_.mtrlMod ==0)
      T1 = (vds - vgs_eff - paramPtr->egidl ) / T0;
    else
      T1 = (vds - vgs_eff - paramPtr->egidl + paramPtr->vfbsd) / T0;

    if ((paramPtr->agidl <= 0.0) || (paramPtr->bgidl <= 0.0)
                  || (T1 <= 0.0) || (paramPtr->cgidl <= 0.0) || (vbd > 0.0))
    {
      Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
    }
    else
    {
      dT1_dVd = 1.0 / T0;
      dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
      T2 = paramPtr->bgidl / T1;
      if (T2 < 100.0)
      {
        Igidl = paramPtr->agidl * paramPtr->weffCJ * T1 * exp(-T2);
        T3 = Igidl * (1.0 + T2) / T1;
        Ggidld = T3 * dT1_dVd;
        Ggidlg = T3 * dT1_dVg;
      }
      else
      {
        Igidl = paramPtr->agidl * paramPtr->weffCJ * 3.720075976e-44;
        Ggidld = Igidl * dT1_dVd;
        Ggidlg = Igidl * dT1_dVg;
        Igidl *= T1;
      }

      T4 = vbd * vbd;
      T5 = -vbd * T4;
      T6 = paramPtr->cgidl + T5;
      T7 = T5 / T6;
      T8 = 3.0 * paramPtr->cgidl * T4 / T6 / T6;
      Ggidld = Ggidld * T7 + Igidl * T8;
      Ggidlg = Ggidlg * T7;
      Ggidlb = -Igidl * T8;
      Igidl *= T7;
    }
    ggidld = Ggidld;
    ggidlg = Ggidlg;
    ggidlb = Ggidlb;

    // Calculate GISL current

    if (model_.mtrlMod == 0)
    {
      T1 = (-vds - vgd_eff - paramPtr->egisl ) / T0;
    }
    else
    {
      T1 = (-vds - vgd_eff - paramPtr->egisl + paramPtr->vfbsd ) / T0;
    }

    if ((paramPtr->agisl <= 0.0) || (paramPtr->bgisl <= 0.0)
      || (T1 <= 0.0) || (paramPtr->cgisl <= 0.0) || (vbs > 0.0))
    {
      Igisl = Ggisls = Ggislg = Ggislb = 0.0;
    }
    else
    {
      dT1_dVd = 1.0 / T0;
      dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
      T2 = paramPtr->bgisl / T1;
      if (T2 < 100.0)
      {
        Igisl = paramPtr->agisl * paramPtr->weffCJ * T1 * exp(-T2);
        T3 = Igisl * (1.0 + T2) / T1;
        Ggisls = T3 * dT1_dVd;
        Ggislg = T3 * dT1_dVg;
      }
      else
      {
        Igisl = paramPtr->agisl * paramPtr->weffCJ * 3.720075976e-44;
        Ggisls = Igisl * dT1_dVd;
        Ggislg = Igisl * dT1_dVg;
        Igisl *= T1;
      }

      T4 = vbs * vbs;
      T5 = -vbs * T4;
      T6 = paramPtr->cgisl + T5;
      T7 = T5 / T6;
      T8 = 3.0 * paramPtr->cgisl * T4 / T6 / T6;
      Ggisls = Ggisls * T7 + Igisl * T8;
      Ggislg = Ggislg * T7;
      Ggislb = -Igisl * T8;
      Igisl *= T7;
    }
    ggisls = Ggisls;
    ggislg = Ggislg;
    ggislb = Ggislb;
  }
  else
  {
    // v4.7 New Gidl/GISL model

    // GISL
    if (model_.mtrlMod == 0)
      T1 = (-vds - paramPtr->rgisl * vgd_eff - paramPtr->egisl) / T0;
    else
      T1 = (-vds - paramPtr->rgisl * vgd_eff - paramPtr->egisl + paramPtr->vfbsd) / T0;

    if ( (paramPtr->agisl <= 0.0) ||
         (paramPtr->bgisl <= 0.0) || (T1 <= 0.0) ||
         (paramPtr->cgisl < 0.0) )
      Igisl = Ggisls = Ggislg = Ggislb = 0.0;
    else
    {
      dT1_dVd = 1 / T0;
      dT1_dVg = - paramPtr->rgisl * dT1_dVd * dvgd_eff_dvg;
      T2 = paramPtr->bgisl / T1;
      if (T2 < CONSTEXPL_THRESHOLD)
      {
        Igisl = paramPtr->weffCJ * paramPtr->agisl * T1 * exp(-T2);
        T3 = Igisl / T1 * (T2 + 1);
        Ggisls = T3 * dT1_dVd;
        Ggislg = T3 * dT1_dVg;
      }
      else
      {
        T3 = paramPtr->weffCJ * paramPtr->agisl * CONSTMIN_EXPL;
        Igisl = T3 * T1;
        Ggisls = T3 * dT1_dVd;
        Ggislg = T3 * dT1_dVg;
      }
      T4 = vbs - paramPtr->fgisl;
      //--chetan dabhi solution for clamping T4-
      if(T4 > model_.gidlclamp)
        T4=model_.gidlclamp;

      if (T4 == 0)
        T5 = CONSTEXPL_THRESHOLD;
      else
        T5 = paramPtr->kgisl / T4;
      if (T5 < CONSTEXPL_THRESHOLD)
      {
        T6 = exp(T5);
        Ggislb = -Igisl * T6 * T5 / T4;
      }
      else
      {
        T6 = CONSTMAX_EXPL;
        Ggislb = 0.0;
      }
      Ggisls *= T6;
      Ggislg *= T6;
      Igisl *= T6;
    }
    ggisls = Ggisls;
    ggislg = Ggislg;
    ggislb = Ggislb;
    // End of GISL

    // GIDL
    if (model_.mtrlMod == 0)
      T1 = (vds - paramPtr->rgidl * vgs_eff - paramPtr->egidl) / T0;
    else
      T1 = (vds - paramPtr->rgidl * vgs_eff - paramPtr->egidl + paramPtr->vfbsd) / T0;

    if ( (paramPtr->agidl <= 0.0) ||
         (paramPtr->bgidl <= 0.0) || (T1 <= 0.0) ||
         (paramPtr->cgidl < 0.0) )
      Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
    else
    {
      dT1_dVd = 1 / T0;
      dT1_dVg = - paramPtr->rgidl * dT1_dVd * dvgs_eff_dvg;
      T2 = paramPtr->bgidl / T1;
      if (T2 < CONSTEXPL_THRESHOLD)
      {
        Igidl = paramPtr->weffCJ * paramPtr->agidl * T1 * exp(-T2);
        T3 = Igidl / T1 * (T2 + 1);
        Ggidld = T3 * dT1_dVd;
        Ggidlg = T3 * dT1_dVg;
      }
      else
      {
        T3 = paramPtr->weffCJ * paramPtr->agidl * CONSTMIN_EXPL;
        Igidl = T3 * T1;
        Ggidld = T3 * dT1_dVd;
        Ggidlg = T3 * dT1_dVg;
      }
      T4 = vbd - paramPtr->fgidl;
      //--chetan dabhi solution for clamping T4-
      if(T4 > model_.gidlclamp)
        T4=model_.gidlclamp;

      if (T4 == 0)
        T5 = CONSTEXPL_THRESHOLD;
      else
        T5 = paramPtr->kgidl / T4;
      if (T5 < CONSTEXPL_THRESHOLD)
      {
        T6 = exp(T5);
        Ggidlb = -Igidl * T6 * T5 / T4;
      }
      else
      {
        T6 = CONSTMAX_EXPL;
        Ggidlb = 0.0;
      }
      Ggidld *= T6;
      Ggidlg *= T6;
      Igidl *= T6;
    }
    ggidld = Ggidld;
    ggidlg = Ggidlg;
    ggidlb = Ggidlb;
    // End of new GIDL
  }
  // End of Gidl


  // Calculate gate tunneling current
  if ((model_.igcMod != 0) || (model_.igbMod != 0))
  {
    Vfb = vfbzb;
    V3 = Vfb - Vgs_eff + Vbseff - CONSTDELTA_3;
    if (Vfb <= 0.0)
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
    else
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
    T1 = 0.5 * (1.0 + V3 / T0);
    Vfbeff = Vfb - 0.5 * (V3 + T0);
    dVfbeff_dVg = T1 * dVgs_eff_dVg;
    dVfbeff_dVb = -T1; // WDLiu: -No surprise? No. -Good!

    Voxacc = Vfb - Vfbeff;
    dVoxacc_dVg = -dVfbeff_dVg;
    dVoxacc_dVb = -dVfbeff_dVb;
    if (Voxacc < 0.0) // WDLiu: Avoiding numerical instability.
    {
      Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;
    }

    T0 = 0.5 * paramPtr->k1ox;
    T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
    if (paramPtr->k1ox == 0.0)
    {
      Voxdepinv = dVoxdepinv_dVg = dVoxdepinv_dVd = dVoxdepinv_dVb = 0.0;
    }
    else if (T3 < 0.0)
    {
      Voxdepinv = -T3;
      dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg + dVgsteff_dVg;
      dVoxdepinv_dVd = dVgsteff_dVd;
      dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
    }
    else
    {
      T1 = sqrt(T0 * T0 + T3);
      T2 = T0 / T1;
      Voxdepinv = paramPtr->k1ox * (T1 - T0);
      dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
      dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
    }

    Voxdepinv += Vgsteff;
    dVoxdepinv_dVg += dVgsteff_dVg;
    dVoxdepinv_dVd += dVgsteff_dVd;
    dVoxdepinv_dVb += dVgsteff_dVb;
  }

  if(model_.tempMod < 2)
  {
    tmp = Vtm;
  }
  else // model_.tempMod = 2, 3
  {
    tmp = Vtm0;
  }

  if (model_.igcMod)
  {
    T0 = tmp * paramPtr->nigc;
    if(model_.igcMod == 1)
    {
      VxNVt = (Vgs_eff - model_.dtype * vth0) / T0;
      if (VxNVt > CONSTEXP_THRESHOLD)
      {
        Vaux = Vgs_eff - model_.dtype * vth0;
        dVaux_dVg = dVgs_eff_dVg;
        dVaux_dVd = 0.0;
        dVaux_dVb = 0.0;
      }
    }
    else if (model_.igcMod == 2)
    {
      VxNVt = (Vgs_eff - von) / T0;
      if (VxNVt > CONSTEXP_THRESHOLD)
      {
        Vaux = Vgs_eff - von;
        dVaux_dVg = dVgs_eff_dVg;
        dVaux_dVd = -dVth_dVd;
        dVaux_dVb = -dVth_dVb;
      }
    }

    if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
    }
    else if ((VxNVt >= -CONSTEXP_THRESHOLD) && (VxNVt <= CONSTEXP_THRESHOLD))
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
      if(model_.igcMod == 1)
      {
        dVaux_dVd = 0.0;
        dVaux_dVb = 0.0;
      }
      else if (model_.igcMod == 2)
      {
        dVaux_dVd = -dVaux_dVg * dVth_dVd; // Synopsis 08/30/2013 modify
        dVaux_dVb = -dVaux_dVg * dVth_dVb; // Synopsis 08/30/2013 modify
      }
      dVaux_dVg *= dVgs_eff_dVg;
    }

    T2 = Vgs_eff * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
    dT2_dVd = Vgs_eff * dVaux_dVd;
    dT2_dVb = Vgs_eff * dVaux_dVb;

    T11 = paramPtr->Aechvb;
    T12 = paramPtr->Bechvb;
    T3 = paramPtr->aigc * paramPtr->cigc
       - paramPtr->bigc;
    T4 = paramPtr->bigc * paramPtr->cigc;
    T5 = T12 * (paramPtr->aigc + T3 * Voxdepinv
       - T4 * Voxdepinv * Voxdepinv);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
      dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
      dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
      dT6_dVg *= dVoxdepinv_dVg;
    }

    Igc = T11 * T2 * T6;
    dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
    dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

    if (model_.pigcdGiven)
    {
      Pigcd = paramPtr->pigcd;
      dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
    }
    else
    {
      T11 = -paramPtr->Bechvb;
      T12 = Vgsteff + 1.0e-20;
      T13 = T11 / T12 / T12;
      T14 = -T13 / T12;
      Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12);
      dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg
                    - 3.0 * Vdseff / T12));
      dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
      dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
    }

    T7 = -Pigcd * Vdseff; // bugfix
    dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
    dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
    dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
    dT7_dVg *= dVgsteff_dVg;
    // dT7_dVb *= dVbseff_dVb;   /* Synopsis, 2013/08/30 */
    T8 = T7 * T7 + 2.0e-4;
    dT8_dVg = 2.0 * T7;
    dT8_dVd = dT8_dVg * dT7_dVd;
    dT8_dVb = dT8_dVg * dT7_dVb;
    dT8_dVg *= dT7_dVg;

    if (T7 > CONSTEXP_THRESHOLD)
    {
      T9 = CONSTMAX_EXP;
      dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
    }
    else if (T7 < -CONSTEXP_THRESHOLD)
    {
      T9 = CONSTMIN_EXP;
      dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
    }
    else
    {
      T9 = exp(T7);
      dT9_dVg = T9 * dT7_dVg;
      dT9_dVd = T9 * dT7_dVd;
      dT9_dVb = T9 * dT7_dVb;
    }

    T0 = T8 * T8;
    T1 = T9 - 1.0 + 1.0e-4;
    T10 = (T1 - T7) / T8;
    dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
    dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
    dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

    Igcs = Igc * T10;
    dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
    dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
    dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

    T1 = T9 - 1.0 - 1.0e-4;
    T10 = (T7 * T9 - T1) / T8;
    dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg
             - T10 * dT8_dVg) / T8;
    dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd
             - T10 * dT8_dVd) / T8;
    dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb
             - T10 * dT8_dVb) / T8;
    Igcd = Igc * T10;
    dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
    dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
    dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

    //    Igcs = Igcs;
    gIgcsg = dIgcs_dVg;
    gIgcsd = dIgcs_dVd;
    gIgcsb =  dIgcs_dVb * dVbseff_dVb;
    //    Igcd = Igcd;
    gIgcdg = dIgcd_dVg;
    gIgcdd = dIgcd_dVd;
    gIgcdb = dIgcd_dVb * dVbseff_dVb;

    T0 = vgs - (paramPtr->vfbsd + paramPtr->vfbsdoff);
    vgs_eff = sqrt(T0 * T0 + 1.0e-4);
    dvgs_eff_dvg = T0 / vgs_eff;

    T2 = vgs * vgs_eff;
    dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
    T11 = paramPtr->AechvbEdgeS;
    T12 = paramPtr->BechvbEdge;
    T3 = paramPtr->aigs * paramPtr->cigs
       - paramPtr->bigs;
    T4 = paramPtr->bigs * paramPtr->cigs;
    T5 = T12 * (paramPtr->aigs + T3 * vgs_eff
       - T4 * vgs_eff * vgs_eff);
    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
                * dvgs_eff_dvg;
    }
    Igs = T11 * T2 * T6;
    dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgs_dVs = -dIgs_dVg;

    T0 = vgd - (paramPtr->vfbsd + paramPtr->vfbsdoff);
    vgd_eff = sqrt(T0 * T0 + 1.0e-4);
    dvgd_eff_dvg = T0 / vgd_eff;

    T2 = vgd * vgd_eff;
    dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
    T11 = paramPtr->AechvbEdgeD;
    T3 = paramPtr->aigd * paramPtr->cigd
      - paramPtr->bigd;
    T4 = paramPtr->bigd * paramPtr->cigd;
    T5 = T12 * (paramPtr->aigd + T3 * vgd_eff
       - T4 * vgd_eff * vgd_eff);
    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff) * dvgd_eff_dvg;
    }
    Igd = T11 * T2 * T6;
    dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgd_dVd = -dIgd_dVg;

    //Igs = Igs;
    gIgsg = dIgs_dVg;
    gIgss = dIgs_dVs;
    //Igd = Igd;
    gIgdg = dIgd_dVg;
    gIgdd = dIgd_dVd;

  }
  else
  {
    Igcs = gIgcsg = gIgcsd = gIgcsb = 0.0;
    Igcd = gIgcdg = gIgcdd = gIgcdb = 0.0;
    Igs = gIgsg = gIgss = 0.0;
    Igd = gIgdg = gIgdd = 0.0;
  }

  if (model_.igbMod)
  {
    T0 = tmp * paramPtr->nigbacc;
    T1 = -Vgs_eff + Vbseff + Vfb;
    VxNVt = T1 / T0;
    if (VxNVt > CONSTEXP_THRESHOLD)
    {
      Vaux = T1;
      dVaux_dVg = -dVgs_eff_dVg;
      dVaux_dVb = 1.0;
    }
    else if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVb = 0.0;
    }
    else
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt);
      dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
    }

    T2 = (Vgs_eff - Vbseff) * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
    dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

    T11 = 4.97232e-7 * paramPtr->weff
            * paramPtr->leff * paramPtr->ToxRatio;
    T12 = -7.45669e11 * toxe;
    T3 = paramPtr->aigbacc * paramPtr->cigbacc
           - paramPtr->bigbacc;
    T4 = paramPtr->bigbacc * paramPtr->cigbacc;
    T5 = T12 * (paramPtr->aigbacc + T3 * Voxacc
         - T4 * Voxacc * Voxacc);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
      dT6_dVb = dT6_dVg * dVoxacc_dVb;
      dT6_dVg *= dVoxacc_dVg;
    }

    Igbacc = T11 * T2 * T6;
    dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);


    T0 = tmp * paramPtr->nigbinv;
    T1 = Voxdepinv - paramPtr->eigbinv;
    VxNVt = T1 / T0;
    if (VxNVt > CONSTEXP_THRESHOLD)
    {
      Vaux = T1;
      dVaux_dVg = dVoxdepinv_dVg;
      dVaux_dVd = dVoxdepinv_dVd;
      dVaux_dVb = dVoxdepinv_dVb;
    }
    else if (VxNVt < -CONSTEXP_THRESHOLD)
    {
      Vaux = T0 * log(1.0 + CONSTMIN_EXP);
      dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
    }
    else
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
      dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
      dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
      dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
      dVaux_dVg *= dVoxdepinv_dVg;
    }

    T2 = (Vgs_eff - Vbseff) * Vaux;
    dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
    dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
    dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

    T11 *= 0.75610;
    T12 *= 1.31724;
    T3 = paramPtr->aigbinv * paramPtr->cigbinv
       - paramPtr->bigbinv;
    T4 = paramPtr->bigbinv * paramPtr->cigbinv;
    T5 = T12 * (paramPtr->aigbinv + T3 * Voxdepinv
       - T4 * Voxdepinv * Voxdepinv);

    if (T5 > CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMAX_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else if (T5 < -CONSTEXP_THRESHOLD)
    {
      T6 = CONSTMIN_EXP;
      dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
    }
    else
    {
      T6 = exp(T5);
      dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
      dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
      dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
      dT6_dVg *= dVoxdepinv_dVg;
    }

    Igbinv = T11 * T2 * T6;
    dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
    dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
    dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

    Igb = Igbinv + Igbacc;
    gIgbg = dIgbinv_dVg + dIgbacc_dVg;
    gIgbd = dIgbinv_dVd;
    gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
  }
  else
  {
    Igb = gIgbg = gIgbd = gIgbs = gIgbb = 0.0;
  } // End of Gate current

  if (nf != 1.0)
  {
    cdrain *= nf;
    gds *= nf;
    gm *= nf;
    gmbs *= nf;
    IdovVds *= nf;

    gbbs *= nf;
    gbgs *= nf;
    gbds *= nf;
    csub *= nf;

    Igidl *= nf;
    ggidld *= nf;
    ggidlg *= nf;
    ggidlb *= nf;

    Igisl *= nf;
    ggisls *= nf;
    ggislg *= nf;
    ggislb *= nf;

    Igcs *= nf;
    gIgcsg *= nf;
    gIgcsd *= nf;
    gIgcsb *= nf;
    Igcd *= nf;
    gIgcdg *= nf;
    gIgcdd *= nf;
    gIgcdb *= nf;

    Igs *= nf;
    gIgsg *= nf;
    gIgss *= nf;
    Igd *= nf;
    gIgdg *= nf;
    gIgdd *= nf;

    Igb *= nf;
    gIgbg *= nf;
    gIgbd *= nf;
    gIgbb *= nf;
  }

  ggidls = -(ggidld + ggidlg + ggidlb);
  ggisld = -(ggisls + ggislg + ggislb);
  gIgbs = -(gIgbg + gIgbd + gIgbb);
  gIgcss = -(gIgcsg + gIgcsd + gIgcsb);
  gIgcds = -(gIgcdg + gIgcdd + gIgcdb);
  cd = cdrain;

  // Calculations for noise analysis

  if (model_.tnoiMod == 0)
  {
    Abulk = Abulk0 * paramPtr->abulkCVfactor;
    Vdsat = Vgsteff / Abulk;
    T0 = Vdsat - Vds - CONSTDELTA_4;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * Vdsat);
    if (T0 >= 0.0)
    {
      Vdseff = Vdsat - 0.5 * (T0 + T1);
    }
    else
    {
      T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
      T4 = 1.0 - T3;
      T5 = Vdsat * T3 / (T1 - T0);
      Vdseff = Vdsat * T4;
    }
    if (Vds == 0.0)
    {
      Vdseff = 0.0;
    }

    T0 = Abulk * Vdseff;
    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
    T2 = Vdseff / T1;
    T3 = T0 * T2;
    qinv = Coxeff * paramPtr->weffCV * nf
                    * paramPtr->leffCV
                    * (Vgsteff - 0.5 * T0 + Abulk * T3);
  }
  else if (model_.tnoiMod == 2)
  {
    noiGd0 = nf * beta * Vgsteff / (1.0 + gche * Rds);
  }

  // C-V begins

  if ((model_.xpart < 0) || (!ChargeComputationNeeded))
  {
    qgate  = qdrn = qsrc = qbulk = 0.0;
    cggb = cgsb = cgdb = 0.0;
    cdgb = cdsb = cddb = 0.0;
    cbgb = cbsb = cbdb = 0.0;
    csgb = cssb = csdb = 0.0;
    cgbb = csbb = cdbb = cbbb = 0.0;
    cqdb = cqsb = cqgb = cqbb = 0.0;
    gtau = 0.0;
  }
  else if (model_.capMod == 0)
  {
    if (Vbseff < 0.0)
    {
      VbseffCV = Vbs;
      dVbseffCV_dVb = 1.0;
    }
    else
    {
      VbseffCV = paramPtr->phi - Phis;
      dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb;
    }

    Vfb = paramPtr->vfbcv;
    Vth = Vfb + paramPtr->phi + paramPtr->k1ox * sqrtPhis;
    Vgst = Vgs_eff - Vth;
    dVth_dVb = paramPtr->k1ox * dsqrtPhis_dVb * dVbseff_dVb;
    dVgst_dVb = -dVth_dVb;
    dVgst_dVg = dVgs_eff_dVg;

    CoxWL = model_.coxe * paramPtr->weffCV * paramPtr->leffCV * nf;
    Arg1 = Vgs_eff - VbseffCV - Vfb;

    if (Arg1 <= 0.0)
    {
      qgate = CoxWL * Arg1;
      qbulk = -qgate;
      qdrn = 0.0;

      cggb = CoxWL * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -CoxWL * dVgs_eff_dVg;
      cbdb = 0.0;
      cbsb = -cgsb;
    } // Arg1 <= 0.0, end of accumulation
    else if (Vgst <= 0.0)
    {
      T1 = 0.5 * paramPtr->k1ox;
      T2 = sqrt(T1 * T1 + Arg1);
      qgate = CoxWL * paramPtr->k1ox * (T2 - T1);
      qbulk = -qgate;
      qdrn = 0.0;

      T0 = CoxWL * T1 / T2;
      cggb = T0 * dVgs_eff_dVg;
      cgdb = 0.0;
      cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);

      cdgb = 0.0;
      cddb = 0.0;
      cdsb = 0.0;

      cbgb = -cggb;
      cbdb = 0.0;
      cbsb = -cgsb;
    } // Vgst <= 0.0, end of depletion
    else
    {
      One_Third_CoxWL = CoxWL / 3.0;
      Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb * dVbseff_dVb;
      dVdsat_dVg = 1.0 / AbulkCV;
      Vdsat = Vgst * dVdsat_dVg;
      dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb) * dVdsat_dVg;

      if (model_.xpart > 0.5)
      {
        // 0/100 Charge partition model
        if (Vdsat <= Vds)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb
              - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.0;

          cggb = One_Third_CoxWL * (3.0
            - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = 0.0;
          cddb = 0.0;
          cdsb = 0.0;

          cbgb = -(cggb
            - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          T7 = 2.0 * Vds - T1 - 3.0 * T3;
          T8 = T3 - T1 - 2.0 * Vds;
          qgate = CoxWL * (Vgs_eff - Vfb
                - paramPtr->phi - 0.5 * (Vds - T3));
          T10 = T4 * T8;
          qdrn = T4 * T7;
          qbulk = -(qgate + qdrn + T10);

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                          * dVgs_eff_dVg;
          T11 = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + T11 + cgdb);
          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
          T7 = T9 * T7;
          T8 = T9 * T8;
          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
          cdgb = (T7 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
          cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
          cdsb = -(cdgb + T12 + cddb);

          T9 = 2.0 * T4 * (1.0 + T5);
          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
              * dVgs_eff_dVg;
          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
          T12 = T4 * (2.0 * T2 + T5 - 1.0);
          T0 = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T10);
          cbdb = -(cgdb + cddb + T12);
          cbsb = -(cgsb + cdsb + T0);
        }
      }
      else if (model_.xpart < 0.5)
      {   // 40/60 Charge partition model
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.4 * T2;

          cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          T3 = 0.4 * Two_Third_CoxWL;
          cdgb = -T3 * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = T3 * dVth_dVb; cdsb = -(T4 + cdgb);

          cbgb = -(cggb - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi
          - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
              * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb
              + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
             + 1.2 * Vds * Vds;
          T8 = T2 / T1;
          T7 = Vds - T1 - T8 * T6;
          qdrn = T4 * T7;
          T7 *= T9;
          tmp = T8 / T1;
          tmp1 = T4 * (2.0 - 4.0 * tmp * T6
               + T8 * (16.0 * Vdsat - 6.0 * Vds));

          cdgb = (T7 * dAlphaz_dVg - tmp1
            * dVdsat_dVg) * dVgs_eff_dVg;
          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
          cddb = T4 * (2.0 - (1.0 / (3.0 * T1
            * T1) + 2.0 * tmp) * T6 + T8
            * (6.0 * Vdsat - 2.4 * Vds));
          cdsb = -(cdgb + T10 + cddb);

          T7 = 2.0 * (T1 + T3);
          qbulk = -(qgate - T4 * T7);
          T7 *= T9;
          T0 = 4.0 * T4 * (1.0 - T5);
          T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg
              - cdgb;
          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5) - cddb;
          tmp = -(T10 + T11 + T12);

          cbgb = -(cggb + cdgb + T12);
          cbdb = -(cgdb + cddb + T10);
          cbsb = -(cgsb + cdsb + tmp);
        }
      }
      else
      {   // 50/50 partitioning
        if (Vds >= Vdsat)
        {   // saturation region
          T1 = Vdsat / 3.0;
          qgate = CoxWL * (Vgs_eff - Vfb
                - paramPtr->phi - T1);
          T2 = -Two_Third_CoxWL * Vgst;
          qbulk = -(qgate + T2);
          qdrn = 0.5 * T2;

          cggb = One_Third_CoxWL * (3.0
                          - dVdsat_dVg) * dVgs_eff_dVg;
          T2 = -One_Third_CoxWL * dVdsat_dVb;
          cgsb = -(cggb + T2);
          cgdb = 0.0;

          cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
          cddb = 0.0;
          T4 = One_Third_CoxWL * dVth_dVb;
          cdsb = -(T4 + cdgb);

          cbgb = -(cggb
                            - Two_Third_CoxWL * dVgs_eff_dVg);
          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
          cbsb = -(cbgb + T3);
          cbdb = 0.0;
        }
        else
        {   // linear region
          Alphaz = Vgst / Vdsat;
          T1 = 2.0 * Vdsat - Vds;
          T2 = Vds / (3.0 * T1);
          T3 = T2 * Vds;
          T9 = 0.25 * CoxWL;
          T4 = T9 * Alphaz;
          qgate = CoxWL * (Vgs_eff - Vfb - paramPtr->phi
          - 0.5 * (Vds - T3));

          T5 = T3 / T1;
          cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
              * dVgs_eff_dVg;
          tmp = -CoxWL * T5 * dVdsat_dVb;
          cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
          cgsb = -(cggb + cgdb + tmp);

          T6 = 1.0 / Vdsat;
          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

          T7 = T1 + T3;
          qdrn = -T4 * T7;
          qbulk = - (qgate + qdrn + qdrn);
          T7 *= T9;
          T0 = T4 * (2.0 * T5 - 2.0);

          cdgb = (T0 * dVdsat_dVg - T7 * dAlphaz_dVg) * dVgs_eff_dVg;
          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
          cddb = T4 * (1.0 - 2.0 * T2 - T5);
          cdsb = -(cdgb + T12 + cddb);

          cbgb = -(cggb + 2.0 * cdgb);
          cbdb = -(cgdb + 2.0 * cddb);
          cbsb = -(cgsb + 2.0 * cdsb);
        } // end of linear region
      } // end of 50/50 partition
    } // end of inversion
  } // end of capMod=0
  else
  {
    if (Vbseff < 0.0)
    {
      VbseffCV = Vbseff;
      dVbseffCV_dVb = 1.0;
    }
    else
    {
      VbseffCV = paramPtr->phi - Phis;
      dVbseffCV_dVb = -dPhis_dVb;
    }

    CoxWL = model_.coxe * paramPtr->weffCV
          * paramPtr->leffCV * nf;

    if (model_.cvchargeMod == 0)
    {
      // Seperate VgsteffCV with noff and voffcv
      noff = n * paramPtr->noff;
      dnoff_dVd = paramPtr->noff * dn_dVd;
      dnoff_dVb = paramPtr->noff * dn_dVb;
      T0 = Vtm * noff;
      voffcv = paramPtr->voffcv;
      VgstNVt = (Vgst - voffcv) / T0;

      if (VgstNVt > CONSTEXP_THRESHOLD)
      {
        Vgsteff = Vgst - voffcv;
        dVgsteff_dVg = dVgs_eff_dVg;
        dVgsteff_dVd = -dVth_dVd;
        dVgsteff_dVb = -dVth_dVb;
      }
      else if (VgstNVt < -CONSTEXP_THRESHOLD)
      {
        Vgsteff = T0 * log(1.0 + CONSTMIN_EXP);
        dVgsteff_dVg = 0.0;
        dVgsteff_dVd = Vgsteff / noff;
        dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
        dVgsteff_dVd *= dnoff_dVd;
      }
      else
      {
        ExpVgst = exp(VgstNVt);
        Vgsteff = T0 * log(1.0 + ExpVgst);
        dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
        dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                        / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
        dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                        / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
        dVgsteff_dVg *= dVgs_eff_dVg;
      }
      // End of VgsteffCV for cvchargeMod = 0
    }
    else
    {
      T0 = n * Vtm;
      T1 = paramPtr->mstarcv * Vgst;
      T2 = T1 / T0;
      if (T2 > CONSTEXP_THRESHOLD)
      {
		      T10 = T1;
		      dT10_dVg = paramPtr->mstarcv * dVgs_eff_dVg;
		      dT10_dVd = -dVth_dVd * paramPtr->mstarcv;
		      dT10_dVb = -dVth_dVb * paramPtr->mstarcv;
      }
      else if (T2 < -CONSTEXP_THRESHOLD)
      {
        T10 = Vtm * log(1.0 + CONSTMIN_EXP);
        dT10_dVg = 0.0;
        dT10_dVd = T10 * dn_dVd;
        dT10_dVb = T10 * dn_dVb;
        T10 *= n;
      }
      else
      {
        ExpVgst = exp(T2);
        T3 = Vtm * log(1.0 + ExpVgst);
        T10 = n * T3;
        dT10_dVg = paramPtr->mstarcv * ExpVgst / (1.0 + ExpVgst);
        dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
        dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
        dT10_dVg *= dVgs_eff_dVg;
      }

      T1 = paramPtr->voffcbncv - (1.0 - paramPtr->mstarcv) * Vgst;
      T2 = T1 / T0;
      if (T2 < -CONSTEXP_THRESHOLD)
      {
        T3 = model_.coxe * CONSTMIN_EXP / paramPtr->cdep0;
        T9 = paramPtr->mstarcv + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
      }
      else if (T2 > CONSTEXP_THRESHOLD)
      {
        T3 = model_.coxe * CONSTMAX_EXP / paramPtr->cdep0;
        T9 = paramPtr->mstarcv + T3 * n;
        dT9_dVg = 0.0;
        dT9_dVd = dn_dVd * T3;
        dT9_dVb = dn_dVb * T3;
      }
      else
      {
        ExpVgst = exp(T2);
        T3 = model_.coxe / paramPtr->cdep0;
        T4 = T3 * ExpVgst;
        T5 = T1 * T4 / T0;
        T9 = paramPtr->mstarcv + n * T4;
        dT9_dVg = T3 * (paramPtr->mstarcv - 1.0) * ExpVgst / Vtm;
        dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
        dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
        dT9_dVg *= dVgs_eff_dVg;
      }

      Vgsteff = T10 / T9;
      T11 = T9 * T9;
      dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
      dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
      dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
      // End of VgsteffCV for cvchargeMod = 1
    }



    if (model_.capMod == 1)
    {
      Vfb = vfbzb;
      V3 = Vfb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (Vfb <= 0.0)
      {
        T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * Vfb);
      }
      else
      {
        T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * Vfb);
      }

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = Vfb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;
      Qac0 = CoxWL * (Vfbeff - Vfb);
      dQac0_dVg = CoxWL * dVfbeff_dVg;
      dQac0_dVb = CoxWL * dVfbeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWL;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWL * T0 / T1;
      }

      Qsub0 = CoxWL * paramPtr->k1ox * (T1 - T0);

      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb
               + dVgsteff_dVb);

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = Vgsteff / AbulkCV;

      T0 = VdsatCV - Vds - CONSTDELTA_4;
      dT0_dVg = 1.0 / AbulkCV;
      dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * VdsatCV);
      dT1_dVg = (T0 + CONSTDELTA_4 + CONSTDELTA_4) / T1;
      dT1_dVd = -T0 / T1;
      dT1_dVb = dT1_dVg * dT0_dVb;
      dT1_dVg *= dT0_dVg;
      if (T0 >= 0.0)
      {
        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
      }
      else
      {
        T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
        T4 = 1.0 - T3;
        T5 = VdsatCV * T3 / (T1 - T0);
        VdseffCV = VdsatCV * T4;
        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
      }

      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
      T2 = VdseffCV / T1;
      T3 = T0 * T2;

      T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
      T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
      T6 = 12.0 * T2 * T2 * Vgsteff;

      qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
      Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
      Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
         + Cgg1 * dVgsteff_dVb;
      Cgg1 *= dVgsteff_dVg;

      T7 = 1.0 - AbulkCV;
                qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
      T4 = -T7 * (T4 - 1.0);
      T5 = -T7 * T5;
      T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
      Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
      Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
      Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
         + Cbg1 * dVgsteff_dVb;
      Cbg1 *= dVgsteff_dVg;

      if (model_.xpart > 0.5)
      {   // 0/100 Charge petition model
        T1 = T1 + T1;
        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
           - T0 * T0 / T1);
        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else if (model_.xpart < 0.5)
      { // 40/60 Charge petition model
        T1 = T1 / 12.0;
        T2 = 0.5 * CoxWL / (T1 * T1);
        T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
           * (Vgsteff - 4.0 * T0 / 3.0))
           - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T2 * T3;
        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
           + 0.4 * T0 * T0;
        T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
           * Vgsteff - 8.0 * T0 / 3.0)
           + 2.0 * T0 * T0 / 3.0);
        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
        Csg = (T4 + T5 * dVdseffCV_dVg);
        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Csg * dVgsteff_dVb;
        Csg *= dVgsteff_dVg;
      }
      else
      { // 50/50 Charge petition model
        qsrc = -0.5 * (qgate + qbulk);
        Csg = -0.5 * (Cgg1 + Cbg1);
        Csb = -0.5 * (Cgb1 + Cbb1);
        Csd = -0.5 * (Cgd1 + Cbd1);
      }

      qgate += Qac0 + Qsub0;
      qbulk -= (Qac0 + Qsub0);
                qdrn = -(qgate + qbulk + qsrc);

      Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
      Cgd = dQsub0_dVd + Cgd1;
      Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
        + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
    }
    // Charge-Thickness capMod (CTM) begins
    else if (model_.capMod == 2)
    {
      V3 = vfbzb - Vgs_eff + VbseffCV - CONSTDELTA_3;
      if (vfbzb <= 0.0)
          T0 = sqrt(V3 * V3 - 4.0 * CONSTDELTA_3 * vfbzb);
      else
          T0 = sqrt(V3 * V3 + 4.0 * CONSTDELTA_3 * vfbzb);

      T1 = 0.5 * (1.0 + V3 / T0);
      Vfbeff = vfbzb - 0.5 * (V3 + T0);
      dVfbeff_dVg = T1 * dVgs_eff_dVg;
      dVfbeff_dVb = -T1 * dVbseffCV_dVb;

      Cox = coxp;
      Tox = 1.0e8 * toxp;
      T0 = (Vgs_eff - VbseffCV - vfbzb) / Tox;
      dT0_dVg = dVgs_eff_dVg / Tox;
      dT0_dVb = -dVbseffCV_dVb / Tox;

      tmp = T0 * paramPtr->acde;
      if ((-CONSTEXP_THRESHOLD < tmp) && (tmp < CONSTEXP_THRESHOLD))
      {
        Tcen = paramPtr->ldeb * exp(tmp);
        dTcen_dVg = paramPtr->acde * Tcen;
        dTcen_dVb = dTcen_dVg * dT0_dVb;
        dTcen_dVg *= dT0_dVg;
      }
      else if (tmp <= -CONSTEXP_THRESHOLD)
      {
        Tcen = paramPtr->ldeb * CONSTMIN_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }
      else
      {
        Tcen = paramPtr->ldeb * CONSTMAX_EXP;
        dTcen_dVg = dTcen_dVb = 0.0;
      }

      LINK = 1.0e-3 * toxp;
      V3 = paramPtr->ldeb - Tcen - LINK;
      V4 = sqrt(V3 * V3 + 4.0 * LINK * paramPtr->ldeb);
      Tcen = paramPtr->ldeb - 0.5 * (V3 + V4);
      T1 = 0.5 * (1.0 + V3 / V4);
      dTcen_dVg *= T1;
      dTcen_dVb *= T1;

      Ccen = epssub / Tcen;
      T2 = Cox / (Cox + Ccen);
      Coxeff = T2 * Ccen;
      T3 = -Ccen / Tcen;
      dCoxeff_dVg = T2 * T2 * T3;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / model_.coxe;

      Qac0 = CoxWLcen * (Vfbeff - vfbzb);
      QovCox = Qac0 / Coxeff;
      dQac0_dVg = CoxWLcen * dVfbeff_dVg
                + QovCox * dCoxeff_dVg;
      dQac0_dVb = CoxWLcen * dVfbeff_dVb
                + QovCox * dCoxeff_dVb;

      T0 = 0.5 * paramPtr->k1ox;
      T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
      if (paramPtr->k1ox == 0.0)
      {
        T1 = 0.0;
        T2 = 0.0;
      }
      else if (T3 < 0.0)
      {
        T1 = T0 + T3 / paramPtr->k1ox;
        T2 = CoxWLcen;
      }
      else
      {
        T1 = sqrt(T0 * T0 + T3);
        T2 = CoxWLcen * T0 / T1;
      }

      Qsub0 = CoxWLcen * paramPtr->k1ox * (T1 - T0);
      QovCox = Qsub0 / Coxeff;
      dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                 + QovCox * dCoxeff_dVg;
      dQsub0_dVd = -T2 * dVgsteff_dVd;
      dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                 + QovCox * dCoxeff_dVb;

      // Gate-bias dependent delta Phis begins
      if (paramPtr->k1ox <= 0.0)
      {   Denomi = 0.25 * paramPtr->moin * Vtm;
                      T0 = 0.5 * paramPtr->sqrtPhi;
      }
      else
      {   Denomi = paramPtr->moin * Vtm
           * paramPtr->k1ox * paramPtr->k1ox;
          T0 = paramPtr->k1ox * paramPtr->sqrtPhi;
      }
      T1 = 2.0 * T0 + Vgsteff;

      DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
      dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
      // End of delta Phis

      // VgDP = Vgsteff - DeltaPhi
      T0 = Vgsteff - DeltaPhi - 0.001;
      dT0_dVg = 1.0 - dDeltaPhi_dVg;
      T1 = sqrt(T0 * T0 + Vgsteff * 0.004);
      VgDP = 0.5 * (T0 + T1);
      dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);

      Tox += Tox; // WDLiu: Tcen reevaluated below due to different Vgsteff
      T0 = (Vgsteff + vtfbphi2) / Tox;
      tmp = exp(model_.bdos * 0.7 * log(T0));
      T1 = 1.0 + tmp;
      T2 = model_.bdos * 0.7 * tmp / (T0 * Tox);
      Tcen = model_.ados * 1.9e-9 / T1;
      dTcen_dVg = -Tcen * T2 / T1;
      dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
      dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
      dTcen_dVg *= dVgsteff_dVg;

      Ccen = epssub / Tcen;
      T0 = Cox / (Cox + Ccen);
      Coxeff = T0 * Ccen;
      T1 = -Ccen / Tcen;
      dCoxeff_dVg = T0 * T0 * T1;
      dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
      dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
      dCoxeff_dVg *= dTcen_dVg;
      CoxWLcen = CoxWL * Coxeff / model_.coxe;

      AbulkCV = Abulk0 * paramPtr->abulkCVfactor;
      dAbulkCV_dVb = paramPtr->abulkCVfactor * dAbulk0_dVb;
      VdsatCV = VgDP / AbulkCV;

      T0 = VdsatCV - Vds - CONSTDELTA_4;
      dT0_dVg = dVgDP_dVg / AbulkCV;
      dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
      T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_4 * VdsatCV);
      dT1_dVg = (T0 + CONSTDELTA_4 + CONSTDELTA_4) / T1;
      dT1_dVd = -T0 / T1;
      dT1_dVb = dT1_dVg * dT0_dVb;
      dT1_dVg *= dT0_dVg;
      if (T0 >= 0.0)
      {
        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
      }
      else
      {
        T3 = (CONSTDELTA_4 + CONSTDELTA_4) / (T1 - T0);
        T4 = 1.0 - T3;
        T5 = VdsatCV * T3 / (T1 - T0);
        VdseffCV = VdsatCV * T4;
        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
      }

      if (Vds == 0.0)
      {
        VdseffCV = 0.0;
        dVdseffCV_dVg = 0.0;
        dVdseffCV_dVb = 0.0;
      }

      T0 = AbulkCV * VdseffCV;
      T1 = VgDP;
      T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
      T3 = T0 / T2;
      T4 = 1.0 - 12.0 * T3 * T3;
      T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
      T6 = T5 * VdseffCV / AbulkCV;

      qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
      QovCox = qgate / Coxeff;
      Cgg1 = CoxWLcen * (T4 * dVgDP_dVg
           + T5 * dVdseffCV_dVg);
      Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1
           * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
           + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


      T7 = 1.0 - AbulkCV;
      T8 = T2 * T2;
      T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
      T10 = T9 * dVgDP_dVg;
      T11 = -T7 * T5 / AbulkCV;
      T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

      qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
      QovCox = qbulk / Coxeff;
      Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
      Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
           * dVgsteff_dVd + QovCox * dCoxeff_dVd;
      Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
           + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
      Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

      if (model_.xpart > 0.5)
      {
        // 0/100 partition
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0
               - 0.5 * T0 * T0 / T2);
        QovCox = qsrc / Coxeff;
        T2 += T2;
        T3 = T2 * T2;
        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
        T5 = T7 * AbulkCV;
        T6 = T7 * VdseffCV;

        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
        Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
            + QovCox * dCoxeff_dVd;
        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
            + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else if (model_.xpart < 0.5)
      {
        // 40/60 partition
        T2 = T2 / 12.0;
        T3 = 0.5 * CoxWLcen / (T2 * T2);
                    T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0
                       * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T3 * T4;
        QovCox = qsrc / Coxeff;
        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
        T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0
           * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
        T7 = T6 * VdseffCV / AbulkCV;

        Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg;
        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd
            + QovCox * dCoxeff_dVd;
        Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb
            + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb;
        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
      }
      else
      {
        // 50/50 partition
        qsrc = -0.5 * qgate;
        Csg = -0.5 * Cgg1;
        Csd = -0.5 * Cgd1;
        Csb = -0.5 * Cgb1;
      }

      qgate += Qac0 + Qsub0 - qbulk;
      qbulk -= (Qac0 + Qsub0);
      qdrn = -(qgate + qbulk + qsrc);

      Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
      Cbd = Cbd1 - dQsub0_dVd;
      Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

      Cgg = Cgg1 - Cbg;
      Cgd = Cgd1 - Cbd;
      Cgb = Cgb1 - Cbb;

      Cgb *= dVbseff_dVb;
      Cbb *= dVbseff_dVb;
      Csb *= dVbseff_dVb;

      cggb = Cgg;
      cgsb = -(Cgg + Cgd + Cgb);
      cgdb = Cgd;
      cdgb = -(Cgg + Cbg + Csg);
      cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                        + Csg + Csd + Csb);
      cddb = -(Cgd + Cbd + Csd);
      cbgb = Cbg;
      cbsb = -(Cbg + Cbd + Cbb);
      cbdb = Cbd;
    }  // End of CTM
  }

  csgb = - cggb - cdgb - cbgb;
  csdb = - cgdb - cddb - cbdb;
  cssb = - cgsb - cdsb - cbsb;
  cgbb = - cgdb - cggb - cgsb;
  cdbb = - cddb - cdgb - cdsb;
  cbbb = - cbgb - cbdb - cbsb;
  csbb = - cgbb - cdbb - cbbb;
  // These three lines are commented out because in SPICE they're
  // here->BSIM4qgate = qgate
  // and used only for bypass.  It's pointless to do this here:
  //  qgate = qgate;
  //  qbulk = qbulk;
  //  qdrn = qdrn;
  // This line in spice is actually:
  // here->BSIM4qsrc = -(qgate + qbulk + qdrn);
  // and that saved version is *never used* --- let's not overwrite
  // our regular qsrc this way.  Don't think it matters, but it's not
  // what spice is doing here.
  //  qsrc = -(qgate + qbulk + qdrn);

  // NQS begins
  if ((trnqsMod) || (acnqsMod))
  {
    qchqs = qcheq = -(qbulk + qgate);
    cqgb = -(cggb + cbgb);
    cqdb = -(cgdb + cbdb);
    cqsb = -(cgsb + cbsb);
    cqbb = -(cqgb + cqdb + cqsb);

    CoxWL = model_.coxe * paramPtr->weffCV * nf
          * paramPtr->leffCV;
    T1 = gcrg / CoxWL; // 1 / tau
    gtau = T1 * ScalingFactor;

    if (acnqsMod)
    {
      taunet = 1.0 / T1;
    }

  }

  // Calculate junction C-V
  if (ChargeComputationNeeded)
  {
    czbd = model_.DunitAreaTempJctCap * Adeff; // bug fix
    czbs = model_.SunitAreaTempJctCap * Aseff;
    czbdsw = model_.DunitLengthSidewallTempJctCap * Pdeff;
    czbdswg = model_.DunitLengthGateSidewallTempJctCap
            * paramPtr->weffCJ * nf;
    czbssw = model_.SunitLengthSidewallTempJctCap * Pseff;
    czbsswg = model_.SunitLengthGateSidewallTempJctCap
            * paramPtr->weffCJ * nf;

    MJS = model_.SbulkJctBotGradingCoeff;
    MJSWS = model_.SbulkJctSideGradingCoeff;
    MJSWGS = model_.SbulkJctGateSideGradingCoeff;

    MJD = model_.DbulkJctBotGradingCoeff;
    MJSWD = model_.DbulkJctSideGradingCoeff;
    MJSWGD = model_.DbulkJctGateSideGradingCoeff;

    // Source Bulk Junction
    if (vbs_jct == 0.0)
    {
      qbs = 0.0;
      capbs = czbs + czbssw + czbsswg;
    }
    else if (vbs_jct < 0.0)
    {
      if (czbs > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBS;
        if (MJS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJS * log(arg));
        }
        qbs = model_.PhiBS * czbs * (1.0 - arg * sarg) / (1.0 - MJS);
        capbs = czbs * sarg;
      }
      else
      {
        qbs = 0.0;
        capbs = 0.0;
      }
      if (czbssw > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBSWS;
        if (MJSWS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWS * log(arg));
        }
        qbs += model_.PhiBSWS * czbssw
           * (1.0 - arg * sarg) / (1.0 - MJSWS);
        capbs += czbssw * sarg;
      }
      if (czbsswg > 0.0)
      {
        arg = 1.0 - vbs_jct / model_.PhiBSWGS;
        if (MJSWGS == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWGS * log(arg));
        }
        qbs += model_.PhiBSWGS * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWGS);
        capbs += czbsswg * sarg;
      }
    }
    else
    {
      T0 = czbs + czbssw + czbsswg;
      T1 = vbs_jct * (czbs * MJS / model_.PhiBS + czbssw * MJSWS
           / model_.PhiBSWS + czbsswg * MJSWGS / model_.PhiBSWGS);

      qbs = vbs_jct * (T0 + 0.5 * T1);
      capbs = T0 + T1;
    }

    // Drain Bulk Junction
    if (vbd_jct == 0.0)
    {
      qbd = 0.0;
      capbd = czbd + czbdsw + czbdswg;
    }
    else if (vbd_jct < 0.0)
    {
      if (czbd > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBD;
        if (MJD == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJD * log(arg));
        }
        qbd = model_.PhiBD* czbd * (1.0 - arg * sarg) / (1.0 - MJD);
        capbd = czbd * sarg;
      }
      else
      {
        qbd = 0.0;
        capbd = 0.0;
      }
      if (czbdsw > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBSWD;
        if (MJSWD == 0.5)
        {
          sarg = 1.0 / sqrt(arg);
        }
        else
        {
          sarg = exp(-MJSWD * log(arg));
        }
        qbd += model_.PhiBSWD * czbdsw
         * (1.0 - arg * sarg) / (1.0 - MJSWD);
        capbd += czbdsw * sarg;
      }
      if (czbdswg > 0.0)
      {
        arg = 1.0 - vbd_jct / model_.PhiBSWGD;
        if (MJSWGD == 0.5)
            sarg = 1.0 / sqrt(arg);
        else
            sarg = exp(-MJSWGD * log(arg));
        qbd += model_.PhiBSWGD * czbdswg
           * (1.0 - arg * sarg) / (1.0 - MJSWGD);
        capbd += czbdswg * sarg;
      }
    }
    else
    {
      T0 = czbd + czbdsw + czbdswg;
      T1 = vbd_jct * (czbd * MJD / model_.PhiBD + czbdsw * MJSWD
         / model_.PhiBSWD + czbdswg * MJSWGD / model_.PhiBSWGD);
      qbd = vbd_jct * (T0 + 0.5 * T1);
      capbd = T0 + T1;
    }
  } // ChargeComputation

  if (rgateMod == 3)
  {
    vgdx = vgmd;
    vgsx = vgms;
  }
  else  // For rgateMod == 0, 1 and 2
  {
    vgdx = vgd;
    vgsx = vgs;
  }

  // gate resistor model currents.  It is not necessary to calculate these
  // directly in spice3f5, but it is necessary in Xyce.
  Igate = IgateMid = 0.0;
  if(rgateMod == 1)
  {
    Igate = grgeltd * (Vgegp);
  }
  else if(rgateMod == 2)
  {
    Igate = (gcrg) * (Vgegp);
  }
  else if(rgateMod == 3)
  {
    Igate = grgeltd * (Vgegm);
    IgateMid = gcrg * (Vgmgp);
  }

  if (model_.capMod == 0)
  {
    cgdo = paramPtr->cgdo;
    qgdo = paramPtr->cgdo * vgdx;
    cgso = paramPtr->cgso;
    qgso = paramPtr->cgso * vgsx;
  }
  else // For both capMod == 1 and 2
  {
    T0 = vgdx + CONSTDELTA_1;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
    T2 = 0.5 * (T0 - T1);

    T3 = paramPtr->weffCV * paramPtr->cgdl;
    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappad);
    cgdo = paramPtr->cgdo + T3 - T3 * (1.0 - 1.0 / T4)
           * (0.5 - 0.5 * T0 / T1);
    qgdo = (paramPtr->cgdo + T3) * vgdx - T3 * (T2
           + 0.5 * paramPtr->ckappad * (T4 - 1.0));

    T0 = vgsx + CONSTDELTA_1;
    T1 = sqrt(T0 * T0 + 4.0 * CONSTDELTA_1);
    T2 = 0.5 * (T0 - T1);
    T3 = paramPtr->weffCV * paramPtr->cgsl;
    T4 = sqrt(1.0 - 4.0 * T2 / paramPtr->ckappas);
    cgso = paramPtr->cgso + T3 - T3 * (1.0 - 1.0 / T4)
         * (0.5 - 0.5 * T0 / T1);
    qgso = (paramPtr->cgso + T3) * vgsx - T3 * (T2
         + 0.5 * paramPtr->ckappas * (T4 - 1.0));
  }

  if (nf != 1.0)
  {
    cgdo *= nf;
    cgso *= nf;
    qgdo *= nf;
    qgso *= nf;
  }
  // This silliness unnecessary, conversion from spice's
  // here->BSIM4cgdo = cgdo;
  // that stuff only needed for bypass
  //  cgdo = cgdo;
  //  qgdo = qgdo;
  //  cgso = cgso;
  //  qgso = qgso;

  setupCapacitors_oldDAE();
  setupCapacitors_newDAE ();

  // Setting up a few currents for the RHS load:
  if (model_.rdsMod == 1)
  {
    Idrain = gdtot * Vddp;
    Isource = gstot * Vssp;
  }
  else
  {
    Idrain = drainConductance * Vddp;
    Isource = sourceConductance * Vssp;
  }

  // More terms that Spice leaves out because of its formulation, but which
  // Xyce absolutely needs in the RHS.
  if (model_.rbodyMod != 0)
  {
    Idbb = grbdb * Vdbb;
    Idbbp = grbpd * Vdbbp;
    Isbb = grbsb * Vsbb;
    Isbbp = grbps * Vsbbp;
    Ibpb = grbpb * Vbpb;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : setupNoiseSources4p82_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
void Instance::setupNoiseSources4p82_ (Xyce::Analysis::NoiseData & noiseData)
{
  int numSources=NUMNOIZ;
  noiseData.numSources = numSources;
  noiseData.resize(numSources);

  noiseData.deviceName = getName().getEncodedName();

  // Note: the letter suffixes (e.g., rd) are used by the DNO() and DNI()
  // operators for .PRINT NOISE
  noiseData.noiseNames[RDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rd");  // noise due to rd
  noiseData.noiseNames[RSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rs");  // noise due to rs
  noiseData.noiseNames[RGNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rg");  // noise due to rg

  noiseData.noiseNames[RBPSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbps");
  noiseData.noiseNames[RBPDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbpd");
  noiseData.noiseNames[RBPBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbpb");
  noiseData.noiseNames[RBSBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbsb");
  noiseData.noiseNames[RBDBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_rbdb");

  noiseData.noiseNames[IDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_id");
  noiseData.noiseNames[FLNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_fl");
  noiseData.noiseNames[IGSNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igs");
  noiseData.noiseNames[IGDNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igd");
  noiseData.noiseNames[IGBNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_igb");

  noiseData.noiseNames[CORLNOIZ] = "noise_" + getName().getEncodedName()+ std::string("_corl");

  // RD thermal:
  noiseData.li_Pos[RDNOIZ] = li_DrainPrime;
  noiseData.li_Neg[RDNOIZ] = li_Drain;

  // RS thermal:
  noiseData.li_Pos[RSNOIZ] = li_SourcePrime;
  noiseData.li_Neg[RSNOIZ] = li_Source;

  // RG thermal:
  //       gNodePrime, gNodeExt,
  // Xyce: li_GatePrime,
  noiseData.li_Pos[RGNOIZ] = li_GatePrime;
  noiseData.li_Neg[RGNOIZ] = li_GateExt;

  int bodymode = 5;
  if (rbodyMod == 2)
  {
    if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    {
      bodymode = 1;
    }
    else if( (!model_.rbsbx0Given && !model_.rbsby0Given) || (!model_.rbdbx0Given && !model_.rbdby0Given) )
    {
      bodymode = 3;
    }
  }

  if (rbodyMod)
  {
    if(bodymode == 5) // all 5 noise sources (RBPS, RBPD, RBPB, RBSB, RBDB)
    {
      // RBPS nodes and conductance: bNodePrime, sbNode, grbps
      //                       Xyce:  li_BodyPrime, li_SourceBody
      noiseData.li_Pos[RBPSNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPSNOIZ] = li_SourceBody;

      // RBPD nodes and conductance: bNodePrime, dbNode, grbpd
      //                       Xyce:  li_BodyPrime, li_DrainBody
      noiseData.li_Pos[RBPDNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPDNOIZ] = li_DrainBody;

      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      // RBSB nodes and conductance: bNode, sbNode, grbsb
      //                       Xyce:  li_Body, li_SourceBody
      noiseData.li_Pos[RBSBNOIZ] = li_Body;
      noiseData.li_Neg[RBSBNOIZ] = li_SourceBody;

      // RBDB nodes and conductance: bNode, dbNode, grbdb
      //                       Xyce:  li_Body, li_DrainBody
      noiseData.li_Pos[RBDBNOIZ] = li_Body;
      noiseData.li_Neg[RBDBNOIZ] = li_DrainBody;

    }
    if(bodymode == 3) // 3 out of 5 noise sources (RBPS, RBPD, RBPB) (RBSB and RBDB excluded)
    {
      // RBPS nodes and conductance: bNodePrime, sbNode, grbps
      //                       Xyce:  li_BodyPrime, li_SourceBody
      noiseData.li_Pos[RBPSNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPSNOIZ] = li_SourceBody;

      // RBPD nodes and conductance: bNodePrime, dbNode, grbpd
      //                       Xyce:  li_BodyPrime, li_DrainBody
      noiseData.li_Pos[RBPDNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPDNOIZ] = li_DrainBody;

      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      noiseData.li_Pos[RBSBNOIZ] = -1;
      noiseData.li_Neg[RBSBNOIZ] = -1;
      noiseData.li_Pos[RBDBNOIZ] = -1;
      noiseData.li_Neg[RBDBNOIZ] = -1;
    }
    if(bodymode == 1) // 1 out of 5 noise sources (only RBPB)
    {
      // RBPB nodes and conductance: bNodePrime, bNode, grbpb
      //                       Xyce:  li_BodyPrime, li_Body
      noiseData.li_Pos[RBPBNOIZ] = li_BodyPrime;
      noiseData.li_Neg[RBPBNOIZ] = li_Body;

      noiseData.li_Pos[RBPSNOIZ] = -1;
      noiseData.li_Neg[RBPSNOIZ] = -1;
      noiseData.li_Pos[RBPDNOIZ] = -1;
      noiseData.li_Neg[RBPDNOIZ] = -1;
      noiseData.li_Pos[RBSBNOIZ] = -1;
      noiseData.li_Neg[RBSBNOIZ] = -1;
      noiseData.li_Pos[RBDBNOIZ] = -1;
      noiseData.li_Neg[RBDBNOIZ] = -1;
    }
  }
  else
  {
    noiseData.li_Pos[RBPSNOIZ] = -1;
    noiseData.li_Neg[RBPSNOIZ] = -1;
    noiseData.li_Pos[RBPDNOIZ] = -1;
    noiseData.li_Neg[RBPDNOIZ] = -1;
    noiseData.li_Pos[RBPBNOIZ] = -1;
    noiseData.li_Neg[RBPBNOIZ] = -1;
    noiseData.li_Pos[RBSBNOIZ] = -1;
    noiseData.li_Neg[RBSBNOIZ] = -1;
    noiseData.li_Pos[RBDBNOIZ] = -1;
    noiseData.li_Neg[RBDBNOIZ] = -1;
  }

  if (model_.tnoiMod == 2)
  {
    if (mode >= 0)
    {
      noiseData.li_Pos[CORLNOIZ] = li_DrainPrime;
      noiseData.li_Neg[CORLNOIZ] = li_SourcePrime;
      noiseData.li_PosCorl[CORLNOIZ] = li_GatePrime;
      noiseData.li_NegCorl[CORLNOIZ] = li_SourcePrime;
    }
    else
    {
      noiseData.li_Pos[CORLNOIZ] = li_SourcePrime;
      noiseData.li_Neg[CORLNOIZ] = li_DrainPrime;
      noiseData.li_PosCorl[CORLNOIZ] = li_GatePrime;
      noiseData.li_NegCorl[CORLNOIZ] = li_DrainPrime;
    }
  }
  else
  {
    noiseData.li_Pos[CORLNOIZ] = -1;
    noiseData.li_Neg[CORLNOIZ] = -1;
    noiseData.li_PosCorl[CORLNOIZ] = -1;
    noiseData.li_NegCorl[CORLNOIZ] = -1;
  }

  //  dNodePrime, sNodePrime,
  //  Xyce:  li_DrainPrime, li_SourcePrime
  noiseData.li_Pos[IDNOIZ] = li_DrainPrime;
  noiseData.li_Neg[IDNOIZ] = li_SourcePrime;

  //dNodePrime, sNodePrime,
  //  Xyce:  li_DrainPrime, li_SourcePrime
  noiseData.li_Pos[FLNOIZ] = li_DrainPrime;
  noiseData.li_Neg[FLNOIZ] = li_SourcePrime;

  //gNodePrime, sNodePrime,
  //  Xyce:  li_GatePrime, li_SourcePrime
  noiseData.li_Pos[IGSNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGSNOIZ] = li_SourcePrime;

  //gNodePrime, dNodePrime,
  //  Xyce:  li_GatePrime, li_DrainPrime
  noiseData.li_Pos[IGDNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGDNOIZ] = li_DrainPrime;

  //gNodePrime, bNodePrime,
  //  Xyce:  li_GatePrime, li_BodyPrime;
  noiseData.li_Pos[IGBNOIZ] = li_GatePrime;
  noiseData.li_Neg[IGBNOIZ] = li_BodyPrime;
}

//-----------------------------------------------------------------------------
// Function      : getNoiseSources4p82_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
void Instance::getNoiseSources4p82_ (Xyce::Analysis::NoiseData & noiseData)
{
  double tmp = 0.0;
  double T1 = 0.0;
  double T2 = 0.0;
  double T3 = 0.0;
  double T4 = 0.0;
  double T5 = 0.0;
  double T6 = 0.0;
  double T7 = 0.0;
  double T8 = 0.0;
  double T9 = 0.0;
  double T10 = 0.0;
  double T11 = 0.0;
  double Ssi = 0.0;
  double Swi = 0.0;
  double npart_beta = 0.0;
  double npart_theta = 0.0;
  double igsquare = 0.0;

  // tnoiMod=2 (v4.7)
  double eta = 0.0;
  double Leff = 0.0;
  double Lvsat = 0.0;
  double gamma = 0.0;
  double delta = 0.0;
  double epsilon = 0.0;
  double GammaGd0 = 0.0;
  double npart_c = 0.0;
  double sigrat = 0.0;
  double C0 = 0.0;
  double omega = 0.0;
  double ctnoi = 0.0;
  double tau = 0.0;

  if (model_.tnoiMod == 0)
  {
    if (model_.rdsMod == 0)
    {
      gspr = sourceConductance;
      gdpr = drainConductance;

      if (grdsw > 0.0)
        tmp = 1.0 / grdsw; /* tmp used below */
      else
        tmp = 0.0;
    }
    else
    {
      gspr = gstot;
      gdpr = gdtot;
      tmp = 0.0;
    }
  }
  else if (model_.tnoiMod == 1)
  {
    T5 = Vgsteff_forNoise / EsatL;
    T5 *= T5;
    npart_beta = model_.rnoia * (1.0 + T5 * model_.tnoia * paramPtr->leff);
    npart_theta = model_.rnoib * (1.0 + T5 * model_.tnoib * paramPtr->leff);
    if (npart_theta > 0.9)
      npart_theta = 0.9;
    if (npart_theta > 0.9 * npart_beta)
      npart_theta = 0.9 * npart_beta;

    if (model_.rdsMod == 0)
    {
      gspr = sourceConductance;
      gdpr = drainConductance;
    }
    else
    {
      gspr = gstot;
      gdpr = gdtot;
    }

    if (vds >= 0.0)
    {
      gspr = gspr * (1.0 + npart_theta * npart_theta * gspr / IdovVds);  /* bugfix */
    }
    else
    {
      gdpr = gdpr * (1.0 + npart_theta * npart_theta * gdpr / IdovVds);
    }
  }
  else
  {
    // tnoiMod=2 (v4.7)
    if (model_.rdsMod == 0)
    {
      gspr = sourceConductance;
      gdpr = drainConductance;
    }
    else
    {
      gspr = gstot;
      gdpr = gdtot;
    }
  }

  devSupport.noiseSupport(noiseData.noiseDens[RDNOIZ],
                          noiseData.lnNoiseDens[RDNOIZ], THERMNOISE,
                          //dNodePrime, dNode,
                          gdpr*numberParallel,temp);

  devSupport.noiseSupport(noiseData.noiseDens[RSNOIZ],
                          noiseData.lnNoiseDens[RSNOIZ],  THERMNOISE,
                          //sNodePrime, sNode,
                          gspr*numberParallel,temp);

  if ((rgateMod == 1) || (rgateMod == 2))
  {
    devSupport.noiseSupport(noiseData.noiseDens[RGNOIZ],
                            noiseData.lnNoiseDens[RGNOIZ],  THERMNOISE,
                            //gNodePrime, gNodeExt,
                            grgeltd*numberParallel,temp);
  }
  else if (rgateMod == 3)
  {
    devSupport.noiseSupport(noiseData.noiseDens[RGNOIZ],
                            noiseData.lnNoiseDens[RGNOIZ],  THERMNOISE,
                            //gNodeMid, gNodeExt,
                            grgeltd*numberParallel,temp);
  }
  else
  {
    noiseData.noiseDens[RGNOIZ] = 0.0;
    noiseData.lnNoiseDens[RGNOIZ] = std::log(std::max(noiseData.noiseDens[RGNOIZ], N_MINLOG));
  }

  int bodymode = 5;
  if (rbodyMod == 2)
  {
    if( ( !model_.rbps0Given) || ( !model_.rbpd0Given) )
    {
      bodymode = 1;
    }
    else if( (!model_.rbsbx0Given && !model_.rbsby0Given) || (!model_.rbdbx0Given && !model_.rbdby0Given) )
    {
      bodymode = 3;
    }
  }

  if (rbodyMod)
  {
    if(bodymode == 5)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPSNOIZ],
                              noiseData.lnNoiseDens[RBPSNOIZ],  THERMNOISE,
                              //bNodePrime, sbNode,
                              grbps*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPDNOIZ],
                              noiseData.lnNoiseDens[RBPDNOIZ],  THERMNOISE,
                              //bNodePrime, dbNode,
                              grbpd*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
                              noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
                              //bNodePrime, bNode,
                              grbpb*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBSBNOIZ],
                              noiseData.lnNoiseDens[RBSBNOIZ],  THERMNOISE,
                              //bNode, sbNode,
                              grbsb*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBDBNOIZ],
                              noiseData.lnNoiseDens[RBDBNOIZ],  THERMNOISE,
                              //bNode, dbNode,
                              grbdb*numberParallel,temp);
    }
    if(bodymode == 3)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPSNOIZ],
                              noiseData.lnNoiseDens[RBPSNOIZ],  THERMNOISE,
                              //bNodePrime, sbNode,
                              grbps*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPDNOIZ],
                              noiseData.lnNoiseDens[RBPDNOIZ],  THERMNOISE,
                              //bNodePrime, dbNode,
                              grbpd*numberParallel,temp);

      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
                              noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
                              //bNodePrime, bNode,
                              grbpb*numberParallel,temp);

      noiseData.noiseDens[RBSBNOIZ] = 0.0;
      noiseData.noiseDens[RBDBNOIZ] = 0.0;
      noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));
    }
    if(bodymode == 1)
    {
      devSupport.noiseSupport(noiseData.noiseDens[RBPBNOIZ],
                              noiseData.lnNoiseDens[RBPBNOIZ],  THERMNOISE,
                              //bNodePrime, bNode,
                              grbpb*numberParallel,temp);

      noiseData.noiseDens[RBPSNOIZ] = 0.0;
      noiseData.noiseDens[RBPDNOIZ] = 0.0;

      noiseData.noiseDens[RBSBNOIZ] = 0.0;
      noiseData.noiseDens[RBDBNOIZ] = 0.0;
      noiseData.lnNoiseDens[RBPSNOIZ] = std::log(std::max(noiseData.noiseDens[RBPSNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBPDNOIZ] = std::log(std::max(noiseData.noiseDens[RBPDNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
      noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));
    }
  }
  else
  {
    noiseData.noiseDens[RBPSNOIZ] = 0.0;
    noiseData.noiseDens[RBPDNOIZ] = 0.0;
    noiseData.noiseDens[RBPBNOIZ] = 0.0;
    noiseData.noiseDens[RBSBNOIZ] = 0.0;
    noiseData.noiseDens[RBDBNOIZ] = 0.0;
    noiseData.lnNoiseDens[RBPSNOIZ] = std::log(std::max(noiseData.noiseDens[RBPSNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBPDNOIZ] = std::log(std::max(noiseData.noiseDens[RBPDNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBPBNOIZ] = std::log(std::max(noiseData.noiseDens[RBPBNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBSBNOIZ] = std::log(std::max(noiseData.noiseDens[RBSBNOIZ], N_MINLOG));
    noiseData.lnNoiseDens[RBDBNOIZ] = std::log(std::max(noiseData.noiseDens[RBDBNOIZ], N_MINLOG));
  }

  if (model_.tnoiMod == 2)
  {
    eta = 1.0 - Vdseff_forNoise * AbovVgst2Vtm;
    T0 = 1.0 - eta;
    T1 = 1.0 + eta;
    T2 = T1 + 2.0 * Abulk_forNoise * model_.vtm / Vgsteff_forNoise;
    Leff = paramPtr->leff;
    Lvsat = Leff * (1.0 + Vdseff_forNoise / EsatL);
    T6 = Leff / Lvsat;
    gamma = T6 * (0.5 * T1 + T0 * T0 / (6.0 * T2));
    T3 = T2 * T2;
    T4 = T0 * T0;
    T5 = T3 * T3;
    delta = (T1 / T3 - (5.0 * T1 + T2) * T4 / (15.0 * T5) + T4 * T4 / (9.0 * T5 * T2)) / (6.0 * T6 * T6 * T6);
    T7 = T0 / T2;
    epsilon = (T7 - T7 * T7 * T7 / 3.0) / (6.0 * T6);

    T8 = Vgsteff_forNoise / EsatL;
    T8 *= T8;
    if (model_.versionDouble <= 4.80)
      {
        npart_c = model_.rnoic * (1.0 + T8 * model_.tnoic * Leff);
        ctnoi = epsilon / sqrt(gamma * delta) * (2.5316 * npart_c);

        npart_beta = model_.rnoia * (1.0 + T8 * model_.tnoia * Leff);
        npart_theta = model_.rnoib * (1.0 + T8 * model_.tnoib * Leff);
        gamma = gamma * (3.0 * npart_beta * npart_beta);
        delta = delta * (3.75 * npart_theta * npart_theta);

        GammaGd0 = gamma * noiGd0;
        C0 = Coxeff * paramPtr->weffCV * nf * paramPtr->leffCV;
        T0 = C0 / noiGd0;
        sigrat = T0 * sqrt(delta / gamma);
      }
      else
      {
        npart_c = model_.rnoic * (1.0 + T8
                                       * model_.tnoic * Leff);
        // Limits added for rnoia, rnoib, rnoic, tnoia, tnoib and tnoic in BSIM4.8.1
        T9 = gamma * delta ;
        if (T9 > 0)
          ctnoi   = epsilon / sqrt( gamma * delta) * (2.5316 * npart_c);
        else
          ctnoi   = 1.0 ;
        if (ctnoi > 1)
          ctnoi=1;
        if (ctnoi < 0)
          ctnoi=0;

        npart_beta = model_.rnoia * (1.0 + T8
                                          * model_.tnoia * Leff);
        npart_theta = model_.rnoib * (1.0 + T8
                                           * model_.tnoib * Leff);
        gamma = gamma * (3.0 * npart_beta * npart_beta);
        delta = delta * (3.75 * npart_theta * npart_theta);

        GammaGd0 = gamma * noiGd0;
        C0 = Coxeff * paramPtr->weffCV * nf * paramPtr->leffCV;
        T0 = C0 / noiGd0;

        if (gamma > 0 && delta > 0)
          sigrat = T0 * sqrt(delta / gamma);
        else
          sigrat = 0.0;
      }
  }

  switch(model_.tnoiMod)
  {
    case 0:
      if (model_.versionDouble <= 4.80)
      {
        T0 = ueff * fabs(qinv);
        T1 = T0 * tmp + paramPtr->leff * paramPtr->leff;

        devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ],
                                noiseData.lnNoiseDens[IDNOIZ], THERMNOISE,
                                //dNodePrime, sNodePrime,
                                (T0 / T1) * model_.ntnoi*numberParallel,temp);
      }
      else
      {
        T0 = ueff * fabs(qinv);
        T1 = T0 * tmp + paramPtr->leff * paramPtr->leff;

        devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ],
                                noiseData.lnNoiseDens[IDNOIZ], THERMNOISE,
                                //dNodePrime, sNodePrime,
                                (T0 / T1) * model_.ntnoi*numberParallel,temp);
        noiseData.noiseDens[CORLNOIZ] = 0.0;
        noiseData.lnNoiseDens[CORLNOIZ] = std::log(std::max(noiseData.noiseDens[CORLNOIZ],N_MINLOG));
      }
      break;
    case 1:
      if (model_.versionDouble <= 4.80)
      {
        T0 = gm + gmbs + gds;
        T0 *= T0;
        igsquare = npart_theta * npart_theta * T0 / IdovVds;
        T1 = npart_beta * (gm + gmbs) + gds;
        T2 = T1 * T1 / IdovVds;

        devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ],
                                noiseData.lnNoiseDens[IDNOIZ], THERMNOISE,
                                //dNodePrime, sNodePrime,
                                (T2 - igsquare)*numberParallel,temp);
      }
      else
      {
        T0 = gm + gmbs + gds;
        T0 *= T0;
        igsquare = npart_theta * npart_theta * T0 / IdovVds;
        T1 = npart_beta * (gm + gmbs) + gds;
        T2 = T1 * T1 / IdovVds;

        devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ],
                                noiseData.lnNoiseDens[IDNOIZ], THERMNOISE,
                                //dNodePrime, sNodePrime,
                                (T2 - igsquare)*numberParallel,temp);
        noiseData.noiseDens[CORLNOIZ] = 0.0;
        noiseData.lnNoiseDens[CORLNOIZ] = std::log(std::max(noiseData.noiseDens[CORLNOIZ],N_MINLOG));
      }
      break;
    case 2:
      T2 = GammaGd0;
      T3 = ctnoi * ctnoi;
      T4 = 1.0 - T3;
      devSupport.noiseSupport(noiseData.noiseDens[IDNOIZ],
                              noiseData.lnNoiseDens[IDNOIZ], THERMNOISE,
                              //dNodePrime, sNodePrime,
                              T2 * T4*numberParallel, temp);

      // Evaluate output noise due to two correlated noise sources
      omega = 2.0 * M_PI * noiseData.freq;
      T5 = omega * sigrat;
      T6 = T5 * T5;
      T7 = T6 / (1.0 + T6);
      noiseData.T0 = sqrt(T2 * T3);
      double T1 = sqrt(T2 * T7);
      noiseData.T2 = T1 * cos(0.5 * M_PI);
      noiseData.T3 = T1 * sin(0.5 * M_PI);

      if (mode >= 0)
      {
        devSupport.noiseSupport(noiseData.noiseDens[CORLNOIZ],
                                noiseData.lnNoiseDens[CORLNOIZ], THERMNOISE,
                                // dNodePrime, sNodePrime, T2 * T3,
                                // gNodePrime, sNodePrime, T2 * T7,
                                // 0.5 * M_PI,
                                1.0*numberParallel, temp);
      }
      else
      {
        devSupport.noiseSupport(noiseData.noiseDens[CORLNOIZ],
                                noiseData.lnNoiseDens[CORLNOIZ], THERMNOISE,
                                // sNodePrime, dNodePrime, T2 * T3,
                                // gNodePrime, dNodePrime, T2 * T7,
                                // 0.5 * M_PI,
                                1.0*numberParallel, temp);
      }
  }

  switch(model_.fnoiMod)
  {
    case 0:
      noiseData.noiseDens[FLNOIZ] = numberParallel * model_.kf * exp(model_.af * std::log(std::max(fabs(cd), N_MINLOG)))
        / (pow(noiseData.freq, model_.ef) * paramPtr->leff * paramPtr->leff * model_.coxe);

      break;
    case 1:
      double VdsLocal = vds;
      if (VdsLocal < 0.0)
      {
        VdsLocal = -VdsLocal;
      }

      Ssi = Eval1ovFNoise(VdsLocal, noiseData.freq, temp);
      T10 = model_.oxideTrapDensityA * CONSTboltz * temp;
      T11 = paramPtr->weff * nf * paramPtr->leff * pow(noiseData.freq, model_.ef) * 1.0e10 * nstar * nstar;
      Swi = T10 / T11 * cd * cd;
      T1 = Swi + Ssi;

      if (T1 > 0.0)
      {
        noiseData.noiseDens[FLNOIZ] = numberParallel*(Ssi * Swi) / T1;
      }
      else
      {
        noiseData.noiseDens[FLNOIZ] = 0.0;
      }
      break;
  }

  noiseData.lnNoiseDens[FLNOIZ] = std::log(std::max(noiseData.noiseDens[FLNOIZ], N_MINLOG));

  if(mode >= 0)
  {  /* bugfix  */
    devSupport.noiseSupport(noiseData.noiseDens[IGSNOIZ],
                            noiseData.lnNoiseDens[IGSNOIZ],  SHOTNOISE,
                            //gNodePrime, sNodePrime,
                            (Igs + Igcs)*numberParallel,temp);

    devSupport.noiseSupport(noiseData.noiseDens[IGDNOIZ],
                            noiseData.lnNoiseDens[IGDNOIZ],  SHOTNOISE,
                            //gNodePrime, dNodePrime,
                            (Igd + Igcd)*numberParallel,temp);
  }
  else
  {
    devSupport.noiseSupport(noiseData.noiseDens[IGSNOIZ],
                            noiseData.lnNoiseDens[IGSNOIZ],  SHOTNOISE,
                            //gNodePrime, sNodePrime,
                            (Igs + Igcd),temp);

    devSupport.noiseSupport(noiseData.noiseDens[IGDNOIZ],
                            noiseData.lnNoiseDens[IGDNOIZ],  SHOTNOISE,
                            //gNodePrime, dNodePrime,
                            (Igd + Igcs)*numberParallel,temp);
  }
  devSupport.noiseSupport(noiseData.noiseDens[IGBNOIZ],
                          noiseData.lnNoiseDens[IGBNOIZ],  SHOTNOISE,
                          //gNodePrime, bNodePrime,
                          Igb*numberParallel,temp);
}

//-----------------------------------------------------------------------------
// Function      : Instance::RdsEndIso4p82_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
int Instance::RdsEndIso4p82_
  (double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
   double nuEnd, int rgeo, int Type, double & Rend)
{
  std::string msg="";
	if (Type == 1)
	{
    switch(rgeo)
    {
      case 1:
      case 2:
      case 5:
          if (nuEnd == 0.0)
              Rend = 0.0;
          else
              Rend = Rsh * DMCG / (Weffcj * nuEnd);
          break;
      case 3:
      case 4:
      case 6:
          if ((DMCG + DMCI) == 0.0)
            msg = "(DMCG + DMCI) can not be equal to zero\n";
          if ((nuEnd == 0.0) || ((DMCG + DMCI) == 0.0))
            Rend = 0.0;
          else
            Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
          break;
      default:
          UserWarning(*this) << "Specified RGEO not matched\n";
    }
	}
	else
	{
    switch(rgeo)
    {
      case 1:
      case 3:
      case 7:
        if (nuEnd == 0.0)
          Rend = 0.0;
        else
          Rend = Rsh * DMCG / (Weffcj * nuEnd);
        break;
      case 2:
      case 4:
      case 8:
        if ((DMCG + DMCI) == 0.0)
          msg = "(DMCG + DMCI) can not be equal to zero\n";
        if ((nuEnd == 0.0) || ((DMCG + DMCI) == 0.0))
          Rend = 0.0;
        else
          Rend = Rsh * Weffcj / (3.0 * nuEnd * (DMCG + DMCI));
        break;
      default:
        UserWarning(*this) << "Specified RGEO not matched\n";
            }
	}
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams4p82_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 14 Sep 2022
//-----------------------------------------------------------------------------
bool Model::processParams4p82_ ()
{
  std::string msg;

  if (SbulkJctPotential < 0.1)
  {
    SbulkJctPotential = 0.1;
    UserWarning(*this) << "Given pbs is less than 0.1. Pbs is set to 0.1";
  }
  if (SsidewallJctPotential < 0.1)
  {
    SsidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbsws is less than 0.1. Pbsws is set to 0.1";
  }
  if (SGatesidewallJctPotential < 0.1)
  {
    SGatesidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbswgs is less than 0.1. Pbswgs is set to 0.1";
  }

  if (DbulkJctPotential < 0.1)
  {
    DbulkJctPotential = 0.1;
     UserWarning(*this) << "Given pbd is less than 0.1. Pbd is set to 0.1";
  }
  if (DsidewallJctPotential < 0.1)
  {
    DsidewallJctPotential = 0.1;
     UserWarning(*this) << "Given pbswd is less than 0.1. Pbswd is set to 0.1";
  }
  if (DGatesidewallJctPotential < 0.1)
  {
    DGatesidewallJctPotential = 0.1;
    UserWarning(*this) << "Given pbswgd is less than 0.1. Pbswgd is set to 0.1";
  }

  // If there are any time dependent parameters, set their values at for
  // the current time.
  // calculate dependent (ie computed) params
  if (!given("TOXP") )
    toxp = toxe;
  if (!given("TOXM") )
    toxm = toxe;
  if (!given("DSUB") )
    dsub = drout;

  if (!given("VTH0") )
    vth0 = (dtype == CONSTNMOS) ? 0.7 : -0.7;

  if (!given("VDDEOT"))
    vddeot= (dtype == CONSTNMOS)?1.5:-1.5;
  if (!given("EU"))
    eu =(dtype == CONSTNMOS) ? 1.67 : 1.0;;
  if (!given("UCS"))
    ucs =(dtype == CONSTNMOS) ? 1.67 : 1.0;;
  if (versionDouble <= 4.80)
  {
    if (!given("UA"))
      ua =(mobMod == 2) ? 1.0E-15 : 1.0E-9; // UNIT M/V
    if (!given("UC"))
      uc = (mobMod == 1) ? -0.0465 : -0.0465E-9;
    if (!given("UC1"))
      uc1 =(mobMod == 1) ? -0.056 : -0.056E-9;
  }
  else
  {
    if (!given("UA"))
      ua = ((mobMod == 2 || mobMod == 6)) ? 1.0e-15 : 1.0e-9; /* unit m/V */
    if (!given("UC"))
      uc = (mobMod == 1 || mobMod == 5) ? -0.0465 : -0.0465e-9;
    if (!given("UC1"))
      uc1 = (mobMod == 1 || mobMod == 5) ? -0.056 : -0.056e-9;
  }
  if (!given("U0"))
    u0 = (dtype == CONSTNMOS) ? 0.067 : 0.025;

  // NOTE!  4.8.2 has different default for fgidl than all previous versions,
  // and so we have to override default from addPar here.
  if (!given("FGIDL"))
      fgidl=1.0;

  if (!given("AGISL"))
    agisl=agidl;
  if (!given("BGISL"))
    bgisl=bgidl;
  if (!given("CGISL"))
    cgisl=cgidl;
  if (!given("EGISL"))
    egisl=egidl;
  if (!given("RGISL"))
    rgisl=rgidl;
  if (!given("KGISL"))
    kgisl=kgidl;
  if (!given("FGISL"))
    fgisl=fgidl;

  if (!given("AIGC"))
    aigc =(dtype == CONSTNMOS) ? 1.36E-2 : 9.80E-3;
  if (!given("BIGC"))
    bigc =(dtype == CONSTNMOS) ? 1.71E-3 : 7.59E-4;
  if (!given("CIGC"))
    cigc = (dtype == CONSTNMOS) ? 0.075 : 0.03;
  if (given("AIGSD"))
  {
    aigs = aigd = aigsd;
  }
  else
  {
    aigsd = (dtype == CONSTNMOS) ? 1.36E-2 : 9.80E-3;
    if (!given("AIGS"))
      aigs = aigsd;
    if (!given("AIGD"))
      aigd = aigsd;
  }

  if (given("BIGSD"))
  {
    bigs = bigd = bigsd;
  }
  else
  {
    bigsd = (dtype == CONSTNMOS) ? 1.71E-3 : 7.59E-4;
    if (!given("BIGS"))
      bigs = bigsd;
    if (!given("BIGD"))
      bigd = bigsd;
  }
  if (given("CIGSD"))
  {
    cigs=cigd=cigsd;
  }
  else
  {
    cigsd = (dtype == CONSTNMOS) ? 0.075 : 0.03;
    if (!given("CIGS"))
      cigs=cigsd;
    if (!given("CIGD"))
      cigd=cigsd;
  }
  if (!given("IJTHDFWD"))
    ijthdfwd = ijthsfwd;
  if (!given("IJTHDREV"))
    ijthdrev = ijthsrev;
  if (!given("XJBVD"))
    xjbvd = xjbvs;
  if (!given("BVD"))
    bvd = bvs;
  if (!given("CKAPPAD"))
    ckappad = ckappas;
  if (!given("DMCI"))
    dmci = dmcg;

  // Length dependence
  if (!given("LAGISL"))
    lagisl = lagidl;
  if (!given("LBGISL"))
    lbgisl = lbgidl;
  if (!given("LCGISL"))
    lcgisl = lcgidl;
  if (!given("LEGISL"))
    legisl = legidl;
  if (!given("LRGISL"))
    lrgisl = lrgidl;
  if (!given("LKGISL"))
    lkgisl = lkgidl;
  if (!given("LFGISL"))
    lfgisl = lfgidl;

  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
    laigs = laigd = laigsd;
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
    lbigs = lbigd = lbigsd;
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
    lcigs = lcigd = lcigsd;

  // Width dependence
  if (!given("WAGISL"))
    wagisl = wagidl;
  if (!given("WBGISL"))
    wbgisl = wbgidl;
  if (!given("WCGISL"))
    wcgisl = wcgidl;
  if (!given("WEGISL"))
    wegisl = wegidl;
  if (!given("WRGISL"))
    wrgisl = wrgidl;
  if (!given("WKGISL"))
    wkgisl = wkgidl;
  if (!given("WFGISL"))
    wfgisl = wfgidl;

  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
    waigs = waigd = waigsd;
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
    wbigs = wbigd = wbigsd;
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
    wcigs = wcigd = wcigsd;

  // Cross-term depdendence
  if (!given("PAGISL"))
    pagisl = pagidl;
  if (!given("PBGISL"))
    pbgisl = pbgidl;
  if (!given("PCGISL"))
    pcgisl = pcgidl;
  if (!given("PEGISL"))
    pegisl = pegidl;
  if (!given("PRGISL"))
    prgisl = prgidl;
  if (!given("PKGISL"))
    pkgisl = pkgidl;
  if (!given("PFGISL"))
    pfgisl = pfgidl;

  if (!(!given("AIGSD") && (given("AIGS") || given("AIGD"))))
    paigs = paigd = paigsd;
  if (!(!given("BIGSD") && (given("BIGS") || given("BIGD"))))
    pbigs = pbigd = pbigsd;
  if (!(!given("CIGSD") && (given("CIGS") || given("CIGD"))))
    pcigs = pcigd = pcigsd;

  if (!given("LLC"))
    Llc = Ll;
  if (!given("LWC"))
    Lwc = Lw;
  if (!given("LWLC"))
    Lwlc = Lwl;
  if (!given("WLC"))
    Wlc = Wl;
  if (!given("WWC"))
    Wwc = Ww;
  if (!given("WWLC"))
    Wwlc = Wwl;
  if (!given("DWC"))
    dwc = Wint;
  if (!given("DLC"))
    dlc = Lint;


  if (!given("DLCIG"))
    dlcig = Lint;
  if (!given("DLCIGD"))
  {
    if (!given("DLCIG"))
      dlcigd = Lint;
    else
      dlcigd = dlcig;
  }
  if (!given("DWJ"))
    dwj = dwc;

  if (!given("JSD"))
    DjctSatCurDensity=SjctSatCurDensity;
  if (!given("JSWD"))
    DjctSidewallSatCurDensity=SjctSidewallSatCurDensity;
  if (!given("JSWGD"))
    DjctGateSidewallSatCurDensity=SjctGateSidewallSatCurDensity;
  if (!given("PBD"))
    DbulkJctPotential=SbulkJctPotential;
  if (!given("NJD"))
    DjctEmissionCoeff = SjctEmissionCoeff;
  if (!given("XTID"))
    DjctTempExponent = SjctTempExponent;
  if (!given("MJD"))
    DbulkJctBotGradingCoeff = SbulkJctBotGradingCoeff;
  if (!given("MJSWD"))
    DbulkJctSideGradingCoeff = SbulkJctSideGradingCoeff;
  if (!given("MJSWGS"))
    SbulkJctGateSideGradingCoeff = SbulkJctSideGradingCoeff;
  if (!given("MJSWGD"))
    DbulkJctGateSideGradingCoeff = SbulkJctGateSideGradingCoeff;
  if (!given("PBSWD"))
    DsidewallJctPotential = SsidewallJctPotential;
  if (!given("PBSWGS"))
    SGatesidewallJctPotential = SsidewallJctPotential;
  if (!given("PBSWGD"))
    DGatesidewallJctPotential = DsidewallJctPotential;
  if (!given("CJD"))
    DunitAreaJctCap=SunitAreaJctCap;
  if (!given("CJSWD"))
    DunitLengthSidewallJctCap = SunitLengthSidewallJctCap;
  if (!given("CJSWGS"))
    SunitLengthGateSidewallJctCap = SunitLengthSidewallJctCap ;
  if (!given("CJSWGD"))
    DunitLengthGateSidewallJctCap = SunitLengthGateSidewallJctCap;

  if (!given("JTSD"))
    jtsd = jtss;
  if (!given("JTSSWD"))
    jtsswd = jtssws;
  if (!given("JTSSWGD"))
    jtsswgd = jtsswgs;

  if (!given("NJTSD"))
  {
    if (given("NJTS"))
      njtsd =  njts;
  }
  if (!given("NJTSSWD"))
  {
    if (given("NJTSSW"))
      njtsswd =  njtssw;
  }
  if (!given("NJTSSWGD"))
  {
    if (given("NJTSSWG"))
      njtsswgd =  njtsswg;
  }

  if (!given("XTSD"))
    xtsd = xtss;
  if (!given("XTSSWD"))
    xtsswd = xtssws;
  if (!given("XTSSWGD"))
    xtsswgd = xtsswgs;

  if (!given("TNJTSD"))
  {
    if (given("TNJTS"))
      tnjtsd =  tnjts;
  }
  if (!given("TNJTSSWD"))
  {
    if (given("TNJTSSW"))
      tnjtsswd =  tnjtssw;
  }
  if (!given("TNJTSSWGD"))
  {
    if (given("TNJTSSWG"))
      tnjtsswgd =  tnjtsswg;
  }

  if (!given("VTSD"))
    vtsd = vtss;
  if (!given("VTSSWD"))
    vtsswd = vtssws;
  if (!given("VTSSWGD"))
    vtsswgd = vtsswgs;

  if (!given("NOIA"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityA = 6.25e41;
    else
      oxideTrapDensityA= 6.188e40;
  }
  if (!given("NOIB"))
  {
    if (dtype == CONSTNMOS)
      oxideTrapDensityB = 3.125e26;
    else
      oxideTrapDensityB = 1.5e25;
  }
  if (!given("NOIC"))
  {
     oxideTrapDensityC = 8.75e9;
  }

  // We have changed model parameters (maybe) and so all size dependent params
  // we may have stored may be invalid.  Clear them
  clearTemperatureData();

  return true;
}

} // namespace MOSFET_B4
} // namespace Device
} // namespace Xyce
