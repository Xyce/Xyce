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
// Purpose        : Defines units and catagory definitions for documentation purposes
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Units_h
#define Xyce_N_DEV_Units_h

namespace Xyce {
namespace Device {
namespace Units {

enum ParameterUnit
  {
    LEVEL_1 = 1,
    LEVEL_2,
    LEVEL_3,
    LEVEL_4,
    LEVEL_5,
    LEVEL_6,
    LEVEL_7,
    LEVEL_8,
    LEVEL_9 = 9,
    LEVEL_10,
    LEVEL_11,
    LEVEL_12,
    LEVEL_13,
    LEVEL_14,
    LEVEL_15,
    LEVEL_16,
    LEVEL_17,
    LEVEL_18,
    LEVEL_19,
    LEVEL_20,
    LEVEL_21,
    LEVEL_22,
    LEVEL_23,
    LEVEL_24,
    LEVEL_25,
    LEVEL_26,
    LEVEL_27,
    LEVEL_28,
    LEVEL_29,
    LEVEL_30,
    U_MAXLEVEL,
    U_INVALID,
    STANDARD,              // Unit and description are taken from standard list
    U_NONE,                // Actual unit specifications start here
    U_UNKNOWN,             // unknown
    U_AMP,                 // A
    U_AMPCMM1,             // A/cm
    U_AMPCMM2,             // A/cm**2
    U_AMPMM1,              // A/m
    U_AMPMM2,              // A/m**2
    U_AMPMM3,              // A/m**3
    U_AMPVM2,              // A/V**2
    U_AMPVM3,              // A/V**3
    U_AMPVM3M,             // Am/V**3
    U_AMPSVM1METERM1,      // As/(Vm)
    U_COULOMB,             // C
    U_CM,                  // cm
    U_CM2,                 // cm**2
    U_CM2VM1SM1,           // cm**2/(Vs)
    U_CMM2,                // 1/cm**2
    U_CMM2VM1SM1,          // 1/(cm**2*V*s)
    U_CMM3,                // 1/cm**3
    U_CMSM1,               // cm/s
    U_MCMM3,               // m/cm**3
    U_M2CMM3,              // m**2/cm**3
    U_CM3SM1,              // cm**3/s
    U_CM6SM1,              // cm**6/s
    U_DEGREE,              // degree
    U_DEGC,                // degree C
    U_DEGK,                // degree K
    U_DEGKM1,              // 1/degree K
    U_DEGCM1,              // 1/degree C
    U_DEGCM2,              // 1/(degree C)**2
    U_PERCENTDEGCM1,       // %/degree C
    U_EV,                  // eV
    U_EVDEGKM1,            // eV/(degree K)
    U_FARAD,               // F
    U_FARADM,              // F*m
    U_FARADMM1,            // F/m
    U_FARADMM2,            // F/m**2
    U_FHGMHSMVM1,          // (F/g)**(1/2)s/mV
    U_FHGMHSMMVM1,         // (F/g)**(1/2)sm/mV
    U_FHGMHSM2MVM1,        // (F/g)**(1/2)sm**2/mV
    U_FS2HGMHMM1,          // (Fs**2/g)**(1/2)/m
    U_FS2HGMHMM1VM1,       // (Fs**2/g)**(1/2)/(Vm)
    U_FVM1MM2,             // F/(V*m**2)
    U_FVM1MM1,             // F/(V*m)
    U_FVM1,                // F/V
    U_GCMM3,               // g/cm^3
    U_HENRY,               // Henry
    U_HMM1,                // Henry/m
    U_HOUR,                // Hour
    U_HZ,                  // Hertz
    U_JKM1,                // J/K
    U_KKGM1JM1,            // K/(Kg*J)
    U_JMM3KM1,             // K/(Kg*J)
    U_KM1,                 // 1/K
    U_KGMM3,               // Kg/m**3
    U_LOGIC,               // True/False
    U_MCMM2VM1SM1,         // m/(cm**2*V*s)
    U_M2CMM2VM1SM1,        // m**2/(cm**2*V*s)
    U_THERMAL,             // unknown thermal
    U_METER,               // m
    U_METERM1,             // m**-1
    U_METERM2,             // m**-2
    U_METERM3,             // m**-3
    U_METER2,              // m**2
    U_METER3,              // m**3
    U_MEXPLL,              // m**(LLN)
    U_MEXPLW,              // m**(LWN)
    U_MEXPLLLW,            // m**(LLN+LWN)
    U_MEXPWL,              // m**(WLN)
    U_MEXPWW,              // m**(WWN)
    U_MEXPWLWW,            // m**(WLN+WWN)
    U_MHVMH,               // m**(1/2)/V**(1/2)
    U_MOM1,                // m/ohm
    U_MSM1,                // m/s
    U_MM3SM1,              // 1/m**3/s
    U_MVMH,                // m/V**(1/2)
    U_M2VMH,               // m**2/V**(1/2)
    U_M3VMH,               // m**2/V**(1/2)
    U_MVM1,                // m/V
    U_M2VM1,               // m**2/V
    U_M2VM1SM1,            // m**2/(V*sec)
    U_M3VM1,               // m**3/V
    U_MVM2,                // m/V**2
    U_MVM2DEGCM1,          // m/(V**2*degree C)
    U_M2OM1,               // m**2/ohm
    U_M2VM2DEGCM1,         // m**2/(V**2*degree C)
    U_M3VM2DEGCM1,         // m**3/(V**2*degree C)
    U_M2SM1,               // m**2/s
    U_M3SM1,               // m**3/s
    U_M2VM2,               // m**2/V**2
    U_M3VM2,               // m**3/V**2
    U_M4VM2,               // m**4/V**2
    U_MOLAR,               // mol/L
    U_OHM,                 // Ohm
    U_OHMM2,               // Ohm*m**2
    U_OHMMICRON,           // Ohm*micron
    U_OHMMICRONM,          // Ohm*micron*m
    U_OHMMICRONM2,         // Ohm*micron*m**2
    U_OHMM,                // Ohm*m
    U_OHMMM1,              // Ohm/m
    U_OHMMM1SM1,           // Ohm/m/s
    U_OHMM1,               // 1/Ohm
    U_OHMM1MM1,            // 1/(Ohm*m)
    U_OHMM1MM2,            // 1/(Ohm*m**2)
    U_OHMPV,               // Ohm/volt
    U_OSQM1,               // Ohm/square
    U_PERUNIT,             // per unit (used in power grid simulations)
    U_RAD,                 // rads
    U_RADPS,               // rads/sec
    U_SECOND,              // s
    U_SECM1,               // 1/s
    U_SQUARES,             // # of squares
    U_VKM1,                // V/K
    U_VHM,                 // V**(1/2)*m
    U_VHM2,                // V**(1/2)*m**2
    U_VHM3,                // V**(1/2)*m**3
    U_VM,                  // V*m
    U_VM2,                 // V*m**2
    U_VM3,                 // V*m**3
    U_VMM1,                // V/m
    U_VMM2,                // V/m**2
    U_VMM3,                // V/m**3
    U_VMMH,                // V/m**(1/2)
    U_VOLT,                // V
    U_VOLT3,               // V**3
    U_VOLTH,               // V**(1/2)
    U_VOLTMH,              // V**(-1/2)
    U_VOLTM1,              // 1/V
    U_VOLTM2,              // 1/V**2
    U_VOLTSAMPM1METERM1    // Vs/(Am)
  };

enum ParameterCategory
  {
    CAT_INVALID,
    CAT_UNKNOWN,
    CAT_NONE,

    // Start categories

    CAT_AC,
    CAT_BASIC,
    CAT_BIN,
    CAT_CAP,
    CAT_CONTROL,
    CAT_CURRENT,
    CAT_DC,
    CAT_DEPENDENCY,
    CAT_DOPING,
    CAT_FLICKER,
    CAT_GEOMETRY,
    CAT_INITIAL,
    CAT_NQS,
    CAT_MATERIAL,
    CAT_RADP,
    CAT_RES,
    CAT_PROCESS,
    CAT_RF,
    CAT_RAD,
    CAT_TEMP,
    CAT_TUNNEL,
    CAT_VBI,
    CAT_VOLT,
    CAT_ASYMRDS,
    CAT_ASYMDDS,
    CAT_IMPACT,
    CAT_GDLEAKAGE,
    CAT_STRESS,
    CAT_WELL,

    CAT_STATIC,
    CAT_DYNAMIC,
    CAT_CARRIER,
    CAT_OUTPUT,
    CAT_PULSE,

    CAT_SCALING,
    CAT_BOUNDARYCONDITIONS,

    CAT_MAX,                // End of Categories

    CAT_MASK     = 0xFF,    // Mask for extracting category from flags

    // Start Flags

    UNDOCUMENTED = 0x100,
    DEPRECATED   = 0x200,
    POSITIONAL   = 0x400,
    COMPATIBILITY = 0x800
  };


struct UnitInfo
{
  ParameterUnit         Unit;
  ParameterUnit         UnitM;
  const char *          description;
  const char *          doc;
};


struct StdDescription
{
  const char *          Name;
  ParameterUnit         Unit;
  ParameterCategory     Category;
  const char *          Description;
};

extern UnitInfo unitTable[];
extern StdDescription descriptionTable[];
extern size_t unitTableSize;
extern size_t descriptionTableSize;

} // namespace Units
} // namespace Device
} // namespace Xyce

namespace N_DEV_Units = Xyce::Device;

using namespace Xyce::Device::Units;

#endif // Xyce_N_DEV_Units_h
