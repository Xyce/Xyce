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

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 2/6/2013
//
//
//
//
//-------------------------------------------------------------------------

#include <cstdlib>

#include <N_DEV_Units.h>

// ---------- Static data ------------
namespace Xyce {
namespace Device {
namespace Units {

UnitInfo unitTable[] =
{
  U_NONE,            U_METER,          "",                       "--",
#ifdef Xyce_DEBUG_UNITS
  U_UNKNOWN,         U_INVALID,        "unknown",                "unknown",
#else
  U_UNKNOWN,         U_INVALID,        "---",                    "---",
#endif
  U_AMP,             U_INVALID,        "A",                      "A",
  U_AMPMM1,          U_AMP,            "A/m",                    "A/m",
  U_AMPMM2,          U_AMPMM1,         "A/m**2",                 "A/m$^{2}$",
  U_AMPMM3,          U_AMPMM2,         "A/m**3",                 "A/m$^{3}$",
  U_AMPVM2,          U_INVALID,        "A/V**2",                 "A/V$^{2}$",
  U_AMPVM3,          U_INVALID,        "A/V**3",                 "A/V$^{3}$",
  U_AMPVM3M,         U_INVALID,        "Am/V**3",                 "Am/V$^{3}$",
  U_AMPSVM1METERM1,  U_INVALID,        "As/(Vm)",                "As/(Vm)",
  U_COULOMB,         U_INVALID,        "C",                      "C",
  U_CM,              U_INVALID,        "cm",                     "cm",
  U_CM2,             U_INVALID,        "cm**2",                  "cm$^{2}$",
  U_CM2VM1SM1,       U_INVALID,        "cm**2/(V*s)",            "cm$^{2}$/(Vs)",
  U_CMM2,            U_INVALID,        "1/cm**2",                "cm$^{-2}$",
  U_CMM2VM1SM1,      U_MCMM2VM1SM1,    "1/(cm**2*V*s)",          "1/(Vcm$^{2}$s)",
  U_CMM3,            U_MCMM3,          "1/cm**3",                "cm$^{-3}$",
  U_MCMM3,           U_M2CMM3,         "m/cm**3",                "m/cm$^{3}$",
  U_M2CMM3,          U_INVALID,        "m**2/cm**3",             "m$^{2}$/cm$^{3}$",
  U_CM3SM1,          U_INVALID,        "cm**3/s",                "cm$^{3}$/s",
  U_CM6SM1,          U_INVALID,        "cm**6/s",                "cm$^{6}$/s",
  U_CM3SM1,          U_INVALID,        "cm**3/s",                "cm$^{3}$/s",
  U_DEGREE,          U_INVALID,        "degree",                 "degree",
  U_DEGC,            U_INVALID,        "degree C",               "$^\\circ$C",
  U_DEGK,            U_INVALID,        "degree K",               "K",
  U_DEGKM1,          U_INVALID,        "1/degree K",             "1/K",
  U_DEGCM1,          U_INVALID,        "1/(degree C)",           "$^\\circ$C$^{-1}$",
  U_DEGCM2,          U_INVALID,        "1/(degree C)**2",        "$^\\circ$C$^{-2}$",
  U_PERCENTDEGCM1,   U_INVALID,        "%/(degree C)",           "\\%$/^\\circ$C",
  U_EV,              U_INVALID,        "eV",                     "eV",
  U_EVDEGKM1,        U_INVALID,        "eV/K",                   "eV/K",
  U_FARAD,           U_FARADM,         "F",                      "F",
  U_FARADM,          U_INVALID,        "F*m",                    "Fm",
  U_FARADMM1,        U_FARAD,          "F/m",                    "F/m",
  U_FARADMM2,        U_FARADMM1,       "F/m**2",                 "F/m$^{2}$",
  U_FHGMHSMVM1,      U_FHGMHSMMVM1,    "(F/g)**(1/2)*s/mV",      "(F/g)$^{1/2}$s/mV",
  U_FHGMHSMMVM1,     U_FHGMHSM2MVM1,   "(F/g)**(1/2)*s*m/mV",    "(F/g)$^{1/2}$sm/mV",
  U_FHGMHSM2MVM1,    U_INVALID,        "(F/g)**(1/2)*s*m**2/mV", "(F/g)$^{1/2}$sm$^{2}$/mV",
  U_FS2HGMHMM1,      U_INVALID,        "(Fs**2/g)**(1/2)/m",     "(Fs$^2$/g)$^{1/2}$/m",
  U_FS2HGMHMM1VM1,   U_INVALID,        "(Fs**2/g)**(1/2)/mV",     "(Fs$^2$/g)$^{1/2}$/mV",
  U_FVM1MM2,         U_FVM1MM1,        "F/(V*m**2)",             "F/(Vm$^{2}$)",
  U_FVM1MM1,         U_FVM1,           "F/(V*m)",                "F/(Vm)",
  U_FVM1,            U_INVALID,        "F/V",                    "F/V",
  U_GCMM3,           U_INVALID,        "g/cm**3",                "g/$\\mbox{cm}^3$",
  U_HENRY,           U_INVALID,        "H",                      "henry",
  U_HMM1,            U_INVALID,        "H/m",                    "Hm$^{-1}$", 
  U_HOUR,            U_INVALID,        "hour",                   "hour",
  U_HZ,              U_INVALID,        "hz",                     "Hz",
  U_JKM1,            U_INVALID,        "Joule/(degree K)",       "J/K",
  U_KKGM1JM1,        U_INVALID,        "(degree K)/(kg*Joule)",  "K/(J-kg)",
  U_JMM3KM1,         U_INVALID,        "(Joule)/(m**3*degree K)",  "$J/(\\mbox{m}^3{}$K)",
  U_KM1,             U_INVALID,        "1/(degree K)",           "K$^{-1}$",
  U_KGMM3,           U_INVALID,        "kg/m**3",                "kg/$\\mbox{m}^3$",
  U_LOGIC,           U_INVALID,        "logical (T/F)",          "logical (T/F)",
  U_MCMM2VM1SM1,     U_M2CMM2VM1SM1,   "m/(cm**2*V*s)",          "m/(Vcm$^{2}$s)",
  U_M2CMM2VM1SM1,    U_INVALID,        "m**2/(cm**2*V*s)",       "m$^{2}$/(Vcm$^{2}$s)",
#ifdef Xyce_DEBUG_UNITS
  U_THERMAL,         U_INVALID,        "unknown thermal",        "unknown thermal",
#else
  U_THERMAL,         U_INVALID,        "---",                    "---",
#endif
  U_METER,           U_METER2,         "m",                      "m",
  U_METERM1,         U_NONE,           "1/m",                    "m$^{-1}$",
  U_METERM2,         U_METERM1,        "1/m**2",                 "m$^{-2}$",
  U_METERM3,         U_METERM2,        "1/m**3",                 "m$^{-3}$",
  U_METER2,          U_METER3,         "m**2",                   "m$^{2}$",
  U_METER3,          U_INVALID,        "m**3",                   "m$^{3}$",
  U_MEXPLL,          U_INVALID,        "m**LLN",                 "m$^{LLN}$",
  U_MEXPLW,          U_INVALID,        "m**LWN",                 "m$^{LWN}$",
  U_MEXPLLLW,        U_INVALID,        "m**(LLN+LWN)",           "m$^{LLN+LWN}$",
  U_MEXPWL,          U_INVALID,        "m**WLN",                 "m$^{WLN}$",
  U_MEXPWW,          U_INVALID,        "m**WWN",                 "m$^{WWN}$",
  U_MEXPWLWW,        U_INVALID,        "m**(WLN+WWN)",           "m$^{WLN+WWN}$",
  U_MHVMH,           U_INVALID,        "(m/s)**(1/2)",           "m$^{1/2}$/V$^{1/2}$",
  U_MOM1,            U_M2OM1,          "m/ohm",                  "m/$\\mathsf{\\Omega}$",
  U_MSM1,            U_M2SM1,          "m/s",                    "m/s",
  U_MM3SM1,          U_INVALID,        "1/m**3/s",               "m$^{-3}$s$^{-1}$",
  U_MVMH,            U_M2VMH,          "m/V**(1/2)",             "m/V$^{1/2}$",
  U_M2VMH,           U_M3VMH,          "m**2/V**(1/2)",          "m$^{2}$/V$^{1/2}$",
  U_M3VMH,           U_INVALID,        "m**3/V**(1/2)",          "m$^{3}$/V$^{1/2}$",
  U_MVM1,            U_M2VM1,          "m/V",                    "m/V",
  U_M2VM1,           U_M3VM1,          "m**2/V",                 "m$^{2}$/V",
  U_M2VM1SM1,        U_INVALID,        "m**2/(Vs)",              "m$^{2}$/(Vs)",
  U_M3VM1,           U_INVALID,        "m**3/V",                 "m$^{3}$/V",
  U_MVM2,            U_M2VM2,          "m/V**2",                 "m/V$^{2}$",
  U_MVM2DEGCM1,      U_M2VM2DEGCM1,    "m/(degree C * V**2)",    "m/($^\\circ$CV$^{2}$)",
  U_M2OM1,           U_INVALID,        "m**2/ohm",               "m$^{2}$/$\\mathsf{\\Omega}$",
  U_M2VM2DEGCM1,     U_M3VM2DEGCM1,    "m**2/(degree C * V**2)", "m$^{2}$/($^\\circ$CV$^{2}$)",
  U_M3VM2DEGCM1,     U_INVALID,        "m**3/(degree C * V**2)", "m$^{3}$/($^\\circ$CV$^{2}$)",
  U_M2SM1,           U_M3SM1,          "m**2/s",                 "m$^{2}$/s",
  U_M3SM1,           U_INVALID,        "m**3/s",                 "m$^{3}$/s",
  U_M2VM2,           U_M3VM2,          "(m/V)**2",               "m$^{2}$/V$^{2}$",
  U_M3VM2,           U_M4VM2,          "m**3/V**2",              "m$^{3}$/V$^{2}$",
  U_M4VM2,           U_INVALID,        "m**4/V**2",              "m$^{4}$/V$^{2}$",
  U_MOLAR,           U_INVALID,        "mol/L",                  "mol/L",
  U_OHM,             U_INVALID,        "ohm",                    "$\\mathsf{\\Omega}$",
  U_OHMM2,           U_INVALID,        "ohm*m**2",               "$\\mathsf{\\Omega}$ m$^{2}$",
  U_OHMMICRON,       U_OHMMICRONM,     "ohm*micron",             "$\\mathsf{\\Omega}$ $\\mu$m",
  U_OHMMICRONM,      U_OHMMICRONM2,    "ohm*micron*m",           "$\\mathsf{\\Omega}$ $\\mu$m m",
  U_OHMMICRONM2,     U_INVALID,        "ohm*micron*m**2",        "$\\mathsf{\\Omega}$ $\\mu$m m$^{2}$",
  U_OHMM,            U_INVALID,        "ohm*m",                  "$\\mathsf{\\Omega}$ m",
  U_OHMMM1,          U_OHM,            "ohm/m",                  "$\\mathsf{\\Omega}/$m",
  U_OHMMM1SM1,       U_INVALID,        "Ohm/m**2/s",             "$\\mathsf{\\Omega}$ m$^{-1}$ s$^{-1}$",
  U_OHMM1,           U_MOM1,           "1/ohm",                  "$\\mathsf{\\Omega}^{-1}$",
  U_OHMM1MM1,        U_MOM1,           "1/(ohm*m)",              "$\\mathsf{\\Omega}^{-1}$ m$^{-1}$",
  U_OHMM1MM2,        U_INVALID,        "1/(ohm*m**2)",           "$\\mathsf{\\Omega}^{-1}$ m$^{-2}$",
  U_OHMPV,           U_INVALID,        "ohm/volt",               "$\\mathsf{\\Omega}$ V$^{-1}$",
  U_OSQM1,           U_INVALID,        "ohm/square",             "$\\mathsf{\\Omega}/\\Box$",
  U_PERUNIT,         U_INVALID,        "per unit",               "per unit",
  U_RAD,             U_INVALID,        "rad",                    "rad",
  U_RADPS,           U_INVALID,        "rad/s",                  "rad/s",
  U_SECOND,          U_INVALID,        "s",                      "s",
  U_SECM1,           U_MSM1,           "1/s",                    "s$^{-1}$",
  U_SQUARES,         U_INVALID,        "squares",                "$\\Box$",
  U_VKM1,            U_INVALID,        "V/K",                    "V/K",
  U_VHM,             U_VHM2,           "V**(1/2)*m",             "V$^{1/2}$m",
  U_VHM2,            U_VHM3,           "V**(1/2)*m**2",          "V$^{1/2}$m$^{2}$",
  U_VHM3,            U_INVALID,        "V**(1/2)*m**3",          "V$^{1/2}$m$^{3}$",
  U_VM,              U_VM2,            "V*m",                    "Vm",
  U_VM2,             U_VM3,            "V*m**2",                 "Vm$^{2}$",
  U_VM3,             U_INVALID,        "V*m**3",                 "Vm$^{3}$",
  U_VMX,             U_VMX,            "V*m**X",                 "Vm$^{X}$",
  U_VMM1,            U_VOLT,           "V/m",                    "Vm$^{-1}$",
  U_VMM2,            U_VOLT,           "V/m**2",                 "Vm$^{-2}$",
  U_VMM3,            U_VOLT,           "V/m**3",                 "Vm$^{-3}$",
  U_VMMH,            U_INVALID,        "V/m**(1/2)",             "V/m$^{1/2}$",
  U_VOLT,            U_VM,             "V",                      "V",
  U_VOLT3,           U_INVALID,        "V**3",                   "V$^3$",
  U_VOLTH,           U_VHM,            "V**(1/2)",               "V$^{1/2}$",
  U_VOLTMH,          U_MVMH,           "V**(-1/2)",              "V$^{-1/2}$",
  U_VOLTM1,          U_MVM1,           "V**(-1)",                "V$^{-1}$",
  U_VOLTM2,          U_MVM2,           "V**(-2)",                "V$^{-2}$",
  U_VOLTSAMPM1METERM1, U_INVALID,      "Vs/(Am)",           "Vs/(Am)"
//  U_INVALID, U_INVALID, "", ""
};

size_t unitTableSize = sizeof(Units::unitTable)/sizeof(Units::unitTable[0]);

StdDescription descriptionTable[] =
{
  // General parameters
  "IC",        U_VOLT,    CAT_INITIAL, "Initial voltage drop across device",

  // Current parameters
  "DEV_I",     U_AMP,     CAT_CURRENT, "Device current",
  "DEV_IB",    U_AMP,     CAT_CURRENT, "Bulk current",
  "DEV_ID",    U_AMP,     CAT_CURRENT, "Drain current",
  "DEV_IG",    U_AMP,     CAT_CURRENT, "Gate current",
  "DEV_IS",    U_AMP,     CAT_CURRENT, "Source current",
  "DEV_IE",    U_AMP,     CAT_CURRENT, "External body current",

  // Temperature parameters
  "TEMP",      U_DEGC,    CAT_TEMP,    "Device temperature",
  "TNOM",      U_DEGC,    CAT_TEMP,    "Nominal device temperature",
  "TC1",       U_DEGCM1,  CAT_TEMP,    "Linear temperature coefficient",
  "TC2",       U_DEGCM2,  CAT_TEMP,    "Quadratic temperature coefficient",

  // Radiation parameters
  // JCV Note: 10/28/2013 This really belongs with the pulse parameters as far
  // as a user is concerned, so I changed the category from CAT_RAD to CAT_RADP.
  "GRORDER",   U_NONE,    CAT_RADP,     "Order of magnitude of radiation generation rate",

  // Pulsedata parameters
  "FUNCTIONTYPE", U_NONE, CAT_RADP,    "Internal function type for pulse (1-6) or user defined (7)",
  "PULSEDATA", U_NONE,    CAT_RADP,    "Expression (in braces \\{\\}) or filename (using \\{tablefile(\"filename\")\\} syntax) of (Time,Value) pairs to use instead of internal radiation source (functiontype = 7)",
  "PULSEV1",   U_NONE,    CAT_RADP,    "Initial value of radiation pulse (functiontype = 1)",
  "PULSEV2",   U_NONE,    CAT_RADP,    "Maximum value of radiation pulse (functiontype = 1)",
  "PULSETD",   U_SECOND,  CAT_RADP,    "Time delay of radiation pulse (functiontype = 1)",
  "PULSETR",   U_SECOND,  CAT_RADP,    "Rise time of radiation pulse (functiontype = 1)",
  "PULSETF",   U_SECOND,  CAT_RADP,    "Fall time of radiation pulse (functiontype = 1)",
  "PULSEPW",   U_SECOND,  CAT_RADP,    "Width of radiation pulse (functiontype = 1)",
  "PULSEPER",  U_SECOND,  CAT_RADP,    "Period of radiation pulse (functiontype = 1)",
  "SINVO",     U_NONE,    CAT_RADP,    "Offset of sinusoidal radiation source (functiontype = 2)",
  "SINVA",     U_NONE,    CAT_RADP,    "Amplitude of sinusoidal radiation source (functiontype = 2)",
  "SINFREQ",   U_HZ,      CAT_RADP,    "Frequency of sinusoidal radiation source (functiontype = 2)",
  "SINTD",     U_SECOND,  CAT_RADP,    "Delay of sinusoidal radiation source (functiontype = 2)",
  "SINTHETA",  U_SECM1,   CAT_RADP,    "Inverse of time constant for exponential envelope of sinusoidal"
                                       " radiation source (functiontype = 2)",
  "EXPV1",     U_NONE,    CAT_RADP,    "Initial value of exponential cusp radiation source (functiontype = 3)",
  "EXPV2",     U_NONE,    CAT_RADP,    "Height of exponential cusp radiation source (functiontype = 3)",
  "EXPTD1",    U_SECOND,  CAT_RADP,    "Rise delay time of exponential cusp radiation source (functiontype = 3)",
  "EXPTAU1",   U_SECOND,  CAT_RADP,    "Rise time constant of exponential cusp radiation source (functiontype = 3)",
  "EXPTD2",    U_SECOND,  CAT_RADP,    "Fall delay time of exponential cusp radiation source (functiontype = 3)",
  "EXPTAU2",   U_SECOND,  CAT_RADP,    "Fall time constant of exponential cusp radiation source (functiontype = 3)",
  "SFFMVO",    U_NONE,    CAT_RADP,    "Offset of single frequency FM radiation source (functiontype = 4)",
  "SFFMVA",    U_NONE,    CAT_RADP,    "Amplitude of single frequency FM radiation source (functiontype = 4)",
  "SFFMFC",    U_HZ,      CAT_RADP,    "Carrier frequency of single frequency FM radiation source (functiontype = 4)",
  "SFFMMDI",   U_NONE,    CAT_RADP,    "Modulation index of single frequency FM radiation source (functiontype = 4)",
  "SFFMFS",    U_HZ,      CAT_RADP,    "Signal frequency of single frequency FM radiation source (functiontype = 4)",
  "PWLT1",     U_SECOND,  CAT_RADP,    "Time of point 1 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLT2",     U_SECOND,  CAT_RADP,    "Time of point 2 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLT3",     U_SECOND,  CAT_RADP,    "Time of point 3 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLT4",     U_SECOND,  CAT_RADP,    "Time of point 4 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLT5",     U_SECOND,  CAT_RADP,    "Time of point 5 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLV1",     U_NONE,    CAT_RADP,    "Value of point 1 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLV2",     U_NONE,    CAT_RADP,    "Value of point 2 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLV3",     U_NONE,    CAT_RADP,    "Value of point 3 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLV4",     U_NONE,    CAT_RADP,    "Value of point 4 of piecewise linear curve radiation source (functiontype = 5)",
  "PWLV5",     U_NONE,    CAT_RADP,    "Value of point 5 of piecewise linear curve radiation source (functiontype = 5)",
  "LINET1",    U_SECOND,  CAT_RADP,    "Time of point 1 of linear radiation source (functiontype = 6)",
  "LINEV1",    U_NONE,    CAT_RADP,    "Value of point 1 of linear radiation source (functiontype = 6)",
  "LINEK",     U_SECM1,   CAT_RADP,    "Slope of point 1 of linear radiation source (functiontype = 6)"
//  NULL,        U_NONE,    CAT_NONE,    ""
};

size_t descriptionTableSize = sizeof(Units::descriptionTable)/sizeof(Units::descriptionTable[0]);

} // namespace Units
} // namespace Device
} // namespace Xyce
