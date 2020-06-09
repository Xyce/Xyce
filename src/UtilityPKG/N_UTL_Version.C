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
// Filename       : N_UTL_Version.C
//
// Purpose        : set version string
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  :
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>
#include <timestamp.h>

#include <N_UTL_FeatureTest.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : Version::getFullVersionString
// Purpose       : get full banner string for Xyce version
// Special Notes : version should be properly formatted in configure.ac
//               : -Dmacros used to get version and timestamp inserted
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
std::string Version::getFullVersionString()
{
  const std::string tmpVer( VERSION );
  std::string version;

  // create developement version string
  if( tmpVer[ 0 ] == 'D' || tmpVer[ 0 ] == 'd' )
  {
    version += "DEVELOPMENT-";

    // add the timestamp

    std::ostringstream ver("");

    ver << XYCEBUILDTIMESTAMP;

    version += std::string( ver.str() );
  }

  // create release version string
  else if( tmpVer[ 0 ] == 'R' || tmpVer[ 0 ] == 'r' )
  {
    // prepend release phase if necessary
    if( tmpVer[ 2 ] != ':' )
    {
      version += "(" + tmpVer.substr( 2, 1 ) + ")";
    }

    // add the release major-minor-rev number
    int i = tmpVer.find_last_of( ":" );
    version += "Release " + tmpVer.substr( i + 1, tmpVer.length() - i );
  }

  version = getXyceName()+" "+version+getBuildVariant();

  return version;
}


//-----------------------------------------------------------------------------
// Function      : Version::getShortVersionString
// Purpose       : get the maj-min-rev number for Xyce version
// Special Notes : version should be properly formatted in configure.ac
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
std::string Version::getShortVersionString()
{
  const std::string tmpVer( VERSION );

  // get position of the major-minor-rev number
  int i = tmpVer.find_last_of( ":" );

  return tmpVer.substr( i + 1, tmpVer.length() - i );
}


//-----------------------------------------------------------------------------
// Function      : Version::getBuildVariant
// Purpose       : Return a short string indicating special builds
// Special Notes : Current special builds are:
//               :  (blank) = full ECI build
//               :  OS      = open source build
//               :  norad   = no rad models, but not open source
// Scope         :
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string Version::getBuildVariant()
{
  std::string variant;
#ifdef Xyce_RAD_MODELS
#ifndef Xyce_ATHENA
  variant="-noATHENA";
#else
  variant="";
#endif // Xyce_ATHENA
#ifndef Xyce_NONFREE_MODELS
  variant+="-nononfree";
#endif // Xyce_NONFREE_MODELS
#else // Xyce_RAD_MODELS
#ifdef Xyce_NONFREE_MODELS
  variant="";
#else
  variant="-opensource";
#endif // Xyce_NONFREE_MODELS
#endif // Xyce_RAD_MODELS

#ifdef Xyce_MODSPEC_MODELS
  variant+="-modspec";
#endif // Xyce_MODSPEC_MODELS

#ifdef Xyce_Dakota
  variant+="-dakota";
#endif // Xyce_Dakota

  return variant;
}



//-----------------------------------------------------------------------------
// Function      : Version::getCapabilities
// Purpose       : Return a string indicating compiled-in features
// Special Notes :
// Scope         :
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string Version::getCapabilities()
{
  std::string capabilities="";

#ifdef Xyce_PARALLEL_MPI
  capabilities += "Parallel with MPI\n";
#else
  capabilities += "Serial\n";
#endif

#ifdef Xyce_RAD_MODELS
  capabilities += "Radiation models\n";
#endif

#ifdef Xyce_NONFREE_MODELS
  capabilities += "Non-GPL device models\n";
#endif

#ifdef Xyce_MODSPEC_MODELS
  capabilities += "ModSpec enabled\n";
#endif

#ifdef Xyce_USE_FFT
  capabilities += "FFT";
#ifdef Xyce_USE_INTEL_FFT
  capabilities += "(Intel FFT)\n";
#else
  capabilities += "(FFTW)\n";
#endif
#endif

#ifdef Xyce_USE_HDF5
  capabilities += "HDF5\n";
#endif

#ifdef Xyce_USE_CURL
  capabilities += "Metrics reporting via cURL\n";
#ifdef Xyce_TRACKING_URL
  capabilities += "Reporting to URL " + std::string(Xyce_TRACKING_URL) + "\n";
#endif
#endif


#ifdef Xyce_REACTION_PARSER
  capabilities += "Reaction parser\n";
#endif

#ifdef Xyce_ATHENA
  capabilities += "ATHENA\n";
#endif

#ifdef Xyce_Dakota
  capabilities += "Dakota direct linkage\n";
#endif

#if __cplusplus>=201103L
  capabilities += "Build compiler is C++11 compliant\n";
#endif

#if __cplusplus>=201402L
  capabilities += "Build compiler is C++14 compliant\n";
#endif

#if Xyce_STOKHOS_ENABLE
  capabilities += "Stokhos enabled\n";
#endif

#ifdef Xyce_AMESOS2
  capabilities += "Amesos2 (Basker and Klu2) enabled\n";
#endif

  if (VERBOSE_TIME)
    capabilities += "Verbose output - time integrator\n";

  if (VERBOSE_LINEAR)
    capabilities += "Verbose output - linear solver\n";

  if (VERBOSE_NONLINEAR)
    capabilities += "Verbose output - nonlinear solver\n";

  return capabilities;
}


//-----------------------------------------------------------------------------
// Function      : Version::getLicense
// Purpose       : Return a string indicating license
// Special Notes :
// Scope         :
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string Version::getLicense()
{
  std::string License="";

#ifdef Xyce_RAD_MODELS
  License += "\n EXPORT CONTROLLED SOFTWARE\n";
  License += "\n";
  License += " Copyright 2002 National Technology & Engineering Solutions of Sandia,\n";
  License += " LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS,\n";
  License += " there is a non-exclusive license for use of this work by or on behalf\n";
  License += " of the U.S. Government.  Export of this data may require a license\n";
  License += " from the United States Government.\n";
#else

#ifdef Xyce_NONFREE_MODELS
  License += "\n THIS IS NOT AN OPEN SOURCE BINARY\n";
  License += "\n";
  License += " The EKV3 model is not an open source model.  It was developed by the\n";
  License += " EKV Team of the Electronics Laboratory-TUC of the Technical University\n";
  License += " of Crete.  TUC has granted Sandia National Laboratories a\n";
  License += " non-exclusive, royalty-free license to use, copy, modify, implement and\n";
  License += " distribute the EKV3 Model Code in the form of compiled, executable\n";
  License += " code.  Sandia National Laboratories is not licensed to distribute this\n";
  License += " code in source form.  Documentation of the EKV3 model is available on\n";
  License += " the Xyce internal web site:\n";
  License += " http://info.sandia.gov/xyce/\n";
  License += "\n";
  License += " More information about the EKV model (including contact information\n";
  License += " for the EKV Model Team) is available at the EKV web site:\n";
  License += " http://ekv.epfl.ch/.\n\n";
  License += "\n";
  License += " All components of this version of Xyce OTHER THAN the EKV3 model are\n";
  License += " released under the GNU Public License.\n";
#endif

  License +="\n Xyce(TM) Parallel Electrical Simulator\n";
  License +=" Copyright 2002 National Technology & Engineering Solutions of Sandia,\n";
  License +=" LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the\n";
  License +=" U.S. Government retains certain rights in this software.\n";
  License +="\n";
  License +=" This program is free software: you can redistribute it and/or modify\n";
  License +=" it under the terms of the GNU General Public License as published by\n";
  License +=" the Free Software Foundation, either version 3 of the License, or\n";
  License +=" (at your option) any later version.\n";
  License +="\n";
  License +=" This program is distributed in the hope that it will be useful,\n";
  License +=" but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  License +=" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
  License +=" GNU General Public License for more details.\n";
  License +="\n";
  License +=" You should have received a copy of the GNU General Public License\n";
  License +=" along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
  License +="\n";
  License +="\n";

#ifdef Xyce_NONFREE_MODELS
  License +="To obtain source code for components of Xyce EXCLUDING the EKV model,\n";
  License +="see the Xyce web site, http://xyce.sandia.gov/.\n";
#endif

#endif

  return License;
}


//-----------------------------------------------------------------------------
// Function      : Version::getXyceName
// Purpose       : Return a string with the copyrighted name of this build
// Special Notes : "XyceRad" means it is the OUO/ECI version, "Xyce" for others
// Scope         :
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string Version::getXyceName()
{
  std::string XyceName="";
#if defined(Xyce_RAD_MODELS) || defined(Xyce_ATHENA)
  XyceName="XyceRad";
#else
  XyceName="Xyce";
#endif
  return XyceName;
}

} // namespace Util
} // namespace Xyce
