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
  capabilities += "Non-Free device models\n";
#endif

#ifdef Xyce_MODSPEC_MODELS
  capabilities += "ModSpec enabled\n";
#endif

#ifdef Xyce_ADMS_MODELS
  capabilities += "Verilog-A derived (ADMS) models\n";
  #ifdef Xyce_ADMS_SENSITIVITIES
    capabilities += "Analytic sensitivities in ADMS models\n";
  #endif
#endif


#ifdef Xyce_USE_FFT
  capabilities += "FFT ";
#ifdef Xyce_USE_INTEL_FFT
  capabilities += "(Intel FFT)\n";
#else
  capabilities += "(FFTW)\n";
#endif
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

#ifdef Xyce_ROL
  capabilities += "ROL enabled\n";
#endif

#if __cplusplus>=201402L
  capabilities += "Build compiler is C++14 compliant\n";
#endif

#ifdef Xyce_STOKHOS_ENABLE
  capabilities += "Stokhos enabled\n";
#endif

#ifdef Xyce_AMESOS2
  #if defined(Xyce_AMESOS2_KLU2) || defined(Xyce_AMESOS2_BASKER)
    capabilities += "Amesos2 (Basker and KLU2) enabled\n";
  #elif defined(Xyce_AMESOS2_KLU2) 
    capabilities += "Amesos2 (KLU2) enabled\n";
  #elif defined(Xyce_AMESOS2_BASKER)
    capabilities += "Amesos2 (Basker) enabled\n";
  #endif
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
  License += " Copyright 2002-2023 National Technology & Engineering Solutions of\n";
  License += " Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n";
  License += " NTESS, there is a non-exclusive license for use of this work by or\n";
  License += " on behalf of the U.S. Government.  Export of this data may require\n";
  License += " a license from the United States Government.\n";

#elif defined Xyce_NONFREE_MODELS
  License += "\n Copyright 2020-2023 National Technology & Engineering Solutions of\n";
  License += " Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n";
  License += " NTESS, there is a non-exclusive license for use of this work by or\n";
  License += " on behalf of the U.S. Government.  Export of this data may require\n";
  License += " a license from the United States Government.\n";
  License += "\n";
  License += " NOTICE:\n";
  License += " For five (5) years from 6/17/2020, the United States Government is\n";
  License += " granted for itself and others acting on its behalf a paid-up,\n";
  License += " nonexclusive, irrevocable worldwide license in this data to reproduce,\n";
  License += " prepare derivative works, and perform publicly and display publicly,\n";
  License += " by or on behalf of the Government.  There is provision for the\n";
  License += " possible extension of the term of this license.  Subsequent to that\n";
  License += " period or any extension granted, the United States Government is\n";
  License += " granted for itself and others acting on its behalf a paid-up,\n";
  License += " nonexclusive, irrevocable worldwide license in this data to reproduce,\n";
  License += " prepare derivative works, distribute copies to the public, perform\n";
  License += " publicly and display publicly, and to permit others to do so.  The\n";
  License += " specific term of the license can be identified by inquiry made to\n";
  License += " National Technology and Engineering Solutions of Sandia, LLC or DOE.\n";
  License += "\n";
  License += " NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT\n";
  License += " OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF\n";
  License += " SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS\n";
  License += " OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY,\n";
  License += " COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR\n";
  License += " PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE\n";
  License += " PRIVATELY OWNED RIGHTS.\n";
  License += "\n";
  License += " Any licensee of \"XyceNF\" has the obligation and responsibility to abide\n";
  License += " by the applicable export control laws, regulations, and general\n";
  License += " prohibitions relating to the export of technical data.  Failure to\n";
  License += " obtain an export control license or other authority from the\n";
  License += " Government may result in criminal liability under U.S. laws.\n";
#else

  License +="\n Xyce(TM) Parallel Electrical Simulator\n";
  License +=" Copyright 2002-2023 National Technology & Engineering Solutions of\n";
  License +=" Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with\n";
  License +=" NTESS, the U.S. Government retains certain rights in this software.\n";
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
#endif

  return License;
}


//-----------------------------------------------------------------------------
// Function      : Version::getXyceName
// Purpose       : Return a string with the copyrighted name of this build
// Special Notes : "XyceRad" means it is the OUO/ECI variant, "XyceNF" if it
//                 includes the "non-free" models (but not the ECI models),
//                 and "Xyce" for the open source variant
// Scope         :
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string Version::getXyceName()
{
  std::string XyceName="";
#if defined(Xyce_RAD_MODELS) || defined(Xyce_ATHENA)
  XyceName="XyceRad";
#elif defined Xyce_NONFREE_MODELS
  XyceName="XyceNF";
#else
  XyceName="Xyce";
#endif
  return XyceName;
}

} // namespace Util
} // namespace Xyce
