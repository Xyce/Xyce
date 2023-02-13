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
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Util {

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class. If set to true then 1A=1e-18 rather than 1.
bool useHspiceUnits = false;

//-----------------------------------------------------------------------------
// Description   : Adjustment of data values when recognizable scale factor is
//                 present (e.g. 2.5n ==> 2.5E-9).  This function is used with
//                 the scaling factors defined for Xyce.
// Special Notes :
// Creator       : Alan Lundin
// Creation Date :
//-----------------------------------------------------------------------------
double Value(const std::string &tmpStr)
{
   char tmp[3];
   double value = atof(tmpStr.c_str());
   int j = tmpStr.find_first_not_of("0123456789.-+eE", 0);
   if (j == tmpStr.npos) return value;
   switch (tmpStr[j])
   {
      case 'T' : case 't' :
         return value*1.0e12;
      case 'G' : case 'g' :
         return value*1.0e9;
      case 'K' : case 'k' :
         return value*1.0e3;
      case 'M' :
      case 'm' :
         tmp[0] = tolower (tmpStr[j]);
         if (tmpStr.size() > j+2)
         {
           tmp[1] = tolower (tmpStr[j+1]);
           tmp[2] = tolower (tmpStr[j+2]);
           if (tmp[1]  == 'i' && tmp[2] == 'l')
             return value*25.4e-6;
           else if (tmp[1] == 'e' && tmp[2] == 'g')
             return value*1.0e6;
         }
         return value*1.0e-3;
      case 'x' : case 'X' :
         return value*1.0e6;
      case 'u' : case 'U' :
         return value*1.0e-6;
      case 'n' : case 'N' :
         return value*1.0e-9;
      case 'p' : case 'P' :
         return value*1.0e-12;
      case 'f' : case 'F' :
         return value*1.0e-15;
      case 'a' : case 'A' :
	 if (Util::useHspiceUnits) {return value*1.0e-18;}
      default :
         return value;
   }
   return 0;
}

//-----------------------------------------------------------------------------
// Description   : Adjustment of data values when a recognizable scale factor is
//                 present (e.g. 2.5n ==> 2.5E-9).  This version uses the 
//                 scaling factors defined in the IBIS standard which are 
//                 different from those defined in Xyce
// Special Notes : The IBIS scaling factors are CASE-SENSITVE.  They are
//                 case-insensitive in Xyce.
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/26/2018
//-----------------------------------------------------------------------------
double IBISValue(const std::string &tmpStr)
{
   char tmp[3];
   double value = atof(tmpStr.c_str());
   int j = tmpStr.find_first_not_of("0123456789.-+eE", 0);
   if (j == tmpStr.npos) return value;

   // the characters in the case statements below are CASE-SENSITIVE
   switch (tmpStr[j])
   {
      case 'T' : 
         return value*1.0e12;
      case 'G' :
         return value*1.0e9;
      case 'k' :
         return value*1.0e3;
      case 'M' :
         return value*1.0e6;
      case 'm' :
         return value*1.0e-3;
      case 'u' :
         return value*1.0e-6;
      case 'n' :
         return value*1.0e-9;
      case 'p' :
         return value*1.0e-12;
      case 'f' :
         return value*1.0e-15;
      default :
         return value;
   }
   return 0;
}


bool isValue(const std::string & tmpStr)
{
  int stringPos = 0;
  int i = 0;
  int stringSize = tmpStr.size();
  
  static const char *units[] = {"C", "V", "VOLT", "VDC", "A", "AMP", "AMPERE", "F", 
                          "FARAD", "HENRY", "HY", "IL", "EG", "H", "HZ", 
                          "HERTZ", "OHM", "SECOND", "S", "METER", "M",
                          "MEG", "MIL", NULL};

  // Check for leading + or - sign.
  char ch ( tmpStr[stringPos] );
  if ( ch == '+' || ch == '-' )
    ++stringPos;

  if ( stringPos == stringSize )
    return false; // string ended too soon to be a numeric value.

  ch = tmpStr[stringPos];
  if ( (!isdigit(ch)) && ch != '.' )
    return false;

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  if ( ch == '.' )
  {
    // Must have a digit before or after the decimal.
    if ( stringPos + 1 >= stringSize )
    {
      // Decimal is at end of string, must have digit before decimal.
      if ( stringPos == 0 ) 
        return false;
      else if ( !isdigit(tmpStr[stringPos-1]) )
        return false;
    }
    else if ( stringPos == 0 )
    {
      // Decimal is at beginning of string, must have digit after decimal.
      if ( stringPos + 1 >= stringSize )
        return false;
      else if ( !isdigit(tmpStr[stringPos+1]) )
        return false;
    }
    else if ( !isdigit(tmpStr[stringPos-1]) &&
              !isdigit(tmpStr[stringPos+1]) )
      return false;

    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  while (isdigit(ch))
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];
  }

  // Check for exponent.
  if ( ch == 'E' || ch == 'e' )
  {
    ++stringPos;
    if ( stringPos == stringSize ) return true; // all is well.
    ch = tmpStr[stringPos];

    // Look for exponent sign.
    if ( ch == '+' || ch == '-' )
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = tmpStr[stringPos];
    }

    while (isdigit(ch))
    {
      ++stringPos;
      if ( stringPos == stringSize ) return true; // all is well.
      ch = tmpStr[stringPos];
    }
  }

  if (tmpStr[stringPos] == 'T' || tmpStr[stringPos] == 't' ||
      tmpStr[stringPos] == 'G' || tmpStr[stringPos] == 'g' ||
      tmpStr[stringPos] == 'X' || tmpStr[stringPos] == 'x' ||
      tmpStr[stringPos] == 'K' || tmpStr[stringPos] == 'k' ||
      tmpStr[stringPos] == 'U' || tmpStr[stringPos] == 'u' ||
      tmpStr[stringPos] == 'N' || tmpStr[stringPos] == 'n' ||
      tmpStr[stringPos] == 'P' || tmpStr[stringPos] == 'p' ||
      tmpStr[stringPos] == 'F' || tmpStr[stringPos] == 'f' ||
      (Util::useHspiceUnits && (tmpStr[stringPos] == 'A' || tmpStr[stringPos] == 'a')) )
    ++stringPos; 

  if (tmpStr[stringPos] == 'M' || tmpStr[stringPos] == 'm') 
  {
     char tmp[3];
     tmp[0] = tolower (tmpStr[stringPos]);
     ++stringPos;

    if (stringSize >= stringPos +2)
    {
           tmp[1] = tolower (tmpStr[stringPos]);
           tmp[2] = tolower (tmpStr[stringPos + 1]);
           if ((tmp[1]  == 'i' && tmp[2] == 'l') || (tmp[1] == 'e' && tmp[2] == 'g'))
             stringPos = stringPos + 2;
    }

  }

  if (stringPos == stringSize)
    return true;

  ExtendedString u(tmpStr.substr(stringPos, stringSize-stringPos));
  u = tmpStr.substr(stringPos, stringSize-stringPos);
  u.toUpper();
  i = 0;
  while (units[i] != NULL)
  {
    if (u == units[i++])
      return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : isInt
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool isInt(const std::string & tmpStr)
{
  int j;

  if (tmpStr.empty())
    return false;

  if (tmpStr[0] == '-' || tmpStr[0] == '+')
    j = tmpStr.find_first_not_of("0123456789", 1);
  else
    j = tmpStr.find_first_not_of("0123456789");

  if (j == (int)std::string::npos)
    return true;

  // But there's one case where we *could* still be an int.  That would be
  // if we've got .[0]* at the current point.

  if (tmpStr[j] == '.')
  {
    std::string::size_type i = tmpStr.find_first_not_of("0",j+1);

    // If we find nothing but 0 after the ., we are still an int.
    if (i == std::string::npos)
     return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : isTableFileKeyword
// Purpose       : Test is string is a tablefile keyword 
// Special Notes :
// Scope         : public
// Creator       : R. Schiek
// Creation Date : 09/25/18
//-----------------------------------------------------------------------------
bool isTableFileKeyword(const std::string & tmpStr)
{
  size_t matchingPosUC = tmpStr.find( "TABLEFILE" );
  size_t matchingPosLC = tmpStr.find( "tablefile" );
  if (((matchingPosUC == std::string::npos) && (matchingPosLC == 0)) || 
      ((matchingPosUC == 0) && (matchingPosLC == std::string::npos)) )
  {
    return true;
  }
  return false;
}

} // namespace Util
} // namespace Xyce
