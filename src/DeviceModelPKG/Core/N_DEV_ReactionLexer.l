/*-------------------------------------------------------------------------
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
*/


%{
#include <Xyce_config.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <map>
// these here because Bison 2.4 no longer puts the pre-prologue into the
// header file... grrrr.
#include <vector>
#include <N_DEV_Specie.h>
#include <N_DEV_ReactionNetwork.h>
#include <N_UTL_Math.h>
namespace Xyce {
namespace Device {
class ReactionLexer;
}}
// now we can include this without errors from bison 2.4
#include "N_DEV_ReactionParser.hxx"
// BLEAH!   This is here DUPLICATED from ReactionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
/* CAREFUL watch continuations! */
#undef YY_DECL
#define YY_DECL \
  int Xyce::Device::ReactionLexer::getToken(XyceDevice::ReactionParser::semantic_type *lvalp, \
                                            XyceDevice::location *llocp, \
                                            std::map<std::string, int> &theSpeciesIDs)

  // YECH!  Work around very stupid way that multiple parsers/lexers are 
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#include <FlexLexer.h>
#include <N_DEV_ReactionLexer.h>
  // undo that kludge.  Note that because of this stupidity, if the 
  // "%name-prefix" is changed, this line needs to be changed, too.
#define yylex XyceDevicelex
//#include <ReactionLexer.h>
#include "location.hh"
#include "position.hh"

%}
%option noyywrap
%option yyclass="Xyce::Device::ReactionLexer"
%s speciedef
%s reactdef
%s sourcedef
%s ics_def
%s enhancementdef

DIGIT [0-9]
ID [a-zA-Z][a-zA-Z0-9_]*

%{
#define YY_USER_ACTION llocp->columns(YYLeng());
%}

%%
%{
  llocp->step();
%}

[\-]?{DIGIT}+ {
        lvalp->ival=strtol(YYText(),0,10);
        return XyceDevice::ReactionParser::token::RXN_INT;
        }

[\-]?{DIGIT}*"."{DIGIT}* {
        lvalp->dval=strtod(YYText(),0);
        return XyceDevice::ReactionParser::token::RXN_FLOAT;
        }

[\-]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3} {
        lvalp->dval=strtod(YYText(), 0);
        return XyceDevice::ReactionParser::token::RXN_FLOAT;
       }

<INITIAL>species {
  BEGIN(speciedef);
  return (XyceDevice::ReactionParser::token::SPECIES_BEGIN);
  }

<INITIAL>reactions {
  BEGIN(reactdef);
  return (XyceDevice::ReactionParser::token::REACTIONS_BEGIN);
  }

<INITIAL>enhancements {
  BEGIN(enhancementdef);
  return (XyceDevice::ReactionParser::token::ENHANCEMENTS_BEGIN);
  }

<INITIAL>sources {
  BEGIN(sourcedef);
  return (XyceDevice::ReactionParser::token::SOURCES_BEGIN);
  }

<INITIAL>initial_conditions {
  BEGIN(ics_def);
  return (XyceDevice::ReactionParser::token::ICS_BEGIN);
  }

end {
  BEGIN(INITIAL);
  return (XyceDevice::ReactionParser::token::KEY_END);
  }

electron_capture  { return(XyceDevice::ReactionParser::token::KEY_ECAPTURE);}
electron_emission { return(XyceDevice::ReactionParser::token::KEY_EEMISSION);}
fdelectron_emission { return(XyceDevice::ReactionParser::token::KEY_FDEEMISSION);}
hole_capture  { return(XyceDevice::ReactionParser::token::KEY_HCAPTURE);}
hole_emission { return(XyceDevice::ReactionParser::token::KEY_HEMISSION);}
fdhole_emission { return(XyceDevice::ReactionParser::token::KEY_FDHEMISSION);}
simple   { return(XyceDevice::ReactionParser::token::KEY_SIMPLE);}
complex  { return(XyceDevice::ReactionParser::token::KEY_COMPLEX);}
complexmultiplier  { return(XyceDevice::ReactionParser::token::KEY_COMPLEXMULTIPLIER);}
decomplex  { return(XyceDevice::ReactionParser::token::KEY_DECOMPLEX);}
capture  { return(XyceDevice::ReactionParser::token::KEY_CAPTURE);}
emission { return(XyceDevice::ReactionParser::token::KEY_EMISSION);}
bourgoin_corbett_hole {return(XyceDevice::ReactionParser::token::KEY_BCHOLE);}
bourgoin_corbett_electron {return(XyceDevice::ReactionParser::token::KEY_BCELECTRON);}
bourgoin_corbett {return(XyceDevice::ReactionParser::token::KEY_BC);}
const { return(XyceDevice::ReactionParser::token::KEY_CONST);}


<speciedef>{ID}   {
     lvalp->sval=new std::string(YYText());
    return XyceDevice::ReactionParser::token::SPECIE;
     }

<reactdef,sourcedef,ics_def,enhancementdef>{ID}   {
  lvalp->sval=new std::string(YYText());
     if (theSpeciesIDs.find(YYText()) != theSpeciesIDs.end())
     { 
       return XyceDevice::ReactionParser::token::SPECIE;
     }
     else
     {
       return XyceDevice::ReactionParser::token::IDENTIFIER;
     }
}

<sourcedef>"{"[^}]*"}" {
  lvalp->sval=new std::string(YYText()); 
  (lvalp->sval)->erase(0,1);// lop off open brace
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop of close brase
  return XyceDevice::ReactionParser::token::XYCE_EXPRESSION;
}

[+;,()]   {
     return((YYText())[0]);
    }


"->" {
     return(XyceDevice::ReactionParser::token::ARROW);
    }
 

"#"[^\n\r]*     /* eat up one-line comments */

[ \t]+        llocp->step();  /* eat up whitespace */
(\n+|(\r\n)+|\r+)         llocp->lines(YYLeng());

.           std::cout << "Unrecognized character: " <<  YYText() << std::endl;

%%
