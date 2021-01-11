/*-------------------------------------------------------------------------
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
*/

%{

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <complex>

#include <ExpressionType.h>
#include <ast.h>
#include <newExpression.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Util {

template <typename ScalarT>
class tableArgs
{
  public:
    std::string * keywordPtr;
    Teuchos::RCP<astNode<ScalarT> > node;
};

class ExpressionLexer;
}}
#include "ExpressionParser.hxx"

#undef YY_DECL
#define YY_DECL \
  int Xyce::Util::ExpressionLexer::getToken(XyceExpression::ExpressionParser::semantic_type *lvalp, \
                                            XyceExpression::location *llocp)

#undef yylex
#include "ExpressionLexer.h"
#define yylex XyceExpressionlex


#include "location.hh"

%}
%option noyywrap
%option yyclass="Xyce::Util::ExpressionLexer"
%s voltagenodedef
%s currentdevdef
%s internaldevdef
%s noisedevdef
   /* %s polydef */
%s tablefilename
%s polyargs
%s polyvoltagenodedef
%s polycurrentdevdef
%s sparamdef
%s yparamdef
%s zparamdef



 /* 
    The bison file is officially case-insensitive.  So converting the
    various strings to upper case is not necessary for it to work by 
    itself.

    However, Xyce mostly deals with netlist strings by converting to upper case.
    So, when the Xyce source code outside of expression library interacts with it, 
    Bison will not have converted the expression string or any subordinate substrings
    to upper or lower case on its own.  They will be the original case, 
    whatever that was in the netlist. 

    I previously set up the expression library so that it would up-case the entire 
    expression string prior to parsing.  However, this doesn't work now because 
    the tableOp class now supports reading data from a user-specified file, 
    and this needs to work on systems that are case sensitive like Linux.  So,
    up-casing the whole expression string won't work anymore.

    Instead, in this file, certain requisite tokens are up-cased via "toUpper".  
    These include voltage node names, parameter names and function names.  
    Most other stuff (including filenames especially) is left in its original 
    case.
 */

 /* ID is used for a variety of random strings, including function and 
    parameter names.  So, it cannot have characters in it that could 
    also be interpretted as operators, such as '!', '+' and '-'. This 
    distinguishes ID from DEV, below. */

 /* DEV is used for device instance names, particularly inside of current vars.  I(DEV) */
 /* NODE is used for voltage node names, particularly inside of voltage vars  V(NODE) and V(NODE,NODE) */

 /* Both DEV and NODE include the colon, to handle Xyce subcircuit syntax */
 /* Other simulators use other delimiters for subcircuit names, so might */
 /* want to consider supporting delimiters other than the colon. */
DIGIT [0-9]
ID [a-zA-Z_`@#\$][a-zA-Z0-9_:`@#\.\$]*

 /* this is to support table files */
FILENAME [a-zA-Z0-9_:`@#\.\$\/]*

 /* to support bug 1034, the DEV needs to recognize some weird characters, including '+' and '-' */
DEV [a-zA-Z][a-zA-Z0-9_:!`@#\.\+\-]*

 /* to support bug 1034, the NODE needs to recognize a lot of weird characters: */
 /* from bug 1034: */
 /* ` ~ ! @ # $ % ^ & - _ + [ ] | \ < > . ? / */
 /* unfortunately, some of these are special regex characters, so they need to be */
 /* handled with a preceding backslash.  Probably more elegant regex, such as a 
    range of ASCI characters could be used for NODE, but this seems to work for now */
NODE [a-zA-Z0-9_:\$\-`~!@#%&_\+|<>\?\.\\/\^\*\[\]]* 

VOP V[RIMP]?
IOP I[RIMP]?
NOP N[RIMP]?

  /* RF parameters */
SOP S[RIMP]?
YOP Y[RIMP]?
ZOP Z[RIMP]?

EXPONENT [eE][\-\+]?{DIGIT}{1,3}

SUFFIX [tTgGkKmMxXuUnNpPfF]

MEG [mM][eE][gG]
MIL [mM][iI][lL]

UNIT [sSvVaAhHfF]
HZ   [Hh][Zz]

TABLE "fast"?"table"

LEAD I[SDGBEC1-9]
%{
#define YY_USER_ACTION llocp->columns(YYLeng());
%}

%%
%{
  llocp->step();
%}

 /* ------------------------------------------------------------------------*/
 /* real numbers with suffixes */
<INITIAL>{DIGIT}+"."?{SUFFIX}?{UNIT}?     |
<INITIAL>{DIGIT}*"."{DIGIT}+{SUFFIX}?{UNIT}?   |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{UNIT}? |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{UNIT}? |
 /* real numbers with 3-character suffixes; meg */
<INITIAL>{DIGIT}+"."?{MEG}{UNIT}? |
<INITIAL>{DIGIT}*"."{DIGIT}+{MEG}{UNIT}? |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MEG}{UNIT}? |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MEG}{UNIT}? |
 /* real numbers with 3-character suffixes; mil */
<INITIAL>{DIGIT}+"."?{MIL}{UNIT}? |
<INITIAL>{DIGIT}*"."{DIGIT}+{MIL}{UNIT}? |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MIL}{UNIT}? |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MIL}{UNIT}? |
 /* real numbers dealing with Hz */
<INITIAL>{DIGIT}+"."?{SUFFIX}?{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}+{SUFFIX}?{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{HZ} |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{HZ} |
 /* real numbers with 3-character suffixes; meg */
<INITIAL>{DIGIT}+"."?{MEG}{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}*{MEG}{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MEG}{HZ} |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MEG}{HZ} |
 /* real numbers with 3-character suffixes; mil */
<INITIAL>{DIGIT}+"."?{MIL}{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}*{MIL}{HZ} |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MIL}{HZ} |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{1,3}{MIL}{HZ} {
     lvalp->build<double>()=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }


 /* ------------------------------------------------------------------------*/
 /* imag numbers with suffixes */
<INITIAL>{DIGIT}+"."?{SUFFIX}?{UNIT}?[ \t]*"J"     |
<INITIAL>{DIGIT}*"."{DIGIT}+{SUFFIX}?{UNIT}?[ \t]*"J"   |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{UNIT}?[ \t]*"J" |
 /* imag numbers with 3-character suffixes; meg */
<INITIAL>{DIGIT}+"."?{MEG}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+{MEG}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MEG}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MEG}{UNIT}?[ \t]*"J" |
 /* imag numbers with 3-character suffixes; mil */
<INITIAL>{DIGIT}+"."?{MIL}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+{MIL}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MIL}{UNIT}?[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MIL}{UNIT}?[ \t]*"J" |
 /* imag numbers dealing with Hz */
<INITIAL>{DIGIT}+"."?{SUFFIX}?{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+{SUFFIX}?{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{HZ}[ \t]*"J" |
 /* imag numbers with 3-character suffixes; meg */
<INITIAL>{DIGIT}+"."?{MEG}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}*{MEG}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MEG}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{MEG}{HZ}[ \t]*"J" |
 /* imag numbers with 3-character suffixes; mil */
<INITIAL>{DIGIT}+"."?{MIL}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}*{MIL}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}*"."{DIGIT}+[eE][\-\+]?{DIGIT}{1,3}{MIL}{HZ}[ \t]*"J" |
<INITIAL>{DIGIT}+"."?[eE][\-\+]?{DIGIT}{1,3}{1,3}{MIL}{HZ}[ \t]*"J" {
     lvalp->build<double>()=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_IMAG_NUM;
  }

 /* ----------------------------------------------------------*/
 /* single character ops */

  /* In newer versions of Bison (3.5 and later) it is advantageous to not have literals */
<INITIAL>"-" { return XyceExpression::ExpressionParser::token::TOK_MINUS; }
<INITIAL>"+" { return XyceExpression::ExpressionParser::token::TOK_PLUS; }
<INITIAL>"*" { return XyceExpression::ExpressionParser::token::TOK_MULTIPLY; }
<INITIAL>"/" { return XyceExpression::ExpressionParser::token::TOK_DIVIDE; }
<INITIAL>"=" { return XyceExpression::ExpressionParser::token::TOK_EQUALS; }
<INITIAL>"%" { return XyceExpression::ExpressionParser::token::TOK_MOD; }
<INITIAL>"," { return XyceExpression::ExpressionParser::token::TOK_COMMA; }
<INITIAL>";" { return XyceExpression::ExpressionParser::token::TOK_SEMICOLON; }
<INITIAL>":" { return XyceExpression::ExpressionParser::token::TOK_COLON; }
<INITIAL>"?" { return XyceExpression::ExpressionParser::token::TOK_QUESTION_MARK; }


<INITIAL>")" { return XyceExpression::ExpressionParser::token::TOK_RIGHTPAREN; }
<INITIAL>"(" { return XyceExpression::ExpressionParser::token::TOK_LEFTPAREN; }
<INITIAL>"}" { return XyceExpression::ExpressionParser::token::TOK_RIGHTCURLYBRACE; }
<INITIAL>"{" { return XyceExpression::ExpressionParser::token::TOK_LEFTCURLYBRACE; }

 /* logic ops */
<INITIAL>"~"    { return XyceExpression::ExpressionParser::token::TOK_NOT; } /* Xyce */
<INITIAL>"|"    { return XyceExpression::ExpressionParser::token::TOK_OR; }  /* Xyce */
<INITIAL>"||"   { return XyceExpression::ExpressionParser::token::TOK_OR; }    /* Hspice */
<INITIAL>"&"    { return XyceExpression::ExpressionParser::token::TOK_AND; } /* Xyce */
<INITIAL>"&&"   { return XyceExpression::ExpressionParser::token::TOK_AND; }  /* Hspice */
 /* <INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_XOR; } */ /* Xyce: caret (^) = XOR; Hspice caret (^) = power */
<INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_CARET; } /* Xyce: caret (^) = XOR; Hspice caret (^) = power */

 /* comparison ops */
<INITIAL>">"     { return XyceExpression::ExpressionParser::token::TOK_GT; }
<INITIAL>"<"     { return XyceExpression::ExpressionParser::token::TOK_LT; }
<INITIAL>"<>"    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"!="    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"=="    { return XyceExpression::ExpressionParser::token::TOK_EQ; }
<INITIAL>">="    { return XyceExpression::ExpressionParser::token::TOK_GE; }
<INITIAL>"<="    { return XyceExpression::ExpressionParser::token::TOK_LE; }

 /* invalid ops */
<INITIAL>"#"    { return XyceExpression::ExpressionParser::token::TOK_INVALID_OP; }
<INITIAL>"=#"    { return XyceExpression::ExpressionParser::token::TOK_INVALID_OP; }
<INITIAL>"#="    { return XyceExpression::ExpressionParser::token::TOK_INVALID_OP; }

 /* variants on voltage node ops */

<INITIAL>{VOP}    |
<INITIAL>"VDB"    |
<INITIAL>{LEAD}   |
<INITIAL>{IOP}    |
<INITIAL>"IDB"    |
<INITIAL>{NOP}    |
<INITIAL>"NDB"    |
<INITIAL>{SOP}    |
<INITIAL>"SDB"    |
<INITIAL>{YOP}    |
<INITIAL>"YDB"    |
<INITIAL>{ZOP}    |
<INITIAL>"ZDB"    |
<INITIAL>"DNO"    |
<INITIAL>"DNI"    {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"exp"   { return XyceExpression::ExpressionParser::token::TOK_EXP_CONST; }

<INITIAL>{VOP}[ \t]*"("     |
<INITIAL>"VDB"[ \t]*"("     {
  BEGIN(voltagenodedef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGE;
}

 /* single node name by itself */
<voltagenodedef>{NODE}")" {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* first node in a comma-separated list of two nodes */
<voltagenodedef>{NODE}"," {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

<INITIAL>{LEAD}[ \t]*"("    {
  BEGIN(currentdevdef);
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_LEAD;
} 

<INITIAL>{IOP}[ \t]*"("    |
<INITIAL>"IDB"[ \t]*"("    {
  BEGIN(currentdevdef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_CURRENT;
} 

 /* single device name by itself */
<currentdevdef>{DEV}")" {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_CURRENTDEV;
}

<INITIAL>{NOP}[ \t]*"("    |
<INITIAL>"NDB"[ \t]*"("    {
  BEGIN(internaldevdef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_INTERNAL_VAR;
}

  /* S-parameters */
<INITIAL>{SOP}[ \t]*"("     |
<INITIAL>"SDB"[ \t]*"("     {
  BEGIN(sparamdef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM;
}

 /* second index in a comma-separated list of two sparam indices */
<sparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two sparam indices */
<sparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

  /* Y-parameters */
<INITIAL>{YOP}[ \t]*"("     |
<INITIAL>"YDB"[ \t]*"("     {
  BEGIN(yparamdef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM;
}

 /* second index in a comma-separated list of two yparam indices */
<yparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two yparam indices */
<yparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

  /* Z-parameters */
<INITIAL>{ZOP}[ \t]*"("     |
<INITIAL>"ZDB"[ \t]*"("     {
  BEGIN(zparamdef);
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM;
}

 /* second index in a comma-separated list of two zparam indices */
<zparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two zparam indices */
<zparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->build<int>()=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* single device name by itself */
 /* this needs to be a "NODE" macro rather than "DEV" b/c the N() can be used to obtain voltage nodes */
<internaldevdef>{NODE}")" {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_INTERNAL_VAR_DEV;
}

<INITIAL>"DNO"[ \t]*"("    |
<INITIAL>"DNI"[ \t]*"("    {
  BEGIN(noisedevdef);
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_NOISE;
}

 /* single device name by itself */
<noisedevdef>{DEV}")" {
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
}

 /* first device in a comma-separated list of two device */
<noisedevdef>{DEV}"," { 
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
}

<INITIAL>"ONOISE"     {
  lvalp->build<std::string *>()=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_ONOISE;
}

<INITIAL>"INOISE"     {
  lvalp->build<std::string *>()=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_INOISE;
}

 /* single device name by itself */
 /* <noisedevdef>{ID}")" {
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
  } */

<INITIAL>"P"[ \t]*"("    | 
<INITIAL>"W"[ \t]*"("    { 
  BEGIN(currentdevdef);
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_POWER; 
}

 /* built in functions */
<INITIAL>"Ph"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_PHASE; }
<INITIAL>"R"[ \t]*"("     { return XyceExpression::ExpressionParser::token::TOK_REAL; }
<INITIAL>"Re"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_REAL; }
 /* <INITIAL>"Im"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_IMAG; } */   /* cannot support this and also IM for imaginary currents */
<INITIAL>"Img"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_IMAG; }

<INITIAL>"max"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_MAX; }
<INITIAL>"min"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_MIN; }

<INITIAL>"exp"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_EXP; }
<INITIAL>"log"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_LOG; }  /* Hspice version of ln. Xyce version of log=log10 */

<INITIAL>"sin"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_SIN; }
<INITIAL>"cos"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_COS; }
<INITIAL>"M"[ \t]*"("     { return XyceExpression::ExpressionParser::token::TOK_MABS; }
<INITIAL>"abs"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_ABS; }
<INITIAL>"acos"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_ACOS; }
<INITIAL>"acosh"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_ACOSH; }
<INITIAL>"asin"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_ASIN; }
<INITIAL>"asinh"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_ASINH; }
<INITIAL>"arctan"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_ATAN; }
<INITIAL>"atan"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_ATAN; }
<INITIAL>"atan2"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_ATAN2; }
<INITIAL>"atanh"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_ATANH; }
<INITIAL>"cosh"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_COSH; }
<INITIAL>"ln"[ \t]*"("     { return XyceExpression::ExpressionParser::token::TOK_LN; }     /* Xyce version of ln */
<INITIAL>"log10"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_LOG10; }
<INITIAL>"sinh"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_SINH; }
<INITIAL>"tan"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_TAN; }
<INITIAL>"tanh"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_TANH; }

<INITIAL>"ceil"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_CEIL; }
<INITIAL>"floor"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_FLOOR; }
<INITIAL>"limit"[ \t]*"("  { return XyceExpression::ExpressionParser::token::TOK_LIMIT; }

<INITIAL>"sdt"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_SDT; }
<INITIAL>"ddt"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_DDT; }
<INITIAL>"ddx"[ \t]*"("   { return XyceExpression::ExpressionParser::token::TOK_DDX; }

<INITIAL>"if"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_IF; }
<INITIAL>"int"[ \t]*"(" { return XyceExpression::ExpressionParser::token::TOK_INT; }

<INITIAL>"pwr"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_POWFUNC; }
<INITIAL>"pow"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_POWFUNC; }
<INITIAL>"**"   {return XyceExpression::ExpressionParser::token::TOK_POW; }
  /* <INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_POW; } */ /* Xyce: caret (^) = XOR; Hspice caret (^) = power */
<INITIAL>"pwrs"[ \t]*"("  {return XyceExpression::ExpressionParser::token::TOK_POWSFUNC; }

<INITIAL>"time" {return XyceExpression::ExpressionParser::token::TOK_TIME; }
<INITIAL>"temp" {return XyceExpression::ExpressionParser::token::TOK_TEMP; }
<INITIAL>"temper" {return XyceExpression::ExpressionParser::token::TOK_TEMP; }
<INITIAL>"vt" {return XyceExpression::ExpressionParser::token::TOK_VT; }
<INITIAL>"freq" {return XyceExpression::ExpressionParser::token::TOK_FREQ; }
<INITIAL>"hertz" {return XyceExpression::ExpressionParser::token::TOK_FREQ; }
<INITIAL>"gmin" {return XyceExpression::ExpressionParser::token::TOK_GMIN; }
<INITIAL>"pi" {return XyceExpression::ExpressionParser::token::TOK_PI; }
<INITIAL>"constctok" {return XyceExpression::ExpressionParser::token::TOK_C_TO_K; }
<INITIAL>"ctok" {return XyceExpression::ExpressionParser::token::TOK_C_TO_K; }

 /* random sampling related functions */
<INITIAL>"rand"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_RAND; }
<INITIAL>"unif"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_UNIF; }
<INITIAL>"aunif"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_AUNIF; }
<INITIAL>"agauss"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_AGAUSS; }
<INITIAL>"gauss"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_GAUSS; }

<INITIAL>"sgn"[ \t]*"("    {return XyceExpression::ExpressionParser::token::TOK_SGN; }
<INITIAL>"sign"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_SIGN; }
<INITIAL>"stp"[ \t]*"("    {return XyceExpression::ExpressionParser::token::TOK_STP; }
<INITIAL>"sqrt"[ \t]*"("   {return XyceExpression::ExpressionParser::token::TOK_SQRT; }

<INITIAL>"spline"[ \t]*"("  |
<INITIAL>"akima"[ \t]*"("   |
<INITIAL>"cubic"[ \t]*"("   |
<INITIAL>"wodicka"[ \t]*"(" |
<INITIAL>"bli"[ \t]*"("     |
<INITIAL>TABLE[ \t]*"("     {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_TABLE; 
}

<INITIAL>"spline"[ \t]*"{"  |
<INITIAL>"akima"[ \t]*"{"   |
<INITIAL>"cubic"[ \t]*"{"   |
<INITIAL>"wodicka"[ \t]*"{" |
<INITIAL>"bli"[ \t]*"{"     |
<INITIAL>TABLE[ \t]*"{"     {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_TABLE2; 
}

<INITIAL>"spline"[ \t]*"("[ \t]*"{"  |
<INITIAL>"akima"[ \t]*"("[ \t]*"{"   |
<INITIAL>"cubic"[ \t]*"("[ \t]*"{"   |
<INITIAL>"wodicka"[ \t]*"("[ \t]*"{" |
<INITIAL>"bli"[ \t]*"("[ \t]*"{"     |
<INITIAL>TABLE[ \t]*"("[ \t]*"{"     {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_TABLE3; 
} 

<INITIAL>"spline"[ \t]*"("[ \t]*"\""      |
<INITIAL>"akima"[ \t]*"("[ \t]*"\""       |
<INITIAL>"cubic"[ \t]*"("[ \t]*"\""       |
<INITIAL>"wodicka"[ \t]*"("[ \t]*"\""     |
<INITIAL>"bli"[ \t]*"("[ \t]*"\""         |
<INITIAL>"table"[ \t]*"("[ \t]*"\""       |
<INITIAL>"fasttable"[ \t]*"("[ \t]*"\""   |
<INITIAL>"splinefile"[ \t]*"("[ \t]*"\""  |
<INITIAL>"akimafile"[ \t]*"("[ \t]*"\""   |
<INITIAL>"cubicfile"[ \t]*"("[ \t]*"\""   |
<INITIAL>"wodickafile"[ \t]*"("[ \t]*"\"" |
<INITIAL>"blifile"[ \t]*"("[ \t]*"\""     |
<INITIAL>"fasttablefile"[ \t]*"("[ \t]*"\""  |
<INITIAL>"tablefile"[ \t]*"("[ \t]*"\""   {
    BEGIN(tablefilename);
    std::string tmp = std::string(YYText());
    Xyce::Util::toUpper(tmp);
    lvalp->build<std::string *>()=new std::string( tmp );
    return XyceExpression::ExpressionParser::token::TOK_TABLEFILE; 
  }

<tablefilename>{FILENAME}+ { 
  lvalp->build<std::string *>()=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_FILENAME; 
}

<tablefilename>"\""[ \t]*")" {
    BEGIN(INITIAL);
    return XyceExpression::ExpressionParser::token::TOK_TABLEFILE_END; 
  }

<INITIAL>"schedule"[ \t]*"("  {return XyceExpression::ExpressionParser::token::TOK_SCHEDULE; }
<INITIAL>"uramp"[ \t]*"("  {return XyceExpression::ExpressionParser::token::TOK_URAMP; }

<INITIAL>"poly"  { return XyceExpression::ExpressionParser::token::TOK_BADPOLY; } 

 /*<INITIAL>"poly"[ \t]*"("    {
    BEGIN(polydef);
    return XyceExpression::ExpressionParser::token::TOK_POLY; 
  } 

<polydef>{DIGIT}+")" {
    std::string tmp(YYText());
    tmp.erase(tmp.length()-1);
    lvalp->build<double>()=Xyce::Util::Value(tmp);
    BEGIN(polyargs);
    return XyceExpression::ExpressionParser::token::TOK_POLYSIZE; 
  } */

<INITIAL>"poly"[ \t]*"("[ \t]*{DIGIT}+")" {
    std::string tmp(YYText());
    tmp.erase(std::remove_if(tmp.begin(),tmp.end(), ::isspace), tmp.end());
    tmp.erase(tmp.length()-1);

    std::size_t found = tmp.find_first_of("(");
    if (found!=std::string::npos) { tmp.erase(tmp.begin(),tmp.begin()+found+1); }

    lvalp->build<double>()=Xyce::Util::Value(tmp);
    BEGIN(polyargs);
    return XyceExpression::ExpressionParser::token::TOK_POLYSIZE; 
  } 

 /* numbers with suffixes */
<polyargs>[-+]?{DIGIT}+{SUFFIX}?{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."{DIGIT}*{SUFFIX}?{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}{SUFFIX}?{UNIT}? |
 /* numbers with 3-character suffixes; meg and mil */
<polyargs>[-+]?{DIGIT}+{MEG}{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."{DIGIT}*{MEG}{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}{MEG}{UNIT}? |
<polyargs>[-+]?{DIGIT}+{MIL}{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."{DIGIT}*{MIL}{UNIT}? |
<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}{MIL}{UNIT}? {
     lvalp->build<double>()=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>"V""("     {
  BEGIN(polyvoltagenodedef); 
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_POLY_VOLTAGE;
}

 /* single node name by itself */
<polyvoltagenodedef>{NODE}")" {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(polyargs);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* first node in a comma-separated list of two nodes */
<polyvoltagenodedef>{NODE}"," {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

<polyargs>"I""("    {
  BEGIN(polycurrentdevdef);
  std::string tmp = std::string(YYText());
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_POLY_CURRENT;
}

 /* single device name by itself */
<polycurrentdevdef>{DEV}")" {
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  BEGIN(polyargs);
  return XyceExpression::ExpressionParser::token::TOK_CURRENTDEV;
  }

<polyargs>")" { return XyceExpression::ExpressionParser::token::TOK_RIGHTPAREN; }
<polyargs>"(" { return XyceExpression::ExpressionParser::token::TOK_LEFTPAREN; }
<polyargs>"}" { return XyceExpression::ExpressionParser::token::TOK_RIGHTCURLYBRACE; }
<polyargs>"{" { return XyceExpression::ExpressionParser::token::TOK_LEFTCURLYBRACE; }


 /* letters */
<polyargs>{ID}+ { 
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_POLY_WORD; 
}

<polyargs>[ \t]+  { llocp->step(); }   /* eat up whitespace */

<polyargs>.  {
  BEGIN(INITIAL);
  lvalp->build<std::string *>()=new std::string(YYText());
  yyless(0);
}

 /* a bunch of spice source functions.  Doing them the right way this time */
<INITIAL>"spice_pulse"[ \t]*"(" {return XyceExpression::ExpressionParser::token::TOK_SPICE_PULSE; }
<INITIAL>"spice_sin"[ \t]*"(" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SIN; }
<INITIAL>"spice_exp"[ \t]*"(" {return XyceExpression::ExpressionParser::token::TOK_SPICE_EXP; }
<INITIAL>"spice_sffm"[ \t]*"(" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SFFM; }

 /* letters */
<INITIAL>{ID}+ { 
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp );
  return XyceExpression::ExpressionParser::token::TOK_WORD; 
}

 /* functions */
<INITIAL>{ID}+[ \t]*"("  { 
  std::string tmp = std::string(YYText());
  Xyce::Util::toUpper(tmp);
  lvalp->build<std::string *>()=new std::string( tmp.substr(0, tmp.length()-1) );
  return XyceExpression::ExpressionParser::token::TOK_FUNCTION; 
}

 /* in order for # to be an allowed character in parameter and function names, this one-line comment rule had to go */
 /* "#"[^\n\r]* */   /* eat up one-line comments */
[ \t]+                     { llocp->step(); }   /* eat up whitespace */
(\n+|(\r\n)+|\r+)          { llocp->lines(YYLeng()); }

.  { 
    std::vector<std::string> errStr(1,std::string("Unrecognized character: " + std::string(YYText()) + " in expression " + expressionString )) ; 
    yyerror(errStr);
   }

%%
