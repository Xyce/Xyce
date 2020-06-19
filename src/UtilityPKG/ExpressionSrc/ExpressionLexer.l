
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

// these here because Bison 2.4 no longer puts the pre-prologue into the
// header file... grrrr.
namespace Xyce {
namespace Util {
class ExpressionLexer;
}}
// now we can include this without errors from bison 2.4
#include "ExpressionParser.hxx"
// BLEAH!   This is here DUPLICATED from ExpressionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
//#undef yyFlexLexer
/* CAREFUL watch continuations! */
#undef YY_DECL
#define YY_DECL \
  int Xyce::Util::ExpressionLexer::getToken(XyceExpression::ExpressionParser::semantic_type *lvalp, \
                                            XyceExpression::location *llocp)

  // YECH!  Work around very stupid way that multiple parsers/lexers are
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#include "ExpressionLexer.h"
  // undo that kludge.  Note that because of this stupidity, if the
  // "%name-prefix" is changed, this line needs to be changed, too.
#define yylex XyceExpressionlex


#include "location.hh"

%}
%option noyywrap
%option yyclass="Xyce::Util::ExpressionLexer"
%s voltagenodedef
%s currentdevdef
%s internaldevdef
%s noisedevdef
%s polydef
%s polyargs
%s polyvoltagenodedef
%s polycurrentdevdef
%s sparamdef
%s yparamdef
%s zparamdef

 /* ID is used for a variety of random strings, including function and 
    parameter names.  So, it cannot have characters in it that could 
    also be interpretted as operators, such as '+' and '-'.   This 
    distinguishes ID from DEV, below. */

 /* DEV is used for device instance names, particularly inside of current vars.  I(DEV) */
 /* NODE is used for voltage node names, particularly inside of voltage vars  V(NODE) and V(NODE,NODE) */

 /* Both DEV and NODE include the colon, to handle Xyce subcircuit syntax */
 /* Other simulators use other delimiters for subcircuit names, so might */
 /* want to consider supporting delimiters other than the colon. */
DIGIT [0-9]
ID [a-zA-Z!_`@#][a-zA-Z0-9_:!`@#]*

 /* to support bug 1034, the DEV needs to recognize some weird characters, including '+' and '-' */
DEV [a-zA-Z][a-zA-Z0-9_:!`@#\+\-]*

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

LEAD I[SDGBEC1-9]
%{
#define YY_USER_ACTION llocp->columns(YYLeng());
%}

%%
%{
  llocp->step();
%}

 /* ------------------------------------------------------------------------*/
 /* numbers with suffixes */
<INITIAL>{DIGIT}+[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

 /* numbers with 3-character suffixes; meg and mil */
<INITIAL>{DIGIT}+[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}+[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }


 /* dealing with Hz */
<INITIAL>{DIGIT}+[tTgGkKmMxXuUnNpPfF]?[Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[tTgGkKmMxXuUnNpPfF]?[Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[tTgGkKmMxXuUnNpPfF]?[Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

 /* numbers with 3-character suffixes; meg and mil */
<INITIAL>{DIGIT}+[mM][eE][gG][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[mM][eE][gG][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][eE][gG][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}+[mM][iI][lL][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."{DIGIT}*[mM][iI][lL][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

<INITIAL>{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][iI][lL][Hh][Zz] {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_NUM;
  }

 /* ----------------------------------------------------------*/
 /* single character ops */
<INITIAL>"+" |
"-" |
"*" |
"/" |
"=" |
"%" |
"," |
";" |
":" |
"?" |
"J" { return((YYText())[0]); }
  /* |
"{" |
"}" |
"(" |
")"      */
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
<INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_XOR; } /* Xyce: caret (^) = XOR; Hspice caret (^) = power */

 /* comparison ops */
<INITIAL>">"     { return XyceExpression::ExpressionParser::token::TOK_GT; }
<INITIAL>"<"     { return XyceExpression::ExpressionParser::token::TOK_LT; }
<INITIAL>"<>"    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"!="    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"=="    { return XyceExpression::ExpressionParser::token::TOK_EQ; }
<INITIAL>">="    { return XyceExpression::ExpressionParser::token::TOK_GE; }
<INITIAL>"<="    { return XyceExpression::ExpressionParser::token::TOK_LE; }

 /* variants on voltage node ops */

<INITIAL>{VOP}    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{VOP}[ \t]*"("     {
  BEGIN(voltagenodedef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"VDB"    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"VDB"[ \t]*"("     {
  BEGIN(voltagenodedef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single node name by itself */
<voltagenodedef>{NODE}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* first node in a comma-separated list of two nodes */
<voltagenodedef>{NODE}"," {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off comma
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

<INITIAL>{LEAD}     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{LEAD}[ \t]*"("    {
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_LEAD;
} 

<INITIAL>{IOP}     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
} 

<INITIAL>{IOP}[ \t]*"("    {
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
} 


<INITIAL>"IDB"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
} 

<INITIAL>"IDB"[ \t]*"("    {
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
} 

 /* single device name by itself */
<currentdevdef>{DEV}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_CURRENTDEV;
}

<INITIAL>{NOP}     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{NOP}[ \t]*"("    {
  BEGIN(internaldevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}


<INITIAL>"NDB"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"NDB"[ \t]*"("    {
  BEGIN(internaldevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

  /* S-parameters */
<INITIAL>{SOP}    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{SOP}[ \t]*"("     {
  BEGIN(sparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"SDB"    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"SDB"[ \t]*"("     {
  BEGIN(sparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* second index in a comma-separated list of two sparam indices */
<sparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two sparam indices */
<sparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

  /* Y-parameters */
<INITIAL>{YOP}    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{YOP}[ \t]*"("     {
  BEGIN(yparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"YDB"    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"YDB"[ \t]*"("     {
  BEGIN(yparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* second index in a comma-separated list of two yparam indices */
<yparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two yparam indices */
<yparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

  /* Z-parameters */
<INITIAL>{ZOP}    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>{ZOP}[ \t]*"("     {
  BEGIN(zparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"ZDB"    {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"ZDB"[ \t]*"("     {
  BEGIN(zparamdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* second index in a comma-separated list of two zparam indices */
<zparamdef>{DIGIT}+")" {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* first index in a comma-separated list of two zparam indices */
<zparamdef>{DIGIT}+"," {
  std::string tmp(YYText());
  tmp.erase(tmp.length()-1);
  lvalp->ival=Xyce::Util::Value(tmp);
  return XyceExpression::ExpressionParser::token::TOK_RFPARAM_ARG;
}

 /* single device name by itself */
 /* this needs to be a "NODE" macro rather than "DEV" b/c the N() can be used to obtain voltage nodes */
<internaldevdef>{NODE}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_INTERNAL_DEV_VAR;
}

<INITIAL>"DNO"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"DNO"[ \t]*"("    {
  BEGIN(noisedevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single device name by itself */
<noisedevdef>{DEV}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
}

 /* first device in a comma-separated list of two device */
<noisedevdef>{DEV}"," { 
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off comma
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
}

<INITIAL>"DNI"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"DNI"[ \t]*"("    {
  BEGIN(noisedevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

<INITIAL>"ONOISE"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_ONOISE;
}

<INITIAL>"INOISE"     {
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_INOISE;
}

 /* single device name by itself */
 /* <noisedevdef>{ID}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_NOISE_DEV_VAR;
  } */

<INITIAL>"P"[ \t]*"("    { 
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_POWER; 
}


<INITIAL>"W"[ \t]*"("    { 
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_POWER; 
}

 /* built in functions */
<INITIAL>"Ph"[ \t]*"("    { return XyceExpression::ExpressionParser::token::TOK_PHASE; }
<INITIAL>"R"[ \t]*"("     { return XyceExpression::ExpressionParser::token::TOK_REAL; }
<INITIAL>"Re""("    { return XyceExpression::ExpressionParser::token::TOK_REAL; }
 /* <INITIAL>"Im"    { return XyceExpression::ExpressionParser::token::TOK_IMAG; } */   /* cannot support this and also IM for imaginary currents */
<INITIAL>"Img"   { return XyceExpression::ExpressionParser::token::TOK_IMAG; }

<INITIAL>"max"   { return XyceExpression::ExpressionParser::token::TOK_MAX; }
<INITIAL>"min"   { return XyceExpression::ExpressionParser::token::TOK_MIN; }

<INITIAL>"exp"   { return XyceExpression::ExpressionParser::token::TOK_EXP; }
<INITIAL>"log"   { return XyceExpression::ExpressionParser::token::TOK_LOG; }  /* Hspice version of ln. Xyce version of log=log10 */

<INITIAL>"sin"   { return XyceExpression::ExpressionParser::token::TOK_SIN; }
<INITIAL>"cos"   { return XyceExpression::ExpressionParser::token::TOK_COS; }
<INITIAL>"M"[ \t]*"("     { return XyceExpression::ExpressionParser::token::TOK_MABS; }
<INITIAL>"abs"   { return XyceExpression::ExpressionParser::token::TOK_ABS; }
<INITIAL>"acos"  { return XyceExpression::ExpressionParser::token::TOK_ACOS; }
<INITIAL>"acosh" { return XyceExpression::ExpressionParser::token::TOK_ACOSH; }
<INITIAL>"asin"  { return XyceExpression::ExpressionParser::token::TOK_ASIN; }
<INITIAL>"asinh" { return XyceExpression::ExpressionParser::token::TOK_ASINH; }
<INITIAL>"arctan"  { return XyceExpression::ExpressionParser::token::TOK_ATAN; }
<INITIAL>"atan"  { return XyceExpression::ExpressionParser::token::TOK_ATAN; }
<INITIAL>"atan2"  { return XyceExpression::ExpressionParser::token::TOK_ATAN2; }
<INITIAL>"atanh" { return XyceExpression::ExpressionParser::token::TOK_ATANH; }
<INITIAL>"cosh"  { return XyceExpression::ExpressionParser::token::TOK_COSH; }
<INITIAL>"ln"    { return XyceExpression::ExpressionParser::token::TOK_LN; }     /* Xyce version of ln */
<INITIAL>"log10" { return XyceExpression::ExpressionParser::token::TOK_LOG10; }
<INITIAL>"sinh"  { return XyceExpression::ExpressionParser::token::TOK_SINH; }
<INITIAL>"tan"   { return XyceExpression::ExpressionParser::token::TOK_TAN; }
<INITIAL>"tanh"  { return XyceExpression::ExpressionParser::token::TOK_TANH; }

<INITIAL>"ceil"  { return XyceExpression::ExpressionParser::token::TOK_CEIL; }
<INITIAL>"floor" { return XyceExpression::ExpressionParser::token::TOK_FLOOR; }
<INITIAL>"limit" { return XyceExpression::ExpressionParser::token::TOK_LIMIT; }

<INITIAL>"sdt"  { return XyceExpression::ExpressionParser::token::TOK_SDT; }
<INITIAL>"ddt"  { return XyceExpression::ExpressionParser::token::TOK_DDT; }
<INITIAL>"ddx"  { return XyceExpression::ExpressionParser::token::TOK_DDX; }

<INITIAL>"if"   { return XyceExpression::ExpressionParser::token::TOK_IF; }
<INITIAL>"int"  { return XyceExpression::ExpressionParser::token::TOK_INT; }

<INITIAL>"pwr"  {return XyceExpression::ExpressionParser::token::TOK_POWFUNC; }
<INITIAL>"pow"  {return XyceExpression::ExpressionParser::token::TOK_POWFUNC; }
<INITIAL>"**"   {return XyceExpression::ExpressionParser::token::TOK_POW; }
  /* <INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_POW; } */ /* Xyce: caret (^) = XOR; Hspice caret (^) = power */
<INITIAL>"pwrs" {return XyceExpression::ExpressionParser::token::TOK_POWSFUNC; }

<INITIAL>"time" {return XyceExpression::ExpressionParser::token::TOK_TIME; }
<INITIAL>"temp" {return XyceExpression::ExpressionParser::token::TOK_TEMP; }
<INITIAL>"vt" {return XyceExpression::ExpressionParser::token::TOK_VT; }
<INITIAL>"freq" {return XyceExpression::ExpressionParser::token::TOK_FREQ; }
<INITIAL>"gmin" {return XyceExpression::ExpressionParser::token::TOK_GMIN; }
<INITIAL>"pi" {return XyceExpression::ExpressionParser::token::TOK_PI; }
<INITIAL>"constctok" {return XyceExpression::ExpressionParser::token::TOK_C_TO_K; }
<INITIAL>"ctok" {return XyceExpression::ExpressionParser::token::TOK_C_TO_K; }

 /* random sampling related functions */
<INITIAL>"rand"  {return XyceExpression::ExpressionParser::token::TOK_RAND; }
<INITIAL>"unif"  {return XyceExpression::ExpressionParser::token::TOK_UNIF; }
<INITIAL>"aunif"  {return XyceExpression::ExpressionParser::token::TOK_AUNIF; }
<INITIAL>"agauss"  {return XyceExpression::ExpressionParser::token::TOK_AGAUSS; }
<INITIAL>"gauss"  {return XyceExpression::ExpressionParser::token::TOK_GAUSS; }

<INITIAL>"sgn"   {return XyceExpression::ExpressionParser::token::TOK_SGN; }
<INITIAL>"sign"  {return XyceExpression::ExpressionParser::token::TOK_SIGN; }
<INITIAL>"stp"   {return XyceExpression::ExpressionParser::token::TOK_STP; }
<INITIAL>"sqrt"  {return XyceExpression::ExpressionParser::token::TOK_SQRT; }
<INITIAL>"table" {return XyceExpression::ExpressionParser::token::TOK_TABLE; }
<INITIAL>"schedule" {return XyceExpression::ExpressionParser::token::TOK_SCHEDULE; }
<INITIAL>"uramp" {return XyceExpression::ExpressionParser::token::TOK_URAMP; }

<INITIAL>"poly""("  {
    BEGIN(polydef);
    return XyceExpression::ExpressionParser::token::TOK_POLY; 
  } 

<polydef>{DIGIT}+")" {
    std::string tmp(YYText());
    tmp.erase(tmp.length()-1);
    lvalp->dval=Xyce::Util::Value(tmp);
    BEGIN(polyargs);
    return XyceExpression::ExpressionParser::token::TOK_POLYSIZE; 
  }

 /* numbers with suffixes */
<polyargs>[-+]?{DIGIT}+[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."{DIGIT}*[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[tTgGkKmMxXuUnNpPfF]?[sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

 /* numbers with 3-character suffixes; meg and mil */
<polyargs>[-+]?{DIGIT}+[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."{DIGIT}*[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][eE][gG][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}+[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."{DIGIT}*[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>[-+]?{DIGIT}*"."?{DIGIT}*[eE][\-\+]?{DIGIT}{1,3}[mM][iI][lL][sSvVaA]? {
     lvalp->dval=Xyce::Util::Value(std::string(YYText()));
     return XyceExpression::ExpressionParser::token::TOK_POLYCOEF;
  }

<polyargs>"V""("     {
  BEGIN(polyvoltagenodedef); 
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single node name by itself */
<polyvoltagenodedef>{NODE}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(polyargs);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* first node in a comma-separated list of two nodes */
<polyvoltagenodedef>{NODE}"," {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off comma
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

<polyargs>"I""("    {
  BEGIN(polycurrentdevdef);
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off open paren
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single device name by itself */
<polycurrentdevdef>{DEV}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(polyargs);
  return XyceExpression::ExpressionParser::token::TOK_CURRENTDEV;
  }

<polyargs>")" { return XyceExpression::ExpressionParser::token::TOK_RIGHTPAREN; }
<polyargs>"(" { return XyceExpression::ExpressionParser::token::TOK_LEFTPAREN; }
<polyargs>"}" { return XyceExpression::ExpressionParser::token::TOK_RIGHTCURLYBRACE; }
<polyargs>"{" { return XyceExpression::ExpressionParser::token::TOK_LEFTCURLYBRACE; }

<polyargs>[ \t]+  { llocp->step(); }   /* eat up whitespace */

<polyargs>.  {
  BEGIN(INITIAL);
  lvalp->sval=new std::string(YYText());
  yyless(0);
}


 /* a bunch of spice source functions.  Doing them the right way this time */
<INITIAL>"spice_pulse" {return XyceExpression::ExpressionParser::token::TOK_SPICE_PULSE; }
<INITIAL>"spice_sin" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SIN; }
<INITIAL>"spice_exp" {return XyceExpression::ExpressionParser::token::TOK_SPICE_EXP; }
<INITIAL>"spice_sffm" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SFFM; }

 /* letters */
{ID}+ { lvalp->sval=new std::string(YYText());  return XyceExpression::ExpressionParser::token::TOK_WORD; }

 /* in order for # to be an allowed character in parameter and function names, this one-line comment rule had to go */
 /* "#"[^\n\r]* */   /* eat up one-line comments */
[ \t]+                     { llocp->step(); }   /* eat up whitespace */
(\n+|(\r\n)+|\r+)          { llocp->lines(YYLeng()); }

.  { std::cout << "Unrecognized character: " <<  YYText() << std::endl; }

%%
