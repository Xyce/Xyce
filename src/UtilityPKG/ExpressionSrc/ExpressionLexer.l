
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

 /* ID is used for device instance names, particularly inside of current vars */
 /* NODE is used for voltage node names, particularly inside of voltage vars */
 /* Both ID and NODE include the colon, to handle Xyce subcircuit syntax */
 /* Other simulators use other delimiters for subcircuit names, so might */
 /* want to consider supporting delimiters other than the colon. */
DIGIT [0-9]
ID [a-zA-Z][a-zA-Z0-9_:]*
NODE [a-zA-Z0-9_:]*

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

 /* ----------------------------------------------------------*/
 /* single character ops */
<INITIAL>"+" |
"-" |
"*" |
"/" |
"=" |
"%" |
"|" |
"," |
";" |
":" |
"?" |
"J" |
"{" |
"}" |
"(" |
")"     { return((YYText())[0]); }

 /* logic ops */
<INITIAL>"~"    { return XyceExpression::ExpressionParser::token::TOK_NOT; } /* Xyce */
<INITIAL>"|"    { return XyceExpression::ExpressionParser::token::TOK_OR; }  /* Xyce */
<INITIAL>"||"   { return XyceExpression::ExpressionParser::token::TOK_OR; }  /* Hspice */
<INITIAL>"&"    { return XyceExpression::ExpressionParser::token::TOK_AND; } /* Xyce */
<INITIAL>"&&"   { return XyceExpression::ExpressionParser::token::TOK_AND; } /* Hspice */
<INITIAL>"^"    {return XyceExpression::ExpressionParser::token::TOK_XOR; } /* Xyce: caret (^) = XOR; Hspice caret (^) = power */

 /* comparison ops */
<INITIAL>">"     { return XyceExpression::ExpressionParser::token::TOK_GT; }
<INITIAL>"<"     { return XyceExpression::ExpressionParser::token::TOK_LT; }
<INITIAL>"<>"    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"!="    { return XyceExpression::ExpressionParser::token::TOK_NE; }
<INITIAL>"=="    { return XyceExpression::ExpressionParser::token::TOK_EQ; }
<INITIAL>">="    { return XyceExpression::ExpressionParser::token::TOK_GE; }
<INITIAL>"<="    { return XyceExpression::ExpressionParser::token::TOK_LE; }

<INITIAL>"V"     {
  BEGIN(voltagenodedef);
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single node name by itself */
<voltagenodedef>"("{NODE}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase(0,1);// lop off open paren
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* first node in a comma-separated list of two nodes */
<voltagenodedef>"("{NODE}"," {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase(0,1);// lop off open paren
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off comma
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* second node in a comma-separated list of two nodes. This one calls the next context (expression) */
<voltagenodedef>{NODE}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_VOLTAGENODE;
}

 /* followed by unrecognized character; In this case "V" is probably a parameter, rather than a current variable */
 /* yyless is used here so that the unrecognized character will be put back on the stack, and processed */
 /* by the rules under the INITIAL state */
<voltagenodedef>. { BEGIN(INITIAL); yyless(0); }

<INITIAL>"I"     {
  BEGIN(currentdevdef);
  lvalp->sval=new std::string(YYText());
  return XyceExpression::ExpressionParser::token::TOK_WORD;
}

 /* single device name by itself */
<currentdevdef>"("{ID}")" {
  lvalp->sval=new std::string(YYText());
  (lvalp->sval)->erase(0,1);// lop off open paren
  (lvalp->sval)->erase((lvalp->sval)->length()-1); // lop off close paren
  BEGIN(INITIAL);
  return XyceExpression::ExpressionParser::token::TOK_CURRENTDEV;
  }

 /* followed by unrecognized character; In this case "I" is probably a parameter, rather than a current variable */
 /* yyless is used here so that the unrecognized character will be put back on the stack, and processed */
 /* by the rules under the INITIAL state */
<currentdevdef>. { BEGIN(INITIAL); yyless(0); }

 /* built in functions */
<INITIAL>"P"     { return XyceExpression::ExpressionParser::token::TOK_PHASE; }

<INITIAL>"R"     { return XyceExpression::ExpressionParser::token::TOK_REAL; }
<INITIAL>"Re"    { return XyceExpression::ExpressionParser::token::TOK_REAL; }
<INITIAL>"Im"    { return XyceExpression::ExpressionParser::token::TOK_IMAG; }
<INITIAL>"Img"   { return XyceExpression::ExpressionParser::token::TOK_IMAG; }

<INITIAL>"max"   { return XyceExpression::ExpressionParser::token::TOK_MAX; }
<INITIAL>"min"   { return XyceExpression::ExpressionParser::token::TOK_MIN; }

<INITIAL>"exp"   { return XyceExpression::ExpressionParser::token::TOK_EXP; }
<INITIAL>"log"   { return XyceExpression::ExpressionParser::token::TOK_LOG; }  /* Hspice version of ln. Xyce version of log=log10 */

<INITIAL>"sin"   { return XyceExpression::ExpressionParser::token::TOK_SIN; }
<INITIAL>"cos"   { return XyceExpression::ExpressionParser::token::TOK_COS; }
<INITIAL>"M"     { return XyceExpression::ExpressionParser::token::TOK_ABS; }
<INITIAL>"abs"   { return XyceExpression::ExpressionParser::token::TOK_ABS; }
<INITIAL>"acos"  { return XyceExpression::ExpressionParser::token::TOK_ACOS; }
<INITIAL>"acosh" { return XyceExpression::ExpressionParser::token::TOK_ACOSH; }
<INITIAL>"asin"  { return XyceExpression::ExpressionParser::token::TOK_ASIN; }
<INITIAL>"asinh" { return XyceExpression::ExpressionParser::token::TOK_ASINH; }
<INITIAL>"atan"  { return XyceExpression::ExpressionParser::token::TOK_ATAN; }
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
<INITIAL>"pi" {return XyceExpression::ExpressionParser::token::TOK_PI; }

 /* random sampling related functions */
<INITIAL>"rand"  {return XyceExpression::ExpressionParser::token::TOK_RAND; }
<INITIAL>"agauss"  {return XyceExpression::ExpressionParser::token::TOK_AGAUSS; }
<INITIAL>"gauss"  {return XyceExpression::ExpressionParser::token::TOK_GAUSS; }

<INITIAL>"sgn"   {return XyceExpression::ExpressionParser::token::TOK_SGN; }
<INITIAL>"sign"  {return XyceExpression::ExpressionParser::token::TOK_SIGN; }
<INITIAL>"stp"   {return XyceExpression::ExpressionParser::token::TOK_STP; }
<INITIAL>"sqrt"  {return XyceExpression::ExpressionParser::token::TOK_SQRT; }
<INITIAL>"table" {return XyceExpression::ExpressionParser::token::TOK_TABLE; }
<INITIAL>"uramp" {return XyceExpression::ExpressionParser::token::TOK_URAMP; }

 /* a bunch of spice source functions.  Doing them the right way this time */
<INITIAL>"spice_pulse" {return XyceExpression::ExpressionParser::token::TOK_SPICE_PULSE; }
<INITIAL>"spice_sin" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SIN; }
<INITIAL>"spice_exp" {return XyceExpression::ExpressionParser::token::TOK_SPICE_EXP; }
<INITIAL>"spice_sffm" {return XyceExpression::ExpressionParser::token::TOK_SPICE_SFFM; }

 /* letters */
{ID}+ { lvalp->sval=new std::string(YYText());  return XyceExpression::ExpressionParser::token::TOK_WORD; }

"#"[^\n\r]*     /* eat up one-line comments */
[ \t]+                     { llocp->step(); }   /* eat up whitespace */
(\n+|(\r\n)+|\r+)          { llocp->lines(YYLeng()); }

.  { std::cout << "Unrecognized character: " <<  YYText() << std::endl; }

%%
