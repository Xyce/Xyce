
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <string>

#include "newNetlist.h"
#include "netlistData.h"

//-------------------------------------------------------------------------------
// Netlist Lexer/Parser header stuff
//
// Grrrr.  Stupid bison 2.4 stopped putting the pre-prologue into the header.
// need this forward declaration
namespace Xyce {
namespace Util {
class NetlistLexer;
}}

#include "NetlistParser.hxx"
// BLEAH!   This is here DUPLICATED from NetlistParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
  /* CAREFUL watch continuations! */
#define YY_DECL \
int NetlistLexer::getToken(NetlistParser::semantic_type *lvalp,  \
                            location *llocp)

  // YECH!  Work around very stupid way that multiple parsers/lexers are
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#undef yyFlexLexer 
#define yyFlexLexer netFlexLexer
#include <FlexLexer.h>
#include "NetlistLexer.h"
  // if we actually *used* yylex anywhere here it would be necessary to
  // undo that kludge.  Note that because of this stupidity, if the
  // "%name-prefix" is changed, this line needs to be changed, too.
  // BUT we don't actually use yylex anywhere in this file, so let's
  // leave yylex undefined.  If later it turns out that this *becomes*
  // necessary, uncomment the next line.
  //  #define yylex XyceDevicelex
//-------------------------------------------------------------------------------

namespace Xyce {
namespace Util {
void lexAndParseNetlist(std::string & netlistString, netlistData & nd)
{

#if 0
  std::ifstream netFile(netlistString.c_str(),std::ios::in);
  Xyce::Util::NetlistLexer netLexer(netlistString, &netFile);

  // play around with the lexer by itself:
  int tok;
  XyceExpression::NetlistParser::semantic_type lval;
  XyceExpression::location lloc;

  while ( (tok = netLexer.getToken( &lval, &lloc)) )
  {
    std::cout << "tok = " << tok ;
   
    if (tok ==  XyceExpression::NetlistParser::token::TOK_EXPSTR)
    {
      std::cout << "\tTOK_EXPSTR";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_DOTFUNC_BEGIN)
    {
      std::cout << "\tTOK_DOTFUNC_BEGIN";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_FUNCARG)
    {
      std::cout << "\tTOK_FUNCARG";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_FUNCNAME)
    {
      std::cout << "\tTOK_FUNCNAME";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }

    else if (tok ==  XyceExpression::NetlistParser::token::TOK_PARAM_BEGIN)
    {
      std::cout << "\tTOK_PARAM_BEGIN";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_PARAMNAME)
    {
      std::cout << "\tTOK_PARAMNAME";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_GLOBAL_PARAM_BEGIN)
    {
      std::cout << "\tTOK_GLOBAL_PARAM_BEGIN";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_GLOBAL_PARAMNAME)
    {
      std::cout << "\tTOK_GLOBAL_PARAMNAME";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else if (tok ==  XyceExpression::NetlistParser::token::TOK_NAME)
    {
      std::cout << "\tTOK_NAME";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }
    else 
    {
      std::cout << "\tTOK mystery";
      std::cout << "\tstring = \"" << *(lval.sval) << "\"";
    }

    std::cout <<std::endl;
  }
#endif

  std::cout << "Testing and evaluating netlist: " << netlistString << std::endl;
  std::ifstream netlistFile(netlistString.c_str(),std::ios::in);
  Xyce::Util::NetlistLexer netlistLexer(netlistString,&netlistFile);

  XyceExpression::NetlistParser netlistParser(&netlistLexer,nd);

  if (netlistParser.parse() != 0)
  {
    std::cout << "netlist OOPS" <<std::endl;
  }
  else
  {
    std::cout << "Successful netlist Parse!" <<std::endl;
  }

  return;
}

}
}
