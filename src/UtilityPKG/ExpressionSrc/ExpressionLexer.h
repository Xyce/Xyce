
#ifndef ExpressionLexer_H
#define ExpressionLexer_H

#include <iostream>
#include <string>

namespace Xyce {
namespace Util {

class ExpressionLexer: public yyFlexLexer
//class ExpressionLexer: public expFlexLexer
{
public:
  ExpressionLexer(
    //const NetlistLocation &     netlist_location,
    const std::string &         expression_filename,
    std::istream *              input = 0,
    std::ostream *              output = 0)
    : yyFlexLexer(input,output),
      //netlistLocation_(netlist_location),
      expressionFilename_(expression_filename)
  {}

  virtual ~ExpressionLexer()
  {}

  int getToken(XyceExpression::ExpressionParser::semantic_type *lvalp, XyceExpression::location *llocp);

  //const NetlistLocation         netlistLocation_;
  const std::string &           expressionFilename_;
};

}
}

#endif
