
#ifndef ExpressionLexer_H
#define ExpressionLexer_H

#include <iostream>
#include <string>

namespace Xyce {
namespace Util {

class ExpressionLexer: public yyFlexLexer
{
public:
  ExpressionLexer(
    const std::string & exprStr,
    std::istream *              input = 0,
    std::ostream *              output = 0)
    : yyFlexLexer(input,output),
      expressionString(exprStr)
  {}

  virtual ~ExpressionLexer()
  {}

  int getToken(XyceExpression::ExpressionParser::semantic_type *lvalp, XyceExpression::location *llocp);

  const std::string expressionString;
};

}
}

#endif
