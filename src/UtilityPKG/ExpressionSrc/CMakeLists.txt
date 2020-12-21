# class source list
target_sources(XyceLib PRIVATE
  ast.C
  expressionGroup.C
  newExpression.C
  deviceExpressionGroup.C
  mainXyceExpressionGroup.C
  outputsXyceExpressionGroup.C
  xyceExpressionGroup.C)

target_include_directories(XyceLib PUBLIC
 "${CMAKE_CURRENT_SOURCE_DIR}"
 "${CMAKE_CURRENT_BINARY_DIR}"
 "${FLEX_INCLUDE_DIR}" )
