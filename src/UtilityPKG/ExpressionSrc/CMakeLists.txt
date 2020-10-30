# class source list
target_sources(XyceStatic PUBLIC
  ast.C
  expressionGroup.C
  newExpression.C
  deviceExpressionGroup.C
  mainXyceExpressionGroup.C
  outputsXyceExpressionGroup.C
  xyceExpressionGroup.C)

target_include_directories(XyceStatic PUBLIC
 "${CMAKE_CURRENT_SOURCE_DIR}"
 "${CMAKE_CURRENT_BINARY_DIR}"
 "${FLEX_INCLUDE_DIR}" )
