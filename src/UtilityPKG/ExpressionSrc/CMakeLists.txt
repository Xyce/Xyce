# class source list
target_sources(XyceLib PRIVATE
  ast.C
  expressionGroup.C
  newExpression.C
  deviceExpressionGroup.C
  mainXyceExpressionGroup.C
  outputsXyceExpressionGroup.C)

target_include_directories(XyceLib PUBLIC
 $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> 
 $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> 
 $<BUILD_INTERFACE:${FLEX_INCLUDE_DIR}> )
