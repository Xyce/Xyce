# class source list
set ( Expression_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/expressionGroup.C
  ${CMAKE_CURRENT_SOURCE_DIR}/newExpression.C
  ${CMAKE_CURRENT_SOURCE_DIR}/deviceExpressionGroup.C
  ${CMAKE_CURRENT_SOURCE_DIR}/mainXyceExpressionGroup.C
  ${CMAKE_CURRENT_SOURCE_DIR}/outputsXyceExpressionGroup.C
  ${CMAKE_CURRENT_SOURCE_DIR}/xyceExpressionGroup.C)

set ( XYCE_SOURCES ${XYCE_SOURCES} ${Expression_SOURCES} PARENT_SCOPE)

include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "${FLEX_INCLUDE_DIR}" )

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
      PARENT_SCOPE )
