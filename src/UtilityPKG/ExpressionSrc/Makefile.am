
AM_CPPFLAGS = @Xyce_INCS@

noinst_LTLIBRARIES = libExpressionSrc.la

# these are files built by bison and flex
BUILT_SOURCES= ExpressionParser.cxx \
            ExpressionParser.hxx \
            ExpressionLexer.cxx \
            location.hh 


CLEANFILES= $(BUILT_SOURCES)

libExpressionSrc_la_SOURCES = \
  expressionGroup.C \
  newExpression.C \
  mainXyceExpressionGroup.C \
  outputsXyceExpressionGroup.C \
  xyceExpressionGroup.C \
  ExpressionLexer.h \
  ExpressionType.h \
  ast.h \
  ast_spice_src.h \
  astbinary.h \
  astcomp.h \
  astfuncs.h \
  expressionGroup.h \
  newExpression.h \
  mainXyceExpressionGroup.h \
  outputsXyceExpressionGroup.h \
  xyceExpressionGroup.h  \
	expressionParamTypes.h

# these can't be in "SOURCES" above or automake screams about YACC.
EXTRA_DIST=ExpressionLexer.l \
  ExpressionParser.yxx

libExpressionSrc_la_SOURCES += $(BUILT_SOURCES)


#We have to do it this way, because automake appears to assume it can use
# yacc and bison interchangably, and we can't.

ExpressionParser.cxx: $(srcdir)/ExpressionParser.yxx
	$(BISON) --graph -o$@ $(srcdir)/ExpressionParser.yxx

ExpressionLexer.cxx: $(srcdir)/ExpressionLexer.l
	$(LEX)  -i -+ --prefix=exp  -o$@ $(srcdir)/ExpressionLexer.l

ExpressionParser.hxx: ExpressionParser.cxx
location.hh: ExpressionParser.cxx

