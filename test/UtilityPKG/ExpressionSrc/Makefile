#

CXX = /usr/bin/clang++
FLEX = /opt/local/bin/flex
BISON = /opt/local/bin/bison
LIBDIRTRI=-L/Users/erkeite/XYCE/Xyce_Lib/intel64_osx_clang_trilinos12.12_withStokhos_withBasker/lib
LIBPATH = -L/opt/local/lib  $(LIBDIRTRI)

CXXFLAGS_AST =  -std=c++11  -g  -I/opt/local/include/ $(INCLUDE)
INCLUDE = -I/Users/erkeite/XYCE/Xyce_Lib/intel64_osx_clang_trilinos12.12_withStokhos_withBasker/include  -I../ExpressionSrc  -I./  -I../ -I../../ErrorHandlingPKG/ -I../../ParallelDistPKG/ -I../../IOInterfacePKG/
CXXFLAGS =  -std=c++11 -DUSE_TYPE_DOUBLE  -g -I/opt/local/include/ $(INCLUDE)

#    clang++ Test.C -I/opt/local/include -L/opt/local/lib -lgtest -lgtest_main -pthread

RCPLIBS =  -lteuchoscore -lkokkoscore -lteuchosparameterlist -lteuchoscomm -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosnumerics -lteuchosremainder

all:  test

test:  astUnitTest parserUnitTest   astCopyConstUnitTest  parserCopyUnitTest

# main program
astUnitTest: astUnitTest.C ast.o  ../ExpressionSrc/ast.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h  N_UTL_BreakPoint.o N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o  N_UTL_NetlistLocation.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o
	$(CXX) $(CXXFLAGS_AST) $(LIBPATH) -o astUnitTest astUnitTest.C ast.o N_UTL_BreakPoint.o N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o  N_UTL_NetlistLocation.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o  -lm -lgtest -lgtest_main  $(RCPLIBS)

astCopyConstUnitTest: astCopyConstUnitTest.C ast.o  ../ExpressionSrc/ast.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../N_UTL_Interpolators.h ../N_UTL_CheckIfvalidFile.h  N_UTL_BreakPoint.o N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o  N_UTL_NetlistLocation.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o
	$(CXX) $(CXXFLAGS_AST) $(LIBPATH) -o astCopyConstUnitTest astCopyConstUnitTest.C ast.o  N_UTL_BreakPoint.o N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o  N_UTL_NetlistLocation.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o -lm -lgtest -lgtest_main $(RCPLIBS)

# expression lexer/parser
ExpressionLexer.cxx: ../ExpressionSrc/ExpressionLexer.l  ../ExpressionSrc/ast.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../ExpressionSrc/ExpressionType.h ../N_UTL_ExtendedString.h ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h
	$(FLEX) -i -+ --full -p -p  --prefix=exp  -o ExpressionLexer.cxx ../ExpressionSrc/ExpressionLexer.l

ExpressionParser.cxx: ../ExpressionSrc/ExpressionParser.yxx  ../ExpressionSrc/ast.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../ExpressionSrc/ExpressionType.h ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h 
	$(BISON) --graph -v -o$@ ../ExpressionSrc/ExpressionParser.yxx

ExpressionLexer.o: ExpressionLexer.cxx
	$(CXX) $(CXXFLAGS) -c $<

ast.o:  ../ExpressionSrc/ast.C  ../ExpressionSrc/ast.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h  ../N_UTL_CheckIfValidFile.h
	$(CXX) $(CXXFLAGS) -c $<

ExpressionParser.o: ExpressionParser.cxx
	$(CXX) $(CXXFLAGS) -c $<

expressionGroup.o: ../ExpressionSrc/expressionGroup.C ../ExpressionSrc/expressionGroup.h 
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_BreakPoint.o: ../N_UTL_BreakPoint.C ../N_UTL_BreakPoint.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_ExtendedString.o: ../N_UTL_ExtendedString.C ../N_UTL_ExtendedString.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_NoCase.o: ../N_UTL_NoCase.C ../N_UTL_NoCase.h
	$(CXX) $(CXXFLAGS) -c $<

N_ERH_Message.o: ../../ErrorHandlingPKG/N_ERH_Message.C ../../ErrorHandlingPKG/N_ERH_Message.h
	$(CXX) $(CXXFLAGS) -c $<

N_ERH_Messenger.o: ../../ErrorHandlingPKG/N_ERH_Messenger.C ../../ErrorHandlingPKG/N_ERH_Messenger.h
	$(CXX) $(CXXFLAGS) -c $<

N_ERH_ErrorMgr.o: ../../ErrorHandlingPKG/N_ERH_ErrorMgr.C ../../ErrorHandlingPKG/N_ERH_ErrorMgr.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_NetlistLocation.o: ../N_UTL_NetlistLocation.C ../N_UTL_NetlistLocation.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_Marshal.o: ../N_UTL_Marshal.C ../N_UTL_Marshal.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_LogStream.o: ../N_UTL_LogStream.C ../N_UTL_LogStream.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_ReportHandler.o: ../N_UTL_ReportHandler.C ../N_UTL_ReportHandler.h
	$(CXX) $(CXXFLAGS) -c $<

N_UTL_CheckIfValidFile.o: ../N_UTL_CheckIfValidFile.C ../N_UTL_CheckIfValidFile.h
	$(CXX) $(CXXFLAGS) -c $<

newExpression.o: ../ExpressionSrc/newExpression.C ../ExpressionSrc/newExpression.h ../ExpressionSrc/ExpressionLexer.h  ../ExpressionSrc/ast.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../ExpressionSrc/ExpressionType.h  ../N_UTL_BreakPoint.h  ExpressionParser.o ExpressionLexer.o ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h
	$(CXX) $(CXXFLAGS) -c $<

parserUnitTest: parserUnitTest.C ast.o  ../ExpressionSrc/ast.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../ExpressionSrc/newExpression.h newExpression.o ExpressionParser.o  ExpressionLexer.o expressionGroup.o N_UTL_BreakPoint.o N_UTL_ExtendedString.o N_UTL_NoCase.o N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o  N_UTL_NetlistLocation.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o  ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h ../../IOInterfacePKG/N_IO_fwd.h
	$(CXX) $(CXXFLAGS) $(LIBPATH) -o parserUnitTest parserUnitTest.C ast.o newExpression.o ExpressionLexer.o ExpressionParser.o  expressionGroup.o N_UTL_BreakPoint.o N_UTL_ExtendedString.o N_UTL_NoCase.o  N_ERH_ErrorMgr.o N_ERH_Message.o N_ERH_Messenger.o   N_UTL_NetlistLocation.o N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o -lm -lgtest -lgtest_main  $(RCPLIBS)


parserCopyUnitTest: parserCopyUnitTest.C ast.o  ../ExpressionSrc/ast.h ../ExpressionSrc/astbinary.h ../ExpressionSrc/astfuncs.h ../ExpressionSrc/astcomp.h ../ExpressionSrc/ast_spice_src.h ../ExpressionSrc/ast_random.h ../ExpressionSrc/astRandEnum.h ../ExpressionSrc/newExpression.h newExpression.o ExpressionParser.o  ExpressionLexer.o expressionGroup.o  N_UTL_BreakPoint.o N_UTL_ExtendedString.o N_UTL_NoCase.o N_UTL_NetlistLocation.o N_ERH_ErrorMgr.o  N_ERH_Message.o N_ERH_Messenger.o   N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o  ../N_UTL_Interpolators.h ../N_UTL_CheckIfValidFile.h ../../IOInterfacePKG/N_IO_fwd.h
	$(CXX) $(CXXFLAGS) $(LIBPATH) -o parserCopyUnitTest parserCopyUnitTest.C ast.o newExpression.o ExpressionLexer.o ExpressionParser.o  expressionGroup.o  N_UTL_BreakPoint.o N_UTL_ExtendedString.o N_UTL_NoCase.o N_UTL_NetlistLocation.o N_ERH_ErrorMgr.o N_ERH_Message.o  N_ERH_Messenger.o  N_UTL_Marshal.o N_UTL_LogStream.o N_UTL_ReportHandler.o N_UTL_CheckIfValidFile.o -lm -lgtest -lgtest_main  $(RCPLIBS)

clean:
	rm -f astUnitTest  parserUnitTest parserCopyUnitTest astCopyConstUnitTest  *.o \
		NetlistLexer.cxx NetlistParser.cxx  	stack.hh  position.hh location.hh  \
		ExpressionLexer.cxx  ExpressionParser.cxx  ExpressionParser.hxx ExpressionParser.dot

