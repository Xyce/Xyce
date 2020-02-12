
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <chrono>
#include <ctime>

#include "newNetlist.h"
#include "netlistData.h"
#include "newExpression.h"
#include "prototypeExpressionGroup.h"
#include "ast.h"

#define ALLOCATE_LOOP(NAME)   \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it )    \
  {    \
    it->second.exp_ = Teuchos::rcp(new Xyce::Util::newExpression(it->second.value_,grp));    \
    if ( !(it->second.exp_->lexAndParseExpression ()) )    \
    {    \
      std::cout << "OOPS. problem parsing " << it->first << " = " << it->second.value_ <<std::endl;    \
      exit(0);    \
    }    \
  }

//-------------------------------------------------------------------------------
// Function      : allocateExpressions
//
// Purpose       : This will go thru all the netlist containers,
//                 allocate expressions for each of them, and then lex/parse
//                 the initial expression string.  This function does not attempt
//                 to resolve unrecognized symbols.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 10/2?/2019
//-------------------------------------------------------------------------------
//void allocateExpressions(Teuchos::RCP<Xyce::Util::prototypeExpressionGroup> & grp, Xyce::Util::netlistData & nd)
void allocateExpressions(Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, Xyce::Util::netlistData & nd)
{
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;

  std::unordered_map<std::string, Xyce::Util::resistorData > & resistors = nd.getResistors();
  std::unordered_map<std::string, Xyce::Util::capacitorData > & capacitors = nd.getCapacitors();
  std::unordered_map<std::string, Xyce::Util::vsrcData > & vsrcs = nd.getVsrcs();
  std::unordered_map<std::string, Xyce::Util::bsrcData > & bsrcs = nd.getBsrcs();

  // functions:
  for ( auto it = functions.begin(); it != functions.end(); ++it )
  {
    it->second.exp_ = Teuchos::rcp(new Xyce::Util::newExpression(it->second.value_,grp));
    it->second.exp_->setFunctionArgStringVec ( it->second.args_ );
    if ( !(it->second.exp_->lexAndParseExpression ()) )
    {
      std::cout << "OOPS. problem parsing " << it->first << " = " << it->second.value_ <<std::endl;
      exit(0);
    }
  }

  ALLOCATE_LOOP(params)
  ALLOCATE_LOOP(global_params)
  ALLOCATE_LOOP(resistors)
  ALLOCATE_LOOP(capacitors)
  ALLOCATE_LOOP(vsrcs)
  ALLOCATE_LOOP(bsrcs)

  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "allocateExpressions timing: ";
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
}

#define RESOLVE_LOOP(NAME)   \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it ) \
  { \
    it->second.exp_->resolveExpression(); \
  }


#define RESOLVE_LOOP2(NAME)   \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it ) \
  { \
    it->second.exp_->setupVariousAstArrays(); \
  }

//-------------------------------------------------------------------------------
// Function      : resolveExpressions
//
// Purpose       : This function does more setup to the expressions (already
//                 allocated) in the netlist container.  It attempts to resolve
//                 unrecognized symbols.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 10/2?/2019
//-------------------------------------------------------------------------------
void resolveExpressions(Xyce::Util::netlistData & nd)
{

  // find unresolved params, globals and funcs.  Start with funcs
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;

  std::unordered_map<std::string, Xyce::Util::resistorData > & resistors = nd.getResistors();
  std::unordered_map<std::string, Xyce::Util::capacitorData > & capacitors = nd.getCapacitors();
  std::unordered_map<std::string, Xyce::Util::vsrcData > & vsrcs = nd.getVsrcs();
  std::unordered_map<std::string, Xyce::Util::bsrcData > & bsrcs = nd.getBsrcs();

  {
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  RESOLVE_LOOP(params)
  RESOLVE_LOOP(global_params)
  RESOLVE_LOOP(functions)

  RESOLVE_LOOP(resistors)
  RESOLVE_LOOP(capacitors)
  RESOLVE_LOOP(vsrcs)
  RESOLVE_LOOP(bsrcs)
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << std::fixed << std::setprecision(2) << "RESOLVE_LOOP timing CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms,  "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
  }

// new:
  {
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
#if 0
  // this is expensive, as it is currently written.  Need a better solution.
  RESOLVE_LOOP2(params)

  RESOLVE_LOOP2(global_params)

  RESOLVE_LOOP2(functions)

  RESOLVE_LOOP2(resistors)
  RESOLVE_LOOP2(capacitors)
  RESOLVE_LOOP2(vsrcs)
  RESOLVE_LOOP2(bsrcs)
#endif
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << std::fixed << std::setprecision(2) << "RESOLVE_LOOP2 timing CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
  }

}

#define DERIVATIVE_SETUP_LOOP(NAME)   \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it ) \
  { \
    it->second.exp_->setupDerivatives(); \
  }

//-------------------------------------------------------------------------------
// Function      : setupDerivatives
//
// Purpose       : This function does more setup to the expressions (already
//                 allocated) in the netlist container.  It attempts to resolve
//                 unrecognized symbols.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 10/2?/2019
//-------------------------------------------------------------------------------
void setupDerivatives (Xyce::Util::netlistData & nd)
{
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  // find unresolved params, globals and funcs.  Start with funcs
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;

  std::unordered_map<std::string, Xyce::Util::resistorData > & resistors = nd.getResistors();
  std::unordered_map<std::string, Xyce::Util::capacitorData > & capacitors = nd.getCapacitors();
  std::unordered_map<std::string, Xyce::Util::vsrcData > & vsrcs = nd.getVsrcs();
  std::unordered_map<std::string, Xyce::Util::bsrcData > & bsrcs = nd.getBsrcs();

#if 0
  DERIVATIVE_SETUP_LOOP(params)
  DERIVATIVE_SETUP_LOOP(global_params)
  DERIVATIVE_SETUP_LOOP(functions)

  DERIVATIVE_SETUP_LOOP(resistors)
  DERIVATIVE_SETUP_LOOP(capacitors)
  DERIVATIVE_SETUP_LOOP(vsrcs)
  DERIVATIVE_SETUP_LOOP(bsrcs)
#endif

  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "setupDerivatives timing: "; 
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
}

#define OUTPUT_LOOP(NAME,OUTNAME) \
  os << OUTNAME << std::endl; \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it ) \
  { \
    Xyce::Util::newExpression & newExp = *(it->second.exp_); \
    usedType result; \
    std::vector<usedType> derivs; \
    newExp.evaluateFunction(result); \
    newExp.evaluate(result,derivs); \
    it->second.output(); \
    newExp.dumpParseTree(os); \
 \
    os << it->first << ": result = " << result << " ["; \
    for (int ii=0;ii<derivs.size();ii++) { os << " " << derivs[ii]; } \
    os << " ]" <<std::endl; \
 \
    os << "---------------------------------------------------------------------------------" << std::endl; \
  }

#define TERSE_OUTPUT_LOOP(NAME,OUTNAME) \
  os << OUTNAME << std::endl; \
  for ( auto it = NAME.begin(); it != NAME.end(); ++it ) \
  { \
    Xyce::Util::newExpression & newExp = *(it->second.exp_); \
    usedType result; \
    std::vector<usedType> derivs; \
    newExp.evaluate(result,derivs); \
 \
    os << it->first << ": result = " << result << " ["; \
    for (int ii=0;ii<derivs.size();ii++) { os << " " << derivs[ii]; } \
    os << " ]" <<std::endl; \
 \
    os << "-------------------------------------------------------------------------------" << std::endl; \
  }

//-------------------------------------------------------------------------------
// Function      : outputExpressions
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/1/2019
//-------------------------------------------------------------------------------
void outputExpressions (Xyce::Util::netlistData & nd, std::ostream & os)
{

  std::unordered_map<std::string, Xyce::Util::nodeData > & nodes = nd.getNodes() ;
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;

  std::unordered_map<std::string, Xyce::Util::resistorData > & resistors = nd.getResistors();
  std::unordered_map<std::string, Xyce::Util::capacitorData > & capacitors = nd.getCapacitors();
  std::unordered_map<std::string, Xyce::Util::vsrcData > & vsrcs = nd.getVsrcs();
  std::unordered_map<std::string, Xyce::Util::bsrcData > & bsrcs = nd.getBsrcs();

  os << "-------------------------------------------------------------------------------" << std::endl;

  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  std::unordered_map<std::string, std::string> & options = nd.getOptions();
  if (options.find("verbose") != options.end())
  {
    if (!nodes.empty())
    {
      os << "nodes:" << std::endl;
      for ( auto it = nodes.begin(); it != nodes.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    OUTPUT_LOOP(params,"params:")
    OUTPUT_LOOP(global_params,"global_params:")
    OUTPUT_LOOP(functions,"functions:")
    OUTPUT_LOOP(resistors,"resistors:")
    OUTPUT_LOOP(capacitors,"capacitors:")
    OUTPUT_LOOP(vsrcs,"vsrcs:")
    OUTPUT_LOOP(bsrcs,"bsrcs:")
  }
  else
  {
    if (!nodes.empty())
    {
      os << "nodes:" << std::endl;
      for ( auto it = nodes.begin(); it != nodes.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    TERSE_OUTPUT_LOOP(global_params,"global_params:")
  }

  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "outputExpressions timing: ";
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
}

//-------------------------------------------------------------------------------
// Function      : outputExpressionCopies
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/1/2019
//-------------------------------------------------------------------------------
void outputExpressionCopies (Xyce::Util::netlistData & nd, std::ostream & os)
{

  std::unordered_map<std::string, Xyce::Util::nodeData > & nodes = nd.getNodes() ;
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;

  std::unordered_map<std::string, Xyce::Util::resistorData > & resistors = nd.getResistors();
  std::unordered_map<std::string, Xyce::Util::capacitorData > & capacitors = nd.getCapacitors();
  std::unordered_map<std::string, Xyce::Util::vsrcData > & vsrcs = nd.getVsrcs();
  std::unordered_map<std::string, Xyce::Util::bsrcData > & bsrcs = nd.getBsrcs();


  // make and evaluate the copies.  Do it twice to check the timings.
  {
  std::unordered_map<std::string, Xyce::Util::nodeData > nodesCopy = nodes;
  std::unordered_map<std::string, Xyce::Util::paramData > paramsCopy = params;
  std::unordered_map<std::string, Xyce::Util::globalParamData > global_paramsCopy = global_params;
  std::unordered_map<std::string, Xyce::Util::functionData > functionsCopy = functions;

  std::unordered_map<std::string, Xyce::Util::resistorData > resistorsCopy = resistors;
  std::unordered_map<std::string, Xyce::Util::capacitorData > capacitorsCopy = capacitors;
  std::unordered_map<std::string, Xyce::Util::vsrcData > vsrcsCopy = vsrcs;
  std::unordered_map<std::string, Xyce::Util::bsrcData > bsrcsCopy = bsrcs;

  os << "-------------------------------------------------------------------------------" << std::endl;

  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  std::unordered_map<std::string, std::string> & options = nd.getOptions();
  if (options.find("verbose") != options.end())
  {
    if (!nodesCopy.empty())
    {
      os << "nodesCopy:" << std::endl;
      for ( auto it = nodesCopy.begin(); it != nodesCopy.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    OUTPUT_LOOP(paramsCopy,"paramsCopy:")
    OUTPUT_LOOP(global_paramsCopy,"global_paramsCopy:")
    OUTPUT_LOOP(functionsCopy,"functionsCopy:")
    OUTPUT_LOOP(resistorsCopy,"resistorsCopy:")
    OUTPUT_LOOP(capacitorsCopy,"capacitorsCopy:")
    OUTPUT_LOOP(vsrcsCopy,"vsrcsCopy:")
    OUTPUT_LOOP(bsrcsCopy,"bsrcsCopy:")
  }
  else
  {
    if (!nodesCopy.empty())
    {
      os << "nodesCopy:" << std::endl;
      for ( auto it = nodesCopy.begin(); it != nodesCopy.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    TERSE_OUTPUT_LOOP(global_paramsCopy,"global_paramsCopy:")
  }

  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "outputExpressionCopies timing: ";
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
  }

  {
  std::unordered_map<std::string, Xyce::Util::nodeData > nodesCopy = nodes;
  std::unordered_map<std::string, Xyce::Util::paramData > paramsCopy = params;
  std::unordered_map<std::string, Xyce::Util::globalParamData > global_paramsCopy = global_params;
  std::unordered_map<std::string, Xyce::Util::functionData > functionsCopy = functions;

  std::unordered_map<std::string, Xyce::Util::resistorData > resistorsCopy = resistors;
  std::unordered_map<std::string, Xyce::Util::capacitorData > capacitorsCopy = capacitors;
  std::unordered_map<std::string, Xyce::Util::vsrcData > vsrcsCopy = vsrcs;
  std::unordered_map<std::string, Xyce::Util::bsrcData > bsrcsCopy = bsrcs;

  os << "-------------------------------------------------------------------------------" << std::endl;

  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  std::unordered_map<std::string, std::string> & options = nd.getOptions();
  if (options.find("verbose") != options.end())
  {
    if (!nodesCopy.empty())
    {
      os << "nodesCopy:" << std::endl;
      for ( auto it = nodesCopy.begin(); it != nodesCopy.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    OUTPUT_LOOP(paramsCopy,"paramsCopy:")
    OUTPUT_LOOP(global_paramsCopy,"global_paramsCopy:")
    OUTPUT_LOOP(functionsCopy,"functionsCopy:")
    OUTPUT_LOOP(resistorsCopy,"resistorsCopy:")
    OUTPUT_LOOP(capacitorsCopy,"capacitorsCopy:")
    OUTPUT_LOOP(vsrcsCopy,"vsrcsCopy:")
    OUTPUT_LOOP(bsrcsCopy,"bsrcsCopy:")
  }
  else
  {
    if (!nodesCopy.empty())
    {
      os << "nodesCopy:" << std::endl;
      for ( auto it = nodesCopy.begin(); it != nodesCopy.end(); ++it )
      {
        os << "key = " << it->first
           << "  value = " << it->second.value_
           << "  index = " << it->second.index_
           <<std::endl;
      }
      os << std::endl;
      os << "---------------------------------------------------------------------------------" << std::endl;
    }

    TERSE_OUTPUT_LOOP(global_paramsCopy,"global_paramsCopy:")
  }

  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "outputExpressionCopies2 timing: ";
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";
  }
}

template <typename ScalarT>
std::string basicTypeString()
{
  return std::string("");
}

template <>
std::string basicTypeString<double>()
{
  return std::string("double");
}

template <>
std::string basicTypeString<std::complex<double>>()
{
  return std::string("std::complex<double>");
}


template <typename ScalarT>
std::string assignmentTypeString()
{
  return std::string("");
}

template <>
std::string assignmentTypeString<double>()
{
  return std::string("");
}

template <>
std::string assignmentTypeString<std::complex<double>>()
{
  return std::string("std::complex<double>");
}

//-------------------------------------------------------------------------------
// Function      : codeGen
//
// Purpose       : This function emits a compilable C++ program, that can be
//                 used to compare to the results produced by the expression
//                 library.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/6/2019
//-------------------------------------------------------------------------------
void codeGen (Xyce::Util::netlistData & nd, std::ostream & os)
{
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd.getParams() ;
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd.getGlobalParams() ;
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd.getFunctions() ;

  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  os << "//---------------------------------------------------------------------------------\n";
  os << "// test program \n";
  os << "// File generated at " << std::ctime(&end_time);
  os << "//---------------------------------------------------------------------------------\n";

  os << "#include <iostream>\n";
  os << "#include <cmath>\n";
  os << "#include <complex>\n";
  os << "#include <string>\n";
  os << std::endl;

  os << "#include <Sacado_No_Kokkos.hpp>\n\n";
  os << "typedef Sacado::Fad::DFad<" << basicTypeString<usedType>() << "> fadType;\n\n";

  // ----------------------------------------------------------------------------
  // output class definition (which contains all the data)
  // ----------------------------------------------------------------------------
  os << "template <typename ScalarT>\n";
  os << "class netlistExpressions {\n\n";
  os << "public:\n";

  // ----------------------------------------------------------------------------
  // constructor
  // ----------------------------------------------------------------------------
  os << "  // constructor \n";
  os << "  netlistExpressions ()\n";
  os << "  {\n";
  os << "    // params:\n";
  for ( auto it = params.begin(); it != params.end(); ++it )
  {
    Xyce::Util::newExpression & newExp = *(it->second.exp_);
    os << "    " << it->first << " = ";
    newExp.codeGen(os);
  }

  os << std::endl;
  os << "    // global_params:\n";
  for ( auto it = global_params.begin(); it != global_params.end(); ++it )
  {
    Xyce::Util::newExpression & newExp = *(it->second.exp_);
    os << "    " << it->first << " = ";
    newExp.codeGen(os);
  }
  os << "  };\n\n";
  // ----------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  os << "  // params:\n";
  for ( auto it = params.begin(); it != params.end(); ++it )
  {
    os << "  ScalarT " << it->first << ";\n";
  }
  os << std::endl;
  os << "  // global_params:\n";
  for ( auto it = global_params.begin(); it != global_params.end(); ++it )
  {
    os << "  ScalarT " << it->first << ";\n";
  }
  os << std::endl;
  // ----------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  os << "  // function declarations:\n";
  for ( auto it = functions.begin(); it != functions.end(); ++it )
  {
    Xyce::Util::newExpression & newExp = *(it->second.exp_);
    os << "  ScalarT " << it->first << "(";
    for (int ii=0;ii<it->second.args_.size();ii++)
    {
      os << "ScalarT ";
      os << it->second.args_[ii];
      if (ii < (it->second.args_.size())-1 ) { os << ","; }
    }
    os << ") \n";
    os << "  {\n";
    os <<  "    return ";
    newExp.codeGen(os);
    os << "  }\n\n";
  }
  // ----------------------------------------------------------------------------
  os << "};\n" <<std::endl;

  // ----------------------------------------------------------------------------
  // output main function
  // ----------------------------------------------------------------------------
  os << "int main(int iargs, char *cargs[])\n";
  os << "{\n";

  os << "  netlistExpressions<fadType> ne;\n";

  os << std::endl << "  // param evaluations:\n";
  for ( auto it = params.begin(); it != params.end(); ++it )
  {
    os << "  {" <<std::endl;
    os << "    std::cout << \"" << it->first << " = \" << ne." << it->first << "<<std::endl;\n";
    os << "  }\n\n";
  }

  os << std::endl << "  // global param evaluations:\n";
  for ( auto it = global_params.begin(); it != global_params.end(); ++it )
  {
    os << "  {\n";
    os << "    std::cout << \"" << it->first << " = \" << ne." << it->first << "<<std::endl;\n";
    os << "  }\n\n";
  }

  os << std::endl << "  // function calls:\n";
  for ( auto it = functions.begin(); it != functions.end(); ++it )
  {
    int numArgs = it->second.args_.size();

    std::vector<usedType> args(numArgs);
    for (int ii=0;ii<numArgs;ii++)
    {
      args[ii] = ii; // need to use something ...
    }

    //  Each function call output is piped to std output.
    //  Prior to the function call, the
    os << "  {\n";
    os << "    std::cout << \" " << it->first << " = \" <<";

    os << " ne." << it->first << "(";

    for (int ii=0;ii<numArgs;ii++)
    {
      os << assignmentTypeString<usedType>();
      os << args[ii];
      if (ii < (it->second.args_.size())-1 ) { os << ","; }
    }

    os << ") << std::endl;\n";
    os << "  }\n\n";
  }

  os << std::endl;;
  os << "  return 0;\n";
  os << "}" <<std::endl;

}

//-------------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 10/2?/2019
//-------------------------------------------------------------------------------
int main (int iargs, char *cargs[])
{  
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  std::cout << "---------------------------------------------------------------------------------" << std::endl;

// Process netlist.  Note, this part of the prototype (netlist lexer/parser) is unlikely to be retained in Xyce
  Xyce::Util::netlistData nd;
  std::string netlistString("funcs.cir");

  if (iargs > 1) // assume that second argument is a netlist name
  {
    netlistString = std::string(cargs[1]);
  }

  std::cout << "Parsing netlist " << netlistString <<std::endl;

  Xyce::Util::lexAndParseNetlist(netlistString,nd);


  std::unordered_map<std::string, std::string> & options = nd.getOptions();
  if (options.find("verbose") != options.end())
  {
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "netlistData contains:" <<std::endl;
    nd.output();
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
  }

  Teuchos::RCP<Xyce::Util::prototypeExpressionGroup>  expGroup = Teuchos::rcp(new Xyce::Util::prototypeExpressionGroup (nd) );
  Teuchos::RCP<Xyce::Util::baseExpressionGroup>  baseGroup = expGroup;
  allocateExpressions(baseGroup,nd);

  // at this point, every object in the "nd" containers (params, global params, and functions)
  // have an expression set up for them.  Next step is to resolve them, for any that are unresolved.
  resolveExpressions(nd);
  setupDerivatives(nd);

  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "Final expression output" << std::endl;
  outputExpressions(nd,std::cout);
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "Final expression COPIES output" << std::endl;
  outputExpressionCopies(nd,std::cout);
  std::cout << "-------------------------------------------------------------------------------" << std::endl;

  if (options.find("codegen") != options.end())
  {
    // code gen:
    std::string cppFileName = netlistString + ".C";
    std::cout << "new filename = " << cppFileName <<std::endl;
    std::ofstream cppFile;
    cppFile.open(cppFileName);
    codeGen(nd,cppFile);
    cppFile.close();
  }


  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "Total timings: "<<std::endl;
  std::cout << std::fixed << std::setprecision(2) << "CPU: "
            << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms, "
            << "Wall clock: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n";

  return 0;
}
