//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        :  Low-level Expression support
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 06/07/01
//
//
//
//
//-----------------------------------------------------------------------------
/// Originally intended to support the behavioral modeling ("B") source
/// device, the expression package allows arbitrary arithmetic operations
/// to be computed during the run of a netlist.
///
/// The package's use has extended far beyond its original "B" source
/// intent, and is now used to print arbitrary expressions on .print lines,
/// as device model parameters to compute them from other constants, in
/// global function definition, in ".param" or ".global_param" statements, and
/// in many other contexts.
///

#ifndef N_UTL_ExpressionInternals_H
#define N_UTL_ExpressionInternals_H

// ---------- Standard Includes ----------
#include <vector>
#include <string>
#include <list>
#include <iosfwd>

#include <N_UTL_fwd.h>
#include <N_UTL_Interface_Enum_Types.h>
#include <N_UTL_ExpressionSymbolTable.h>

namespace Xyce {
namespace Util {


// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to false then AGAUSS() and GAUSS() will
// just return the mean rather than a random number.  The default is true.
extern bool enableRandomExpression;

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true then logical AND is &&, logical
// OR is || and ^ is a synonym for exponentiation.  The default is false.
extern bool useHspiceMath;

// This function is used to substitute the HSPICE logarithm functions
// for the Xyce ones.
void updateExpFuncArray();

// This function is used to substitute the "HSPICE math" strings for the
// EXPR_OR and EXPR_AND tokens for the Xyce ones.
void updateExpOpsArray();

//-----------------------------------------------------------------------------
// Class         : ExpressionNode
// Purpose       : Represent a node in an expression tree
// Special Notes :
// Scope         : public
// Creator       : David Shirley
// Creation Date : 6/7/2001
//-----------------------------------------------------------------------------
///
/// Expressions are converted to trees of ExpressionNode objects.
/// To evaluate an expression (or perform other operations such as
/// differentiation or substitution), the tree is walked and the operations
/// specified by the node performed.
///
class ExpressionNode
{
  public:
   /// Constructor
   ExpressionNode () :
     type(0),
     constant(0),
     valueIndex(0),
     funcname(""),
     eval_value(0),
     eval_num(0),
     funcnum(0),
     op(0)
   {
    operands.clear();
    state.clear();
    fptr.unary = NULL;
    fptr.binary = NULL;
   };

   /// Destructor
   ~ExpressionNode () {};

  public:
    int type;                                ///< From EXPR_OPS
    std::vector<ExpressionNode *> operands;  ///< list of pointers to operands
    double constant;                         ///< If type is EXPR_CONSTANT,
                                             ///< or if a  table entry then
                                             ///< use this as a cache for the
                                             ///< entry (which is
                                             ///< required to be constant)
                                             ///<
    std::vector<double> state;               ///< state data for operation
                                             ///< such as sdt or ddt
                                             ///<
    int valueIndex;                          ///< If EXPR_VAR, EXPR_ARGUMENT,
                                             ///< or (EXPR_FUNCTION and
                                             ///< EXPR_F_USER), use as
                                             ///< index into vars.  If
                                             ///< EXPR_FUNCTION and
                                             ///< funcnum == EXPR_F_TABLE then
                                             ///< this is used as a cache to
                                             ///< record the section that the
                                             ///< table last evaluated in
                                             ///<
    std::string funcname;                    ///< If EXPR_FUNCTION, name of
                                             ///< function,
                                             ///< If EXPR_VAR, name of var
                                             ///<
    double eval_value;                       ///< Evaluated value at current
                                             ///< eval_num
                                             ///<
    int eval_num;                            ///< Evaluation number
    int funcnum;                             ///< one of EXPR_F_FUNCS
    Op::Operator *        op;                ///< If type is EXPR_F_OP,
                                             ///< the operator that provides
                                             ///< the node's value
                                             ///<
    union
    {
      double (*unary)(double);
      double (*binary)(double, double);
      double (*op)(Op::Operator *);
    } fptr;                                 ///< pointer to the function to be
                                            ///< called if type is EXPR_FUNCTION
                                            ///<
} ;

//-----------------------------------------------------------------------------
// Class         : ExpressionElement
// Purpose       : Used in lexical analysis of input string
// Special Notes :
// Scope         : public
// Creator       : David Shirley
// Creation Date : 6/7/2001
//-----------------------------------------------------------------------------
///
/// The ExpressionElement class is used only in setting up an expression.
/// The input string is tokenized by the lexical analyzer ("tokenize_"),
/// and converted into a list of ExpressionElement objects (the tokens).
///
class ExpressionElement
{
  public:
    /// Constructor
    ExpressionElement () :
      token(0),
      type(0),
      name(""),
      number(0),
      node(NULL)
    {};

    /// Destructor
    ~ExpressionElement () {};

  public:
    int token;                 ///< see enum TOK_LIST
    int type;                  ///< see enum TYP_NODES
    std::string name;          ///< If the type is a string, the string
    double number;             ///< If the type is a number, the value
    ExpressionNode *node;      ///< If type is a node, pointer to the node.
};

// ---------- Enumerations ----------

enum EXPR_OPS
{
  EXPR_PLACEHOLDER,      //  0       Used during parsing.  This could be a
                         //          device or node name
                         //          as in V() or I() or it could just be a
                         //          string encountered in an expression which
                         //          will later be turned into a
                         //          variable or constant.
  EXPR_OR,               //  1       Binary functions are 1 - 15
  EXPR_XOR,              //  2
  EXPR_AND,              //  3
  EXPR_EQUAL,            //  4
  EXPR_NOTEQ,            //  5
  EXPR_GREAT,            //  6
  EXPR_GREATEQ,          //  7
  EXPR_LESS,             //  8
  EXPR_LESSEQ,           //  9
  EXPR_PLUS,             // 10
  EXPR_MINUS,            // 11
  EXPR_TIMES,            // 12
  EXPR_DIVIDE,           // 13
  EXPR_REMAINDER,        // 14
  EXPR_POWER,            // 15
  EXPR_FUNCTION,         // 16       Function, type according to EXPR_F_FUNCS
  EXPR_CONSTANT,         // 17       Constant
  EXPR_VAR               // 18       Variable, type according to EXPR_T_TYPES
};

enum EXPR_F_FUNCS                 // These differentiate between different types
                                  // of EXPR_FUNCTIONs
{
  EXPR_F_ABS,             //  0
  EXPR_F_ACOS,            //  1
  EXPR_F_ACOSH,           //  2
  EXPR_F_ASIN,            //  3
  EXPR_F_ASINH,           //  4
  EXPR_F_ATAN,            //  5
  EXPR_F_ATANH,           //  6
  EXPR_F_COS,             //  7
  EXPR_F_COSH,            //  8
  EXPR_F_DDT,             //  9
  EXPR_F_DDX,             // 10
  EXPR_F_EXP,             // 11
  EXPR_F_IF,              // 12
  EXPR_F_LN,              // 13
  EXPR_F_LOG,             // 14
  EXPR_F_NOT,             // 15
  EXPR_F_RAND,            // 16
  EXPR_F_SDT,             // 17
  EXPR_F_SGN,             // 18
  EXPR_F_SIN,             // 19
  EXPR_F_SINH,            // 20
  EXPR_F_SQRT,            // 21
  EXPR_F_TABLE,           // 22
  EXPR_F_F_TABLE,         // 23
  EXPR_F_R_TABLE,         // 24
  EXPR_F_TAN,             // 25
  EXPR_F_TANH,            // 26
  EXPR_F_UMINUS,          // 27
  EXPR_F_URAMP,           // 28
  EXPR_F_USER,            // 29
  EXPR_F_AGAUSS,          // 30
  EXPR_F_GAUSS,           // 31
  EXPR_F_INT,             // 32
  EXPR_F_SCHEDULE,        // 33
  EXPR_F_OP,              // 34
  EXPR_F_FLOOR,           // 35
  EXPR_F_CEIL             // 36
};

enum EXPR_T_TYPES                 // These differentiate between different types
                                  // of EXPR_VARs
{
  EXPR_T_NODE=10,         // 10
  EXPR_T_STRING,          // 11
  EXPR_T_INSTANCE,        // 12
  EXPR_T_SPECIAL,         // 13
  EXPR_T_VARIABLE,        // 14
  EXPR_T_FUNCTION,         // 15
  EXPR_T_NODAL_COMPUTATION // 16
};

enum TOK_LIST
{
  TOK_END,             //  0
  TOK_OR,              //  1
  TOK_XOR,             //  2
  TOK_AND,             //  3
  TOK_EQUAL,           //  4
  TOK_NOTEQ,           //  5
  TOK_GREAT,           //  6
  TOK_GREATEQ,         //  7
  TOK_LESS,            //  8
  TOK_LESSEQ,          //  9
  TOK_PLUS,            // 10
  TOK_MINUS,           // 11
  TOK_TIMES,           // 12
  TOK_DIVIDE,          // 13
  TOK_REMAINDER,       // 14
  TOK_POWER,           // 15
  TOK_NOT,             // 16
  TOK_UMINUS,          // 17
  TOK_LPAREN,          // 18
  TOK_RPAREN,          // 19
  TOK_VALUE,           // 20
  TOK_COMMA,           // 21
  TOK_QUESTION,        // 22
  TOK_COLON,           // 23
  TOK_LBRACE,          // 24
  TOK_RBRACE,          // 25
  TOK_SPACE            // 26
};

enum TYP_NODES        // These are for expression element types
{
  TYP_NUM,             //  0
  TYP_STRING,          //  1
  TYP_PNODE            //  2
};

enum EXPR_WARNINGS   // These indicate improper operations which are likely
                     // inconsequential
{
  EXPR_WARNING = 100,    // 100
  EXPR_OR_WARNING,       // 101
  EXPR_XOR_WARNING,      // 102
  EXPR_AND_WARNING,      // 103
  EXPR_SUM_WARNING,      // 104
  EXPR_DIFF_WARNING,     // 105
  EXPR_TIMES_WARNING,    // 106
  EXPR_DIVIDE_WARNING,   // 107
  EXPR_POWER_WARNING,    // 108
  EXPR_NOT_WARNING       // 109
};

enum EXPR_ERRORS   // These indicate possible loss of precision
{
  EXPR_ERROR = 200,      // 200
  EXPR_REMAINDER_ERROR,  // 201
  EXPR_POWER_ERROR,      // 202
  EXPR_COSH_ERROR,       // 203
  EXPR_EXP_ERROR,        // 204
  EXPR_SINH_ERROR,       // 205
  EXPR_SQRT_ERROR       // 206
};

enum EXPR_FATALS   // These indicate ill-defined operations, and wrong answers
{
  EXPR_FATAL = 300,      // 300
  EXPR_DIVIDE_FATAL,     // 301
  EXPR_REMAINDER_FATAL,  // 302
  EXPR_LOG_FATAL,        // 303
  EXPR_TABLE_FATAL,      // 304
  EXPR_NODERIV_FATAL,    // 305
  EXPR_ARGUMENT_FATAL,   // 306
  EXPR_FUNCTION_FATAL    // 307
};

// Some miscellaneous definitions:

#define EXPR_HUGE (1.e+50)
#define N_INT_STATE     6

double EXPRor(double arg1, double arg2);
double EXPRxor(double arg1, double arg2);
double EXPRand(double arg1, double arg2);
double EXPRequal(double arg1, double arg2);
double EXPRnoteq(double arg1, double arg2);
double EXPRgreat(double arg1, double arg2);
double EXPRgreateq(double arg1, double arg2);
double EXPRless(double arg1, double arg2);
double EXPRlesseq(double arg1, double arg2);
double EXPRplus(double arg1, double arg2);
double EXPRminus(double arg1, double arg2);
double EXPRtimes(double arg1, double arg2);
double EXPRdivide(double arg1, double arg2);
double EXPRremainder(double arg1, double arg2);
double EXPRpower(double arg1, double arg2);

double EXPRabs (double arg);
double EXPRacos (double arg);
double EXPRacosh (double arg);
double EXPRasin (double arg);
double EXPRasinh (double arg);
double EXPRatan (double arg);
double EXPRatanh (double arg);
double EXPRcos (double arg);
double EXPRcosh (double arg);
double EXPRexp (double arg);
double EXPRln (double arg);
double EXPRlog (double arg);
double EXPRnot (double arg);
double EXPRsgn (double arg);
double EXPRsin (double arg);
double EXPRsinh (double arg);
double EXPRsqrt (double arg);
double EXPRtan (double arg);
double EXPRtanh (double arg);
double EXPRuminus (double arg);
double EXPRuramp (double arg);
double EXPRint (double arg);
double EXPRfloor (double arg);
double EXPRceil (double arg);

//-----------------------------------------------------------------------------
// Class         : ExpressionInternals
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/20/08
//-----------------------------------------------------------------------------
///
/// The ExpressionInternals class is the core of the expression system,
/// and provides all the operations needed to create, manipulate, and evaluate
/// expressions.
///
class ExpressionInternals
{

public:

  /// Constructor of expression from a string
  ExpressionInternals (std::string const & exp = std::string());
  /// Constructor of expression from an ExpressionNode, a previously constructed expression tree
  ExpressionInternals( const ExpressionNode &node );
  /// Copy constructor
  ExpressionInternals (const ExpressionInternals & right);
  //  ExpressionInternals& operator=(const ExpressionInternals& right) ;
  /// Destructor
  ~ExpressionInternals (void);

  bool parsed() const
  {
    return parsed_;
  }

  bool set (std::string const & exp);
  void getSymbolTable (std::vector< ExpressionSymbolTableEntry > & names) const;
  void get_names (int type, std::vector< std::string > & names);
  int get_type (std::string const & var);
  bool make_constant (std::string const & var, double const & val);
  bool make_var (std::string const & var);

  int differentiate();

  bool set_var (const std::string &, const double &);
  bool set_vars (const std::vector< double > &);

  std::string get_expression (void) const;
  std::string get_derivative(std::string const & var);
  int get_num(int type);

  int evaluate (double &result, std::vector< double > &derivs, std::vector< double > &vals);
  int evaluateFunction (double &result, std::vector< double > &vals);

  int evaluate (double &result, std::vector< double > &derivs);
  int evaluateFunction (double &result);

  bool set_sim_time (double time);
  bool set_sim_dt (double dt);
  bool set_temp (double const & temp);
  bool set_sim_freq (double freq);
  void set_accepted_time (double const time);
  double get_break_time (void) {if(time_index == -2) return 0; else return(get_break_time_i());};
  double get_break_time_i (void);
  const std::string & get_input (void);
  int order_names (std::vector< std::string > const & new_names);
  int replace_func (std::string const & func_name,
                    ExpressionInternals & func_def, int numArgs);
  int replace_var (const std::string & var_name, const ExpressionInternals & subexpr);
  int replace_var (const std::string & var_name, Op::Operator *op);
  bool replace_name (const std::string & old_name, const std::string & new_name);
  int getNumDdt();
  void getDdtVals (std::vector<double> &);
  void setDdtDerivs (std::vector<double> &);
  const ExpressionNode *get_tree() const {return tree_;}
  ExpressionNode *get_tree() {return tree_;}
  void set_tree(ExpressionNode *new_root) {tree_ = new_root;}
  int num_vars () {return num_N_+num_I_+num_lead_+num_string_+num_special_+num_var_+num_func_+num_node_computation_;}
  bool isTimeDepedent() const {return timeDependent_;}
  bool isRandomDepedent() const {return randomDependent_;}
  void dumpParseTree();

  static void seedRandom(long seed);

private:

  std::string Input_;            ///< The input string, if this expression
                                 ///< was constructed from a string.
                                 ///<
  bool parsed_;                  ///< true if expression was successfully parsed
  bool differentiated_;          ///< true if expression has been differentiated
                                 ///< with respect to all variables
                                 ///<
  bool ddxProcessed_;            ///< true if expression has been differentiated
                                 ///< with respect to a specific variable
                                 ///<
  int num_N_;                    ///< Number of nodal variables we depend on
  int num_I_;                    ///< Number of instances we depend on,
                                 ///< e.g. lead or branch currents
                                 ///<
  int num_lead_;                 ///< number of lead currents we depend on
  int num_string_;               ///< number of strings in expression
  int num_special_;              ///< number of special variables
  int num_var_;                  ///< number of variables
  int num_func_;                 ///< number of functions
  int num_node_computation_;     ///< number of "nodal computation" terms
                                 ///< such as "IM(a,b)" or "IP(a,b)"
                                 ///<
  double sim_time_;              ///< Current simpulation time
  double sim_dt_;                ///< Current simpulation time step
  int time_index;                ///< index of TIME variable in symbol table
  bool timeDependent_;           ///< true if expression contains DDT or SDT
  bool randomDependent_;         ///< true if expression contains GAUSS, AGAUSS or RAND
  bool breakpointed_;            ///< true if breakpoints have been computed

  int numVars_;                  ///< total number of external variables
                                 ///< needed for evaluation of expression
                                 ///<

  std::vector<int> varTypes_;    ///< array of types of variables
  std::vector<std::string> varValues_;   ///< array of values of variables
  std::string leadDesignator_;   ///< lead designator for current variables
  std::vector<double> var_vals_;  ///< Values of variables, nodes, instances

  ExpressionNode *tree_;                      ///< The expression tree
  std::vector<ExpressionNode *> derivs_;      ///< The derivative parse trees

  std::vector<ExpressionNode *> free_list_;   ///< List of ExpressionNodes to be freed
                                              ///< in destructor
                                              ///<
  std::vector<ExpressionNode *> done_list_;   ///< list used to keep track
                                              ///< of nodes to free
  std::vector<ExpressionElement *> ee_list_;  ///< List of ExpressionElements to be freed
                                              ///< after parsing complete
  std::vector<ExpressionNode *> breaks_;      ///< Expression for computing
                                              ///< breakpoints
                                              ///<

  ExpressionNode *PThead_;                    ///<A second pointer to the
                                              ///<head of parse tree?
                                              ///<

  int Rmode_;
  int ind_replace_;
  std::string Rstring_;
  double Rcval_;
  int curr_magic_;
  int curr_num_;
  bool values_changed_;

  // functions
  int find_num_ (const std::string &);

  void set_nums_ ();
  void create_vars_ ();
  void copy_elements_ (std::list<ExpressionElement *> &to, std::list<ExpressionElement *> *from);
  void copy_element_ (ExpressionElement *to, ExpressionElement *from);
  ExpressionNode *makepnode_ (ExpressionElement *elem);
  ExpressionNode *mkfnode_ (const std::string & fname, int num_args, std::vector<ExpressionNode *> args);
  ExpressionNode *mkfnode_ (const std::string & fname, int num_args, ExpressionNode *n);
  ExpressionNode *mksnode_ (const std::string & name);
  ExpressionNode *mkcon_ (double value);
  ExpressionNode *mkb_ (int type, ExpressionNode *left, ExpressionNode *right);
  ExpressionNode *mkf_(int type, ExpressionNode *arg);
  ExpressionNode *newExpressionNode_ ();
  void deleteExpressionNode_ (ExpressionNode *p);
  ExpressionElement *newExpressionElement_ ();
  void RpTree_ (const ExpressionNode * pt, std::ostringstream & s) const;
  std::string varStr_ (int i) const;
  ExpressionNode * PTcheck_(ExpressionNode *p);
  ExpressionNode * diffDDX_(ExpressionNode *p);
  ExpressionNode * PTdiffDDX_(ExpressionNode *p);
  ExpressionNode * com_expr_ (ExpressionNode *c, ExpressionNode *p);
  ExpressionNode * Differentiate_ (ExpressionNode *arg, int varnum);

  void Rconvert_ (ExpressionNode & node);
  void RcountDDT_ (ExpressionNode & node);
  void RgetDDT_ (ExpressionNode & node, std::vector<double> & vals);
  void RsetDDT_ (ExpressionNode & node, std::vector<double> & vals);
  void convert_to_constant_ (int i, double c_value);
  void convert_to_variable_ (int i);

  // method to support copy constructor:
  ExpressionNode * copy_exprNode_ (ExpressionNode *n);

  // Methods to support breakpoints:
  void simplify_ (ExpressionNode & node);
  void get_breaks_ (ExpressionNode & node);
  bool dependent_ (ExpressionNode & node, int ind);
  bool dependent_other_ (ExpressionNode & node, int ind);
  bool arithmatic_ (ExpressionNode & node);

  // Methods to suupport order_names:
  void Rmap_ (ExpressionNode & node, int mode, std::vector<int> &nmap);
  int EXPRaddDummyString_ (std::string & dummy);

  // Methods to support replace_func:
  void addNode_ (ExpressionNode *n, int ind, const ExpressionNode *f,
                 const ExpressionInternals & func_expr,
                 int na_func,
                 std::vector<ExpressionNode *> operands);
  void Nreplace_ (ExpressionNode *n, const ExpressionNode *f, const ExpressionInternals & func_expr, int na_func, std::vector<ExpressionNode *> operands);
  int Freplace_ (ExpressionNode *n, std::string const & func_name,
                 ExpressionInternals & func_expr,
                 int na_func);
  void RemoveFentry_ (ExpressionNode *n, int new_ind, int old_ind);

  // method to support replace_var
  int Vreplace_ (ExpressionNode *n, const std::string &varName, const ExpressionInternals & subexpr);

  // Methods to support evaluate:
  void EXPReval_ (ExpressionNode & node, double & res, std::vector<double> &vals);
  void clear_eval_num_ (ExpressionNode *n);

  // Miscellaneous utility routines
  void compactLine_(std::string & inputLine, std::string & compactedLine);
  void tokenize_(std::string & inputLine, std::list<ExpressionElement *> & tokenList);
  bool convertPolyToExpr_(std::list<ExpressionElement *> & tokenList);
  void standardizeTable_(std::list<ExpressionElement *> & tokenList);
  void braceToParen_(std::list<ExpressionElement *> & tokenList);
  // Debugging methods

  void dumpParseTree_(const ExpressionNode *tree, int indentLevel=0);
  void indentWithDashes_(int level);

  bool isNameSpecial_(const std::string &name) const;

  // Catch attempts to use operator= at compile time, by
  // defining this null-op as a private member.
  // This doesn't work because parameters are often put inside of
  // STL objects, which rely on operator=.
  //ExpressionInternals& operator=(const ExpressionInternals& right) ;

};

} // namespace Util
} // namespace Xyce

#endif // N_UTL_EXPRESSION_H
