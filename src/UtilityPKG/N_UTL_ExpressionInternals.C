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

//-------------------------------------------------------------------------
//
// Purpose       : Low-level Expression support
//
// Special Notes : 
//
// Creator       : Dave Shirley, PSSI, SNL
//
// Creation Date : 06/07/01
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>
#include <N_DEV_Const.h>

// ---------- Standard Includes ----------

#include <iterator>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionInternals.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_Op.h>
#include <N_UTL_RandomNumbers.h>

namespace Xyce {
namespace Util {

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to false then AGAUSS() and GAUSS() will
// just return the mean rather than a random number.
bool enableRandomExpression = true;

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true then logical AND is &&, logical
// OR is || and ^ is a synonym for exponentiation.
bool useHspiceMath = false;

namespace {

double op_eval(Op::Operator *op)
{
  return (*op)(MPI_COMM_NULL, Util::Op::OpData()).real();
}

}

// Operator precedence table:
enum
{
  G=1,      // Greater
  L,        // Lesser
  E,        // Equal
  R         // Error
};

// This constant is the wrap around for eval_num.  I doesn't matter
// what is is as long as it is used consistently
enum { MAX_EVAL = 1000000 };

// Note: Order of this table must be consistent with TOK_LIST order
// Also, may need to have to add EXPR_F_FUNCS function and/or EXPR_OPS
// in include file
// Also, may have to add to ops[] and other changes in this file.
// Also,, don't forget op_name below.
// Also, add to EXPR* functions, if needed

static char prectable[24][24] = {
//         $  |  ^  &  == != >  >= <  <= +  -  *  /  %  ** ~  u- (  )  v  ,  ?  :
/* $  */ { R, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, R, L, R, L, L},
/* |  */ { G, G, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* ^  */ { G, G, G, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* &  */ { G, G, G, G, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* == */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* != */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* >  */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* >= */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* <  */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* <= */ { G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, L, L, G, L, G, G, G},
/* +  */ { G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, G, L, G, G, G},
/* -  */ { G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, L, L, L, L, G, L, G, G, G},
/* *  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, L, G, L, G, G, G},
/* /  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, L, G, L, G, G, G},
/* %  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, L, G, L, G, G, G},
/* ** */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, L, L, G, L, G, G, G},
/* ~  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, G, L, G, G, G},
/* u- */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, L, G, G, L, G, L, G, G, G},
/* (  */ { R, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, E, L, E, L, R},
/* )  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, R, G, R, G, G, G},
/* v  */ { G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, R, G, G, G},
/* ,  */ { G, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, E, L, E, L, L},
/* ?  */ { R, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, G, L, G, L, L},
/* :  */ { G, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, L, G, L, G, L, G}
} ;

static const char *op_name[24] = {"$","|","^","&","==","!=",">",">=","<","<=","+",
                       "-","*","/","%","**","~","u-","(",")","v",",","?",":"};

static const char *specials = {"{}()+-*/%, <>=!&|^~?:"};
static const char *specialsNoColon = {"{}()+-*/%, <>=!&|^~?"};
// These characters are "special" with respect to determining if we are
// parsing inside of an operator like V() 
static const char *operatorSpecials = {"{}(), "};
// These characters can appear inside of a number 
const std::string numberChars(".0123456789");

const char *toks[27] = {"END", "OR", "XOR", "AND", "EQUAL",
                        "NOTEQ", "GREAT", "GREATEQ", "LESS", "LESSEQ",
                        "PLUS", "MINUS", "TIMES", "DIVIDE", "REMAINDER",
                        "POWER", "NOT", "UMINUS", "LPAREN", "RPAREN",
                        "VALUE", "COMMA", "QUESTION", "COLON",
                        "LBRACE", "RBRACE", "SPACE"};
const  char *typs[3] = {"NUM", "STRING", "PNODE"};
const  char *expr_ops[22] = {"PLACEHOLDER", "OR", "XOR", "AND", "EQUAL",
         "NOTEQ", "GREAT", "GREATEQ", "LESS", "LESSEQ",
         "PLUS", "MINUS", "TIMES", "DIVIDE", "REMAINDER",
         "POWER", "FUNCTION", "CONSTANT", "VAR", "ARGUMENT",
         "TOAST"};

// Binary operation table:

struct op
{
    int number;
    bool commutative;
    std::string name;
    double (*funcptr)(double, double);
};

// There binary operations MUST be in the same order as EXPR_OPS in the
// include file.  This struct is used in RpTree_ to help map a parsed
// expression back into an "expression string".
static op ops[] =
{
  { EXPR_PLACEHOLDER, false, "",   NULL } ,
  { EXPR_OR,          true,  "|",  EXPRor } ,
  { EXPR_XOR,         true,  "^",  EXPRxor } ,
  { EXPR_AND,         true,  "&",  EXPRand } ,
  { EXPR_EQUAL,       true,  "==", EXPRequal } ,
  { EXPR_NOTEQ,       true,  "!=", EXPRnoteq } ,
  { EXPR_GREAT,       false, ">",  EXPRgreat } ,
  { EXPR_GREATEQ,     false, ">=", EXPRgreateq } ,
  { EXPR_LESS,        false, "<",  EXPRless } ,
  { EXPR_LESSEQ,      false, "<=", EXPRlesseq } ,
  { EXPR_PLUS,        true,  "+",  EXPRplus } ,
  { EXPR_MINUS,       false, "-",  EXPRminus } ,
  { EXPR_TIMES,       true,  "*",  EXPRtimes } ,
  { EXPR_DIVIDE,      false, "/",  EXPRdivide } ,
  { EXPR_REMAINDER,   false, "%",  EXPRremainder } ,
  { EXPR_POWER,       false, "**", EXPRpower }
};

#define NUM_OPS (sizeof (ops) / sizeof (struct op))
// This function is used to substitute the "HSPICE math" strings for the
// EXPR_OR and EXPR_AND tokens for the Xyce ones.
void updateExpOpsArray()
{
  for (int i = 0; i < NUM_OPS; ++i)
  {
    if (ops[i].name == "|") { ops[i].name = "||";}
    if (ops[i].name == "&") { ops[i].name = "&&";}
  }
  return;
}

#define BINARY(A) ((A)>=EXPR_OR && (A)<=EXPR_POWER)
#define MODULUS(NUM,LIMIT) ((NUM) - (static_cast<long> ((NUM) / (LIMIT))) * (LIMIT))
#define EXPR_ERROR(A) EXPRerrno = ((A>EXPRerrno) ? A : EXPRerrno)
#define LOGIC(Input,Error) if(Input != 0. && Input != 1.) EXPR_ERROR(Error)

#define PRECISION 15        // digits of precision when coding constants from
                            // expression into a string

// Function table:

static struct func
{
  std::string name;
  int number;
  int num_args;
  double (*funcptr)(double);
} funcs[] =
{
  { "ABS",         EXPR_F_ABS,      1, EXPRabs } ,
  { "ACOS",        EXPR_F_ACOS,     1, EXPRacos } ,
  { "ACOSH",       EXPR_F_ACOSH,    1, EXPRacosh } ,
  { "ASIN",        EXPR_F_ASIN,     1, EXPRasin } ,
  { "ASINH",       EXPR_F_ASINH,    1, EXPRasinh } ,
  { "ATAN",        EXPR_F_ATAN,     1, EXPRatan } ,
  { "ATANH",       EXPR_F_ATANH,    1, EXPRatanh } ,
  { "COS",         EXPR_F_COS,      1, EXPRcos } ,
  { "COSH",        EXPR_F_COSH,     1, EXPRcosh } ,
  { "DDT",         EXPR_F_DDT,      1, NULL } ,
  { "DDX",         EXPR_F_DDX,      2, NULL } ,
  { "EXP",         EXPR_F_EXP,      1, EXPRexp } ,
  { "IF",          EXPR_F_IF,       3, NULL } ,
  { "LN",          EXPR_F_LN,       1, EXPRln } ,
  { "LOG",         EXPR_F_LOG,      1, EXPRlog } ,
  { "~",           EXPR_F_NOT,      1, EXPRnot } ,
  { "RAND",        EXPR_F_RAND,     0, NULL } ,
  { "SDT",         EXPR_F_SDT,      1, NULL } ,
  { "SGN",         EXPR_F_SGN,      1, EXPRsgn } ,
  { "SIN",         EXPR_F_SIN,      1, EXPRsin } ,
  { "SINH",        EXPR_F_SINH,     1, EXPRsinh } ,
  { "SQRT",        EXPR_F_SQRT,     1, EXPRsqrt } ,
  { "TABLE",       EXPR_F_TABLE,    3, NULL } ,
  { "F~TABLE",     EXPR_F_F_TABLE,  1, NULL } ,
  { "R~TABLE",     EXPR_F_R_TABLE,  1, NULL } ,
  { "TAN",         EXPR_F_TAN,      1, EXPRtan } ,
  { "TANH",        EXPR_F_TANH,     1, EXPRtanh } ,
  { "-",           EXPR_F_UMINUS,   1, EXPRuminus } ,
  { "URAMP",       EXPR_F_URAMP,    1, EXPRuramp },
  { "AGAUSS",      EXPR_F_AGAUSS,   3, NULL },
  { "GAUSS",       EXPR_F_GAUSS,    3, NULL },
  { "INT",         EXPR_F_INT,      1, EXPRint },
  { "SCHEDULE",    EXPR_F_SCHEDULE, 2, NULL },
  { "_OP",         EXPR_F_OP,       0, NULL },
  { "FLOOR",       EXPR_F_FLOOR,    1, EXPRfloor },
  { "CEIL",        EXPR_F_CEIL,     1, EXPRceil }
} ;

#define NUM_FUNCS (sizeof (funcs) / sizeof (struct func))

// This function is used to substitute the HSPICE logarithm functions
// for the Xyce ones.
void updateExpFuncArray()
{
  for (int i = 0; i < NUM_FUNCS; ++i)
  {
    if (funcs[i].name == "LOG") { funcs[i].name = "LOG10";}
    if (funcs[i].name == "LN")  { funcs[i].name = "LOG";}
  }
  return;
}

static struct constant {
    const char *name;
    double value;
} constants[] = {
    { "EXP",  exp(1.0) },
    { "PI", M_PI }
} ;

#define NUM_CONSTANTS (sizeof (constants) / sizeof (struct constant))

static std::string specSigs[] = {"TIME", "TEMP", "VT", "FREQ"};
static double Epsilon = 1.e-12;

// static error number for evaluate:
static int EXPRerrno;

static std::string msg;

static double tentative_accepted_time, accepted_time;
static int numDDT;

static Xyce::Util::RandomNumbers *theRandomNumberGenerator=0;


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::ExpressionInternals
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Construct an expression object from an input string
///
/// @param[in] exp  The expression string
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
ExpressionInternals::ExpressionInternals( const std::string & exp )
  : parsed_           (false),
    differentiated_   ( false ),
    ddxProcessed_     (false),
    num_N_            ( 0 ),
    num_I_            ( 0 ),
    num_lead_         ( 0 ),
    num_string_       ( 0 ),
    num_special_      ( 0 ),
    num_var_          ( 0 ),
    num_func_         ( 0 ),
    sim_time_         ( 0.0 ),
    sim_dt_           ( 1.0e-10 ),
    time_index        ( -1 ),
    timeDependent_    ( false ),
    randomDependent_  ( false ),
    breakpointed_     ( false ),
    numVars_          ( 0 )
{
  Input_ = "";
  if (!exp.empty())
    parsed_ = set(exp);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::ExpressionInternals
// Purpose       : Copy Constructor
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/05/04
//-----------------------------------------------------------------------------
///
/// Construct a copy of an expression object (deep copy)
///
/// @param[in] right  The expression object to copy
///
/// @author Dave Shirley, PSSI
/// @date 11/05/04
///
ExpressionInternals::ExpressionInternals(
  const ExpressionInternals & right)
  : Input_             (right.Input_),
    differentiated_    (false),
    ddxProcessed_      (right.ddxProcessed_),
    num_N_             (right.num_N_),
    num_I_             (right.num_I_),
    num_lead_          (right.num_lead_),
    num_string_        (right.num_string_),
    num_special_       (right.num_special_),
    num_var_           (right.num_var_),
    num_func_          (right.num_func_),
    sim_time_          (right.sim_time_),
    sim_dt_            (right.sim_dt_),
    time_index         (right.time_index),
    timeDependent_     (right.timeDependent_),
    randomDependent_   (right.randomDependent_),
    breakpointed_      (right.breakpointed_),
    numVars_           (right.numVars_),
    varTypes_           (right.varTypes_),
    varValues_          (right.varValues_),
    leadDesignator_     (right.leadDesignator_),
    var_vals_          (right.var_vals_),
  curr_magic_        (right.curr_magic_)
{
  std::vector<ExpressionNode *>::iterator free_i;

  done_list_.clear();
  ee_list_.clear();
  breaks_.clear();
  PThead_ = NULL;
  Rmode_ = 0;
  ind_replace_ = 0;
  Rstring_ = "";
  Rcval_ = 0;
  curr_num_ = 1;
  values_changed_ = true;

  tree_ = copy_exprNode_(right.tree_);

  PThead_ = tree_;
  done_list_.clear();

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Calling PTcheck_ from copy constructor" << std::endl;

  if (tree_ != PTcheck_ (PThead_))
  {
    Report::DevelFatal().in("ExpressionInternals::ExpressionInternals") << "ExpressionInternals::ExpressionInternals: Internal error in copy constructor";
  }
  for (free_i = done_list_.begin() ; free_i != done_list_.end() ; ++free_i)
    deleteExpressionNode_(*free_i);

  if (right.differentiated_)
    differentiate();

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::ExpressionInternals
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 02/23/2015
//-----------------------------------------------------------------------------
///
/// Construct an expression object from an expression node
///
/// @param[in] node  The node from which to create the expression
///
/// This method is primarily used to create a subexpression object in the
/// special case where a variable such as GMIN or VT is being replaced with
/// an operator to retrieve those global variables.
///
/// @author Dave Baur, Raytheon
/// @date 02/23/2015
///
ExpressionInternals::ExpressionInternals( const ExpressionNode &node )
  : parsed_           (false),
    differentiated_   ( false ),
    ddxProcessed_     (false),
    num_N_            ( 0 ),
    num_I_            ( 0 ),
    num_lead_         ( 0 ),
    num_string_       ( 0 ),
    num_special_      ( 0 ),
    num_var_          ( 0 ),
    num_func_         ( 0 ),
    sim_time_         ( 0.0 ),
    sim_dt_           ( 1.0e-10 ),
    time_index        ( -1 ),
    timeDependent_    ( false ),
    randomDependent_  ( false ),
    breakpointed_     ( false ),
    numVars_          ( 0 ),
    tree_(new ExpressionNode(node))
{}


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::~ExpressionInternals
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Destroy this expression object.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
ExpressionInternals::~ExpressionInternals ()
{
  std::vector<ExpressionNode *>::iterator free_i;

  for (free_i = free_list_.begin() ; free_i != free_list_.end() ; ++free_i)
  {
    delete *free_i;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set
// Purpose       : Set the value of the expression to a string
// Special Notes : A lot of this routine is in fact based directly on the
//                 stuff done by the SPICE 3F5 routine "inpptree."  That
//                 accounts for a lot of the unreadability.
//                 Further Issues are introduced because far more
//                 functionality was inserted into the parsing than SPICE
//                 ever had.
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/05
//-----------------------------------------------------------------------------
///
/// Parse an expression string and set up internal datastructures of ExpressionInternals
///
/// @param[in] exp  The expression string to parse
/// @return  true on success
///
/// This routine is called by the ExpressionInternals constructor.  At
/// the conclusion of this routine, if the parsing has been
/// successful, the ExpressionInternals object is ready to use.
///
/// Operation of the function is as follows:
///    - The input is scanned for matching parens and extraneous white space is removed
///    - The input is tokenized, creating a list of ExpressionElement objects
///    - The PSpice-style "POLY" format is converted into normal polynomial form
///    - TABLE and SCHEDULE functions are cleaned up and put into standard form
///    - Special processing is done to handle YPDE devices
///    - An expression tree is generated from the token stack
///    - Common subexpressions are cleaned up
///    - ddx functions are processed
///    - The ExpressionElement token list is deleted
///
/// @author Dave Shirley, PSSI
/// @date 03/10/05
///
bool ExpressionInternals::set ( const std::string & exp )
{
  std::string line, compact_line, val_string;
  int last_token, len, i;
  std::string::size_type j, k;
  std::string token;
  ExpressionElement *el, *next, *top;
  std::vector<ExpressionElement *>::iterator ee_i;
  std::list<ExpressionElement *> tok;
  std::list<ExpressionElement *>::iterator tok_i, tok_j;
  std::vector<ExpressionNode *>::iterator free_i;
  std::vector<ExpressionElement *> stack;
  std::vector<ExpressionNode *> node_list;
  ExpressionNode *pn, *lpn, *rpn;
  int sp, st, num_nodes;

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Original input expression: " << exp << std::endl;

  if (!Input_.empty())
  {
    std::vector<ExpressionNode *>::iterator free_i;

    for (free_i = free_list_.begin() ; free_i != free_list_.end() ; ++free_i)
      delete *free_i;

    differentiated_ = false;
    ddxProcessed_ = false;
    breakpointed_ = false;
    time_index = -1;
    sim_time_ = 0;
  }

  varTypes_.clear();
  varValues_.clear();
  leadDesignator_ = "";
  var_vals_.clear();
  tree_ = static_cast<ExpressionNode *> (NULL);
  derivs_.clear();
  free_list_.clear();
  ee_list_.clear();

  Rmode_ = 0;
  ind_replace_ = 0;
  Rstring_ = "";
  Rcval_ = 0;
  curr_magic_ = 0;
  curr_num_ = 1;
  values_changed_ = true;

// First step in parsing expression is to check for extraneous braces, check
// matching of parens, and reduce white space to single spaces to form the
// string 'compact_line'.  The only reason not to completely remove white
// space is for poly, which relys on white space as a delimiter

  line = exp;
  Input_ = exp;
  compact_line = "";

  compactLine_(line,compact_line);

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Compact input line is: \'" << compact_line << "\'" << std::endl;

// Second step is to do the actual parsing into ExpressionElement structures
// which describe the type of each token and the value, if applicable

// TVR Note:  In SPICE, lexical analysis is done token by token by a
// "PTlexer" function that returns a single token at a time, and the
// parse tree building steps just call that when they need a new token.
// We don't do that, because many of the things we try to do involve massaging
// the input expression before parsing it, such as by converting "poly"
// expressions into regular polyomials, and massaging "table" input.  So
// rather than tokenize on the fly as SPICE does, we do it all up front first.
  tokenize_(compact_line, tok);

  if (DEBUG_EXPRESSION) {
    Xyce::dout() << std::endl << "Token List:\n";
    for (tok_i = tok.begin() ; tok_i != tok.end() ; ++tok_i)
    {
      Xyce::dout() << "Token: " << toks[(*tok_i)->token];
      if ((*tok_i)->token == TOK_VALUE)
      {
        Xyce::dout() << ": ";
        switch ((*tok_i)->type)
        {
          case TYP_NUM:
            Xyce::dout() << (*tok_i)->number;
            break;
          case TYP_STRING:
            Xyce::dout() << "\"" << ((*tok_i)->name) << "\"";
            break;
          case TYP_PNODE:
            Xyce::dout() << "Pointer to node";
            break;
        }
      }
      Xyce::dout() << std::endl;
    }
  }

  // Now that we've already tokenized the entire string, we process the
  // tokenized version further.

  // Third step is to handle the special case formats like table {} ()
  // () and poly.  Do poly first since this uses TOK_SPACE as a
  // delimiter.  Then remove the TOK_SPACE elements and do table.


  // TVR:  "Handling" POLY means converting the obscure "POLY" format
  // into an actual polynomial.  e.g., if the input is
  //  POLY(1) V(2) 3 2 1
  // then the tokenized version coming out of this would be the same as if
  // the original input had instead been
  //   3+2*v(2)+1*v(2)*V(2)
  // This completely replaces the token list from the first form with that of
  // the second.

  convertPolyToExpr_(tok);

  // Now we can remove spaces from the token list.  We only needed
  // them because POLY requires space as delimiter, and we rewrote the
  // input if it was POLY by now.
  for (tok_i = tok.begin() ; tok_i != tok.end() ; ++tok_i)
  {
    el = *tok_i;
    if (el->token == TOK_SPACE)
      tok_i = tok.erase(tok_i);
  }

  // Handle TABLE and SCHEDULE.  This is mostly clean-up, and is meant to
  // deal with the fact that SCHEDULE needs an implicit argument made
  // explicit, and that there are multiple variants of the form of TABLE.
  standardizeTable_(tok);

  // Now fake out the expression by replacing all remaining braces with
  // the equivalent parens
  braceToParen_(tok);
  
  // We allow PDE devices to be in expressions as:
  // Y%PDE%device_name  --> after tokenization is Y PDE device_name
  //                        with a TOK_REMINDER after the Y and PDE (i.e. the % sign).
  // YPDE device_name   --> after tokenization is YPDE devcie_name
  // device_name        --> after tokenization is device_name
  // So, look the sequences:
  //  type         Value
  //  TOK_STRING  "YPDE"
  //  OR
  //  TOK_STRING  "Y"
  //  TOK_REMAINDER
  //  TOK_STRING  "PDE"
  //  TOK_REMAINDER
  // and remove them.

  int numTok = tok.size();
  int tokCount = 0;
  for (tok_i = tok.begin() ; tok_i != tok.end() ; ++tok_i)
  {
    if (((*tok_i)->type == TYP_STRING) && ((*tok_i)->name == "YPDE"))
    {
      Xyce::dout() << "Erasing 1" << (*tok_i)->name << std::endl;
      // safe to erase this token
      tok_i = tok.erase(tok_i);
    }
    else
    {
      if ( (((*tok_i)->type == TYP_STRING) && ((*tok_i)->name == "Y")) && ((tokCount + 4) < numTok) )
      {
        // found a "Y" and at least enough room on the tok list that there could be '%' "PDE" '%'
        // so check for that.
        tok_j=tok_i;  // remember position of "Y"
        std::list<ExpressionElement *>::iterator tok_j_plus_1 = ++tok_i;
        std::list<ExpressionElement *>::iterator tok_j_plus_2 = ++tok_i;
        std::list<ExpressionElement *>::iterator tok_j_plus_3 = ++tok_i;
        if( ((*tok_j_plus_1)->token == TOK_REMAINDER) &&
            (((*tok_j_plus_2)->type == TYP_STRING) && ( (*tok_j_plus_2)->name == "PDE" )) &&
            ((*tok_j_plus_3)->token == TOK_REMAINDER) )
        {
          // safe to erase both tokens
          Xyce::dout() << "Erasing 2" << (*tok_j)->name << std::endl;
          Xyce::dout() << "Erasing 3" << (*tok_j_plus_3)->name << std::endl;
          tok_i++;
          tok_i = tok.erase(tok_j, tok_i);
        }
        tokCount+=3;
      }
    }
    tokCount++;
  }

  if (DEBUG_EXPRESSION) {
    Xyce::dout() << std::endl << "After Table/Schedule/Poly, Token List:\n";
    for (tok_i = tok.begin() ; tok_i != tok.end() ; ++tok_i)
    {
      Xyce::dout() << "Token: " << toks[(*tok_i)->token];
      if ((*tok_i)->token == TOK_VALUE)
      {
        Xyce::dout() << ": ";
        switch ((*tok_i)->type)
        {
          case TYP_NUM:
            Xyce::dout() << (*tok_i)->number;
            break;
          case TYP_STRING:
            Xyce::dout() << "\"" << ((*tok_i)->name) << "\"";
            break;
          case TYP_PNODE:
            Xyce::dout() << "Pointer to node";
            break;
        }
      }
      Xyce::dout() << std::endl;
    }
  }

  // Now the tokenizing of the expression is completed, and the various
  // obscure and complex formats have been massaged down to much simpler
  // formats --- everything is now a straight formula, which could include
  // function calls.  We now parse that into a parse tree.

// Fourth step, the actual formation of the tree

// We read from left to right and create a stack of symbols until the next
// symbol has greater precedence than what we have at the top.  At that point,
// the stack is reduced by searching backward
  if (DEBUG_EXPRESSION)
    Xyce::dout() << " Forming parse tree." << std::endl;

  el = newExpressionElement_();
  el->token = TOK_END;
  stack.clear();
  stack.push_back(el);
  tok_i = tok.begin();
  next = *tok_i++;
  sp = 0;

  while (sp > 1 || next->token != TOK_END)
  {
    if (DEBUG_EXPRESSION)
    {
      Xyce::dout() << "  In while loop, sp is " << sp << std::endl;
      Xyce::dout() << "         Next token is " << toks[next->token] << std::endl;
      Xyce::dout() << "         Full stack at top of loop:" << std::endl;
      for (j=0 ; j<=sp ; ++j)
      {
        Xyce::dout() << "Token: " << op_name[stack[j]->token] << std::endl;
        if (stack[j]->token == TOK_VALUE)
        {
          if (stack[j]->type == TYP_NUM)
          {
            Xyce::dout() << "  Numeric element: " << stack[j]->number << std::endl;
          }
          if (stack[j]->type == TYP_STRING)
          {
            Xyce::dout() << "  String element: '" << stack[j]->name << "'" << std::endl;
          }
          if (stack[j]->type == TYP_PNODE)
          {
            Xyce::dout() << "  parsenode, type = " << expr_ops[stack[j]->node->type]
                         << std::endl;
          }
        }
      }
      Xyce::dout() << " Now working back from top of stack to first non-value token" << std::endl;
    }

    stack.resize(sp+1);
    i = sp;
    do
    {
      top = stack[i--];

    } while (top->token == TOK_VALUE);

    if (DEBUG_EXPRESSION)
    {
      Xyce::dout() << " Top has been popped: " << std::endl;
      Xyce::dout() << "Token: " << toks[top->token];
      if (top->token == TOK_VALUE)
      {
        Xyce::dout() << ": ";
        switch (top->type)
        {
          case TYP_NUM:
            Xyce::dout() << top->number;
            break;
          case TYP_STRING:
            Xyce::dout() << "\"" << (top->name) << "\"";
            break;
          case TYP_PNODE:
            Xyce::dout() << "Pointer to node";
            break;
        }
      }
      Xyce::dout() << std::endl;
      Xyce::dout() << " Precedence between " << toks[top->token] << " and " << toks[next->token] << " is " << (int) prectable[top->token][next->token] << std::endl;
    }

    switch (prectable[top->token][next->token])
    {
      case L:
      case E:
// Push the token read
        stack.push_back(next);
        if (DEBUG_EXPRESSION)
        {
          Xyce::dout() << "Token " << toks[next->token] << " has lower or equal precedence to top of stack token " << toks[top->token] << " so pushing next token." << std::endl;
        }
        ++sp;
        next = *tok_i++;
        continue;

      case R:
        goto err;

      case G:
        // TVR note:  The helpful comment below
        // is verbatim from SPICE.  Love the specificity of
        // "try and do some stuff."
        //------------
// Reduce. Make st and sp point to the elts on the
// stack at the end and beginning of the junk to
// reduce, then try and do some stuff. When scanning
// back for a <, ignore VALUES.
        //------------

        // As soon as we hit a token of greater precedence than what
        // we have on the stack, we create a parse tree fragment for
        // what we have so far.  The reduction is from a subexpression
        // string to a single value, even though that value is
        // represented by a parse tree fragment.

        if (DEBUG_EXPRESSION) {
          Xyce::dout() << "Token " << toks[next->token] << " has greater precedence to top of stack token " << toks[top->token] << " so reducing." << std::endl;
          Xyce::dout() << "\nFull stack so far:" << std::endl;
          for (j=0 ; j<=sp ; ++j)
          {
            Xyce::dout() << "Token: " << op_name[stack[j]->token] << std::endl;
            if (stack[j]->token == TOK_VALUE)
            {
              if (stack[j]->type == TYP_NUM)
              {
                Xyce::dout() << "  Numeric element: " << stack[j]->number << std::endl;
              }
              if (stack[j]->type == TYP_STRING)
              {
                Xyce::dout() << "  String element: '" << stack[j]->name << "'" << std::endl;
              }
              if (stack[j]->type == TYP_PNODE)
              {
                Xyce::dout() << "  parsenode, type = " << expr_ops[stack[j]->node->type]
                             << std::endl;
              }
            }
          }
        }

        // Figure out what fragment we're supposed to reduce.  Search backward
        // from the end of the stack.  If the thing under "sp" and the thing
        // prior to that (skipping VALUES) have "L" precedence, then we stop,
        // and everything from sp to the end is our fragment.  Otherwise, move
        // sp back and keep checking until we *do* have an "L" precedence.

        st = sp;
        if (stack[sp]->token == TOK_VALUE)
          sp--;
        while (sp > 0)
        {
          if (stack[sp - 1]->token == TOK_VALUE)
            i = 2;  // No 2 pnodes together...
          else
            i = 1;
          if (prectable[stack[sp - i]->token][stack[sp]->token] == L)
            break;
          else
            sp = sp - i;
        }
        if ((sp>0) && (stack[sp - 1]->token == TOK_VALUE))
          sp--;
        if ((stack[st-1]->token == TOK_COLON) && (sp>0) && (stack[sp-1]->token == TOK_QUESTION))
        {
          sp--;
          if ((sp>0) && (stack[sp - 1]->token == TOK_VALUE))
            sp--;
        }
        // everything between sp and st is to be reduced.

// Now try and see what we can make of this.
// The possibilities are: - node
//              node op node
//              ( node )
//              func ( node )
//              func ( node, node, node, ... )
//              node

        if (DEBUG_EXPRESSION) {
          Xyce::dout() << "\nProcessing stack chunk:\n";
          for (j=sp ; j<=st ; ++j)
          {
            if (stack[j]->token == TOK_VALUE)
            {
              if (stack[j]->type == TYP_NUM)
              {
                Xyce::dout() << "Numeric element: " << stack[j]->number << std::endl;
              }
              if (stack[j]->type == TYP_STRING)
              {
                Xyce::dout() << "String element: '" << stack[j]->name << "'" << std::endl;
              }
              if (stack[j]->type == TYP_PNODE)
              {
                Xyce::dout() << "parsenode, type = " << expr_ops[stack[j]->node->type]
                             << std::endl;
              }
            }
            else
            {
              Xyce::dout() << "Token: " << op_name[stack[j]->token] << std::endl;
            }
          }

          Xyce::dout() << " reducing to a pnode..." << std::endl;
        }

        if (st == sp)
        {
          // the fragment is just one element, make a pnode out of it.
          pn = makepnode_(stack[st]);
          if (pn == NULL)
            goto err;
        }
        else if ((stack[sp]->token == TOK_UMINUS) && (st == sp + 1))
        {
          // The fragment is a unary minus, and the fragment is just the
          // minus token and one more element. First make a node out of the
          // element, then a function node for the unary minus operation
          lpn = makepnode_(stack[st]);
          if (lpn == NULL)
            goto err;
          pn = mkfnode_("-", 1, lpn);
        }
        else if ((stack[sp]->token == TOK_NOT) && (st == sp + 1))
        {
          // This is "~foo", the NOT operation.  Same technique as unary minus
          lpn = makepnode_(stack[st]);
          if (lpn == NULL)
            goto err;
          pn = mkfnode_("~", 1, lpn);
        }
        else if ((stack[sp]->token == TOK_LPAREN) && (stack[st]->token == TOK_RPAREN))
        {
          // something in parentheses, make a node out of it
          pn = makepnode_(stack[sp + 1]);
          if (pn == NULL)
            goto err;
        }
        else if ((stack[sp + 1]->token == TOK_LPAREN) && (stack[st]->token == TOK_RPAREN))
        {
          // we have
          //     <something>(<some number of tokens>)
          // This could be a nodal expression (e.g. V(X)), a function call
          // (including TABLE or SCHEDULE), etc.  Regardless, the first token
          // better be the name (e.g. of string type).  All of these are
          // handled by the "mkfnode_" node "factory", and things like V(x)
          // are just special cases in that routine.
          int offset=1;
          if (stack[sp]->type != TYP_STRING)
            goto err;

          // comma separated arguments.., so each one is 2 tokens.  Count 'em
          // and make a list of pnodes.
          num_nodes = (st-sp-offset)/2;
          node_list.clear();
          for (i=0 ; i<num_nodes ; ++i)
          {
            node_list.push_back(makepnode_(stack[sp + (i+offset)*2]));
            if (node_list[i] == NULL)
              goto err;
          }

          // make the "function" node (which could wind up just being a VAR)
          if (!(pn = mkfnode_(stack[sp]->name, num_nodes, node_list)))
            goto err;
        }
        else if (sp+4<=st &&
                 (stack[sp+1]->token == TOK_QUESTION)
                 && (stack[sp+3]->token==TOK_COLON)
                 && stack[sp]->token==TOK_VALUE
                 && stack[sp+2]->token==TOK_VALUE
                 && stack[sp+4]->token==TOK_VALUE)
        {
          if (DEBUG_EXPRESSION)
          {
            Xyce::dout() << "Found ternary operator of simple type" << std::endl;
          }
          num_nodes = 3;
          node_list.clear();
          node_list.push_back(makepnode_(stack[sp]));
          node_list.push_back(makepnode_(stack[sp+2]));
          node_list.push_back(makepnode_(stack[sp+4]));
          if (!(pn = mkfnode_("IF",num_nodes,node_list)))
            goto err;
        }
        else    // node op node
        {
          // Otherwise this is a simple node op node case, so make pnodes for
          //  nodes, then make a pnode for the operator.
          lpn = makepnode_(stack[sp]);
          rpn = makepnode_(stack[st]);
          if ((lpn == NULL) || (rpn == NULL))
            goto err;
          pn = mkb_(stack[sp + 1]->token, lpn, rpn);
          if (pn == NULL)
            goto err;
        }

        // We now have a parse node for the subexpression we're reducing.
        // Replace whatever we had at the beginning of the subexpression
        // with a special value token of type "pnode" and store our pointer.
        // Note that the entire stack will be resized to sp+1 meaning we've
        // now got our new node right at the end of the stack, effectively
        // having popped off the subexpression.
        stack[sp]->token = TOK_VALUE;
        stack[sp]->type = TYP_PNODE;
        stack[sp]->node = pn;

        if (DEBUG_EXPRESSION) {
          Xyce::dout() << " The stack element was replaced by a pnode whose tree is:"
                       << std::endl;
          dumpParseTree_(pn);
        }
        continue;
    default:
      goto err;
    }
  }

  // The token list has now been reduced to a single element.  Make a
  // tree out of it.  It is entirely possible it is already a tree, in
  // which case the makepnode_ method will just return the pointer to
  // the tree stored in element.

  if (DEBUG_EXPRESSION) {
    Xyce::dout() << " Attempting to form tree from stack[1], which is: " << std::endl;
    if (stack[1])
    {
      Xyce::dout() << "Token: " << toks[stack[1]->token];
      if (stack[1]->token == TOK_VALUE)
      {
        Xyce::dout() << ": ";
        switch (stack[1]->type)
        {
          case TYP_NUM:
            Xyce::dout() << stack[1]->number;
            break;
          case TYP_STRING:
            Xyce::dout() << "\"" << (stack[1]->name) << "\"";
            break;
          case TYP_PNODE:
            Xyce::dout() << "Pointer to node";
            break;
        }
      }
      else
      {
        Xyce::dout() << " NULL" << std::endl;
      }
      Xyce::dout() << std::endl;
    }
  }

  tree_ = makepnode_(stack[1]);
  if (tree_ == NULL)
  {
    goto err;
  }

  if (DEBUG_EXPRESSION) {
    Xyce::dout() << " Tree is made. Before PTcheck_:" << std::endl;
    dumpParseTree_(tree_);
  }

// Fifth step, resolve any remaining placeholder and find and use any
// common subexpressions.  This is tricky because the subexpression
// condensation will free some ExpressionNodes which can only really be
// freed after the operation is done.  This is because the tree must
// be fully traversable during this operation.  The items to be freed
// are collected in done_list_

  PThead_ = tree_;
  done_list_.clear();

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Calling PTcheck_" << std::endl;

  if (tree_ != PTcheck_ (PThead_))
    goto err;
  for (free_i = done_list_.begin() ; free_i != done_list_.end() ; ++free_i)
    deleteExpressionNode_(*free_i);

  if (DEBUG_EXPRESSION) {
    Xyce::dout() << " After PTcheck_:" << std::endl;
    dumpParseTree_(tree_);
  }

// Sixth step, differentiate any ddx() functions and replace the functions
// with the differentiated expression

  set_nums_();
  create_vars_();

  if (!ddxProcessed_ && num_string_ == 0 && num_func_ == 0)
  {
    tree_ = diffDDX_ (PThead_);
    PThead_ = tree_;
  }

// Seventh step, zero vars, record number of each type, and delete all
// ExpressionElements, as they are only used for parsing

  for (ee_i=ee_list_.begin() ; ee_i != ee_list_.end() ; ++ee_i)
  {
    el = *ee_i;
    delete el;
  }

  return true;

err:
  Report::UserError() << "Syntax error in expression " << Input_;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::getSymbolTable
// Purpose       : Returns the entire symbol table for this expression
// Special Notes : Names in varValues_ are not necessarily unique, and
//                 are NOT unique if more than one "INSTANCE" type variable
//                 is a lead current.  IG(M1) and ID(M1) both have the
//                 name "M1" --- the leadDesignator_ and possibly
//                 varTypes_ differ.
//                 Because of this, it is not correct to do a
//                 get_names(XEXP_ALL,names) and then iterate over names,
//                 calling "get_type" on that name.
//                 To avoid this mistake, this function has been created
//                 so that the entire tuple (name, type, lead designator)
//                 is returned for each variable.
//                 
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 08/19/2016
//-----------------------------------------------------------------------------
///
/// Returns a representation of the symbol table for this expression
///
/// @param[out] theSymbolTable  An STL vector of ExpressionSymbolTableEntry
///
/// The internal representation of a "symbol table" is just a set of
/// parallel vectors, varValues_, varTypes_, leadDesignator_, and var_vals.
/// Each of these can be queried by various "get_" functions, but often
/// more information is needed than a single one of these functions provides.
/// In several cases, upstream code needs a more comprehensive view of the
/// symbols required by an expression.  This function provides this view.
///
/// Specifically, this function was created because names in the (misnamed)
/// varValues_ array are not necessarily unique: for lead currents, e.g.
/// "IG(M1)" or "ID(M1)", the name stored is just the instance name (M1),
/// and the other information is spread out between varTypes_ (showing that
/// this is a lead current) and leadDesignator_ (showing that it's the Gate
/// or Drain current).  Attempting to get a list of names, then using those
/// names to query type or lead designator will fail, because the name is
/// not unique.
///
/// @author Tom Russo, SNL
/// @date 08/19/2016
///
void ExpressionInternals::getSymbolTable(
   std::vector< ExpressionSymbolTableEntry > &theSymbolTable)  const
{
  theSymbolTable.clear();

  for (int i = 0; i < numVars_; ++i)
  {
    // the type in varTypes_ is NOT the actual type needed by the users
    // of this function.  It must be converted.  See get_type().
    int var_type=0;
    switch (varTypes_[i])
    {
      case EXPR_T_NODAL_COMPUTATION:
        var_type = XEXP_NODAL_COMPUTATION;
        break;
      case EXPR_T_FUNCTION:
        var_type = XEXP_FUNCTION;
        break;
      case EXPR_T_NODE:
        var_type = XEXP_NODE;
        break;
      case EXPR_T_STRING:
        var_type = XEXP_STRING;
        break;
      case EXPR_T_INSTANCE:
        if (leadDesignator_[i] == ' ')
        {
          std::string s = varValues_[i];
          std::string::size_type pos = s.find_last_of(":");
          if( pos == std::string::npos )
            pos = 0;
          else
            ++pos;
          if (s[pos] == 'V')
            var_type = XEXP_INSTANCE;
          else
            var_type = XEXP_LEAD;
        }
        else
        {
          var_type = XEXP_LEAD;
        }
        break;
      case EXPR_T_SPECIAL:
        var_type = XEXP_SPECIAL;
        break;
      case EXPR_T_VARIABLE:
        var_type = XEXP_VARIABLE;
        break;
    }

    theSymbolTable.push_back(ExpressionSymbolTableEntry(varValues_[i],var_type,leadDesignator_[i]));
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_names
// Purpose       : Returns the names of input quantities by type
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Get a list of names of input quantities, possibly restricting by type
///
/// @param[in] type  Value from enum XEXP_TYPES, which type of name to extract
/// @param[out] names string vector in which to return the names
///
/// If type is XEXP_ALL, all names are returned.  Otherwise, only names
/// of the specific type (XEXP_NODE, XEXP_NODAL_COMPUTATION, XEXP_INSTANCE,
/// etc.) are returned.
///
/// @note It is important to note that in the case of XEXP_LEAD, the names
/// returned might not be unique.  They will not be unique if two different
/// lead currents from the same device are used in a single expression
/// (e.g. IG and ID for a MOSFET), in which case both variables will have
/// the MOSFET instance name returned.  In this case, trying to use the name
/// in one of the other lookup functions (e.g. get_type) will not be correct,
/// and may return the type of the wrong variable.
///
/// @author Dave Shirley, PSSI
///  @date 06/07/01
///
void
ExpressionInternals::get_names(
  int                           type,
  std::vector<std::string> &    names)
{
  names.clear();

  for (int i = 0; i < numVars_; ++i)
  {
    if (type == XEXP_ALL)
    {
      names.push_back(varValues_[i]);
    }
    else
    {
      switch (varTypes_[i])
      {
        case EXPR_T_NODE:
          if (type == XEXP_NODE)
            names.push_back(varValues_[i]);
          break;
        case EXPR_T_NODAL_COMPUTATION:
          if (type == XEXP_NODAL_COMPUTATION)
            names.push_back(varValues_[i]);
          break;
        case EXPR_T_INSTANCE:
          if (leadDesignator_[i] == ' ')
          {
            std::string s = varValues_[i];
            std::string::size_type pos = s.find_last_of(":");
            if( pos == std::string::npos )
              pos = 0;
            else
              ++pos;
            if ((s[pos] == 'V') && type == XEXP_INSTANCE)
              names.push_back(varValues_[i]);
            else if (s[pos] != 'V' && type == XEXP_LEAD)
            {
              names.push_back(varValues_[i]);
            }
          }
          else
          {
            if (type == XEXP_LEAD)
            {
              std::string lead = varValues_[i] + "{" + leadDesignator_[i] + "}";
              names.push_back(lead);
            }
          }
          break;
        case EXPR_T_STRING:
          if (type == XEXP_STRING)
            names.push_back(varValues_[i]);
          break;
        case EXPR_T_SPECIAL:
          if (type == XEXP_SPECIAL)
            names.push_back(varValues_[i]);
          break;
        case EXPR_T_VARIABLE:
          if (type == XEXP_VARIABLE)
            names.push_back(varValues_[i]);
          break;
        case EXPR_T_FUNCTION:
          if (type == XEXP_FUNCTION)
            names.push_back(varValues_[i]);
          break;
        default:
          break;
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_type
// Purpose       : Finds the type of an input quantity name
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Return the type of the named variable
///
/// @param[in] var  The variable name
/// @return The type of the variable
///
/// @note The name is looked up in the list of variable names and the first
/// match is used.  In the case of lead currents, the name returned by
/// get_names is the name of the instance associated with the current,
/// and may not be unique if more than one lead current for the same device
/// is present in the expression.  To guard against this kind of problem,
/// one should use getSymbolTable, which returns name/type/lead designator
/// tuples all at once when trying to process expressions that might include
/// lead currents.
///
/// @author Dave Shirley, PSSI
///  @date 06/07/01
///
int
ExpressionInternals::get_type(
  const std::string &   var)
{
  int index;
  int var_type=-1;

  if ( ( index = find_num_( var ) ) >= 0 )
  {
    switch (varTypes_[index])
    {
      case EXPR_T_NODAL_COMPUTATION:
        var_type = XEXP_NODAL_COMPUTATION;
        break;
      case EXPR_T_FUNCTION:
        var_type = XEXP_FUNCTION;
        break;
      case EXPR_T_NODE:
        var_type = XEXP_NODE;
        break;
      case EXPR_T_STRING:
        var_type = XEXP_STRING;
        break;
      case EXPR_T_INSTANCE:
        if (leadDesignator_[index] == ' ')
        {
          std::string s = varValues_[index];
          std::string::size_type pos = s.find_last_of(":");
          if( pos == std::string::npos )
            pos = 0;
          else
            ++pos;
          if (s[pos] == 'V')
            var_type = XEXP_INSTANCE;
          else
            var_type = XEXP_LEAD;
        }
        else
        {
          var_type = XEXP_LEAD;
        }
        break;
      case EXPR_T_SPECIAL:
        var_type = XEXP_SPECIAL;
        break;
      case EXPR_T_VARIABLE:
        var_type = XEXP_VARIABLE;
        break;
    }
    return var_type;
  }

  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::make_constant
// Purpose       : Convert a 'string' placeholder into a constant
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Convert a 'string' placeholder into a constant
///
/// @param[in] var  The string to replace
/// @param[in] val  The value of the constant to substitute
///
/// All occurances of the variable var in the expression are replaced with
/// the value given, and the tree is recursively simplified:  any nodes that
/// now evaluate to a constant as a result of the replacement are themselves
/// turned into constant nodes.
///
/// @author Dave Shirley, PSSI
///  @date 06/07/01
///
bool ExpressionInternals::make_constant (const std::string & var,
                                      const double & val)
{
  int index;

  values_changed_ = true;
  if ( ( index = find_num_( var ) ) >= 0 )
  {
    if( varTypes_[index] == EXPR_T_STRING )
    {
      num_string_--;
      convert_to_constant_(index, val);
      simplify_ (*tree_);

      return true;
    }
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::make_var
// Purpose       : Convert a 'string' placeholder into a variable
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Convert a string placeholder into a variable
///
/// @param[in] var  variable name
///
/// The parse tree is walked, converting any node of type "EXPR_T_STRING"
/// whose name matches the given name into an EXPR_T_VARIABLE node.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
bool ExpressionInternals::make_var (std::string const & var)
{
  int index;
  bool retval = false;

  values_changed_ = true;
  if ((index = find_num_(var)) >= 0)
  {
    if (varTypes_[index] == EXPR_T_STRING)
    {
      varTypes_[index] = EXPR_T_VARIABLE;
      num_string_--;
      ++num_var_;
      convert_to_variable_(index);

      retval=true;
    }
  }
  set_nums_();
  if (!ddxProcessed_ && num_string_ == 0 && num_func_ == 0)
  {
    tree_ = diffDDX_ (PThead_);
    PThead_ = tree_;
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::differentiate
// Purpose       : Form the analytic derivative trees for all variables
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Form the analytic derivative trees of this expression by differentiating the
/// expression with respect to all variables.
///
/// @return  A positive return value indicates that differentiation has been successfule, and counts the number of variables differentiated with respect to.  A negative value indicates that no differentiation has been performed.  -2 indicates that the expression has no variables.  -1 indicates that the operation has previously been performed, and has been skipped.
///
/// The derivs_ vector contains a list of expression trees that evaluate
/// the derivatives of the main expression with respect to each variable it
/// depends on.
/// 
/// @note This method sets up the expression trees for derivatives, it
/// does not evaluate the derivatives themselves.  If it is
/// successful, this method sets the "differentiated_" member variable
/// to true, so that subsequent calls will not do anything.
///
/// @author Dave Shirley, PSSI
///  @date 06/07/01
///
int ExpressionInternals::differentiate ()
{
  int i;

  values_changed_ = true;
  if (num_string_ > 0 || num_func_ > 0)
    return -2;

  if( !differentiated_ )
  {
    if (!ddxProcessed_ && num_string_ == 0 && num_func_ == 0)
    {
      tree_ = diffDDX_ (PThead_);
      PThead_ = tree_;
    }
    simplify_ (*tree_);
    derivs_.resize(numVars_);

    if (DEBUG_EXPRESSION)
      Xyce::dout() << std::endl;

    for (i=0 ; i<numVars_ ; ++i)
    {

      if (DEBUG_EXPRESSION)
        Xyce::dout() << "Differentiating with respect to: " << varStr_(i) << ":\n" <<
          get_expression() << std::endl;

      derivs_[i] = Differentiate_ (tree_, i);

      if (DEBUG_EXPRESSION) {
        std::ostringstream s("");
        s << std::setprecision(PRECISION);
        RpTree_ (derivs_[i], s);
        Xyce::dout() << "Derivative:\n" << s.str() << std::endl << std::endl;
      }
    }
    differentiated_ = true;
    return numVars_;
  }

  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_var
// Purpose       : Sets the value of an input quantity
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Set the value of a single, named variable
///
/// @param[in] var   The variable name
/// @param[in] val   The value
///
/// @note This function should NOT be used to set the value of individual
/// lead currents such as IG(M1) or ID(M1), for reasons described in
/// comments for get_names and get_type.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
bool ExpressionInternals::set_var ( const std::string & var,
                                 const double & val)
{
  int ind;

  if (DEBUG_EXPRESSION)
  {
    Xyce::dout() << "ExpressionInternals::set_var(" << var << " , " << val << ")" << std::endl;
  }

  if ( ( ind = find_num_( var ) ) >= 0 )
  {
    if (varTypes_[ind] == EXPR_T_INSTANCE ||
        varTypes_[ind] == EXPR_T_NODE ||
        varTypes_[ind] == EXPR_T_SPECIAL ||
        varTypes_[ind] == EXPR_T_VARIABLE)
    {
      if (var_vals_[ind] != val)
      {
        values_changed_ = true;
        var_vals_[ind] = val;
      }
    }
  }

  return values_changed_;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_vars
// Purpose       : Sets the values of all input quantities
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Set the value of all variables in the expression
///
/// @param[in] vals  Vector of variable values
///
/// @note The values in this list must be in the same order as the variable
/// names in the varValues_ array!
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
bool ExpressionInternals::set_vars ( const std::vector<double> & vals )
{
  if (vals.empty()) return values_changed_;

  int jLimit = vals.size();
  for (int i=0, j=0 ; i<numVars_ ; ++i)
  {
    if (var_vals_[i] != vals[j])
    {
      values_changed_ = true;
      var_vals_[i] = vals[j];
    }
    if(++j >= jLimit)
      break;
  }

  return values_changed_;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_expression
// Purpose       : Returns a string of the expression
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Return a string representation of the expression
///
/// @return a string representation of the expression reconstructed from its tree
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
std::string ExpressionInternals::get_expression () const
{
  std::ostringstream s("");
  s << std::setprecision(PRECISION);
  if (tree_ != NULL)
    RpTree_ (tree_, s);

  return s.str();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_derivative
// Purpose       : Returns a string of a derivative
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Return a string representation of the expression for a derivative with respect to one variable
///
/// @param[in] var variable name 
/// @return a string representation of the expression reconstructed from its tree
///
/// @note Recall the caution in get_names, get_type, and other functions:
/// lead current names are not unique in the varValues_ array, and therefore
/// attempting to call this function on a lead current variable can
/// return the wrong one if more than one lead current for a single device
/// is present in the expression (e.g.  ID(M1)+IG(M1))!
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
std::string ExpressionInternals::get_derivative ( std::string const & var )
{
  int index;

  index = find_num_ (var);
  if ( index >= 0)
  {
    if (!differentiated_)
    {
      int diff_retval = differentiate();
      if (diff_retval < 0)
      {
        std::vector< std::string > temparray;
        Report::UserError() << "get_derivative: Unable to differentiate: " << get_expression();
        if (diff_retval == -2)
        {
          if (num_string_>0)
          {
            Xyce::lout() << "Unable to differentiate because of " << num_string_ << " unresolved strings:" << std::endl;
            get_names(XEXP_STRING,temparray);
            for (int i=0;i<temparray.size();i++)
            {
              Xyce::lout() << "    " << temparray[i] << std::endl;
            }
          }
          if (num_func_>0)
          {
            Xyce::lout() << "Unable to differentiate because of " << num_func_ << " unresolved functions:" << std::endl;
            get_names(XEXP_FUNCTION,temparray);
            for (int i=0;i<temparray.size();i++)
            {
              Xyce::lout() << "    " << temparray[i] << std::endl;
            }
          }
        } else {
          Xyce::lout() << "Unable to differentiate because of prior skipped differentiation." << std::endl;
        }
        return NULL;
      }
    }
    std::ostringstream s("");
    s << std::setprecision(PRECISION);
    RpTree_ (derivs_[index], s);
    return s.str();
  }

  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_num
// Purpose       : Returns the number of input quantities of a requested type
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Returns the number of input variables of the given type
///
/// @param[in] type   Type from the XEXP_TYPES enum
/// @return number of input quantities of selected type
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///

int
ExpressionInternals::get_num(
  int           type)
{
  switch (type)
  {
    case XEXP_ALL:
      return numVars_;
    case XEXP_NODE:
      return num_N_;
    case XEXP_INSTANCE:
      return num_I_;
    case XEXP_LEAD:
      return num_lead_;
    case XEXP_STRING:
      return num_string_;
    case XEXP_SPECIAL:
      return num_special_;
    case XEXP_VARIABLE:
      return num_var_;
    case XEXP_FUNCTION:
      return num_func_;
  }
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::evaluate
// Purpose       : Evaluate expression and derivatives using provided input values
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Evaluate expression and derivatives using provided input values
///
/// @param[out] exp_r expression result
/// @param[out] deriv_r vector of derivative results
/// @param[in]  vals  vector of input values
/// @return  error code
///
/// @note the vals array of variable values must be in the same order as the
/// names in the varValues_ array!
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
int ExpressionInternals::evaluate ( double & exp_r,
                                 std::vector<double> & deriv_r,
                                 std::vector<double> & vals )
{
  set_vars(vals);
  return evaluate ( exp_r, deriv_r );
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::evaluateFunction
// Purpose       : Evaluate expression using provided input values.
// Special Notes : This is for cases in which the user does not need
//                 the derivatives.
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 04/14/08
//-----------------------------------------------------------------------------
///
/// Evaluate expression using provided input values
///
/// @param[out] exp_r expression result
/// @param[in]  vals  vector of input values
/// @return  error code
///
/// @note the vals array of variable values must be in the same order as the
/// names in the varValues_ array!
///
/// @author Eric Keiter, SNL
/// @date 04/14/08
///
int ExpressionInternals::evaluateFunction ( double & exp_r,
                                                  std::vector<double> & vals )
{
  set_vars(vals);
  return evaluateFunction ( exp_r );
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Evaluate expression and derivatives using stored input values
///
/// @param[out] exp_r expression result
/// @param[out] deriv_r vector of derivative results
/// @return  error code
///
/// Variable values must previously have been set using calls to set_vars
/// or set_var.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
int ExpressionInternals::evaluate ( double & exp_r,
                                 std::vector<double> & deriv_r)
{
  if (!differentiated_)
  {
    int diff_retval=differentiate();
    if (diff_retval < 0)
    {
      std::vector< std::string > temparray;
      Report::UserError() << "evaluate: Unable to differentiate: " << get_expression();
      if (diff_retval == -2)
      {
        if (num_string_>0)
        {
          Xyce::lout() << "Unable to differentiate because of " << num_string_ << " unresolved strings:" << std::endl;
          get_names(XEXP_STRING,temparray);
          for (int i=0;i<temparray.size();i++)
          {
            Xyce::lout() << "    " << temparray[i] << std::endl;
          }
        }
        if (num_func_>0)
        {
          Xyce::lout() << "Unable to differentiate because of " << num_func_ << " unresolved functions:" << std::endl;
          get_names(XEXP_FUNCTION,temparray);
          for (int i=0;i<temparray.size();i++)
          {
            Xyce::lout() << "    " << temparray[i] << std::endl;
          }
        }
      } else {
        Xyce::lout() << "Unable to differentiate because of prior skipped differentiation." << std::endl;
      }
      EXPR_ERROR (EXPR_NODERIV_FATAL);
      return EXPRerrno;
    }
  }

  (void) evaluateFunction (exp_r);

  if (differentiated_)
  {
    for (int i=0, j=0 ; i<numVars_ ; ++i)
    {
      if ( varTypes_[i] != EXPR_T_SPECIAL )
      {
        if ( deriv_r.size() <= j )
        {
          // prevent null parameter passing
          deriv_r.push_back( 0.0 );
        }
        EXPReval_( *(derivs_[i]), deriv_r[j++], var_vals_ );
      }
    }
  }
  else
  {
    EXPR_ERROR (EXPR_NODERIV_FATAL);
  }

  return EXPRerrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::evaluateFunction
// Purpose       : Evaluate expression using stored input values
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Evaluate expression using stored input values
///
/// @param[out] exp_r expression result
/// @return  error code
///
/// Variable values must previously have been set using calls to set_vars
/// or set_var.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
int ExpressionInternals::evaluateFunction ( double & exp_r )
{
  int i;

  EXPRerrno = 0;
  if (values_changed_)
  {
    values_changed_ = false;
    ++curr_num_;
  }
  if (curr_num_ > MAX_EVAL)
  {
    clear_eval_num_(tree_);
    if (differentiated_)
    {
      for (i=0 ; i<numVars_ ; ++i)
        clear_eval_num_(derivs_[i]);
    }
    curr_num_ = 1;
  }
  if (DEBUG_EXPRESSION)
  {
    Xyce::dout() << "In evaluateFunction with curr_num_ = " << curr_num_ << std::endl
                 << " evaluating expression " << get_expression() << std::endl;
  }

  EXPReval_( *tree_, exp_r, var_vals_ );

  if (DEBUG_EXPRESSION)
  {
    Xyce::dout() << "Evaluated value = " << exp_r << std::endl;
  }

  return EXPRerrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_sim_time
// Purpose       : Set 'time' special variable in expression
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Set the special variable "TIME"
///
/// @param[in] time  The time value to set
///
/// The internal variable "TIME" is special, and should not be set directly
/// using the "set_var" function alone.  This function also sets some internal
/// state variables needed for proper handling of time-dependent expressions.
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///
bool ExpressionInternals::set_sim_time(double time)
{
  if ( time != sim_time_ )
  {
    values_changed_ = true;
  }

  if (set_var(std::string("TIME"), time))
  {
    values_changed_ = true;
  }
  sim_time_ = time;
  tentative_accepted_time = time;

  return values_changed_;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_sim_dt
// Purpose       : Set 'dt' special variable in expression
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// Set the special variable "DT"
///
/// @param[in] time  The time value to set
///
/// The internal variable "DT" is special, and should not be set directly
/// using the "set_var" function alone.  This function also sets some internal
/// state variables needed for proper handling of time-dependent expressions.
///
/// @author Eric Keiter
/// @date 
///
bool ExpressionInternals::set_sim_dt(double dt)
{
  bool retVal=false;
  if ( dt != sim_dt_ )
  {
    retVal = true;
  }

  sim_dt_ = dt;
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_temp
// Purpose       : Set 'temp' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/11/06
//-----------------------------------------------------------------------------
///
/// Set the special variable "TEMP"
///
/// @param[in] temp  The temperature value to set in Celsius
///
/// The internal variable "TEMP" is special, and should not be set directly
/// using the "set_var" function alone.  This function also sets some internal
/// state variables needed for proper handling of time-dependent expressions.
///
/// This function also sets the "VT" variable.
///
/// @note The TEMP variable is stored internally in Kelvin, and the conversion
/// from Celsius to Kelvin is handled in this function.
///
/// @author Dave Shirley, PSSI
/// @date 09/11/06
///
bool ExpressionInternals::set_temp(double const & tempIn)
{
  set_var(std::string("VT"), tempIn*CONSTKoverQ);

  return set_var(std::string("TEMP"), tempIn-CONSTCtoK);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_sim_freq
// Purpose       : Set 'freq' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 09/7/18
//-----------------------------------------------------------------------------
///
/// Set the special variable "FREQ"
///
/// @param[in] freq The frequency value 
///
/// The internal variable "FREQ" is special, and should not be set directly
/// using the "set_var" function alone.  This function also sets some internal
/// state variables needed for proper handling of time-dependent expressions.
///
/// @author Eric Keiter, SNL
/// @date 09/7/18
///
bool ExpressionInternals::set_sim_freq(double freqIn)
{
  return (set_var(std::string("FREQ"), freqIn));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_accepted_time
// Purpose       : Set accepted time for converged solution
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/26/05
//-----------------------------------------------------------------------------
///
/// Set accepted time for converged solution
///
/// Computation of the sdt and ddt function depends on keeping track of
/// the time of converged, accepted simulation results.  set_accepted_time
/// must be called when a time step is accepted by the time integrator.
///
/// @author Dave Shirley, PSSI
/// @date 04/26/05
///
void ExpressionInternals::set_accepted_time(double const time)
{
  if (accepted_time != time)
  {
    values_changed_ = true;
  }
  accepted_time = time;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_break_time
// Purpose       : Returns next breakpoint time
//
// Special Notes : ERK.  The usage of "set_sim_time" is very odd in
//                 this function.  In most cases, sim_time_ is set directly,
//                 and then "set_sim_time_" is called immediately afterwards.
//                 This should probably be refactored.
//
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Returns the next breakpoint time
///
/// @return next breakpoint
///
/// @note The expression library attempts to identify discontinuities in
/// explicitly time-dependent expressions, and report these times on request
/// so that the time integrator can handle them properly.  This is done
/// by constructing a breakpoint expression that must be evaluated at the
/// current simulation time.  
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///
double ExpressionInternals::get_break_time_i()
{
  double bp, bp2, rval, rval_old, delta;
  double old_time, val, t_limit;
  int i;
  std::vector<ExpressionNode *>::iterator break_i, best_break_i;

  if (!breakpointed_)
  {
    if (num_string_ > 0 || num_func_ > 0)
    {
      Report::UserWarning() << "Failed to obtain breakpoint in: " << Input_ << std::endl
                            << "Expression is not fully resolved";
      return 0;
    }
    if (time_index == -1)
    {
      time_index = find_num_("TIME");
      if (time_index == -1)
      {
        time_index = -2;
        return 0;
      }
    }
    if (!differentiated_)
      simplify_ (*tree_);
    breaks_.clear();
    get_breaks_ (*tree_);
    breakpointed_ = true;

    if (DEBUG_EXPRESSION) {
      Xyce::dout() << "Generated " << breaks_.size() << " breakpoint expression";
      if (breaks_.size() != 1)
        Xyce::dout() << "s";
      if (breaks_.size() > 0)
        Xyce::dout() << ":";
      Xyce::dout() << std::endl;
      for (break_i = breaks_.begin() ; break_i != breaks_.end() ; ++break_i)
      {
        std::ostringstream s("");
        s << std::setprecision(5);
        RpTree_ (*break_i, s);
        Xyce::dout() << s.str() << std::endl;
      }
    }
  }
// Now, theoretically we have a list of expressions that will
// provide linear approximations of the breakpoints, based on the
// current time.  Find the best breakpoint and expression.

  rval = -1;
  ++curr_num_;

  old_time = sim_time_;
  if (old_time < 1)
    delta = Epsilon;
  else
    delta = old_time*Epsilon;
//    Xyce::dout() << std::setprecision(15);

  for (break_i = breaks_.begin(), i=0 ; break_i != breaks_.end() ;
        ++break_i, ++i)
  {
    sim_time_ = old_time+delta;
    set_sim_time(sim_time_);
    EXPReval_ (**break_i, bp, var_vals_);
    if (bp != EXPR_HUGE && bp >= 0)
    {
      if (bp < delta)
      {
        sim_time_ += 100*delta;
        set_sim_time(sim_time_);
        ++curr_num_;
        EXPReval_ (**break_i, bp, var_vals_);
        sim_time_ = old_time;
        set_sim_time(sim_time_);
      }
      if (bp >= delta)
      {
        if (rval == -1 ||( rval > bp && bp > 0))
        {
          best_break_i = break_i;
          rval = bp;
        }
      }
    }
  }
  sim_time_ = old_time;
  if (rval == -1)
    return 0;

  sim_time_ = old_time + rval;
  rval += delta;
  bp = 1;
  i = 0;
  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Before refinement, breakpoint = " << rval << " at time = "
                 << old_time << std::endl;

  sim_time_ -= delta;
  set_sim_time(sim_time_);
  (void) evaluateFunction(val);
  EXPReval_ (**best_break_i, bp, var_vals_);
//    Xyce::dout() << "del = " << fabs(sim_time_+bp - (old_time + rval)) << std::endl;
  if (fabs(sim_time_+bp - (old_time + rval)) >= delta)
  {
    while (fabs(bp) > Epsilon && i<20)
    {
      set_sim_time(sim_time_);
      (void) evaluateFunction(val);
      EXPReval_ (**best_break_i, bp, var_vals_);
//        Xyce::dout() << "time = " << sim_time_ << "  bp = " << bp << std::endl;
      sim_time_ += bp;
      ++i;
    }
    rval = sim_time_;
    sim_time_ = old_time;
    set_sim_time(sim_time_);
  }
  else
  {
    rval = old_time + rval;
    sim_time_ = rval;
    set_sim_time(sim_time_);
  }
  if (DEBUG_EXPRESSION)
    Xyce::dout() << "After refinement, breakpoint = " << rval << std::endl;

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_input
// Purpose       : Return expression input string
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Return expression input string
///
/// @return original input string
///
/// @note This is primarily a debugging function, and returns the original
/// string from which this expression was parsed.  Because the expression may
/// later have been simplified or modified (e.g. to handle POLY or TABLE
/// standardization), it is not necessarily an exact representation of the
/// expression as it is actually stored in its internal tree.  The
/// "get_expression" function should be used to view the string representation
/// of the actual tree.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
const std::string & ExpressionInternals::get_input ()
{
  return (Input_);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::order_names
// Purpose       : Put input quantity names in a particular order (used for
//                 replace_func which requires identical ordering for expression
//                 and user defined function
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Reorder input quantity names
///
/// @param[in] new_names   Vector of names in desired order
/// @return error code
///
/// This function reorders the variables as stored in the "symbol table,"
/// which is really four parallel arrays, varValues_, varTypes_,
/// leadDesignator_, and var_vals_.
///
/// Because all of the use of the var_vals_ array (storing the actual
/// values of input quantities) depends on them being provided in the
/// order matching the varValues_ (names) array, and because the users
/// of the expression might want to provide them in in a specific order
/// determined by some other use, this function allows a caller to
/// reorganize the "symbol table" in the desired order.
///
/// This operation is intrusive, as variables are stored in the expression
/// tree by noting their index into the "symbol table".  Reordering the
/// table requires fixing all the indicies in the tree.
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///
int ExpressionInternals::order_names(std::vector<std::string> const & new_names)
{
  std::vector<std::string>::iterator old_i;
  std::vector<std::string>::const_iterator new_i;
  std::string dum;
  int i, j, k, found;
  int istat, n_args;
  std::vector<int> Vmap;

  std::vector<int> t_varTypes_;
  std::vector<std::string> t_varValues_;
  std::string t_leadDesignator_;
  std::vector<double> t_var_vals_;
  std::vector<ExpressionNode *> t_derivs_;
  std::vector<int> nmap;


  if (DEBUG_EXPRESSION) {
    std::vector<std::string>::const_iterator n_i, n_end;
    n_i = new_names.begin();
    n_end = new_names.end();
    Xyce::dout() << "At start of order_names with name list:" << std::endl;
    for ( j=0; n_i != n_end ; ++n_i,j++)
      Xyce::dout() << "   New " << j << ": "  << *n_i << std::endl;
    Xyce::dout() << "Initial expression:" << std::endl
                 << get_expression() << std::endl;
    for (old_i = varValues_.begin(),j=0;old_i != varValues_.end(); ++old_i,j++)
      Xyce::dout() << "   Original " << j << ": " << *old_i << std::endl;
  }

  values_changed_ = true;
  n_args = new_names.size();

  Vmap.resize(n_args);

  found = 0;
  for (new_i=new_names.begin(), j=0; new_i!=new_names.end() ; ++new_i, ++j)
  {
    k = (*new_i).size();
    for (old_i=varValues_.begin(), i=0; old_i!=varValues_.end() ; ++old_i, ++i)
    {
      if ((*new_i)[k-1] == '}')
      {
        if ((*new_i).substr(0,k-3) == *old_i && (*new_i)[k-2] == leadDesignator_[i])
        {
          Vmap[j] = i;
          ++found;
          break;
        }
      }
      else
      {
        if (*new_i == *old_i)
        {
          Vmap[j] = i;
          ++found;
          break;
        }
      }
    }
    if (found == j)
    {
      std::ostringstream s("");
      s << "__EXPRdummyString_ForUserFuncArg::";
      s << j;
      dum = s.str();
      istat = EXPRaddDummyString_ (dum);
      Vmap[j] = istat;
      ++found;
    }
  }
  // Vmap now contains a map of where new names are to be placed into varValues_
  if (numVars_ > n_args)
  {
    std::vector<bool> mapped(numVars_, false);
    for (i=0 ; i<n_args ; ++i)
      mapped[Vmap[i]] = true;
    for (i=0 ; i<numVars_ ; ++i)
    {
      if (!mapped[i])
        Vmap.push_back(i);
    }
  }
  found = 0;
  nmap.resize(numVars_);
  for (i=0 ; i<numVars_ ; ++i)
    nmap[i] = i;
  for (i=0 ; i<numVars_ ; ++i)
  {
    for (j=0 ; j<numVars_ ; ++j)
    {
      if (Vmap[j] == i)
      {
        nmap[Vmap[found++]] = i;
        break;
      }
    }
  }

  t_varTypes_.resize(numVars_);
  t_varValues_.resize(numVars_);
  t_leadDesignator_.resize(numVars_);
  t_var_vals_.resize(numVars_);
  if (differentiated_)
    t_derivs_.resize(numVars_);
  for (i=0 ; i<numVars_ ; ++i)
  {
    t_varTypes_[i] = varTypes_[i];
    t_varValues_[i] = varValues_[i];
    t_leadDesignator_[i] = leadDesignator_[i];
    t_var_vals_[i] = var_vals_[i];
    if (differentiated_)
      t_derivs_[i] = derivs_[i];
  }
  for (i=0 ; i<numVars_ ; ++i)
  {
    varTypes_[nmap[i]] = t_varTypes_[i];
    varValues_[nmap[i]] = t_varValues_[i];
    leadDesignator_[nmap[i]] = t_leadDesignator_[i];
    var_vals_[nmap[i]] = t_var_vals_[i];
    if (differentiated_)
      derivs_[nmap[i]] = t_derivs_[i];
  }

  for (i=0 ; i<2 ; ++i)
  {
    Rmap_(*tree_, i, nmap);
    if (differentiated_)
    {
      for (j=0 ; j<numVars_ ; ++j)
      {
        Rmap_(*(derivs_[j]), i, nmap);
      }
    }
  }

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Final expression:" << std::endl
                 << get_expression() << std::endl;

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::replace_func
// Purpose       : Replace user defined function with its definition in expression
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Replace a user defined function with its definition
///
/// @param[in] func_name The function name to replace
/// @param[in] func_def  The expression tree representing the function definition
/// @param[in] numArgs   The number of arguments taken by the function
/// @return error code
///
/// The function's name is located in the symbol table, and all nodes
/// referencing the function are replaced by the tree of the function
/// definition.  This is a recursive and intrusive operation.
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///
int ExpressionInternals::replace_func (std::string const & func_name,
                                   ExpressionInternals & func_def,
                                   int numArgs)
{
// This initiates the insertion of a previously parsed function into
// an expression.  The routine Freplace_ will search for
// instances of the function and initiate replacement at each of
// these.
  int i;
  int Ferrno;

  if (differentiated_)
  {
    Report::UserError() << "Attempt to do replacement in differentiated expression: " << get_expression();
    Ferrno = -1;
    return Ferrno;
  }

  values_changed_ = true;
  Ferrno = 0;
  // find the name of this function in the symbol table.
  for (i=0 ; i<numVars_ ; ++i)
    if (varTypes_[i] == EXPR_T_FUNCTION && varValues_[i] == func_name)
      break;
  if (i < numVars_)
  {

    Ferrno=Freplace_(tree_, func_name, func_def, numArgs);

    numVars_--;
    varTypes_[i] = varTypes_[numVars_];
    varValues_[i] = varValues_[numVars_];
    leadDesignator_[i] = leadDesignator_[numVars_];
    var_vals_[i] = var_vals_[numVars_];

    RemoveFentry_ (tree_, i, numVars_);

    varTypes_.resize(numVars_);
    varValues_.resize(numVars_);
    leadDesignator_.resize(numVars_);
    var_vals_.resize(numVars_);
  }

  set_nums_();

  if (!ddxProcessed_ && num_string_ == 0 && num_func_ == 0)
  {
    tree_ = diffDDX_ (PThead_);
    PThead_ = tree_;
  }

  return Ferrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::replace_var
// Purpose       : Replace a variable with a previously parsed subexpression
// Special Notes : This is primarily needed to handle a case where a subcircuit
//                 parameter is given a definition in terms of a global param,
//                 in which case it resolves to an expression instead of a value
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 08/10/2010
//-----------------------------------------------------------------------------
///
/// Replace a variable with a previously parsed subexpression
///
/// @param[in] varName The variable name to replace
/// @param[in] subexpr  The expression tree to replace the variable
/// @return error code
///
/// Similar in function to replace_func, this method recursively traverses
/// the expression tree and substitutes the given expression tree for any
/// expression node representing the variable.
///
/// The primary purpose of this method is to handle cases where subcircuit
/// parameters are given in terms of a ".global_param" that is itself an
/// expression.
///
/// @author Tom Russo, SNL
/// @date 08/10/2010
///
int ExpressionInternals::replace_var(
  const std::string &           varName,
  const ExpressionInternals &   subexpr)
{
// This initiates the insertion of a previously parsed function into
// an expression.  The routine Vreplace_ will search for
// instances of the variable and initiate replacement at each of
// these.
  int i;
  int Ferrno;

  if (differentiated_)
  {
    Report::UserError() << "Attempt to do replacement in differentiated expression: " << get_expression();
    Ferrno = -1;
    return Ferrno;
  }

  values_changed_ = true;
  Ferrno = 0;
  // find the name of this function in the symbol table.
  for (i=0 ; i<numVars_ ; ++i)
    if ((varTypes_[i] == EXPR_T_VARIABLE || varTypes_[i] == EXPR_T_STRING)
        && varValues_[i] == varName)
      break;
  if (i < numVars_)
  {

    Ferrno=Vreplace_(tree_, varName, subexpr);

    numVars_--;
    varTypes_[i] = varTypes_[numVars_];
    varValues_[i] = varValues_[numVars_];
    leadDesignator_[i] = leadDesignator_[numVars_];
    var_vals_[i] = var_vals_[numVars_];

    RemoveFentry_ (tree_, i, numVars_);

    varTypes_.resize(numVars_);
    varValues_.resize(numVars_);
    leadDesignator_.resize(numVars_);
    var_vals_.resize(numVars_);
  }

  set_nums_();

  if (!ddxProcessed_ && num_string_ == 0 && num_func_ == 0)
  {
    tree_ = diffDDX_ (PThead_);
    PThead_ = tree_;
  }

  return Ferrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::replace_var
// Purpose       : Replace a variable with an operator call
// Special Notes : This is done by creating a new ExpressionNode referencing
//                 the op, and using the replace_var(name,subexpr) function
//                 to do the real replacement of the name in the current
//                 expression.
//                 
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 02/23/2015
//-----------------------------------------------------------------------------
///
/// Replace a variable with an op
///
/// @param[in] var_name The variable name to replace
/// @param[in] op  The operator to evaluate the variable
/// @return error code
///
/// In the case of expressions involving the GMIN or VT global simulator
/// variables, these values are obtained using "Ops" from the IO Interface
/// package.  Rather than set these values time and again using set_var
/// calls, the variable nodes are replaced with operator call nodes, and
/// when the expression is evaluated the operator is called.
///
/// @author Dave Baur, Raytheon
/// @date 02/23/2015
///
int
ExpressionInternals::replace_var(
  std::string const &           var_name,
  Op::Operator *                op)
{
  Util::ExpressionNode node;
  node.type = Util::EXPR_FUNCTION;
  node.funcnum = Util::EXPR_F_OP;
  node.funcname = var_name;
  node.valueIndex = 0;
  node.op = op;
  node.fptr.op = op_eval;

  if (replace_var(var_name, ExpressionInternals(node)) != 0)
  {
    Report::UserWarning0() << "Problem inserting operator as substitute for " << var_name << " in expression " << Input_;
  }

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::replace_name
// Purpose       : Change the name of an input quantity
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/28/01
//-----------------------------------------------------------------------------
///
/// Rename an input quantity
///
/// @param[in] old_name  The name to replace
/// @param[in] new_name  The name with which to replace it
/// @return  true if old_name found, false if not
///
/// If the new name is zero, the variable is converted to a constant.
///
/// @note This method primarily used to handle replacement of node names
/// resulting from node aliases and subcircuiting hierarchy.
///
/// @author Dave Shirley, PSSI
/// @date 08/28/01
///
bool ExpressionInternals::replace_name ( const std::string & old_name,
                                      const std::string & new_name)
{
  int ind;

  values_changed_ = true;
  if ( (ind = find_num_( old_name )) >= 0)
  {
    if (new_name == "0")
      convert_to_constant_ (ind, 0.);
    else
      varValues_[ind] = new_name;
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::getNumDdt
// Purpose       : Return the number of ddt() calls in the expression
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Return the number of ddt() calls in the expression
///
/// @return number of ddt() functions in the expression
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
int ExpressionInternals::getNumDdt ( )
{
  numDDT = 0;
  RcountDDT_ (*tree_);
  return numDDT;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::getDdtVals
// Purpose       : Return the most recent arguments of ddt() in the expression
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Return the most recent arguments to ddt() in the expression
///
/// @param[out] vals  Values of arguments to each ddt call in the expression
///
/// While ddt() in an expression normally computes time derivatives internally,
/// for expressions used in devices (e.g. a B source or switch), it is
/// desirable to compute these terms using the time integrator.  These devices
/// can therefore call getDdtVals to get the value of the terms that are
/// supposed to be differentiated, load them approprately into the state
/// vector, retrieve their derivatives from the time integrator, and
/// set them directly in the expression using setDdtDerivs calls.
///
/// Doing so bypasses the internal computation of ddt() terms in expressions,
/// and allows derivatives computed to be done at the same order of
/// approximation as the simulator's time integration method.
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
void ExpressionInternals::getDdtVals ( std::vector<double> & vals )
{
  double val;
  numDDT = 0;
  if (values_changed_)
    evaluateFunction(val);
  RgetDDT_ (*tree_, vals);
  if (numDDT != vals.size())
  {
    Report::DevelFatal().in("ExpressionInternals::getDdtVals") << "Length of return vector inconsistent";
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::setDdtDerivs
// Purpose       : Set the evaluated value of the ddt functions
// Special Notes : This is normally done with derivative values from the
//                 time integration package
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Bypass internal ddt evaluation using values from the time integrator
///
/// @param[in] vals  Values of derivatives to substitute for ddt()
///
/// While ddt() in an expression normally computes time derivatives internally,
/// for expressions used in devices (e.g. a B source or switch), it is
/// desirable to compute these terms using the time integrator.  These devices
/// can therefore call getDdtVals to get the value of the terms that are
/// supposed to be differentiated, load them approprately into the state
/// vector, retrieve their derivatives from the time integrator, and
/// set them directly in the expression using setDdtDerivs calls.
///
/// Doing so bypasses the internal computation of ddt() terms in expressions,
/// and allows derivatives computed to be done at the same order of
/// approximation as the simulator's time integration method.
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
void ExpressionInternals::setDdtDerivs ( std::vector<double> & vals )
{
  int curr_num_old;

  numDDT = 0;
  curr_num_old = curr_num_;
  ++curr_num_;
  if (curr_num_ > MAX_EVAL)
    curr_num_ = 1;
  RsetDDT_ (*tree_, vals);
  curr_num_ = curr_num_old;
  values_changed_ = true;

  if (numDDT != vals.size())
  {
    Report::DevelFatal().in("ExpressionInternals::setDdtVals") << "Length of return vector inconsistent";
  }
}
// Private methods:

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::find_num_
// Purpose       : Find the index of an input quantity
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Find the index of an input quantity in the symbol table
///
/// @param[in] var Name of variable
/// @return index of variable in symbol table
///
/// @note See notes for get_names, get_type, and getSymbolTable: the
/// "variable name" stored for lead currents such as IG(M1) or ID(M1)
/// is just the device instance name, and therefore these are not
/// unique if more than one current for the same device is present in
/// the expression.  find_num_ on that name will therefore find the
/// index of the FIRST such instance, and not necessarily the one
/// intended.  Care must be used when using any function that calls
/// find_num_ on a variable name when lead currents are involved.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
int ExpressionInternals::find_num_( const std::string & var )
{
  for (int i = 0; i < numVars_; ++i )
  {
    if (varValues_[i] == var)
    {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::set_nums_
// Purpose       : Set tally of input quantity types
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/12/01
//-----------------------------------------------------------------------------
///
/// Set tally of input quantity types
///
/// This function loops over the symbol table and makes a tally of each type
/// of input quantity.  These are stored in member variables such as num_N_,
/// num_lead_, num_I_, etc.
///
/// @author Dave Shirley, PSSI
/// @date 07/12/01
///
void ExpressionInternals::set_nums_()
{
  num_N_ = num_I_ = num_lead_ = num_string_ = num_special_ = num_var_ = num_func_ = num_node_computation_ = 0;

  for (int i = 0;i<numVars_;++i)
  {
    switch (varTypes_[i])
    {
      case EXPR_T_NODE:
        ++num_N_;
        break;
      case EXPR_T_INSTANCE:
        if (leadDesignator_[i] == ' ')
        {
          std::string s = varValues_[i];
          std::string::size_type pos = s.find_last_of(":");
          if( pos == std::string::npos )
            pos = 0;
          else
            ++pos;
          if (s[pos] == 'V')
            ++num_I_;
          else
            ++num_lead_;
        }
        else
          ++num_lead_;
        break;
      case EXPR_T_STRING:
        ++num_string_;
        break;
      case EXPR_T_SPECIAL:
        ++num_special_;
        break;
      case EXPR_T_VARIABLE:     // not possible at this point
        break;
      case EXPR_T_FUNCTION:
        ++num_func_;
        break;
      case EXPR_T_NODAL_COMPUTATION:
        ++num_node_computation_;
        break;
      default:
        Report::DevelFatal().in("ExpressionInternals::set_nums_") << "Bad variable type found in expression";
        break;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::create_vars_
// Purpose       : Create vector for input quantity values
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/07/01
//-----------------------------------------------------------------------------
///
/// Create vector for input quantity values
///
/// Resizes var_vals_ array and initializes to zero.
///
/// @author Dave Shirley, PSSI
/// @date 06/07/01
///
void ExpressionInternals::create_vars_ ()
{
  var_vals_.resize( numVars_ );

  values_changed_ = true;
  for ( int i = 0 ; i < numVars_; ++i )
    var_vals_[i] = 0;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::copy_elements_
// Purpose       : Copy a list of ExpressionElements
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Copy a list of ExpressionElements
///
/// @param to  the destination list
/// @param[in] from the source list
///
/// @note This implementation is very C like, because the lists are of
/// ExpressionElement pointers.  The entire thing should probably be rewritten
/// in a more C++ style, so that ExpressionElement copy constructors and
/// STL list appends can be used in a more straightforward manner.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::copy_elements_ (std::list<ExpressionElement *> &to,
                                      std::list<ExpressionElement *> *from)
{
  ExpressionElement *el;
  std::list<ExpressionElement *>::iterator exp_i;

  for (exp_i = from->begin() ; exp_i != from->end() ; ++exp_i)
  {
    el = newExpressionElement_();
    copy_element_ (el, *exp_i);
    to.push_back(el);
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::copy_element_
// Purpose       : Copy a single ExpressionElement
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Copy an ExpressionElement
///
/// @param to  Pointer to the the destination element
/// @param[in] Pointer to the source element
///
/// @note This implementation is very C like, because the lists are of
/// ExpressionElement pointers.  The entire thing should probably be rewritten
/// in a more C++ style, so that ExpressionElement copy constructors and
/// STL list appends can be used in a more straightforward manner.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::copy_element_ (ExpressionElement *to, ExpressionElement *from)
{
  to->token = from->token;
  if (from->token == TOK_VALUE)
  {
    to->type = from->type;
    if (from->type == TYP_NUM)
    {
      to->number = from->number;
    }
    else if (from->type == TYP_STRING)
    {
      to->name = from->name;
    }
    else if (from->type == TYP_PNODE)
    {
      to->node = from->node;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::makepnode_
// Purpose       : Create a ExpressionNode from a ExpressionElement
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Create an expression node from an expression element
///
/// @param[in] elem Pointer to an ExpressionElement
/// @return pointer to a new expression node
///
/// @note During initial parsing the input string is converted to a list
/// of tokens of type ExpressionElement.  These must ultimately be processed
/// into tree nodes by this method.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::makepnode_ (ExpressionElement *elem)
{
  if (elem->token != TOK_VALUE)
    return (NULL);

  switch (elem->type)
  {
    case TYP_STRING:
      return (mksnode_(elem->name));

    case TYP_NUM:
      return (mkcon_(elem->number));

    case TYP_PNODE:
      return (elem->node);

    default:
      Report::DevelFatal().in("ExpressionInternals::makepnode_") << "Internal error, unknown type";
      return 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mkfnode_
// Purpose       : Create a function ExpressionNode from a function name and args
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Create a function ExpressionNode from a function name and arguments
///
/// @param[in] fname function name
/// @param[in] num_args number of arguments
/// @param[in] args Vector of ExpressionNode pointers representing args
/// @return ExpressionNode representing the function call
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode *
ExpressionInternals::mkfnode_(
  const std::string &                   fname,
  int                                   num_args,
  std::vector<ExpressionNode *>         args)      // Function node
{
  int i;
  ExpressionNode *p, *expp, *sinpp;
  std::string name;
  ExpressionNode *arg1, *e1, *e2;
  ExpressionNode *newp, *newp2, *newp_if, *newp_if1, *newp_if2;

  if (fname == "PWR" || fname == "POW")
  {
    if (num_args != 2) return NULL;
    return mkb_(EXPR_POWER, args[0], args[1]);
  }
  if (fname == "LOG10")
  {
    if (num_args != 1) return NULL;
    return mkf_(EXPR_F_LOG, args[0]);
  }
  if (fname == "LIMIT")
  {
    if (num_args == 3)
    {
      newp = mkf_(EXPR_F_IF, mkb_(EXPR_LESS, args[0], args[1]));
      newp->operands[1] = args[1];
      newp2 = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, args[0], args[2]));
      newp2->operands[1] = args[2];
      newp2->operands[2] = args[0];
      newp->operands[2] = newp2;
      return newp;
    }
    else if (num_args == 2)
    {
      newp = mkb_(EXPR_PLUS,args[0],args[1]);
      return newp;
    }
    else
      return NULL;
  }
  if (fname == "SIGN")
  {
    if (num_args != 2) return NULL;
    newp = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, args[1], mkcon_(0.0)));
    newp->operands[1]= mkf_(EXPR_F_ABS,args[0]);
    newp2= mkf_(EXPR_F_IF, mkb_(EXPR_LESS, args[1], mkcon_(0.0)));
    newp2->operands[2]=mkcon_(0.0);
    newp2->operands[1]= mkb_(EXPR_TIMES,newp->operands[1],mkcon_(-1.0));
    newp->operands[2]=newp2;
    return newp;
  }
  if (fname == "STP")
  {
    if (num_args != 1) return NULL;
    newp = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, args[0], mkcon_(0.0) ));
    newp->operands[1] = mkcon_(1.0) ;
    newp->operands[2] = mkcon_(0.0) ;
    return newp;
  }
  if (fname == "PWRS")
  {
    if (num_args != 2) return NULL;
    newp = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, args[0], mkcon_(0.0)));
    newp->operands[1] = mkb_(EXPR_POWER, args[0], args[1]);
    newp2 = mkf_(EXPR_F_IF, mkb_(EXPR_LESS, args[0], mkcon_(0.0)));
    newp2->operands[1] = mkf_(EXPR_F_UMINUS, mkb_(EXPR_POWER, mkf_(EXPR_F_UMINUS, args[0]), args[1]));
    newp2->operands[2] = mkcon_(0.0);
    newp->operands[2] = newp2;
    return newp;
  }
  if (fname == "ARCTAN")
  {
    if (num_args != 1) return NULL;
    return mkf_(EXPR_F_ATAN, args[0]);
  }
  if (fname == "ATAN2")
  {
    if (num_args != 2) return NULL;
    return mkf_(EXPR_F_ATAN, mkb_(EXPR_DIVIDE, args[0], args[1]));
  }
  if (fname == "M")
  {
    if (num_args != 1) return NULL;
    return mkf_(EXPR_F_ABS, args[0]);
  }
  if (fname == "MIN")
  {
    if (num_args != 2) return NULL;
    newp = mkf_(EXPR_F_IF, mkb_(EXPR_LESS, args[0], args[1]));
    newp->operands[1] = args[0];
    newp->operands[2] = args[1];
    return newp;
  }
  if (fname == "MAX")
  {
    if (num_args != 2) return NULL;
    newp = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, args[0], args[1]));
    newp->operands[1] = args[0];
    newp->operands[2] = args[1];
    return newp;
  }

  if (fname == "SPICE_PULSE")
  {
//    interpret old style spice pulse function:
//       PULSE (V1, V2, TD, TR, TF, PW, PER)
//         V1 = initial value
//         V2 = pulsed value
//         TD = delay
//         TR = rise time
//         TF = fall time
//         PW = pulse width
//         PER = period
    if (num_args < 5)
    {
      Report::UserError() << "Fewer than 5 args in pulse function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args > 7)
    {
      Report::UserError() << "More than 7 args in pulse function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
// periodic pulse case, include remainder in first argument of table
// if(time >= TD, TABLE((time - TD) % PER, t*, v*), V1)
    if (num_args == 7)
      newp = mkf_(EXPR_F_TABLE, mkb_(EXPR_REMAINDER, mkb_(EXPR_MINUS, mksnode_(std::string("TIME")),
                  args[2]), args[6]));
    else
      newp = mkf_(EXPR_F_TABLE, mkb_(EXPR_MINUS, mksnode_(std::string("TIME")), args[2]));

    newp->operands[1] = mkcon_(0);                           //  0, V1
    newp->operands[2] = args[0];
    newp->operands.push_back(args[3]);                      //  TR, V2
    newp->operands.push_back(args[1]);
    if (num_args >= 6)                                      // if pulse width is given
    {
      newp->operands.push_back(mkb_(EXPR_PLUS, args[3], args[5]));   //  TR + PW, V2
      newp->operands.push_back(args[1]);
      newp->operands.push_back(mkb_(EXPR_PLUS, args[3], mkb_(EXPR_PLUS, args[4],
            args[5])));                                     //  TR + PW + TF, V1
      newp->operands.push_back(args[0]);
    }
    return newp;
  }
  if (fname == "SPICE_SIN")
  {
//    interpret old style spice sin function:
//       SIN (V0, VA, FREQ, TD, THETA)
//         V0 = offset
//         VA = amplitude
//         FREQ = frequency (hz)
//         TD = delay
//         THETA = damping factor
    if (num_args < 3)
    {
      Report::UserError() << "Fewer than 3 args in sin function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args > 5)
    {
      Report::UserError() << "More than 5 args in sin function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args < 5)
    {
      expp = mkcon_(1);
    }
    else
    {
      expp = mkf_ (EXPR_F_EXP, mkb_(EXPR_TIMES, mkb_ (EXPR_MINUS, args[3], mksnode_(std::string("TIME"))), args[4]));
    }
    if (num_args < 4)
    {
      sinpp = mkf_ (EXPR_F_SIN, mkb_ (EXPR_TIMES, mkb_ (EXPR_TIMES, mkcon_(2*M_PI), args[2]),
                mksnode_(std::string("TIME"))));
    }
    else
    {
      sinpp = mkf_ (EXPR_F_SIN, mkb_ (EXPR_TIMES, mkb_ (EXPR_TIMES, mkcon_(2*M_PI), args[2]),
                mkb_ (EXPR_MINUS, mksnode_(std::string("TIME")), args[3]) ));
    }
    if (num_args >= 4)
    {
      newp_if = mkf_(EXPR_F_IF, mkb_(EXPR_GREATEQ, mksnode_(std::string("TIME")), args[3]));
      newp_if->operands[1] = mkb_(EXPR_TIMES, args[1], mkb_(EXPR_TIMES, expp, sinpp));
      newp_if->operands[2] = mkcon_(0.);
    }
    else
    {
      newp_if = mkb_(EXPR_TIMES, args[1],sinpp);
    }

    newp = mkb_ (EXPR_PLUS, args[0], newp_if);
    return newp;
  }
  if (fname == "SPICE_EXP")
  {
//    interpret old style spice exp function:
//       EXP (V1, V2, TD1, TAU1, TD2, TAU2)
//         V1 = initial value
//         V2 = pulsed value
//         TD1 = rise delay time
//         TAU1 = rise time constant
//         TD2 = fall delay time
//         TAU2 = fall time constant
    if (num_args < 6)
    {
      Report::UserError() << "Fewer than 6 args in exp function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args > 6)
    {
      Report::UserError() << "More than 6 args in exp function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    e1 = mkf_ (EXPR_F_EXP, mkb_(EXPR_DIVIDE, mkb_(EXPR_MINUS, args[2], mksnode_(std::string("TIME"))), args[3]));
    e2 = mkf_ (EXPR_F_EXP, mkb_(EXPR_DIVIDE, mkb_(EXPR_MINUS, args[4], mksnode_(std::string("TIME"))), args[5]));
    newp_if2 = mkf_(EXPR_F_IF, mkb_(EXPR_GREATEQ, mksnode_(std::string("TIME")), args[4]));
    newp_if2->operands[1] = mkb_(EXPR_MINUS, e2, e1);
    newp_if2->operands[2] = mkb_(EXPR_MINUS, mkcon_(1), e1);
    newp_if1 = mkf_(EXPR_F_IF, mkb_(EXPR_GREATEQ, mksnode_(std::string("TIME")), args[2]));
    newp_if1->operands[1] = mkb_(EXPR_TIMES, mkb_(EXPR_MINUS, args[1], args[0]), newp_if2);
    newp_if1->operands[2] = mkcon_(0);
    newp = mkb_(EXPR_PLUS, args[0], newp_if1);
    return newp;
  }
  if (fname == "SPICE_SFFM")
  {
//    interpret old style spice sffm function:
//       SFFM (V0, VA, FC, MDI, FS)
//         V0 = offset
//         VA = amplitude
//         FC = carrier frequency
//         MDI = modulation index
//         FS = signal frequency
    if (num_args < 5)
    {
      Report::UserError() << "Fewer than 5 args in sffm function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args > 5)
    {
      Report::UserError() << "More than 5 args in sffm function: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    newp = mkb_(EXPR_PLUS, args[0], mkb_(EXPR_TIMES, args[1], mkf_(EXPR_F_SIN,
              mkb_(EXPR_PLUS, mkb_(EXPR_TIMES, mkcon_(2*M_PI), mkb_(EXPR_TIMES,
                args[2], mksnode_(std::string("TIME")))),mkb_(EXPR_TIMES, args[3], mkf_(EXPR_F_SIN,
                  mkb_(EXPR_TIMES, mkcon_(2*M_PI), mkb_(EXPR_TIMES, args[4],
                    mksnode_(std::string("TIME"))))))))));
    return newp;
  }

  p = newExpressionNode_();
  p->operands.resize(num_args);
  for (i=0 ; i<num_args ; ++i)
    p->operands[i] = args[i];


  if (fname == "VR" || fname == "VI" || fname == "VM" || fname == "VP" || fname == "VDB"
      || fname == "IR" || fname == "II" || fname == "IM" || fname == "IP" || fname == "IDB"
      || fname == "N" || fname == "DNI" || fname == "DNO" || fname == "P" || fname == "W")
  {
    // Handle the special frequency-domain output types.
    // In initial implementation, this will ONLY be accessible from
    // the .print line, these won't work in device expressions!

    // We will simply construct the "pretty" name of the request, and let
    // the consumer of this expression take care of setting it up.

    if (DEBUG_EXPRESSION) {
      Xyce::dout() << "processing function node: " << fname << Util::push
                   << std::endl;
      Xyce::dout() << "args:" << Util::push << std::endl;
      for (i=0; i<num_args; i++)
      {
        Xyce::dout()<< args[i]->funcname << std::endl;
      }
      Xyce::dout() << Util::pop << std::endl;
    }

    name = fname + "(";
    for (i=0; i<num_args; i++)
    {
      if (i>0)
        name += ",";
      name += args[i]->funcname;
    }
    name += ")";

    if (DEBUG_EXPRESSION)
      Xyce::dout() << "Pretty name: " << name << std::endl
                   << Util::pop << std::endl;

    for (i = 0; i < numVars_; ++i)
    {
      if (varTypes_[i] == EXPR_T_NODAL_COMPUTATION && varValues_[i] == name)
        break;
    }
    if (i == numVars_)
    {
      varValues_.resize (numVars_ + 1);
      varTypes_.resize (numVars_ + 1);
      leadDesignator_.resize (numVars_ + 1);
      var_vals_.resize (numVars_ + 1);

      varValues_[i] = name;
      varTypes_[i] = EXPR_T_NODAL_COMPUTATION;
      leadDesignator_[i] = ' ';
      var_vals_[i] = 0;
      ++numVars_;
    }
    p->valueIndex = i;
    p->type = EXPR_VAR;
    p->funcname=name;
  }
  else if (fname == "V")
  {
    // Handle V(A) and V(A,B)
    name = "";
    if (args[0]->type == EXPR_PLACEHOLDER)
    {
      name = args[0]->funcname;
    }
    else
    {
      Report::UserError() << "Badly formed node voltage: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (num_args == 2)

    {
      // Change v(a,b) into v(a) - v(b)
      arg1 = mkfnode_(fname, 1, args);
      p = mkb_(EXPR_MINUS, arg1, mkfnode_(fname, 1, args[1]));
    }
    else
    {
      if (name == "0")
      {
        p->type = EXPR_CONSTANT;
        p->constant = 0;
        return (p);
      }
      for (i = 0; i < numVars_; ++i)
      {
        if (varTypes_[i] == EXPR_T_NODE && varValues_[i] == name)
          break;
      }
      if (i == numVars_)
      {
        varValues_.resize (numVars_ + 1);
        varTypes_.resize (numVars_ + 1);
        leadDesignator_.resize (numVars_ + 1);
        var_vals_.resize (numVars_ + 1);

        varValues_[i] = name;
        varTypes_[i] = EXPR_T_NODE;
        leadDesignator_[i] = ' ';
        var_vals_[i] = 0;
        ++numVars_;
      }
      p->valueIndex = i;
      p->type = EXPR_VAR;
      p->funcname = name;
    }
  }
  else if (fname == "I" || (fname.size() == 2 && fname[0] == 'I' && fname[1] != 'F')) // lead current IC, IB, etc
  {
    name = "";
    if (args[0]->type == EXPR_PLACEHOLDER)
    {
      name = args[0]->funcname;
    }
    else
    {
      Report::UserError() << "Badly formed branch current: " << Input_;
      newp = mkcon_(0);
      return newp;
    }
    if (DEBUG_EXPRESSION)
      Xyce::dout() << "getting a device called: " << name << std::endl;

    for (i = 0; i < numVars_; ++i)
    {
      if (varTypes_[i] == EXPR_T_INSTANCE && varValues_[i] == name
          && (fname.size() == 1 ||
              (fname.size() == 2 && leadDesignator_[i] == fname[1])))
        break;
    }
    if (i == numVars_)
    {
      varValues_.resize (numVars_ + 1);
      varTypes_.resize (numVars_ + 1);
      leadDesignator_.resize (numVars_ + 1);
      var_vals_.resize (numVars_ + 1);

      varValues_[i] = name;
      varTypes_[i] = EXPR_T_INSTANCE;
      if (fname.size() == 2)
        leadDesignator_[i] = fname[1];
      else
        leadDesignator_[i] = ' ';
      var_vals_[i] = 0;
      ++numVars_;
    }
    p->valueIndex = i;
    p->type = EXPR_VAR;
  }
  else
  {
    p->type = EXPR_FUNCTION;
    for (i = 0; i < NUM_FUNCS; ++i)
    {
      if (funcs[i].name == fname)
        break;
    }

    if (i == NUM_FUNCS)
    {
      p->funcname = fname;
      p->funcnum = EXPR_F_USER;
      p->fptr.unary = NULL;
      for (i = 0; i < numVars_; ++i)
      {
        if (varTypes_[i] == EXPR_T_FUNCTION && varValues_[i] == fname)
          break;
      }
      if (i == numVars_)
      {
        varValues_.resize (numVars_ + 1);
        varTypes_.resize (numVars_ + 1);
        leadDesignator_.resize (numVars_ + 1);
        var_vals_.resize (numVars_ + 1);

        varValues_[i] = fname;
        varTypes_[i] = EXPR_T_FUNCTION;
        leadDesignator_[i] = ' ';
        var_vals_[i] = 0;
        ++numVars_;
      }
      p->valueIndex = i;
    }
    else
    {
      p->funcname = funcs[i].name;
      p->funcnum = funcs[i].number;
      if (p->funcnum == EXPR_F_TABLE)
      {
        p->operands.resize(num_args);
        if (num_args%2 == 0)
        {
          Report::UserFatal() << "ExpressionInternals::mkfnode_: Internal: Even number of args in table function in " << Input_;
        }
      }
      else if(p->funcnum == EXPR_F_SCHEDULE)
      {
        p->operands.resize(num_args);
        if( num_args%2 == 0 )
        {
          // we've added an implicit argument "TIME", so if the total is now even, then the user
          // didn't give us a good list of pairs.  Issue a error
          Report::UserFatal() << "ExpressionInternals::mkfnode_: Internal: Odd number of args in schedule function in: " << Input_;
        }
      }
      else
      {
        if (num_args != funcs[i].num_args)
        {

          Report::UserFatal()
            << "ExpressionInternals::mkfnode_: number of arguments " << num_args
            << " not correct for function: "
            << funcs[i].name
            <<  "\nin expression: "
            <<  Input_ << std::endl
            << " Required number is " << funcs[i].num_args;
        }
      }
      if (p->funcnum == EXPR_F_SDT || p->funcnum == EXPR_F_DDT)
      {
        timeDependent_ = true;
        p->state.resize(N_INT_STATE,0.0);
      }
      if (p->funcnum == EXPR_F_RAND || p->funcnum == EXPR_F_GAUSS || p->funcnum == EXPR_F_AGAUSS)
      {
        randomDependent_=true;
      }
      p->fptr.unary = funcs[i].funcptr;
    }
  }

  return (p);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mkfnode_
// Purpose       : Make a function ExpressionNode from a single ExpressionNode *
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Make a function ExpressionNode from a single ExpressionNode  pointer
///
/// @param[in] fname Function name
/// @param[in] num_args number of arguments
/// @param[in] n  ExpressionNode pointer
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::mkfnode_ (const std::string & fname, int num_args, ExpressionNode *n)
{
  std::vector<ExpressionNode *> args;

  args.clear();
  args.push_back(n);
  return (mkfnode_ (fname, num_args, args));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mksnode_
// Purpose       : Make a string ExpressionNode
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Make a string ExpressionNode
///
/// @param[in] name String to add
/// @return ExpressionNode pointer for new node
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::mksnode_ (const std::string & name)   // String node
{
  int i, j;
  ExpressionNode *p;

  p = newExpressionNode_();

  // check if it's something special, like TEMP or Time.
  if (isNameSpecial_(name))
  {
    for (j = 0; j < numVars_; ++j)
    {
      if (varTypes_[j] == EXPR_T_SPECIAL && name == varValues_[j])
        break;
    }
    if (j == numVars_)
    {
      varValues_.resize (numVars_ + 1);
      varTypes_.resize (numVars_ + 1);
      leadDesignator_.resize (numVars_ + 1);
      var_vals_.resize (numVars_ + 1);

      varValues_[j] = name;
      varTypes_[j] = EXPR_T_SPECIAL;
      leadDesignator_[j] = ' ';
      var_vals_[j] = 0;
      ++numVars_;
    }
    p->valueIndex = j;
    p->type = EXPR_VAR;
    p->funcname = name;
    return (p);
  }

  bool isConstant=false;
  for (i = 0; i < NUM_CONSTANTS; ++i)
  {
    if (constants[i].name == name)
    {
      isConstant=true;
      break;
    }
  }

  if (!isConstant)
  {
    // save this to be resolved later
    p->type = EXPR_PLACEHOLDER;
    p->funcname = name;
  }
  else
  {
    p->type = EXPR_CONSTANT;
    p->constant = constants[i].value;
  }

  return (p);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mkcon_
// Purpose       : Make a constant ExpressionNode
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Make a constant ExpressionNode
///
/// @param[in] value  constant value
/// @return ExpressionNode pointer to new node
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::mkcon_ (double value)
{
  ExpressionNode *p;

  p = newExpressionNode_();

  p->type = EXPR_CONSTANT;
  p->constant = value;

  return (p);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mkb_
// Purpose       : Make a binary op ExpressionNode with simplification of
//                 constant subexpressions
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Make a binary op ExpressionNode with simplification of constant
/// subexpressions
///
/// @param[in] type Operator type from the EXPR_OPS enum
/// @param[in] left ExpressionNode pointer for left operand
/// @param[in] right ExpressionNode pointer for right operand
/// @return ExpressionNode pointer to new binary op node
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::mkb_ (int type, ExpressionNode *left, ExpressionNode *right)
{
  ExpressionNode *p;
  int i, num_op;

  if (left == NULL || right == NULL)
    return NULL;

  num_op = -1;
  for (i = 0; i < NUM_OPS; ++i)
    if (ops[i].number == type)
    {
      num_op = i;
      break;
    }
  if (num_op == -1)
  {
    Report::DevelFatal() << "ExpressionInternals::mkb_: Internal: bad type";
  }

  if (right->type == EXPR_CONSTANT && left->type == EXPR_CONSTANT)
  {
    EXPRerrno = 0;
    p = mkcon_((ops[num_op].funcptr)(left->constant,right->constant));

    if (EXPRerrno >= EXPR_FATAL)
    {
      Report::DevelFatal()
        << "ExpressionInternals::mkb_: error in evaluation of constant in '"
        << ops[num_op].name
        << "' in expression: " << Input_;
    }
    else if (EXPRerrno >= EXPR_WARNING)
    {
      Report::UserWarning()
        << "ExpressionInternals::mkb_: error in evaluation of constant in '"
        << ops[num_op].name
        << "' in expression: " << Input_;
    }
    return p;
  }
  switch (type)
  {
    case EXPR_TIMES:
      if (left->type == EXPR_CONSTANT && left->constant == 0)
        return left;
      else if (right->type == EXPR_CONSTANT && right->constant == 0)
        return right;
      else if (left->type == EXPR_CONSTANT && left->constant == 1)
        return right;
      else if (right->type == EXPR_CONSTANT && right->constant == 1)
        return left;
      break;

    case EXPR_DIVIDE:
      if (left->type == EXPR_CONSTANT && left->constant == 0)
        return left;
      else if (right->type == EXPR_CONSTANT && right->constant == 1)
        return left;
      break;

    case EXPR_REMAINDER:
      if (left->type == EXPR_CONSTANT && left->constant == 0)
        return left;
      break;

    case EXPR_PLUS:
      if (left->type == EXPR_CONSTANT && left->constant == 0)
        return right;
      else if (right->type == EXPR_CONSTANT && right->constant == 0)
        return left;
      break;

    case EXPR_MINUS:
      if (right->type == EXPR_CONSTANT && right->constant == 0)
        return left;
      else if (left->type == EXPR_CONSTANT && left->constant == 0)
        return mkf_(EXPR_F_UMINUS, right);
      break;

    case EXPR_POWER:
      if (right->type == EXPR_CONSTANT)
      {
        if (right->constant == 0)
          return (mkcon_(1.0));
        else if (right->constant == 1)
          return (left);
      }
      break;
  }

  p = newExpressionNode_();
  p->operands.resize(2);
  p->type = type;
  p->operands[0] = left;
  p->operands[1] = right;

  p->fptr.binary = ops[num_op].funcptr;
  p->funcname = ops[num_op].name;

  return (p);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::mkf_
// Purpose       : Create a single argument function ExpressionNode of a given
//                 type.  Note: only connects first argument.  Additional args
//                 must be filled in by caller, if needed.
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Create a single argument function ExpressionNode of a given
/// type.  Note: only connects first argument.  Additional args
/// must be filled in by caller, if needed.
///
/// @param[in] type function identifier from the funcs array
/// @param[in] arg ExpressionNode pointer to single argument of function
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::mkf_(int type, ExpressionNode *arg)
{
  ExpressionNode *p;
  int i;
  double constval;

  p = newExpressionNode_();

  for (i = 0; i < NUM_FUNCS; ++i)
    if (funcs[i].number == type)
      break;
  if (i == NUM_FUNCS)
  {
    Report::UserFatal()
      << "ExpressionInternals::mkf_: Internal: unknown function type in:\n"
      << Input_;
  }

  if (funcs[i].num_args == 1 && arg->type == EXPR_CONSTANT)
  {
    constval = ((*funcs[i].funcptr) (arg->constant));
    return (mkcon_(constval));
  }

  p->operands.resize(funcs[i].num_args);
  if (funcs[i].number == EXPR_F_SDT || funcs[i].number == EXPR_F_DDT)
  {
    p->state.resize(N_INT_STATE,0.0);
  }
  p->type = EXPR_FUNCTION;
  p->operands[0] = arg;

  p->funcnum = funcs[i].number;
  p->fptr.unary = funcs[i].funcptr;
  p->funcname = funcs[i].name;

  return (p);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::newExpressionNode_
// Purpose       : Allocate a new ExpressionNode
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
/// Allocate a new ExpressionNode
///
/// @return pointer to new node
///
/// This method creates the new node using the default constructor of the
/// ExpressionNode class and adds it to the list of nodes that must be
/// freed when the ExpressionInternals object is destructed.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::newExpressionNode_ ()
{
  ExpressionNode *p;

  p = new ExpressionNode();
  free_list_.push_back(p);

  return p;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::deleteExpressionNode_
// Purpose       : Deletes a ExpressionNode
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/28/04
//-----------------------------------------------------------------------------
///
/// Deletes an ExpressionNode and removes it from the list of nodes to be
/// deleted at destruction
///
/// @param[in] p ExpressionNode pointer to destroy and forget about
///
/// @author Dave Shirley, PSSI
/// @date 09/28/04
///
void ExpressionInternals::deleteExpressionNode_ (ExpressionNode *p)
{
  std::vector<ExpressionNode *>::iterator free_i;

  if (p == NULL)
    return;

  for (free_i = free_list_.begin() ; free_i != free_list_.end() ; ++free_i)
  {
    if (*free_i == p)
    {
      delete *free_i;
      *free_i = free_list_[free_list_.size()-1];
      free_list_.resize(free_list_.size()-1);
      break;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::newExpressionElement_ ()
// Purpose       : Allocate a new ExpressionElement
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
/// Allocate a new ExpressionElement
///
/// @return pointer to new element
///
/// This method creates the new node using the default constructor of the
/// ExpressionElement class and adds it to the list of nodes that must be
/// freed when we are finished with parsing.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionElement * ExpressionInternals::newExpressionElement_ ()
{
  ExpressionElement *p;

  p = new ExpressionElement ();
  ee_list_.push_back(p);

  return p;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::RpTree_
// Purpose       : Recursive method to add string translation of ExpressionNode
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to create string translation of ExpressionNode
///
/// @param[in] pt Root of tree to translate
/// @param s std::ostringstream into which to dump the string translation
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::RpTree_ (const ExpressionNode *pt, std::ostringstream & s) const
{
  int i;

  switch (pt->type)
  {
    case EXPR_CONSTANT:
      s << pt->constant;
    break;

    case EXPR_VAR:
      s << varStr_(pt->valueIndex);
    break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      if (BINARY(pt->operands[0]->type) && prectable[pt->operands[0]->type][pt->type] == L)
        s << "(";
      RpTree_(pt->operands[0], s);
      if (BINARY(pt->operands[0]->type) && prectable[pt->operands[0]->type][pt->type] == L)
        s << ")";
      if (pt->type != EXPR_PLUS || (pt->type == EXPR_PLUS &&
              (pt->operands[1]->type != EXPR_CONSTANT ||
              (pt->operands[1]->type == EXPR_CONSTANT && pt->operands[1]->constant >=0))))
        s << ops[pt->type].name;
      if (BINARY(pt->operands[1]->type) && prectable[pt->type][pt->operands[1]->type] == G)
        s << "(";
      RpTree_(pt->operands[1], s);
      if (BINARY(pt->operands[1]->type) && prectable[pt->type][pt->operands[1]->type] == G)
        s << ")";
      break;

    case EXPR_FUNCTION:
      s << pt->funcname << "(";
      if (pt->operands.size() > 1)
      {
        for (i=0 ; i<pt->operands.size()-1 ; ++i)
        {
          RpTree_(pt->operands[i], s);
          s << ", ";
        }
      }
      if (pt->operands.size() > 0)
      {
        RpTree_(pt->operands[pt->operands.size()-1], s);
      }
      s << ")";
      break;

    case EXPR_PLACEHOLDER:
      s << "P:" << pt->funcname;
      break;

    default:
      s << "oops";
      break;
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::varStr_
// Purpose       : Return string for an input quantity
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Return string for an input quantity
///
/// @param[in] i index of input quantity
/// @return string representation of input quantity name
///
/// @note Used by RpTree_ and in debugging output.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
std::string ExpressionInternals::varStr_(int i) const
{
  std::ostringstream s("");
  s << std::setprecision(PRECISION);

  if (varTypes_[i] == EXPR_T_NODE)
  {
    s << "V(" << varValues_[i] << ")";
  }
  else if (varTypes_[i] == EXPR_T_NODAL_COMPUTATION)
  {
    s << varValues_[i];
  }
  else if (varTypes_[i] == EXPR_T_INSTANCE)
  {
    s << "I";
    if (leadDesignator_[i] != ' ')
      s << leadDesignator_[i];
    s << "(" << varValues_[i] << ")";
  }
  else if (varTypes_[i] == EXPR_T_SPECIAL ||
            varTypes_[i] == EXPR_T_STRING  || varTypes_[i] == EXPR_T_VARIABLE)
  {
    s << varValues_[i];
  }
  else
  {
    s << "V" << i;
  }

  return s.str();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::PTcheck_
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to walk an expression tree and fix up incomplete set-up
///
/// @param p root of tree to check
/// @return normal return: should be same as p.  error return, NULL
///
/// @note Placeholders are set to be variables, common subexpressions are
/// reduced.  If check is failed, the node is saved in the list of
/// objects to be deleted after the check returns.
///
/// This method is called by the constructor and by the
/// ExpressionInternals::set function.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::PTcheck_(ExpressionNode *p)
{
  int i;
  ExpressionNode *arg0, *arg1, *rval;

  if (p == NULL)
  {
    Report::DevelFatal()
      << "ExpressionInternals::PTcheck_: Null pointer";
  }

  switch (p->type)
  {
    case EXPR_PLACEHOLDER:
      for (i = 0; i < numVars_; ++i)
        if (varTypes_[i] == EXPR_T_STRING && varValues_[i] == p->funcname)
          break;

      if (i == numVars_)
      {
        varValues_.push_back(p->funcname);
        if (isNameSpecial_(p->funcname))
          varTypes_.push_back(EXPR_T_SPECIAL);
        else
          varTypes_.push_back(EXPR_T_STRING);
        leadDesignator_.push_back(' ');
        var_vals_.push_back(0);
        ++numVars_;
      }
      p->type = EXPR_VAR;
      p->valueIndex = i;

    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<p->operands.size() ; ++i)
        p->operands[i] = PTcheck_(p->operands[i]);
      break;

    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      p->operands[0] = PTcheck_(p->operands[0]);
      p->operands[1] = PTcheck_(p->operands[1]);
      break;

    default:
      Report::DevelFatal()
        << "ExpressionInternals::PTcheck_: Internal: bad node type";
  }
  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Calling com_expr_ from PTcheck_" << std::endl;

  rval = com_expr_ (PThead_, p);
  if (rval != p)
  {
    if (DEBUG_EXPRESSION) {
      std::ostringstream s("");
      s << std::setprecision(PRECISION);
      if (tree_ != NULL)
        RpTree_ (p, s);
      Xyce::dout() << "Replacing:\n" << s.str() << std::endl;
      std::ostringstream t("");
      t << std::setprecision(PRECISION);
      if (tree_ != NULL)
        RpTree_ (rval, t);
      Xyce::dout() << "With:\n" << t.str() << std::endl;
    }

    done_list_.push_back(p);
  }
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::diffDDX_
// Purpose       : Differentiate ddx() in expression tree
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/12/07
//-----------------------------------------------------------------------------
///
/// Differentiate ddx() in expression tree
///
/// @param[in] p  Root of expression tree to differentiate
/// @return The new parse tree head after all processing is complete
///
/// This method is called to generate the derivative parse trees needed
/// to evaluate ddx functions that appear in the expression.
///
/// @author Dave Shirley, PSSI
/// @date 09/12/07
///
ExpressionNode * ExpressionInternals::diffDDX_(ExpressionNode *p)
{
  ExpressionNode * rval;
  rval = PTdiffDDX_(p);
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::PTdiffDDX_
// Purpose       : Recursive method to differentiate ddx()
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 07/13/07
//-----------------------------------------------------------------------------
///
/// Recursive method to differentiate ddx()
///
/// @param[in] p Root of tree to differentiate recursively
/// @return replacement root of parse tree
///
/// This method recurses through the parse tree based at p.  If it encounters
/// a node of function type with function number "EXPR_F_DDX", it performs
/// a differentiation of the first operand of the ddx function with respect
/// to the variable in the second, and returns the parse tree of the
/// differentiated function.  This returned tree will end up replacing
/// the original function node that contained the ddx() call.
/// 
/// @author Dave Shirley, PSSI
/// @date 07/13/07
///
ExpressionNode * ExpressionInternals::PTdiffDDX_(ExpressionNode *p)
{
  int i;
  ExpressionNode *rval;

  if (p == NULL)
  {
    Report::DevelFatal()
      << "ExpressionInternals::PTdiffDDX_: Null pointer";
  }
  ddxProcessed_ = true;

  rval = p;
  switch (p->type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<p->operands.size() ; ++i)
        p->operands[i] = PTdiffDDX_(p->operands[i]);
      if (p->funcnum == EXPR_F_DDX)
      {
        if (p->operands[1]->type != EXPR_VAR)
        {
          Report::UserFatal()
            << "ExpressionInternals::Differentiate: Attempt to differentiate by non-variable in ddx()";
        }
        rval = Differentiate_(p->operands[0], p->operands[1]->valueIndex);
      }
      break;

    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      p->operands[0] = PTdiffDDX_(p->operands[0]);
      p->operands[1] = PTdiffDDX_(p->operands[1]);
      break;

    default:
      Report::DevelFatal()
        << "ExpressionInternals::PTdiffDDX_: Internal: bad node type";
  }

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::com_expr_
// Purpose       : Recursive method to locate and fix common subexpressions
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/28/04
//-----------------------------------------------------------------------------
///
/// Recursive method to locate and fix common subexpressions
///
/// @author Dave Shirley, PSSI
/// @date 09/28/04
///
ExpressionNode * ExpressionInternals::com_expr_ (ExpressionNode *c, ExpressionNode *p)
{
  int i;
  ExpressionNode *arg0, *arg1;

  if (c == NULL || p == NULL)
  {
    Report::DevelFatal()
      << "ExpressionInternals::com_expr_: Null pointer";
  }
  if (c == p)
    return c;

  switch (c->type)
  {
    case EXPR_CONSTANT:
      if (c->type == p->type && c->constant == p->constant)
        return c;
      break;

    case EXPR_VAR:
      if (c->type == p->type && c->valueIndex == p->valueIndex)
        return c;
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<c->operands.size() ; ++i)
      {
        arg0 = com_expr_ (c->operands[i], p);
        if (arg0 != NULL)
          return arg0;
      }
      if (c->type == p->type && c->funcnum == p->funcnum &&
          c->operands.size() == p->operands.size() && c->funcname == p->funcname)
      {
        for (i=0 ; i<p->operands.size() ; ++i)
          if (c->operands[i] != p->operands[i])
            break;
        if (i == p->operands.size())
          return c;
      }
      break;

    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      arg0 = com_expr_(c->operands[0], p);
      if (arg0 != NULL)
        return arg0;
      arg1 = com_expr_(c->operands[1], p);
      if (arg1 != NULL)
        return arg1;
      if (c->type == p->type &&
          ((c->operands[0] == p->operands[0] && c->operands[1] == p->operands[1]) ||
          (ops[c->type].commutative &&
            c->operands[0] == p->operands[1] && c->operands[1] == p->operands[0])))
        return c;
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::PTcheck_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::PTcheck_: Internal: bad node type";
  }

  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Differentiate_
// Purpose       : Recursive method to differentiate an ExpressionNode and children
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to differentiate an ExpressionNode and children
///
/// @param[in] arg Pointer to head of ExpressionNode tree
/// @param[in] varnum variable number with respect to which we should differentiate
/// @return expression tree for differentiated expression
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
ExpressionNode * ExpressionInternals::Differentiate_ (ExpressionNode *arg, int varnum)
{
  ExpressionNode *arg1(NULL), *arg2(NULL), *newp(NULL);
  int i;

  switch (arg->type)
  {
    case EXPR_CONSTANT:
      newp = mkcon_(0.0);
      break;

    case EXPR_VAR:
      // Is this the variable we're differentiating wrt?
      if (arg->valueIndex == varnum)
        newp = mkcon_(1.0);
      else
        newp = mkcon_(0.0);
      break;

    case EXPR_PLUS:
    case EXPR_MINUS:
      arg1 = Differentiate_(arg->operands[0], varnum);
      arg2 = Differentiate_(arg->operands[1], varnum);
          newp = mkb_(arg->type, arg1, arg2);
      break;

    case EXPR_TIMES:
      /* d(a * b) = d(a) * b + d(b) * a */
      arg1 = Differentiate_(arg->operands[0], varnum);
      arg2 = Differentiate_(arg->operands[1], varnum);

      newp = mkb_(EXPR_PLUS, mkb_(EXPR_TIMES, arg1, arg->operands[1]),
                  mkb_(EXPR_TIMES, arg->operands[0], arg2));
      break;

    case EXPR_DIVIDE:
      /* d(a / b) = (d(a) * b - d(b) * a) / b^2 */
      arg1 = Differentiate_(arg->operands[0], varnum);
      arg2 = Differentiate_(arg->operands[1], varnum);

      newp = mkb_(EXPR_DIVIDE, mkb_(EXPR_MINUS, mkb_(EXPR_TIMES, arg1,
              arg->operands[1]), mkb_(EXPR_TIMES, arg->operands[0], arg2)),
              mkb_(EXPR_POWER, arg->operands[1], mkcon_(2.0)));
      break;

    case EXPR_REMAINDER:
      /* same as input */
      newp = Differentiate_(arg->operands[0], varnum);
      break;

    case EXPR_POWER:
      // Two cases... If the power is a constant then we're cool.
      // Otherwise we have to be tricky.

      if (arg->operands[1]->type == EXPR_CONSTANT)
      {
        arg1 = Differentiate_(arg->operands[0], varnum);

        newp = mkb_(EXPR_TIMES, mkb_(EXPR_TIMES,
                mkcon_(arg->operands[1]->constant),
                mkb_(EXPR_POWER, arg->operands[0],
                mkcon_(arg->operands[1]->constant - 1))),
                arg1);
      }
      else
      {
        // This is complicated.  f(x) ^ g(x) ->
        // g(x)*f(x)**(g(x)-1)*(df/dx)+f(x)**g(x)*ln(f(x))*(dg/dx)

        arg1 = Differentiate_(arg->operands[0], varnum);
        arg2 = Differentiate_(arg->operands[1], varnum);
        newp = mkb_(EXPR_PLUS,
                  mkb_(EXPR_TIMES,
                    mkb_(EXPR_TIMES, arg->operands[1],
                      mkb_(EXPR_POWER, arg->operands[0],
                        mkb_(EXPR_MINUS, arg->operands[1], mkcon_(1.0)))) , arg1),
                  mkb_(EXPR_TIMES,
                    mkb_(EXPR_TIMES,
                      mkb_(EXPR_POWER, arg->operands[0], arg->operands[1]),
                        mkf_(EXPR_F_LOG, arg->operands[0])), arg2));
      }
      break;

    case EXPR_OR:
    case EXPR_AND:
    case EXPR_XOR:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      newp = mkcon_(0.0);
      break;

    case EXPR_FUNCTION:
      // Many cases.  For single argument functions, set arg1 to the
      // derivative of the function, and arg2 to the derivative of
      // the argument.  The "IF" function is a special   case that is
      // handled differently.

      if (arg->funcnum == EXPR_F_IF)
      {
        arg1 = Differentiate_(arg->operands[1], varnum);
        arg2 = Differentiate_(arg->operands[2], varnum);
        newp = mkf_(EXPR_F_IF, arg->operands[0]);
        newp->operands[1] = arg1;
        newp->operands[2] = arg2;
      }
      else
      {
        switch (arg->funcnum)
        {
          case EXPR_F_ABS:  /* sgn(u) */
            arg1 = mkf_(EXPR_F_SGN, arg->operands[0]);
            break;

          case EXPR_F_SDT:  // sdt(u) = 0.5*[u(new)+u(old)]*dt, so derivative of sdt(u) = 0.5*dt*u(new)'
            arg1 = mkb_(EXPR_TIMES, mkcon_(0.5), mkcon_(sim_dt_));
            break;

          case EXPR_F_DDT:
          case EXPR_F_SGN:
          case EXPR_F_RAND:
          case EXPR_F_AGAUSS:
          case EXPR_F_GAUSS:
            arg1 = mkcon_(0.0);
            break;

          case EXPR_F_INT:
          case EXPR_F_FLOOR:
          case EXPR_F_CEIL:
            arg1 = mkcon_(0.0);
            break;

          case EXPR_F_ACOS:  /* - 1 / sqrt(1 - u^2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(-1.0), mkf_(EXPR_F_SQRT,
                    mkb_(EXPR_MINUS, mkcon_(1.0),
                    mkb_(EXPR_POWER, arg->operands[0], mkcon_(2.0)))));
            break;

          case EXPR_F_ACOSH: /* 1 / sqrt(u^2 - 1) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkf_(EXPR_F_SQRT,
                    mkb_(EXPR_MINUS, mkb_(EXPR_POWER, arg->operands[0],
                    mkcon_(2.0)),
                    mkcon_(1.0))));
            break;

          case EXPR_F_ASIN:  /* 1 / sqrt(1 - u^2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkf_(EXPR_F_SQRT,
                    mkb_(EXPR_MINUS, mkcon_(1.0),
                    mkb_(EXPR_POWER, arg->operands[0], mkcon_(2.0)))));
            break;

          case EXPR_F_ASINH: /* 1 / sqrt(u^2 + 1) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkf_(EXPR_F_SQRT,
                    mkb_(EXPR_PLUS, mkb_(EXPR_POWER, arg->operands[0],
                    mkcon_(2.0)), mkcon_(1.0))));
            break;

          case EXPR_F_ATAN:  /* 1 / (1 + u^2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkb_(EXPR_PLUS,
                    mkb_(EXPR_POWER, arg->operands[0], mkcon_(2.0)),
                    mkcon_(1.0)));
            break;

          case EXPR_F_ATANH: /* 1 / (1 - u^2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkb_(EXPR_MINUS,
                    mkcon_(1.0), mkb_(EXPR_POWER,
                    arg->operands[0], mkcon_(2.0))));
            break;

          case EXPR_F_COS:   /* - sin(u) */
            arg1 = mkf_(EXPR_F_UMINUS, mkf_(EXPR_F_SIN, arg->operands[0]));
            break;

          case EXPR_F_COSH:  /* sinh(u) */
            arg1 = mkf_(EXPR_F_SINH, arg->operands[0]);
            break;

          case EXPR_F_EXP:   /* exp(u) */
            arg1 = mkf_(EXPR_F_EXP, arg->operands[0]);
            break;

          case EXPR_F_LN:    /* 1 / u */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), arg->operands[0]);
            break;

          case EXPR_F_LOG:   /* log(e) / u */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(log10(exp(1.0))), arg->operands[0]);
            break;

          case EXPR_F_SIN:   /* cos(u) */
            arg1 = mkf_(EXPR_F_COS, arg->operands[0]);
            break;

          case EXPR_F_SINH:  /* cosh(u) */
            arg1 = mkf_(EXPR_F_COSH, arg->operands[0]);
            break;

          case EXPR_F_SQRT:  /* 1 / (2 * sqrt(u)) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkb_(EXPR_TIMES,
                    mkcon_(2.0), mkf_(EXPR_F_SQRT, arg->operands[0])));
            break;

          case EXPR_F_TABLE:
            {
              int num_table_pairs=(arg->operands.size()-1)/2;
              int j=1;
              // To differentiate a table, we construct a new table with
              // first and last abscissas at the same place as the original
              // table, but all internal abscissas at the midpoints of
              // intervals in the original table.  Thus we have at most
              // one pair more than the original table did.
              // The fly in the ointment is that if the original table
              // has some pathological cases where there are multiple table
              // elements with the same abscissa, which means we have
              // fewer intervals.
              arg1 = mkf_(EXPR_F_TABLE, arg->operands[0]);
              arg1->operands.resize((num_table_pairs+1)*2+1);
              // First and last table element:
              arg1->operands[j++] = arg->operands[1];
              arg1->operands[j++] = mkcon_(0.0);

              // for each interval in the original table, we make
              // an element at the midpoint for the finite-difference
              // derivative, but only if the interval is not pathological:
              for (i=1 ; i<arg->operands.size()-2 ; i+=2)
              {
                if (arg->operands[i+2] != arg->operands[i])
                {
                  // midpoint = (x1+x2)/2
                  arg1->operands[j++] = mkb_(
                      EXPR_DIVIDE,
                      mkb_(EXPR_PLUS,
                          arg->operands[i],
                          arg->operands[i+2]),
                      mkcon_(2.0));
                  // derivative at midpoint = (y2-y1)/(x2-x1):
                  arg1->operands[j++] =
                    mkb_(EXPR_DIVIDE,
                        mkb_(EXPR_MINUS,
                            arg->operands[i+3],
                            arg->operands[i+1]),
                        mkb_(EXPR_MINUS,
                            arg->operands[i+2],
                            arg->operands[i]));
                }
              }
              arg1->operands[j++] = arg->operands[arg->operands.size()-2];
              arg1->operands[j++] = mkcon_(0.0);
              // Allow for case when some intervals were left out:
              if (j < (num_table_pairs+1)*2+1)
              {
                arg1->operands.resize(j);
              }
            }
            break;

          case EXPR_F_TAN:   /* 1 / (cos(u) ^ 2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkb_(EXPR_POWER,
                    mkf_(EXPR_F_COS, arg->operands[0]), mkcon_(2.0)));
            break;

          case EXPR_F_TANH:  /* 1 / (cosh(u) ^ 2) */
            arg1 = mkb_(EXPR_DIVIDE, mkcon_(1.0), mkb_(EXPR_POWER,
                    mkf_(EXPR_F_COSH, arg->operands[0]), mkcon_(2.0)));
            break;

          case EXPR_F_UMINUS:    /* - 1 ; like a constant (was 0 !) */
            arg1 = mkcon_(-1.0);
            break;

          case EXPR_F_URAMP:
            arg1 = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, arg->operands[0], mkcon_(0.0)));
            arg1->operands[1] = mkcon_(1.0);
            arg1->operands[2] = mkcon_(0.0);
            break;

          case EXPR_F_NOT:  /* Should never happen!! */
            arg1 = mkcon_(0.0);
            break;

          case EXPR_F_USER:
            newp = NULL;
            break;

          default:
            Report::DevelFatal()
              << "ExpressionInternals::Differentiate: Internal: bad function";
            newp = NULL;
            break;
        }

        arg2 = Differentiate_(arg->operands[0], varnum);
        newp = mkb_(EXPR_TIMES, arg1, arg2);
      }

      break;

    default:
      Report::DevelFatal()
        << "ExpressionInternals::Differentiate: Internal: bad node type";

      newp = NULL;
      break;
  }

  return (newp);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Rconvert_
// Purpose       : Recursive method to convert ExpressionNode subtree to replace
//                 a placeholder with a constant or variable
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to convert ExpressionNode subtree to replace
/// a placeholder with a constant or variable
///
/// @note This method takes no arguments, but depends instead on the caller
/// having set several member variables to tell it what to do.
/// These member variables are:
///    - Rmode_: if set to zero, convert to constant, otherwise convert to var
///    - Rstring_: Used only to check that replacement index and name are consistent
///    - Rcval_: if Rmode_ is zero, the constant value to use to replace
///    - ind_replace_: the index of the variable to replace
///    - curr_magic_: a counter used while reordering node indices
/// 
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::Rconvert_ (ExpressionNode & node)
{
  int i;

  switch (node.type)
  {
    case EXPR_CONSTANT:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      Rconvert_ ( *(node.operands[0]));
      Rconvert_ ( *(node.operands[1]));
      break;

    case EXPR_VAR:
    case EXPR_FUNCTION:
      if ((node.type == EXPR_VAR || (node.type == EXPR_FUNCTION && node.funcnum == EXPR_F_USER)) &&
            node.constant != curr_magic_)
      {
        if (node.type == EXPR_VAR && ind_replace_ == node.valueIndex)
        {
          if (((varTypes_[node.valueIndex] == EXPR_T_STRING ||
                varTypes_[node.valueIndex] == EXPR_T_VARIABLE ||
                varTypes_[node.valueIndex] == EXPR_T_FUNCTION) && node.funcname != Rstring_) ||
                (varTypes_[node.valueIndex] == EXPR_T_INSTANCE && varValues_[node.valueIndex] != Rstring_))
          {
            Report::DevelFatal()
              << "ExpressionInternals::Rconvert_: function name inconsistency";
          }
          if (Rmode_ == 0)
          {
            node.type = EXPR_CONSTANT;
            node.constant = Rcval_;
          }
          else
          {
            node.type = EXPR_VAR;
          }
        }
        else
        {
/* decrement the index of any variables with higher numbers */
          if (Rmode_ == 0 && node.valueIndex > ind_replace_)
            node.valueIndex--;
            node.constant = curr_magic_;
        }
      }
      if (node.type == EXPR_VAR || node.type == EXPR_CONSTANT)
        break;
      for (i=0 ; i<node.operands.size() ; ++i)
      {
        Rconvert_ ( *(node.operands[i]));
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::Rconvert_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::Rconvert_: Unknown node type";
      break;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::RcountDDT_
// Purpose       : Recursive method to count number of ddt function calls in the
//                 expression.
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Recursive method to count number of ddt function calls in the
/// expression.
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
void ExpressionInternals::RcountDDT_ (ExpressionNode & node)
{
  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      RcountDDT_ ( *(node.operands[0]));
      RcountDDT_ ( *(node.operands[1]));
      break;

    case EXPR_FUNCTION:
      if (node.funcnum == EXPR_F_DDT)
      {
        ++numDDT;
      }
      else
      {
        for (int i=0 ; i<node.operands.size() ; ++i)
        {
          RcountDDT_ ( *(node.operands[i]));
        }
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::RcountDDT_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::RcountDDT_: Unknown node type";
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::RgetDDT_
// Purpose       : Recursive method to get arguments of ddt function calls in the
//                 expression.
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Recursive method to get arguments of ddt function calls in the
/// expression.
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
void ExpressionInternals::RgetDDT_ (ExpressionNode & node, std::vector<double> & vals)
{
  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      RgetDDT_ ( *(node.operands[0]), vals);
      RgetDDT_ ( *(node.operands[1]), vals);
      break;

    case EXPR_FUNCTION:
      if (node.funcnum == EXPR_F_DDT)
      {
        vals[numDDT] = node.operands[0]->eval_value;
        ++numDDT;
      }
      else
      {
        for (int i=0 ; i<node.operands.size() ; ++i)
        {
          RgetDDT_ ( *(node.operands[i]), vals);
        }
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::RgetDDT: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::RgetDDT: Unknown node type";
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::RsetDDT_
// Purpose       : Recursive method to set the evaluated values of ddt function
//                 calls in the expression.
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/05
//-----------------------------------------------------------------------------
///
/// Recursive method to set the evaluated values of ddt function
/// calls in the expression.
///
/// @author Dave Shirley, PSSI
/// @date 12/16/05
///
void ExpressionInternals::RsetDDT_ (ExpressionNode & node, std::vector<double> & vals)
{
  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      RsetDDT_ ( *(node.operands[0]), vals);
      RsetDDT_ ( *(node.operands[1]), vals);
      break;

    case EXPR_FUNCTION:
      if (node.funcnum == EXPR_F_DDT)
      {
        node.eval_value = vals[numDDT];
        node.eval_num = curr_num_;
        ++numDDT;
      }
      else
      {
        for (int i=0 ; i<node.operands.size() ; ++i)
        {
          RsetDDT_ ( *(node.operands[i]), vals);
        }
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::RsetDDT_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::RsetDDT_: Unknown node type";
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::convert_to_constant_
// Purpose       : Resolve a placeholder or variable into a constant
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Resolve a placeholder or variable into a constant
///
/// @param[in] i variable index
/// @param[in] c_value constant to replace
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::convert_to_constant_ (int i, double c_value)
{
  values_changed_ = true;
  Rcval_ = c_value;
  Rmode_ = 0;
  Rstring_ = varValues_[i];
  ind_replace_ = i;
  curr_magic_ += 1;
  Rconvert_(*tree_);
  numVars_--;
  varTypes_.erase(varTypes_.begin()+i);
  varValues_.erase(varValues_.begin()+i);
  leadDesignator_.erase(leadDesignator_.begin()+i);
  var_vals_.erase(var_vals_.begin()+i);

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::convert_to_variable_
// Purpose       : Resolve a placeholder into a variable
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Resolve a placeholder into a variable
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::convert_to_variable_ (int i)
{
  values_changed_ = true;
  Rmode_ = 1;
  Rstring_ = varValues_[i];
  ind_replace_ = i;
  Rconvert_(*tree_);
  varTypes_[i] = EXPR_T_VARIABLE;
  var_vals_[i] = 0;

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::copy_exprNode_
// Purpose       : Recursive method to copy an expression node and subnodes
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/05/04
//-----------------------------------------------------------------------------
///
/// Recursive method to copy an expression node and subnodes
///
/// @param[in] n root of expression tree to copy
/// @return pointer to new expression tree
///
/// This is a deep copy operation.
///
/// @note Seems ripe to refactor this C-like code into proper C++.
///
/// @author Dave Shirley, PSSI
/// @date 11/05/04
///
ExpressionNode * ExpressionInternals::copy_exprNode_ (ExpressionNode *n)
{
  int i;
  ExpressionNode *e;

  e = newExpressionNode_();
  e->type = n->type;
  e->operands.resize(n->operands.size());
  for (i=0 ; i<n->operands.size() ; ++i)
    e->operands[i] = copy_exprNode_ (n->operands[i]);
  e->constant = n->constant;
  e->state = n->state;
  e->valueIndex = n->valueIndex;
  e->funcname = n->funcname;
  e->eval_value = 0;
  e->eval_num = 0;
  e->funcnum = n->funcnum;
  e->fptr = n->fptr;

  return e;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::simplify_
// Purpose       : Recursive method to condense any constant subexpressions
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to condense any constant subexpressions
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::simplify_ (ExpressionNode & node)
{
  int i;

  values_changed_ = true;
  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      simplify_ ( *(node.operands[0]));
      simplify_ ( *(node.operands[1]));
      if (node.operands[0]->type == EXPR_CONSTANT && node.operands[1]->type == EXPR_CONSTANT)
      {
        node.type = EXPR_CONSTANT;
        node.constant = (node.fptr.binary)(node.operands[0]->constant, node.operands[1]->constant);
        node.operands.clear();
      }
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<node.operands.size() ; ++i)
        simplify_ ( *(node.operands[i]));
      if (node.funcnum == EXPR_F_USER)
        break;
//DNS: could implement conversion to constant for constant valued functions
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::simplify_: placeholder found";

    default:
      Report::DevelFatal()
        <<  "ExpressionInternals::simplify_: Unknown node type";
      break;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::get_breaks_
// Purpose       : Recursive method to generate breakpoint expressions.  These
//                 are expressions which will generate the next breakpoint after
//                 the current simulation time, like table segment boundaries.
//
// Special Notes : Note added by Tom Russo on 11 Feb 2014, while investigating
//                 bug 1827 on Charleston:
//
//                 The "Purpose" above is a misstatement of what this method
//                 actually does.
//
//                 The expressions generated here do *not* generate the next
//                 breakpoint at all, they are effectively used in a Newton's
//                 method to approximate the solution of the problem:
//                 ThisExpression(time)-(value of expression when break occurs)=0
//
//                 This is *ripe* for refactoring, especially when the current
//                 expression should have a closed-form solution for
//                 breakpoints (e.g. pulse or PWL expressions, which are
//                 currently reduced to a general TABLE function at the
//                 input parser level before we even see it).
//
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to generate breakpoint expressions.  These
/// are expressions which will generate the next breakpoint after
/// the current simulation time, like table segment boundaries.
///
/// @note The "Purpose" above is a misstatement of what this method
/// actually does. 
/// The expressions generated here do *not* generate the next
/// breakpoint at all, they are effectively used in a Newton's
/// method to approximate the solution of the problem:
/// ThisExpression(time)-(value of expression when break occurs)=0
/// This is *ripe* for refactoring, especially when the current
/// expression should have a closed-form solution for
/// breakpoints (e.g. pulse or PWL expressions, which are
/// currently reduced to a general TABLE function at the
/// input parser level before we even see it).
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::get_breaks_ (ExpressionNode & node)
{
  int i;
  ExpressionNode *p, *t, *t1, *t2, *b;

  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_POWER:
      get_breaks_ ( *(node.operands[0]));
      get_breaks_ ( *(node.operands[1]));
      break;

      // TVR 11 Feb 2014:  This case makes it most obvious what this
      // method is really doing.  The expression pushed back into "breaks_"
      // when the current expression is "F(t)>=A" is in fact
      //  -(F(t)-A)/F'(t) --- exactly the expression that represents the
      // next iteration of a Newton's root finding method for finding the
      // value of "t" at which F(t)=A.
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      get_breaks_ ( *(node.operands[0]));
      get_breaks_ ( *(node.operands[1]));
      if ((dependent_( *(node.operands[0]), time_index) ||
            dependent_( *(node.operands[1]), time_index)) &&
          (!dependent_other_( *(node.operands[0]), time_index) &&
            !dependent_other_( *(node.operands[1]), time_index)))
      {
        t = mkb_(EXPR_MINUS, node.operands[0], node.operands[1]);
        p = Differentiate_ (t, time_index);
        b = mkf_(EXPR_F_UMINUS, mkb_(EXPR_DIVIDE, t, p));
        breaks_.push_back(b);
      }
      break;

    case EXPR_REMAINDER:
      get_breaks_ ( *(node.operands[0]));
      get_breaks_ ( *(node.operands[1]));
      if (dependent_( *(node.operands[0]), time_index) &&
          !dependent_other_( *(node.operands[0]), time_index))
      {
// in a%b, p is da/dt
        p = Differentiate_ (node.operands[0], time_index);
// t is the target value of a%b in the future.  This can be -abs(b), 0, or abs(b)
        t = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, node.operands[0], mkcon_(0)));
        t->operands[1] = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, p, mkcon_(0)));
        t->operands[1]->operands[1] = mkf_(EXPR_F_ABS, node.operands[1]);
        t->operands[1]->operands[2] = mkcon_(0);
        t->operands[2] = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, p, mkcon_(0)));
        t->operands[2]->operands[1] = mkcon_(0);
        t->operands[2]->operands[2] = mkf_(EXPR_F_UMINUS, mkf_(EXPR_F_ABS,
                                        node.operands[1]));
// finally, b is the predicted time of the crossing, (t-a%b)/(da/dt)
        b = mkb_(EXPR_DIVIDE, mkb_(EXPR_MINUS, t, &node), p);
        breaks_.push_back(b);
      }
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<node.operands.size() ; ++i)
        get_breaks_ ( *(node.operands[i]));
      if (node.funcnum == EXPR_F_USER)
        return;
      if (node.funcnum == EXPR_F_TABLE)
      {
        if (dependent_( *(node.operands[0]), time_index) &&
            !dependent_other_( *(node.operands[0]), time_index))
        {
          p = Differentiate_ (node.operands[0], time_index);
          b = mkf_(EXPR_F_IF, mkb_(EXPR_GREAT, p, mkcon_(0)));
          t1 = mkf_(EXPR_F_F_TABLE, node.operands[0]);
          t2 = mkf_(EXPR_F_R_TABLE, node.operands[0]);
          b->operands[1] = mkb_(EXPR_DIVIDE, t1, p);
          b->operands[2] = mkb_(EXPR_DIVIDE, t2, p);
          t1->operands[0] = node.operands[0];
          t2->operands[0] = node.operands[0];
          for (i=1 ; i<node.operands.size() ; i+=2)
          {
            t1->operands.push_back(node.operands[i]);
            t2->operands.push_back(node.operands[i]);
          }
          breaks_.push_back(b);
        }
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::get_breaks_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::get_breaks_: Unknown node type";
      break;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::dependent_
// Purpose       : Recursive method to determine if a subexpression is dependent
//                 on a particular variable, of valueIndex = ind
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to determine if a subexpression is dependent
/// on a particular variable, of valueIndex = ind
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
bool ExpressionInternals::dependent_ (ExpressionNode & node, int ind)
{
  int i;

  switch (node.type)
  {
    case EXPR_CONSTANT:
      break;

    case EXPR_VAR:
      if (node.valueIndex == ind)
        return true;
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_POWER:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_REMAINDER:
      if (dependent_ ( *(node.operands[0]), ind))
        return true;
      if (dependent_ ( *(node.operands[1]), ind))
        return true;
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<node.operands.size() ; ++i)
      {
        if (dependent_ ( *(node.operands[i]), ind))
          return true;
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        <<  "ExpressionInternals::dependent_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::dependent_: Unknown node type";
      break;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::dependent_other_
// Purpose       : Recursive method to determine if a subexpression is
//                 dependent_ on any variable other than
//                 a particular variable, of valueIndex = ind
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/11/05
//-----------------------------------------------------------------------------
///
/// Recursive method to determine if a subexpression is
/// dependent_ on any variable other than
/// a particular variable, of valueIndex = ind
///
/// @author Dave Shirley, PSSI
/// @date 03/11/05
///
bool ExpressionInternals::dependent_other_ (ExpressionNode & node, int ind)
{
  int i;

  switch (node.type)
  {
    case EXPR_CONSTANT:
      break;

    case EXPR_VAR:
      if (node.valueIndex != ind)
        return true;
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_POWER:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_REMAINDER:
      if (dependent_other_ ( *(node.operands[0]), ind))
        return true;
      if (dependent_other_ ( *(node.operands[1]), ind))
        return true;
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<node.operands.size() ; ++i)
      {
        if (dependent_other_ ( *(node.operands[i]), ind))
          return true;
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::dependent_other_: placeholder found";

    default:
      Report::DevelFatal()
        <<  "ExpressionInternals::dependent_other_: Unknown node type";
      break;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::arithmatic_
// Purpose       : Recursive method to determine if a subexpression is arithmatic_
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to determine if a subexpression is arithmatic_
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
bool ExpressionInternals::arithmatic_ (ExpressionNode & node)
{
  int i;

  switch (node.type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;

    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_POWER:
    case EXPR_REMAINDER:
      if (!arithmatic_( *(node.operands[0])))
        return false;
      if (!arithmatic_( *(node.operands[1])))
        return false;
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<node.operands.size() ; ++i)
      {
        if (!arithmatic_( *(node.operands[i])))
          return false;
      }
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
      return false;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        <<  "ExpressionInternals::arithmatic_: placeholder found";

    default:
      Report::DevelFatal()
        <<   "ExpressionInternals::arithmatic_: Unknown node type";
      break;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Rmap_
// Purpose       : Recursive method to fix indexing in a subtree for order_names
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to fix indexing in a subtree for order_names
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::Rmap_ (ExpressionNode & node, int mode,
                                       std::vector<int> & nmap)
{
  int i;

  switch (node.type)
  {
    case EXPR_CONSTANT:
      break;
    case EXPR_VAR:
      if (mode == 0)
      {
        if (node.valueIndex >=0)
          node.valueIndex = nmap[node.valueIndex]-nmap.size();
      }
      else
      {
        if (node.valueIndex < 0)
          node.valueIndex += nmap.size();
      }
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      Rmap_ ( *(node.operands[0]), mode, nmap);
      Rmap_ ( *(node.operands[1]), mode, nmap);
      break;

    case EXPR_FUNCTION:
      if (node.funcnum == EXPR_F_USER)
      {
        if (mode == 0)
        {
          if (node.valueIndex >=0)
            node.valueIndex = nmap[node.valueIndex]-nmap.size();
        }
        else
        {
          if (node.valueIndex < 0)
            node.valueIndex += nmap.size();
        }
      }
      for (i=0 ; i<node.operands.size() ; ++i)
        Rmap_ ( *(node.operands[i]), mode, nmap);
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        <<  "ExpressionInternals::Rmap_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::Rmap_: Unknown node type";
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPRaddDummyString_
// Purpose       : Add a dummy input quantity entry to a function definition
// Special Notes : Used when a function definition has fewer inputs than the
//                 prototype (e.g.: f(x,y) = y*2)
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Add a dummy input quantity entry to a function definition
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
int ExpressionInternals::EXPRaddDummyString_ (std::string & dummy)
{
  varTypes_.resize(numVars_+1);
  varValues_.resize(numVars_+1);
  leadDesignator_.resize(numVars_+1);
  var_vals_.resize(numVars_+1);
  if (differentiated_)
  {
    derivs_.resize(numVars_+1);
  }

  varTypes_[numVars_] = EXPR_T_STRING;
  varValues_[numVars_] = dummy;
  leadDesignator_[numVars_] = ' ';
  var_vals_[numVars_] = 0;
  if (differentiated_)
  {
    derivs_[numVars_] = NULL;
  }
  ++num_string_;

  return numVars_++;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::addNode_
// Purpose       : Recursive method to add a node during user function resolution
// Special Notes : very complex operation, recursive with Nreplace_
//                 This method *assumes* that the names in the symbol table
//                 for funcExpr have been ordered (by order_names) so that
//                 function arguments are the low numbered ones in varTypes.
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to add a node during user function resolution
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void
ExpressionInternals::addNode_(
  ExpressionNode *                      n,
  int                                   ind,
  const ExpressionNode *                f,
  const ExpressionInternals &           func_expr,
  int                                   na_func,
  std::vector<ExpressionNode *>         operands)
{

// This is called from Nreplace_ when there is a (prospective) node to
// add to the expression tree.  The corresponding node on the
// function definition is first examined to see what action to take.
// If it is a EXPR_VAR that points to a EXPR_T_STRING then we just
// point to the proper node in the expression.  Otherwise, we add a
// node and go back to Nreplace_.

  ExpressionNode *arg;
  int i, arg_num;

  arg = f->operands[ind];
  if (arg->type == EXPR_VAR && func_expr.varTypes_[arg->valueIndex] == EXPR_T_STRING)
  {
    arg_num = 0;

    // Figure out whether the index of the argument we're working with is
    // one of the function arguments of func_expr.  If so, make "n"'s operand
    // the same as the function's argument to which we correspond.
    for (i=0 ; i<arg->valueIndex ; ++i)
      if (func_expr.varTypes_[i] == EXPR_T_STRING)
        ++arg_num;
    if (arg_num < na_func)
    {
      n->operands[ind] = operands[arg_num];
      return;
    }
  }

  // If we get here, we're either not a variable or we are but we're not
  // a function argument.  In that case, create a new node for our
  // nodes operand, and recurse to continue the replacement.
  n->operands[ind] = newExpressionNode_();
  Nreplace_ (n->operands[ind], f->operands[ind], func_expr, na_func, operands);

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Nreplace_
// Purpose       : Recursive method to do work of resolving a user defined function
// Special Notes : very complex operation, recursive with addNode_
//                 What it appears this method is doing is to copy the tree
//                 rooted at "f" into "n," fixing up arguments and symbols
//                 appropriately.  Since the copying depends at every node
//                 on what type of node is being copied, and some of those
//                 nodes might be complex, this has to be recursive.
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to do work of resolving a user defined function
///
/// @param n   Expression node that should be replaced
/// @param f   Expression tree for function that should replace this node
/// @param func_expr   Full expression object for replacement function
/// @param na_func  Number of function arguments
/// @param operands  vector of expression nodes that represent n's operands.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void
ExpressionInternals::Nreplace_(
  ExpressionNode *                      n,
  const ExpressionNode *                f,
  const ExpressionInternals &           func_expr,
  int                                   na_func,
  std::vector<ExpressionNode *>         operands)
{

// This is the (recursive) routine that traverses the function that is
// being inserted and does the work of replicating the function in the
// expression tree.  The function addNode_ does the work of deciding
// which nodes need to be added, and merging in addional vars and
// funcs into the IFparseTree structure of the expression as needed.
// Although it may be more efficient to do this one time up front, it
// is simpler to just do it here.

  int i, arg_num;
 
  // Copy the function's tree head to our current node.
  *n=*f;
  switch (n->type)
  {
    case EXPR_CONSTANT:
      break;


      // If the current node is a binary operator, create nodes for its
      // operands. If the operands are more complex than just a function
      // argument, we recurse here.an
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      addNode_ (n, 0, f, func_expr, na_func, operands);
      addNode_ (n, 1, f, func_expr, na_func, operands);
      break;

      // If we're a variable, first determine if we're a function argument
      // or not.  (this depends on the func_expr's symbol table having been
      // ordered by the function order_names, in which case arguments are
      // always the first items in the table.)  If it is NOT an argument, then
      // we need to make sure the node's valueIndex points to the correct
      // element of the symbol table.
    case EXPR_VAR:
      arg_num = 0;
      for (i=0 ; i<f->valueIndex ; ++i)
        if (func_expr.varTypes_[i] == EXPR_T_STRING)
          ++arg_num;

      // Are we a function argument?
      if (arg_num >= na_func)
      {

        // NO, we're not a function argument.  Are we a variable that already
        // exists in our main expression's symbol table (i.e. the one into
        // which we're replacing the function)?
        i = find_num_(n->funcname);
        if (i >= 0)
        {
          // Found it, just use the index
          n->valueIndex = i;
        }
        else
        {
          // Did NOT find it.  Add the new variable to our existing symbol
          // table.
          varValues_.resize (numVars_ + 1);
          varTypes_.resize (numVars_ + 1);
          leadDesignator_.resize (numVars_ + 1);
          var_vals_.resize (numVars_ + 1);

          varValues_[numVars_] = n->funcname;
          if (isNameSpecial_(n->funcname))
            varTypes_[numVars_] = EXPR_T_SPECIAL;
          else
            varTypes_[numVars_] = EXPR_T_STRING;
          leadDesignator_[numVars_] = ' ';
          var_vals_[numVars_] = 0;
          n->valueIndex = numVars_;

          ++numVars_;
        }
      }
      else
      {
        // It is a function argument.  Copy the operand.
        *n = *(operands[arg_num]);
      }
      break;
    case EXPR_FUNCTION:
      if (n->funcnum == EXPR_F_USER)
      {
// This is where user defined funcs are introduced into the expression
// and the indexing is set for the expression rather than the function
// definition.
// We check whether we already know this function's name in the expression's
// symbol table.  If so, use its index.  Otherwise, add it to our symbol
// table, then use the new index.  Then add nodes for each function argument
// (of the user defineaed function).
        for (i=0 ; i<numVars_ ; ++i)
        {
          if (varTypes_[i] == func_expr.varTypes_[f->valueIndex] &&
              varValues_[i] == func_expr.varValues_[f->valueIndex])
          {
            break;
          }
        }
        
        if (i == numVars_)
        {
          ++numVars_;
          varTypes_.resize(numVars_);
          varValues_.resize(numVars_);
          leadDesignator_.resize(numVars_);
          var_vals_.resize(numVars_);
        }
        varTypes_[i] = func_expr.varTypes_[f->valueIndex];
        varValues_[i] = func_expr.varValues_[f->valueIndex];
        leadDesignator_[i] = func_expr.leadDesignator_[f->valueIndex];
        var_vals_[i] = 0;
        n->valueIndex = i;
      }
      for (i=0 ; i<n->operands.size() ; ++i)
      {
        addNode_ (n, i, f, func_expr, na_func, operands);
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::Nreplace_: placeholder found";

    default:
      Report::UserWarning() << "Unknown node type " << n->type;
      break;    // If an unknown type is encountered just return
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Freplace_
// Purpose       : Recursive method to find instances of a user defined function
//                 that needs to be replaced
//                 func_name is name of function to be replaced
//                 func_expr is the expression internals object with which to
//                  replace it
//                 na_func is the number of arguments the function takes
//
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to find instances of a user defined function
/// that needs to be replaced
///
/// @param n ExpressionNode tree on which to operate
/// @param [in] func_name name of function to be replaced
/// @param[in] func_expr  expression internals object with which to replace it
/// @param[in] na_func  number of arguments the function takes
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
int ExpressionInternals::Freplace_ (ExpressionNode *n,
                                          std::string const & func_name,
                                          ExpressionInternals & func_expr,
                                          int na_func)
{

// This recursive routine traverses the expression tree looking for
// instances of the function to be replaced.  When an instance is
// found, the routine Nreplace_ will do the actual work of the
// replacement.  When Nreplace_ is called, it passes a list of
// pointers to the arguments of the function.

  int i;
  int Ferrno=0; // overall return code for the function
  int errno1=0, errno2=0; // return codes for recursive calls to Freplace_ for binary ops

  switch (n->type)
  {
    // Nothing to do for constants and vars
    case EXPR_CONSTANT:
    case EXPR_VAR:
      break;


    // for binary ops, just recursively call Freplace_ for each argument
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      errno1=Freplace_ (n->operands[0], func_name, func_expr, na_func);
      errno2=Freplace_ (n->operands[1], func_name, func_expr, na_func);
      // return Ferrno=-1 if either recursive function call to Freplace_  returns -1
      if ( (errno1 < 0) || (errno2 < 0) )
      {
        Ferrno=-1;   
      }
      break;

    case EXPR_FUNCTION:
      // if we've found a node that is an instance of the function we're
      // replacing, do the real work.
      if (n->funcnum == EXPR_F_USER)
      {
        if (func_name == n->funcname)
        {
          // error out if we have the wrong number of operands for the
          // user-defined function
          if (n->operands.size() != na_func)
          {
            Ferrno = -1;
            break;
          }
          
          // This replaces the node "n" with the function's tree.
          Nreplace_ (n, func_expr.get_tree(), func_expr, na_func, n->operands);
        }
      }
      // Now we do an freplace on all the arguments, too.
      for (i=0 ; i<n->operands.size() && Ferrno==0; ++i)
        Ferrno=Freplace_ (n->operands[i], func_name, func_expr, na_func);
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::Freplace_: placeholder found";

    default:
      Report::UserWarning() << "Unknown node type " << n->type;
      break;    // If an unknown type is encountered just return
  }

  return Ferrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::Vreplace_
// Purpose       : Recursive method to find instances of a variable
//                 that needs to be replaced with a previously parsed expression
//                 varName is name of variable/string to be replaced
//                 subexpr is the expression internals object with which to
//                  replace it
//
// Special Notes : This function is needed for handling the case where a
//                 subcircuit parameter is passed an expression involving
//                 global params.
//                 It is very similar to, but slightly simpler than, Freplace_
// Scope         :
// Creator       : Tom Russo
// Creation Date : 8/10/2010
//-----------------------------------------------------------------------------
///
/// Recursive method to find instances of a variable
/// that needs to be replaced with a previously parsed expression
///
/// @param n Root of ExpressionNode tree on which to operate
/// @param[in] varName name of variable/string to be replaced
/// @param[in] subexpr is the expression internals object with which to replace it
///
/// @author Tom Russo
/// @date 8/10/2010
///
int ExpressionInternals::Vreplace_(
  ExpressionNode *              n,
  std::string const &           varName,
  const ExpressionInternals &   subexpr)
{

// This recursive routine traverses the expression tree looking for
// instances of the variable to be replaced.  When an instance is
// found, the routine Nreplace_ will do the actual work of the
// replacement.  When Nreplace_ is called, it passes a list of
// pointers to the arguments of the function, which in this case is empty.
// We only need the operands vector because we're reusing the Nreplace_
// method, which expects it.  But that method only uses operands if its
// fourth argument is non-zero.

  std::vector<ExpressionNode *> operands;
  operands.clear();

  int i;
  int Ferrno=0;   // overall return code for the function
  int errno1=0, errno2=0; // return codes for recursive calls to Vreplace_ for binary ops

  switch (n->type)
  {
    // Nothing to do for constants.
    case EXPR_CONSTANT:
      break;
    case EXPR_VAR:
      // if we have found a node that is an instance of the variable we're
      // replacing, do the real work.
      if (varName == varValues_[n->valueIndex])
      {
        // pass dummy operands and 0 arguments.
        Nreplace_(n, subexpr.get_tree(), subexpr, 0, operands);
      }
      break;


      // for binary ops, just recurse for each argument
    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      errno1=Vreplace_ (n->operands[0], varName, subexpr);
      errno2=Vreplace_ (n->operands[1], varName, subexpr);
      // return Ferrno=-1 if either recursive function call to Vreplace_  returns -1
      if ( (errno1 < 0) || (errno2 < 0) )
      {
        Ferrno=-1;   
      }
      break;

    case EXPR_FUNCTION:
      // do a Vreplace on all the arguments of any function we encounter.
      for (i=0 ; i<n->operands.size() && Ferrno==0; ++i)
        Ferrno=Vreplace_ (n->operands[i], varName, subexpr);
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::Freplace_: placeholder found";

    default:
      Report::UserWarning() << "Unknown node type " << n->type;
      break;    // If an unknown type is encountered just return
  }

  return Ferrno;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::RemoveFentry_
// Purpose       : Recursive method to fix indexing associated with resolution
//                 of a user defined function.  Once we've replaced a function
//                 using Freplace, that name is no longer needed in the
//                 symbol table.
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to fix indexing associated with resolution
/// of a user defined function.  Once we've replaced a function
/// using Freplace, that name is no longer needed in the
/// symbol table.
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::RemoveFentry_ (ExpressionNode *n, int new_ind, int old_ind)
{

// This recursive routine traverses the expression tree and modifies
// index numbers into the varValues_ and varTypes_ arrays to reflect
// removal of the user defined function from the expression.

  int i;

  switch (n->type)
  {
    case EXPR_CONSTANT:
      break;
    case EXPR_VAR:
      if (n->valueIndex == old_ind)
        n->valueIndex = new_ind;
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      RemoveFentry_ (n->operands[0], new_ind, old_ind);
      RemoveFentry_ (n->operands[1], new_ind, old_ind);
      break;

    case EXPR_FUNCTION:
      if (n->funcnum == EXPR_F_USER)
      {
        if (n->valueIndex == old_ind)
          n->valueIndex = new_ind;
      }
      for (i=0 ; i<n->operands.size() ; ++i)
        RemoveFentry_ (n->operands[i], new_ind, old_ind);
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::RemoveFentry_: placeholder found";

    default:
      Report::UserWarning() << "Unknown node type " << n->type;
      break;    // If an unknown type is encountered just return
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPReval_
// Purpose       : Recursive method to perform expression evaluation
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/24/04
//-----------------------------------------------------------------------------
///
/// Recursive method to perform expression evaluation
///
/// @author Dave Shirley, PSSI
/// @date 09/24/04
///
void ExpressionInternals::EXPReval_ (ExpressionNode & node, double & res, std::vector<double> &vals)
{
  double r1, r2, x1, x2, y1, y2, frac, condition, raw_value;// result;
  int i, t_len;
  bool interpolate;

  if (node.eval_num == curr_num_)
  {
    res = node.eval_value;
    return;
  }

  switch (node.type)
  {
    case EXPR_CONSTANT:
      res = node.constant;
      break;

    case EXPR_VAR:
      if (node.valueIndex == -1)
      {
      Report::DevelFatal()
        << "ExpressionInternals::EXPReval_: Internal Error: node.valueIndex == -1";
      }
      else
      {
        res = vals[node.valueIndex];
      }
      break;

    case EXPR_FUNCTION:
    // Multi-argument and functions that require data from the
    // simulator are evaluated here
      if (node.funcnum == EXPR_F_IF)
      {
        EXPReval_( *(node.operands[0]), condition, vals);
        if (condition > 0)
          EXPReval_( *(node.operands[1]), res, vals);
        else
          EXPReval_( *(node.operands[2]), res, vals);
      }
      else if (node.funcnum == EXPR_F_SDT)
      {
        double & sdtOld  = node.state[0];
        double & timeOld = node.state[1];
        double & fOld    = node.state[2];
        double & timeNew = node.state[3];
        double & fNew    = node.state[4];
        double & sdtNew  = node.state[5];

        bool rotate = (accepted_time >= timeNew);
        if (rotate)
        {
            sdtOld = sdtNew;
            fOld = fNew;
            timeOld = timeNew;
        }

        // evaluates the value of the function inside the integral (raw_value)
        EXPReval_( *(node.operands[0]), raw_value, vals);
        fNew = raw_value;
        timeNew = sim_time_;

        double dt = timeNew - timeOld;
        double faverage = (fNew + fOld)*0.5;
        double dS = faverage*dt;

        sdtNew = sdtOld + dS;
        res = sdtNew;
      }
      else if (node.funcnum == EXPR_F_AGAUSS)
      {
        if (node.eval_value == 0)
        {

          double mu, abs_var, sigma,sigma_unit;
          EXPReval_( *(node.operands[0]), mu, vals);
          EXPReval_( *(node.operands[1]), abs_var, vals);
          EXPReval_( *(node.operands[2]), sigma_unit, vals);

          sigma = abs_var/sigma_unit;

          if (enableRandomExpression)
          {
            if (theRandomNumberGenerator == 0)
            {
              theRandomNumberGenerator = new Xyce::Util::RandomNumbers(0);
            }
            do
            {
              res = theRandomNumberGenerator->gaussianRandom(mu,sigma);
            } while (abs(res-mu)>abs_var);
          }
          else
          {
            res = mu;
          }
        }
        else
          res = static_cast<double>(node.eval_value);

      }
      else if (node.funcnum == EXPR_F_GAUSS)
      {
        if (node.eval_value == 0)
        {

          double mu, abs_var,rel_var, sigma,sigma_unit;
          EXPReval_( *(node.operands[0]), mu, vals);
          EXPReval_( *(node.operands[1]), rel_var, vals);
          EXPReval_( *(node.operands[2]), sigma_unit, vals);
          abs_var=rel_var*mu;
          sigma = abs_var/sigma_unit;
          if (enableRandomExpression)
          {
            if (theRandomNumberGenerator == 0)
            {
              theRandomNumberGenerator = new Xyce::Util::RandomNumbers(0);
            }

            do
            {
              res = theRandomNumberGenerator->gaussianRandom(mu,sigma);
            } while (abs(res-mu)>abs_var);
          }
          else
          {
            res = mu;
          }
        }
        else
          res = static_cast<double>(node.eval_value);

      }
      else if (node.funcnum == EXPR_F_RAND)
      {
        if (node.eval_value == 0)
        {
          if (theRandomNumberGenerator == 0)
          {
            theRandomNumberGenerator = new Xyce::Util::RandomNumbers(0);
          }
          res = theRandomNumberGenerator->uniformRandom();
        }
        else
          res = static_cast<double>(node.eval_value);
      }
      else if (node.funcnum == EXPR_F_DDT)
      {
        // The definition of the state variables are:
        //  0    Current derivative
        //  1    Time at accepted derivative
        //  2    Value of function at accepted derivative
        //  3    Current time
        //  4    Current value of funtion
        double & ddtOld  = node.state[0];
        double & timeOld = node.state[1];
        double & fOld    = node.state[2];
        double & timeNew = node.state[3];
        double & fNew    = node.state[4];
        double & ddtNew  = node.state[5];

        bool rotate = (accepted_time >= timeNew);
        if (rotate)
        {
          ddtOld = ddtNew;
          fOld = fNew;
          timeOld = timeNew;
        }

        EXPReval_( *(node.operands[0]), raw_value, vals);
        timeNew = sim_time_;
        fNew = raw_value;
        double df = fNew-fOld;
        double dt = timeNew-timeOld;

        if (dt != 0.0)
        {
          ddtNew = df/dt;
        }
        else
        {
          ddtNew = ddtOld;
        }

        res = ddtNew;
      }
      else if (node.funcnum == EXPR_F_OP)
      {
        res = (*node.fptr.op)(node.op);
      }
      else
      {
        EXPReval_( *(node.operands[0]), r1, vals);
        if (node.funcnum == EXPR_F_TABLE)
        {
          t_len = node.operands.size()/2;
          i = node.valueIndex;
          EXPReval_( *(node.operands[2*i+1]), x1, vals);
          while (i >= 0 && i< t_len)
          {
            interpolate = false;
            if (x1 < r1)
            {
              if (i < t_len-1)
              {
                EXPReval_( *(node.operands[2*i+3]), x2, vals);
                if (x2 < x1)
                {
                  EXPR_ERROR (EXPR_TABLE_FATAL);
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  res = y1;
                  break;
                }
                if (x2 > r1)
                {
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  EXPReval_( *(node.operands[2*i+4]), y2, vals);
                  interpolate = true;
                }
                else
                {
                  ++i;
                  x1 = x2;
                }
              }
              else
              {
                EXPReval_( *(node.operands[2*i+2]), y1, vals);
                res = y1;
                break;
              }
            }
            else if (x1 > r1)
            {
              if (i > 0)
              {
                x2 = x1;
                EXPReval_( *(node.operands[2*i-1]), x1, vals);
                if (x2 < x1)
                {
                  EXPR_ERROR (EXPR_TABLE_FATAL);
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  res = y1;
                  break;
                }
                if (x1 < r1)
                {
                  EXPReval_( *(node.operands[2*i]), y1, vals);
                  EXPReval_( *(node.operands[2*i+2]), y2, vals);
                  interpolate = true;
                }
                else
                {
                  i--;
                }
              }
              else
              {
                EXPReval_( *(node.operands[2]), y1, vals);
                res = y1;
                break;
              }
            }
            else if (x1 == r1)
            {
              EXPReval_( *(node.operands[2*i+2]), y1, vals);
              res = y1;
              break;
            }
            else  // The only way to get here is if r1 is a NaN
            {
              EXPR_ERROR ( EXPR_TABLE_FATAL);
              res=EXPR_HUGE;
              break;
            }
            if (interpolate)
            {
              if (x1 == x2)
              {
                res = 0.5 * (y1+y2);
              }
              else
              {
                frac = (r1-x1)/(x2-x1);
                res = frac*y2 + (1-frac)*y1;
              }
              break;
            }
          }
          node.valueIndex = i;
        }
        else if (node.funcnum == EXPR_F_F_TABLE || node.funcnum == EXPR_F_R_TABLE)
        {
          i = node.valueIndex;
          if (i == 0)
            i = 1;
          EXPReval_( *(node.operands[i]), x1, vals);
          while (i > 0 && i<node.operands.size())
          {
            interpolate = false;
            if (x1 < r1)
            {
              if (i < node.operands.size()-1)
              {
                EXPReval_( *(node.operands[i+1]), x2, vals);
                if (x2 > r1)
                {
                  interpolate = true;
                }
                else
                {
                  x1 = x2;
                  ++i;
                }
              }
              else
              {
                if (node.funcnum == EXPR_F_F_TABLE)
                  res = EXPR_HUGE;
                else
                  res = x1-r1;
                break;
              }
            }
            else if (x1 > r1)
            {
              if (i > 1)
              {
                x2 = x1;
                EXPReval_( *(node.operands[i-1]), x1, vals);
                if (x1 < r1)
                {
                  interpolate = true;
                }
                else
                {
                  i--;
                }
              }
              else
              {
                if (node.funcnum == EXPR_F_R_TABLE)
                  res = -EXPR_HUGE;
                else
                  res = x1-r1;
                break;
              }
            }
            else
            {
              if (node.funcnum == EXPR_F_R_TABLE)
              {
                if (i == 1)
                {
                  res = -EXPR_HUGE;
                }
                else
                {
                  EXPReval_( *(node.operands[i-1]), x2, vals);
                  res = x2-x1;
                }
              }
              else
              {
                if (i == node.operands.size()-1)
                {
                  res = EXPR_HUGE;
                }
                else
                {
                  EXPReval_( *(node.operands[i+1]), x2, vals);
                  res = x2-x1;
                }
              }
              break;
            }
            if (interpolate)
            {
              if (node.funcnum == EXPR_F_R_TABLE)
                res = x1-r1;
              else
                res = x2-r1;
              break;
            }
          }
          node.valueIndex = i;
        }
        else if(node.funcnum == EXPR_F_SCHEDULE)
        {
          // a schedule is handled like a table except that there is
          // no interpolation between points.  If the schedule is (t0,
          // dt0, t1, dt1, t2, dt2) then the value is: if time < t0
          // value = 0 if t0 < time < t1 value = dt0 if t1 < time < t2
          // value = dt1 if t2 < time value = dt2
          t_len = node.operands.size()/2;
          i = node.valueIndex;
          // this evaluation gets the first time value from the
          // schedule data and puts the result in x1
          EXPReval_( *(node.operands[2*i+1]), x1, vals);

          // now we search futher time values (the odd values in the
          // array of data) to find values that bracket the current
          // time.
          while (i >= 0 && i< t_len)
          {
            // r1 is evaluated above as the implicit argument to this
            // function. (line 5461) in this case r1 is TIME
            if (x1 < r1)  // test if t0 (or tn) is less than TIME
            {
              if (i < t_len-1)
              {
                // there are still elements in the argument array
                // so evaluate the next time argument (t_(n+1))
                EXPReval_( *(node.operands[2*i+3]), x2, vals);
                if (x2 < x1)
                {
                  // mal-formed table with t_(n+1) < t_(n)
                  // given an error
                  EXPR_ERROR (EXPR_TABLE_FATAL);
                  // return the y values associated with t_(n)
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  res = y1;
                  break;
                }
                if (x2 > r1)
                {
                  // next time t_(n+1) > TIME so we found bracketing values
                  // get the y value associated with t_n and return it.
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  res = y1;
                  break;
                }
                else
                {
                  // didn't find times that bracket TIME yet.  So incrment
                  // our offset and copy t_(n_1) to x1 for the next loop.
                  ++i;
                  x1 = x2;
                }
              }
              else
              {
                // there isn't any more data in the array, so return
                // the last value
                EXPReval_( *(node.operands[2*i+2]), y1, vals);
                res = y1;
                break;
              }
            }
            else if (x1 > r1)    // test if t0 (or tn) is > TIME
            {
              if (i > 0)
              {
                // we may be able to back up the argument array
                // set x2 to our upper time bound
                x2 = x1;
                // get the next lower time bound
                EXPReval_( *(node.operands[2*i-1]), x1, vals);
                if (x2 < x1)
                {
                  // this is a mal-formed table with t_{n+1} < t_{n}
                  // so issue an error
                  EXPR_ERROR (EXPR_TABLE_FATAL);

                  // return value associated with one good point we have
                  EXPReval_( *(node.operands[2*i+2]), y1, vals);
                  res = y1;
                  break;
                }
                if (x1 < r1)     // test if t_n < TIME
                {
                  // found a set of points bracketing the current TIME
                  // so return the value assocaited wtih x1
                  EXPReval_( *(node.operands[2*i]), y1, vals);
                  res = y1;
                  break;
                }
                else
                {
                  // try looking further down the argument array
                  i--;
                }
              }
              else
              {
                // no more data in the argument array so return end result.
                EXPReval_( *(node.operands[2]), y1, vals);
                res = y1;
                break;
              }
            }
            else if (x1 == r1)
            {
              // TIME is equal to the left end of the range, so return the value associated
              // this with this end.
              EXPReval_( *(node.operands[2*i+2]), y1, vals);
              res = y1;
              break;
            }
          }
          node.valueIndex = i;
          // dubugging output
          //Xyce::dout() << "Leaving schedule eval (x1, y1), (x2, y2), (r1, res): (" << x1 << ", "
          //  << y1 << "), (" << x2 << ", " << y2 << "), (" << r1 << ", " << res << ")" << std::std::endl;
        }
        else
        {
          if (*node.fptr.unary == NULL)
          {
            EXPR_ERROR (EXPR_FUNCTION_FATAL);
            res = 0;
          }
          else
          {
            res = (*node.fptr.unary)(r1);
          }
        }
      }
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      EXPReval_( *(node.operands[0]), r1, vals);
      EXPReval_( *(node.operands[1]), r2, vals);
      res = (*node.fptr.binary)(r1, r2);
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        <<  "ExpressionInternals::EXPReval_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::EXPReval_: Internal Error: bad node type";
  }
  node.eval_value = res;
  node.eval_num = curr_num_;

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::clear_eval_num_
// Purpose       : Recursive method to set eval_num = 0 in a tree
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/29/04
//-----------------------------------------------------------------------------
///
/// Recursive method to set eval_num = 0 in a tree
///
/// @author Dave Shirley, PSSI
/// @date 09/29/04
///
void ExpressionInternals::clear_eval_num_ (ExpressionNode *n)
{
  int i;

  switch (n->type)
  {
    case EXPR_CONSTANT:
    case EXPR_VAR:
      n->eval_num = 0;
      break;

    case EXPR_OR:
    case EXPR_XOR:
    case EXPR_AND:
    case EXPR_EQUAL:
    case EXPR_NOTEQ:
    case EXPR_GREAT:
    case EXPR_GREATEQ:
    case EXPR_LESS:
    case EXPR_LESSEQ:
    case EXPR_PLUS:
    case EXPR_MINUS:
    case EXPR_TIMES:
    case EXPR_DIVIDE:
    case EXPR_REMAINDER:
    case EXPR_POWER:
      clear_eval_num_ (n->operands[0]);
      clear_eval_num_ (n->operands[1]);
      break;

    case EXPR_FUNCTION:
      for (i=0 ; i<n->operands.size() ; ++i)
      {
        clear_eval_num_ (n->operands[i]);
      }
      break;

    case EXPR_PLACEHOLDER:
      Report::DevelFatal()
        << "ExpressionInternals::clear_eval_num_: placeholder found";

    default:
      Report::DevelFatal()
        << "ExpressionInternals::clear_eval_num_: Unknown node type";
      break;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::compactLine_
// Purpose       : This routine is cut out of the original "set" routine.
//                 It takes an input line and squeezes out multiple white space,
//                 extraneous braces/parens, and checks parenthesis matching.
// Creator       : Tom Russo
// Creation Date : 08/03/10
//-----------------------------------------------------------------------------
///
/// This routine is cut out of the original "set" routine.
/// It takes an input line and squeezes out multiple white space,
/// extraneous braces/parens, and checks parenthesis matching.
///
/// @author Tom Russo
/// @date 08/03/10
///
void ExpressionInternals::compactLine_(std::string & inputLine,
                                             std::string & compactedLine)
{
  int brace_depth, paren_depth, start_depth, end_depth;
  std::string::iterator line_i;
  char thisChar;
  bool space;

  brace_depth = 0;
  paren_depth = 0;
  start_depth = 0;
  end_depth = 0;
  for (line_i = inputLine.begin() ; line_i != inputLine.end() ; ++line_i)
  {
    thisChar = toupper(*line_i);
    if (thisChar == '{')
      ++brace_depth;
    else if (thisChar == '}')
      brace_depth--;
    else if (thisChar == '(')
      ++paren_depth;
    else if (thisChar == ')')
      paren_depth--;
    if (brace_depth < 0 || paren_depth < 0)
    {
      Report::UserFatal()
        << "ExpressionInternals::ExpressionInternals: incorrect braces/parentheses found in expression: " << Input_;
    }
    if (compactedLine.size() == 0)
    {
      if (thisChar != ' ' && thisChar != '\t')
      {
        if (thisChar != '{')
        {
          compactedLine.push_back(thisChar);
          space = false;
          start_depth = brace_depth;
          end_depth = brace_depth;
        }
      }
    }
    else
    {
      if (thisChar != ' ' && thisChar != '\t' && thisChar != '\n')
      {
        if (thisChar != '}')
        {
          if (brace_depth < start_depth)
          {
            compactedLine.push_back('}');
            space = false;
            start_depth--;
            compactedLine.insert(0,"{");
          }
          compactedLine.push_back(thisChar);
          space = false;
          end_depth = brace_depth;
        }
        else
        {
          if (brace_depth >= start_depth)
          {
            compactedLine.push_back(thisChar);
            space = false;
            end_depth = brace_depth;
          }
        }
      }
      else
      {
        if (!space)
          compactedLine.push_back(' ');
          space = true;
      }
    }
  }
  if (brace_depth != 0 || paren_depth != 0 || start_depth != end_depth)
  {
      Report::UserFatal()
        <<   "ExpressionInternals::ExpressionInternals: unmatched parentheses found in expression: " << Input_;
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::tokenize_
// Purpose       : Lexical Analysis:  Given an input string (which should have
//                 been compacted by compactLine_ before calling this method),
//                 create a token list representing the input.
// Creator       : Tom Russo
// Creation Date : 08/0310
//-----------------------------------------------------------------------------
///
/// Lexical Analysis:  Given an input string (which should have
/// been compacted by compactLine_ before calling this method),
/// create a token list representing the input.
///
/// @author Tom Russo
/// @date 08/0310
///
void ExpressionInternals::tokenize_(std::string &inputLine,
                                          std::list<ExpressionElement *> & tokenList)
{
  int last_token, len;
  ExpressionElement *el;
  // Flag for whether we're inside of an operator like V(a+).  This is used to 
  // allow arithmetic symbols like + to be used in node names and device names
  // inside of operators in expressions.
  bool inOperator;  
  bool new_element;
  int name_len;
  double val;
  std::string::size_type j, k;

  last_token = -1;
  len = inputLine.size();
  el = newExpressionElement_();
  inOperator = false;
  for (int i=0 ; i<len ; ++i)
  {
    if (DEBUG_EXPRESSION)
      Xyce::dout() << "Processing character: " << i << " = '" << inputLine[i] << "'\n";

    new_element = false;
    if ( (inputLine[i] == '?') && (!inOperator) )
    {
      el->token = TOK_QUESTION;
      new_element = true;
    }
    else if (inputLine[i] == ':')
    {
      el->token = TOK_COLON;
      new_element = true;
    }
    else if (inputLine[i] == ',')
    {
      el->token = TOK_COMMA;
      new_element = true;
    }
    else if (inputLine[i] == ' ')
    {
      el->token = TOK_SPACE;
      new_element = true;
    }
    else if ( (inputLine[i] == '-') && (!inOperator) )
    {
      if (last_token == TOK_VALUE || last_token == TOK_RPAREN || last_token == TOK_RBRACE)
          el->token = TOK_MINUS;
      else
          el->token = TOK_UMINUS;
      new_element = true;
    }
    else if ( (inputLine[i] == '+') && (!inOperator) )
    {
      if (last_token == TOK_VALUE || last_token == TOK_RPAREN || last_token == TOK_RBRACE)
      {
        el->token = TOK_PLUS;
        new_element = true;
      }
    }
    else if ( (inputLine[i] == '*') && (!inOperator) )
    {
      if (i<len-1 && inputLine[i+1] == '*')
      {
        el->token = TOK_POWER;
        ++i;
      }
      else
      {
        el->token = TOK_TIMES;
      }
      new_element = true;
    }
    else if ( (inputLine[i] == '<')  && (!inOperator) )
    {
      if (i<len-1 && inputLine[i+1] == '=')
      {
        el->token = TOK_LESSEQ;
        ++i;
      }
      else
      {
        el->token = TOK_LESS;
      }
      new_element = true;
    }
    else if ( (inputLine[i] == '>') && (!inOperator) ) 
    {
      if (i<len-1 && inputLine[i+1] == '=')
      {
        el->token = TOK_GREATEQ;
        ++i;
      }
      else
      {
        el->token = TOK_GREAT;
      }
      new_element = true;
    }
    else if (inputLine[i] == '=')
    {
      if (i<len-1 && inputLine[i+1] == '=')
      {
        el->token = TOK_EQUAL;
        ++i;
      }
      else
      {
        el->token = TOK_EQUAL;
      }
      new_element = true;
    }
    else if ( (inputLine[i] == '!') && (!inOperator) )
    {
      if (i<len-1 && inputLine[i+1] == '=')
      {
        el->token = TOK_NOTEQ;
        ++i;
      }
      else
      {
        Report::UserFatal()
          << "ExpressionInternals::ExpressionInternals: Exclamation point found in expression:\n"
          << Input_ + "\n"
          <<  "Not equal operator is '!='";
      }
      new_element = true;
    }
    else if ( (inputLine[i] == '/') && (!inOperator) )
    {
      el->token = TOK_DIVIDE;
      new_element = true;
     }
    else if ( (inputLine[i] == '%') && (!inOperator) )
    {
      el->token = TOK_REMAINDER;
      new_element = true;
     }
    else if ( (inputLine[i] == '|') && (!inOperator) )
    {
      if (!useHspiceMath)
      {
        // Xyce uses one | for Logical OR
        el->token = TOK_OR;
        new_element = true;
      }
      else if (i<len-1 && inputLine[i+1] == '|')
      {
        // HSPICE uses two || for Logical OR
        el->token = TOK_OR;
        new_element = true;
        ++i;
      }
    }
    else if ( (inputLine[i] == '^') && (!inOperator) )
    {
      // ^ is a synonym for ** in HSPICE.  In Xyce, it is XOR.
      el->token = (useHspiceMath) ? TOK_POWER : TOK_XOR;
      new_element = true;
    }
    else if ( (inputLine[i] == '&') && (!inOperator) )
    {
      if (!useHspiceMath)
      {
        // Xyce uses one & for Logical AND
        el->token = TOK_AND;
        new_element = true;
      }
      else if (i<len-1 && inputLine[i+1] == '&')
      {
        // Xyce uses two && for Logical AND
        el->token = TOK_AND;
        new_element = true;
        ++i;
      }
    }
    else if ( (inputLine[i] == '~') && (!inOperator) )
    {
      el->token = TOK_NOT;
      new_element = true;
    }
    else if (inputLine[i] == '(') 
    {
      el->token = TOK_LPAREN;
      new_element = true;
          }
    else if (inputLine[i] == ')')
    {
      el->token = TOK_RPAREN;
      new_element = true;
      // can't still be in an operator, like V(a), if we see a right paren
      inOperator = false; 
    }
    else if (inputLine[i] == '{')
    {
      el->token = TOK_LBRACE;
      new_element = true;
    }
    else if (inputLine[i] == '}')
    {
      el->token = TOK_RBRACE;
      new_element = true;
    }
    // numberChars is a string with the characters (.0123456789) that are valid in a number.
    // So, this next clause handles numbers when they are outside of operators.
    else if ( (numberChars.find(inputLine[i],0) != std::string::npos) && (!inOperator) )
    {
      j = inputLine.find_first_not_of (".0123456789", i);
      if (j == std::string::npos)
        j = len-1;
      if (inputLine[j] == 'E')
      {
        ++j;
        if (j<len-1 && (inputLine[j] == '-' || inputLine[j] == '+'))
          ++j;
        j = inputLine.find_first_not_of ("0123456789", j);
      }
      k = inputLine.find_first_of (specials, j);
      if (k == std::string::npos)
        k = inputLine.length()-1;
      else
        k = k-1;
      ExtendedString diag = inputLine.substr(i,k-i+1);
      if (diag.isValue())
      {
        i = k;
        val = diag.Value();
      }
      else
      {
        Report::UserFatal()
          << "Illegal format for number (" << diag << ") in expression: "
          << Input_;
      }
      el->token = TOK_VALUE;
      el->type = TYP_NUM;
      el->number = val;
      new_element = true;
    }      
    else
    {
      std::string specialStrings = specialsNoColon;
      // this is the default case.  
      if (inOperator == 0)
      {
        // specials are the characters that are not allowed in either 
        // node names or device names.
        j = inputLine.find_first_of(specialStrings, i);
      }
      else
      {
        // operatorSpecials is the less restrictive set of characters that 
        // are not allowed in node and device names.  This is a smaller set because
        // node names only occur within operators (e.g, V(a)) in expressions.  This also
        // handles devices names that occur within operators in expressions.
        j = inputLine.find_first_of(operatorSpecials, i);
      }

      if (j == std::string::npos)
        j = len;
      while (inputLine[j] == '!')
      {
        if (inputLine[j + 1] != '=')
        {
          if (inOperator == 0)
          {
            j = inputLine.find_first_of(specialStrings, j + 1);
          }
          else
	  {
	    j = inputLine.find_first_of(operatorSpecials, j + 1);
          }
          if (j == std::string::npos)
            j = len;
        }
        else
          break;
      }

      el->token = TOK_VALUE;
      el->type = TYP_STRING;
      el->name = "";
      name_len = j-i;
      el->name = inputLine.substr(i, name_len);

      // determine if the following set of characters, in inputLine, is the start of an 
      // operator, like VR(a)
      if (el->name == "VR" || el->name == "VI" || el->name == "VM" || el->name == "VP" || el->name == "VDB" ||
	  el->name == "IR" || el->name == "II" || el->name == "IM" || el->name == "IP" || el->name == "IDB" ||
          el->name == "DNI" || el->name == "DNO")
      {
        if ( (inputLine.length() > i+name_len) && (inputLine[i+name_len] == '(') ) inOperator = true;
      }
      else if (name_len == 1)
      {
        // operators like V(a) and I(a)
        if (inputLine[i] == 'V' || inputLine[i] == 'I' || inputLine[i] == 'P' || inputLine[i] == 'W' ||
            inputLine[i] == 'N' )
        {
	  if ( (inputLine.length() > i+1) && (inputLine[i+1] == '(') ) inOperator = true;
        }
      }
      else if (name_len == 2)
      {
        // Handle multi-letter branch current operators (but not IF).  Examples are IC(Q1), IB(Q1)
        // and IE(Q1) for the Q device.
        if (inputLine[i] == 'I' && inputLine[i + 1] != 'F')
        {
          if ( (inputLine.length() > i+2) && (inputLine[i+2] == '(') ) inOperator = true;
        }
      }
      i = j-1;
      new_element = true;

    }
    if (new_element)
    {
      if (el->token != TOK_SPACE)
        last_token = el->token;
      tokenList.push_back(el);
      el = newExpressionElement_();
    }
  }
  el->token = TOK_END;
  tokenList.push_back(el);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::convertPolyToExpr_
// Purpose       : This method massages the awkward "POLY" format into
//                 a normal polynomial expression.
// Special Notes : This implementation works on the tokenized output of the
//                 lexical analyzer, and replaces it with the equivalent
//                 token list that would have been created had the user
//                 specified the polynomial without the POLY statement.
// Creator       : Tom Russo
// Creation Date : 08/04/10
//-----------------------------------------------------------------------------
///
/// This method massages the awkward "POLY" format into
/// a normal polynomial expression.
///
/// @author Tom Russo
/// @date 08/04/10
///
bool ExpressionInternals::convertPolyToExpr_(std::list<ExpressionElement *> & tokenList)
{
  std::vector<std::list<ExpressionElement *> *> tok2;
  std::list<ExpressionElement *>::iterator tok_i, tok_j;
  ExpressionElement *el;
  int n_nodes, n_space, pos;
  std::vector<int> node_count;
  std::vector<int>::iterator node_i;
  std::list<ExpressionElement *>  *e_list;
  int i,ii;

  tok_i = tokenList.begin();
  el = *tok_i;
  if (el->token == TOK_VALUE && el->type == TYP_STRING && el->name == "POLY")
  {
    n_nodes = -1;
    el = *++tok_i;
    if (el->token == TOK_SPACE)
      el = *++tok_i;
    if (el->token == TOK_LPAREN)
    {
      el = *++tok_i;
      if (el->token == TOK_VALUE && el->type == TYP_NUM)
      {
        n_nodes = static_cast<int> (el->number);
        el = *++tok_i;
        if (el->token == TOK_RPAREN)
          el = *++tok_i;
        else
          n_nodes = -1;
      }
    }
    if (n_nodes <= 0)
    {
      Report::UserError0() << "Syntax error in number of nodes in expression: " << Input_;
      return false;
    }
    if (el->token == TOK_SPACE)
      ++tok_i;
    tok_j = tok_i;
    n_space = 0;
    tok2.clear();
    while (el->token != TOK_END)
    {
      e_list = new std::list<ExpressionElement *>;
      if (true || n_space < n_nodes)
      {
        el = newExpressionElement_();
        el->token = TOK_LPAREN;
        e_list->push_back(el);
      }
      el = *tok_j;
      if (el->token == TOK_MINUS)
        el->token = TOK_UMINUS;
      while (el->token != TOK_SPACE && el->token != TOK_END)
      {
        e_list->push_back(el);
        el = *++tok_j;
      }
      if (true || n_space < n_nodes)
      {
        el = newExpressionElement_();
        el->token = TOK_RPAREN;
        e_list->push_back(el);
      }
      ++n_space;
      el = *tok_j;
      if (el->token == TOK_SPACE)
        el = *++tok_j;
      tok2.push_back(e_list);
    }
    if (n_space <= n_nodes)
    {
      Report::UserFatal()
        << "ExpressionInternals::ExpressionInternals: Apparently no coefficients in expression:"
        << Input_;
    }
    tokenList.resize(0);
    copy_elements_ (tokenList, tok2[n_nodes]);
    node_count.clear();
    node_count.push_back(-1);
    for (i=n_nodes+1 ; i<n_space ; ++i)
    {
      pos = 0;
      for (node_i = node_count.begin() ; node_i != node_count.end() ; ++node_i)
      {
        ++(*node_i);
        if (*node_i >= n_nodes)
        {
          if ((node_i+1) == node_count.end())
          {
            node_count.push_back(0);
            for (node_i = node_count.begin() ; node_i != node_count.end() ; ++node_i)
              (*node_i) = 0;
            break;
          }
        }
        else
        {
          if (pos > 0)
          {
            for (ii=0 ; ii<pos ; ++ii)
              node_count[ii] = node_count[pos];
          }
          break;
        }
        ++pos;
      }
      pos = 0;
      for (node_i = node_count.begin() ; node_i != node_count.end() ; ++node_i)
      {
        if (pos++ == 0)
        {
          el = newExpressionElement_();
          el->token = TOK_PLUS;
          tokenList.push_back(el);
          copy_elements_ (tokenList, tok2[i]);
        }
        el = newExpressionElement_();
        el->token = TOK_TIMES;
        tokenList.push_back(el);
        copy_elements_ (tokenList, tok2[*node_i]);
      }
    }
    el = newExpressionElement_();
    el->token = TOK_END;
    tokenList.push_back(el);

    for (i=0 ; i<tok2.size() ; ++i)
      delete tok2[i];
  }
  // Here endeth POLY handling.  If we had POLY, then now we have an actual
  // polynomial instead.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::standardizeTable_
// Purpose       : Put all TABLE input into a clean, standard, function-like
//                 format.  There are alternate inputs thanks to PSPICE
//                 compatibility.
// Creator       : Tom Russo
// Creation Date : 07/28/10
//-----------------------------------------------------------------------------
///
/// Put all TABLE input into a clean, standard, function-like
/// format.  There are alternate inputs thanks to PSPICE
/// compatibility.
///
/// @author Tom Russo
/// @date 07/28/10
///
void ExpressionInternals::standardizeTable_(std::list<ExpressionElement *> &tokenList)
{
  std::list<ExpressionElement *>::iterator tok_i;
  ExpressionElement *el;

  // The next section makes sure that what we end up with is a tokenized
  // representation of:
  // TABLE(<expression>,[<pair of values>]*)
  // or
  // SCHEDULE(TIME[,<comma separated pair of values>]*)
  //
  // At issue is that because this feature was added to comply with PSPICE
  // compatibility, there are multiple ways the user might have written a
  // table.
  //.
  // For example, in an E source, one could have:
  //    Efoo a b TABLE <expression>=[(<c.s. pair of values>)]*
  // or  in a B source:
  //    Bjunk a b V={TABLE {<expression>} [(<c.s. pair of values>)]*
  // or the way we really want it:
  //    Bjunk a b V={TABLE(<expression>,[<c.s. pair of values>])*
  //
  // In the first case (the E source), pass 1 of the parser has already
  // helpfully massaged it into the following format:
  //    Efoo a b TABLE {<expression}>=[(<c.s. pair of values>)]*
  // adding the braces around the expression.
  //
  // Note that in the last case, there are no parentheses around the pairs.
  // That appears only to be legal if we've also enclosed the expression
  // in braces.  Furthermore, if we had braces around the expression, we
  // *MUST* have parens around the pairs.  What a mess.
  //
  // So this section attempts to bring all the variants into the common,
  // simpler function-like form.  If there are no braces around the
  // expression, we apparently assume all else is fine, because we do
  // nothing.

  tok_i = tokenList.begin();
  el = *tok_i;
  if (el->token == TOK_VALUE && el->type == TYP_STRING && ((el->name == "TABLE") || (el->name == "SCHEDULE")))
  {
    // for schedules there is an implicit TIME variable that we need
    // to insert.  we do that insertion just after parsing the
    // TOK_LPAREN
    if (el->name=="SCHEDULE")
    {
      std::list<ExpressionElement *>::iterator tokenIteratorCopy( tok_i );
      tokenIteratorCopy++;  // this points us to LPAREN
      tokenIteratorCopy++;  // this is the value after LPAREN
      ExpressionElement * newExpressionElement = newExpressionElement_();
      newExpressionElement->token = TOK_VALUE;
      newExpressionElement->type = TYP_STRING;
      newExpressionElement->name = "TIME";
      tokenList.insert( tokenIteratorCopy, newExpressionElement );
      newExpressionElement = newExpressionElement_();
      newExpressionElement->token = TOK_COMMA;
      tokenList.insert( tokenIteratorCopy, newExpressionElement );
   }

    // advance to the next token. if it is a LBRACE, then we need to
    // reformat this.  if it isn't then we don't do anything
    el = *(++tok_i);
    if (el->token == TOK_LBRACE)
    {
      // Convert the LBRACE to an LPAREN, search forward for the RBRACE, and
      // convert *it* to a comma.  If we don't find it, boom.
      el->token = TOK_LPAREN;
      while (el->token != TOK_RBRACE)
      {
        if (++tok_i == tokenList.end())
          break;
        el = *tok_i;
      }
      if (el->token == TOK_RBRACE)
      {
        el->token = TOK_COMMA;
      }
      else
      {
        Report::UserFatal()
          << "ExpressionInternals::ExpressionInternals: Right brace not found in table expression: "
          << Input_;
      }

      // So now we have TABLE(<expression>,
      // Now get rid of an unnecessary equals sign if there is one
      // (as there might be if this were an E source, for example)
      ++tok_i;
      if ((*tok_i)->token == TOK_EQUAL)
        tok_i = tokenList.erase(tok_i);

      // If we had braces around the expression, then we MUST have parens
      // around the pairs.  Check that here.
      int i = 0;
      if ((*tok_i)->token != TOK_LPAREN)
      {
        Report::UserFatal() 
          << "ExpressionInternals::ExpressionInternals: Left parenthesis missing in:"
          << Input_;
      }

      // So now we know we do have parens around all our pairs.  Remove them.
      // Also sanity check the syntax of the pairs themselves.  Only unary minus
      // and values are allowed in the pairs.
      while ((*tok_i)->token != TOK_END)
      {
        ++i;
        if ((*tok_i)->token == TOK_LPAREN)
        {
          tok_i = tokenList.erase(tok_i);
        }
        else
        {
          Report::UserFatal()
             << "ExpressionInternals::ExpressionInternals: Left parenthesis not found in pair # "
             << i
             << " of expression:"
             << Input_;
        }
        if ((*tok_i)->token == TOK_UMINUS)
          ++tok_i;
        if ((*tok_i)->token != TOK_VALUE)
        {
          Report::UserFatal()
            << "ExpressionInternals::ExpressionInternals: Syntax error in first value of pair # "
            << i
            << " of expression: "
            << Input_;
        }
        if ((*++tok_i)->token != TOK_COMMA)
        {
          Report::UserFatal()
             << "ExpressionInternals::ExpressionInternals: Comma not found in pair # "
             << i
             << " of expression: "
             << Input_;
        }
        if ((*++tok_i)->token == TOK_UMINUS)
          ++tok_i;
        if ((*tok_i)->token != TOK_VALUE)
        {
          Report::UserFatal()
             << "ExpressionInternals::ExpressionInternals: Syntax error in second value of pair # "
             << i
             << " of expression: "
             << Input_;
        }
        if ((*++tok_i)->token != TOK_RPAREN)
        {
          Report::UserFatal()
            << "ExpressionInternals::ExpressionInternals: Syntax error in first value of pair # "
            << i
            << " of expression: "
            << Input_;
        }
        el = *tok_i++;
        if ((*tok_i)->token != TOK_END)
          el->token = TOK_COMMA;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::braceToParen_
// Purpose       : convert any braces that remain inside an expression to
//                 parens
// Special Notes : Must be called *AFTER* standardizeTable_, which depends
//                 on some braces being present in some syntaxes.
// Creator       : Tom Russo
// Creation Date : 01/08/19
//-----------------------------------------------------------------------------
///
/// Convert all braces that remain in the token list after a standardizeTable_
/// call into the equivalent parentheses.
///
/// @author Tom Russo
/// @date 01/08/19
///
void ExpressionInternals::braceToParen_(std::list<ExpressionElement *> &tokenList)
{
  std::list<ExpressionElement *>::iterator tok_i;
  ExpressionElement *el;

  for (tok_i = tokenList.begin(); tok_i != tokenList.end(); ++tok_i)
  {
    el = *tok_i;
    if (el->token == TOK_LBRACE)
      el->token = TOK_LPAREN;
    if (el->token == TOK_RBRACE)
      el->token = TOK_RPAREN;
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::seedRandom
// Purpose       : Public method to initialize random number generator
//                 used by rand(), gauss() and agauss() functions.
// Creator       : Tom Russo
// Creation Date : 03/28/17
//-----------------------------------------------------------------------------
///
/// Public method to initialize random number generator
///
/// @author Tom Russo
/// @date 03/28/17
///
void ExpressionInternals::seedRandom(long seed)
{
  if (enableRandomExpression)
  {
    if (theRandomNumberGenerator == 0)
    {
      theRandomNumberGenerator = new Xyce::Util::RandomNumbers(seed);
    }
    else
    {
      theRandomNumberGenerator->seedRandom(seed);
    }
  }
}
  

// These methods only for debugging and exploring the package behavior
//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::dumpParseTree
// Purpose       : Public method to cause dump of a parse tree
// Creator       : Tom Russo
// Creation Date : 07/28/10
//-----------------------------------------------------------------------------
///
/// Public method to cause dump of a parse tree
///
/// @author Tom Russo
/// @date 07/28/10
///
void ExpressionInternals::dumpParseTree()
{
  dumpParseTree_(tree_);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::dumpParseTree_
// Purpose       : Make a formatted dump of a parse tree.  This is a recursive
//                 method.
// Creator       : Tom Russo
// Creation Date : 07/28/10
//-----------------------------------------------------------------------------
///
/// Make a formatted dump of a parse tree.  This is a recursive
/// method.
///
/// @author Tom Russo
/// @date 07/28/10
///
void ExpressionInternals::dumpParseTree_(const ExpressionNode *tree, int indentLevel)
{
  std::string varTypesStrings[]={"NODE","STRING","INSTANCE","SPECIAL","VARIABLE","FUNCTION","NODAL_COMPUTATION"};

  if (indentLevel == 0)
  {
    Xyce::dout() << " Dumping parse tree: " << std::endl;
  }

  indentWithDashes_(indentLevel);   Xyce::dout() << expr_ops[tree->type] << std::endl;

  if (tree->type == EXPR_CONSTANT)
  {
    indentWithDashes_(indentLevel);   Xyce::dout() << " Value : "
                                        << tree->constant << std::endl;
  }
  else if (tree->type == EXPR_VAR)
  {
    indentWithDashes_(indentLevel);  Xyce::dout() << " Variable:  " << tree->valueIndex
                                          << " Name: "
                                          << varValues_[tree->valueIndex]
                                          << " Type: "
                                          << varTypesStrings[varTypes_[tree->valueIndex]-10]
                                          << " and node constant is "
                                          << tree->constant
                                          << std::endl;
  }
  else
  {
    indentWithDashes_(indentLevel);   Xyce::dout() << " Number of operands: "
                                          << tree->operands.size() << std::endl;
  }

  int nOps=tree->operands.size();
  for (int i=0; i< nOps; i++)
    dumpParseTree_(tree->operands[i],indentLevel+1);
}


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::indentWithDashes__
// Purpose       : Format a line with dashes as indentation, for use by
//                  debugging tool above
// Creator       : Tom Russo
// Creation Date : 07/28/10
//-----------------------------------------------------------------------------
///
/// Format a line with dashes as indentation, for use by
/// debugging tool above
///
/// @author Tom Russo
/// @date 07/28/10
///
void ExpressionInternals::indentWithDashes_(int indentLevel)
{
  for (int i=0;i<indentLevel;i++)
    Xyce::dout() << "-";
}


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::isNameSpecial_
// Purpose       : Return true if given name is one of the "special" names
// Creator       : Tom Russo
// Creation Date : 10/05/16
//-----------------------------------------------------------------------------
///
/// Return true if given name is one of the "special" names
///
/// @param name  string name to check
///
/// @author Tom Russo
/// @date 10/05/16
///
bool ExpressionInternals::isNameSpecial_(const std::string &name) const
{
  bool retval=false;

  if (name[0] == '#')
    retval = true;
  else
  {
    for (int i=0; i< sizeof(specSigs)/sizeof(specSigs[0]); i++)
    {
      if (specSigs[i] == name)
      {
        retval = true;
        break;
      }
    }
  }

  return retval;
}

// basic functions pointed to by function pointers in the expression tree.

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPRor
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRor(double arg1, double arg2)
{
  LOGIC(arg1, EXPR_OR_WARNING);
  LOGIC(arg2, EXPR_OR_WARNING);
  if (arg1 > 0 || arg2 > 0)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRxor(double arg1, double arg2)
{
  LOGIC(arg1, EXPR_XOR_WARNING);
  LOGIC(arg2, EXPR_XOR_WARNING);
  if ( (arg1 > 0 && arg2 <= 0) ||
       (arg1 <= 0 && arg2 > 0) )
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRand(double arg1, double arg2)
{
  LOGIC(arg1, EXPR_AND_WARNING);
  LOGIC(arg2, EXPR_AND_WARNING);
  if (arg1 > 0 && arg2 > 0)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRequal(double arg1, double arg2)
{
  if (arg1 == arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRnoteq(double arg1, double arg2)
{
  if (arg1 != arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRgreat(double arg1, double arg2)
{
  if (arg1 > arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRgreateq(double arg1, double arg2)
{
  if (arg1 >= arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRless(double arg1, double arg2)
{
  if (arg1 < arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRlesseq(double arg1, double arg2)
{
  if (arg1 <= arg2)
    return 1.;
  return 0.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRplus(double arg1, double arg2)
{
  double rval;

  rval = arg1 + arg2;
  if (fabs(rval) > EXPR_HUGE)
    EXPR_ERROR(EXPR_SUM_WARNING);
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRminus(double arg1, double arg2)
{
  double rval;

  rval = arg1 - arg2;
  if (fabs(rval) > EXPR_HUGE)
    EXPR_ERROR(EXPR_DIFF_WARNING);
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRtimes(double arg1, double arg2)
{
  double rval;

  rval = arg1 * arg2;
  if (fabs(rval) > EXPR_HUGE)
    EXPR_ERROR(EXPR_TIMES_WARNING);
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRdivide(double arg1, double arg2)
{
  double rval;

  if (arg2 == 0.0)
  {
      EXPR_ERROR(EXPR_DIVIDE_FATAL);
      return EXPR_HUGE;
  }
  rval = arg1 / arg2;
  if (fabs(rval) > EXPR_HUGE)
    EXPR_ERROR(EXPR_DIVIDE_WARNING);
  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRremainder(double arg1, double arg2)
{
  double res, rval;

  if (arg2 == 0.0)
  {
    EXPR_ERROR(EXPR_REMAINDER_FATAL);
    return EXPR_HUGE;
  }
  res = fabs(arg1/arg2);
  if (sizeof(double) >= 8)
  {
    if (res > 1.e+12)
    {
      EXPR_ERROR(EXPR_REMAINDER_ERROR);
      if (res > 1.e+15)
      {
        EXPR_ERROR(EXPR_REMAINDER_FATAL);
        return 0;
      }
    }
  }
  else
  {
    if (res > 1.e+3)
    {
      EXPR_ERROR(EXPR_REMAINDER_ERROR);
      if (res > 1.e+6)
      {
        EXPR_ERROR(EXPR_REMAINDER_FATAL);
        return 0;
      }
    }
  }
  rval = fabs(arg1) - (static_cast<long> (res))*fabs(arg2);
  if (arg1 < 0)
    rval = -rval;

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRpower(double arg1, double arg2)
{
  int i, j;
  double x(arg1);
  double y(arg2);
  double iy(0.0), rval(0.0);

  if (arg1 == 0)
  {
    return 0.;
  }

  if (arg2 == 0)
  {
    return 1.;
  }
  else if (arg2 < 0)
  {
    x = 1.0/x;
  }

  y = fabs(arg2);
  j = static_cast<int> (y);
  iy = static_cast<double> (j);
  if (y == iy) // this means arg2 is really an int.
  {
    rval = x;
    for (i=1 ; i<j ; ++i)
    {
      rval *= x;
    }
  }
  else
  {
    if (arg1 < 0)
    {
      x *= -1.0;
    }

    rval = pow(x, y);
    if (arg1 < 0)
    {
      EXPR_ERROR(EXPR_POWER_ERROR);
      rval = -rval;
    }
  }

  if (fabs(rval) > EXPR_HUGE)
  {
    EXPR_ERROR(EXPR_POWER_WARNING);
  }

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRabs (double arg)
{
  return fabs(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRacos (double arg)
{
  return acos(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRacosh (double arg)
{
  if (arg < 1.0)
      arg = 1.0;
  return (log(arg + sqrt(arg*arg-1.0)));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRasin (double arg)
{
  return asin(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRasinh (double arg)
{
  return (log(arg + sqrt(arg * arg + 1.0)));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPRint
// Creator       : Keith Santarelli, 01437
// Creation Date : 06/25/08
//-----------------------------------------------------------------------------
double EXPRint (double arg)
{
  int b;
  double c;

  b = arg;
  c = b;
  return c;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPRfloor
// Creator       : Tom Russo
// Creation Date : 12/12/2017
//-----------------------------------------------------------------------------
double EXPRfloor (double arg)
{
  return floor(arg);
}


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::EXPRceil
// Creator       : Tom Russo
// Creation Date : 12/12/2017
//-----------------------------------------------------------------------------
double EXPRceil (double arg)
{
  return ceil(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRatan (double arg)
{
  return atan(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRatanh (double arg)
{
  if (arg < Epsilon -1.0)
    arg = Epsilon - 1.0;
  else if (arg > 1.0 - Epsilon)
    arg = 1.0 - Epsilon;
  return (log((1.0 + arg) / (1.0 - arg)) / 2.0);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRcos (double arg)
{
  return (cos(MODULUS(arg, 2 * M_PI)));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRcosh (double arg)
{
  if (fabs(arg) < log(0.999*EXPR_HUGE))
  {
    return cosh(arg);
  }
  EXPR_ERROR(EXPR_COSH_ERROR);
  return(EXPR_HUGE);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRexp (double arg)
{
  if (arg < log(EXPR_HUGE))
    return exp(arg);
  EXPR_ERROR(EXPR_EXP_ERROR);
  return(EXPR_HUGE);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRln (double arg)
{
  if (arg <= 0)
  {
    EXPR_ERROR(EXPR_LOG_FATAL);
    return(-EXPR_HUGE);
  }
  return log(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRlog (double arg)
{
  if (arg <= 0)
  {
    EXPR_ERROR(EXPR_LOG_FATAL);
    return(-EXPR_HUGE);
  }
  return log10(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRnot (double arg)
{
  LOGIC(arg, EXPR_NOT_WARNING);
  if (arg > 0)
    return 0.;
  return 1.;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRsgn (double arg)
{
  return (arg > 0.0 ? 1.0 : (arg < 0.0 ? -1.0 : 0.0));
}


//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRsin (double arg)
{
  return sin(MODULUS(arg, 2 * M_PI));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRsinh (double arg)
{
  if (fabs(arg) < log(0.999*EXPR_HUGE))
  {
    return sinh(arg);
  }
  EXPR_ERROR(EXPR_SINH_ERROR);
  return(EXPR_HUGE);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRsqrt (double arg)
{
  if (arg < 0.0)
  {
    arg = -arg;
    EXPR_ERROR(EXPR_SQRT_ERROR);
  }
  return sqrt(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRtan (double arg)
{
  return tan(MODULUS(arg, M_PI));
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRtanh (double arg)
{
  if (arg > 20)
    return 1.;
  if (arg < -20)
    return -1.;
  return tanh(arg);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRuminus (double arg)
{
  return -arg;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::
// Creator       : Dave Shirley, PSSI
// Creation Date : 09/09/04
//-----------------------------------------------------------------------------
double EXPRuramp (double arg)
{
  if (arg < 0.0)
    return 0.0;
  else
    return arg;
}


} // namespace Util
} // namespace Xyce
