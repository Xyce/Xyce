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
// Purpose        : Provide solution/state/store/etc. access to expressions
//                  and print statements
//
// Special Notes  :
//
// Creator        : Dave Baur, Raytheon 
//
// Creation Date  : 11/15/2013
//
//
//
//
//-------------------------------------------------------------------------
///
/// Ops are used to calculate a single complex value.  A collection of
/// Ops is generally created during initialization and then at specific
/// times, the list of operations is executed which calculates the
/// values.
///
/// For example, the circuit time of the simulation must be printed at
/// the conclusion of each time step during the results output.  To
/// implement this, a OutputMgrTimeOp operator is created during
/// initialization which implements the operator(comm) function.  The
/// construction of the OutputMgrTimeOp operator takes the OutputMgr
/// reference as an argument which is queired for the curcuit time
/// during the operator(comm) call as output is being performed.
///
/// This implementation allows at wide variety of values to the
/// calculated while the caller only needs to provide the communicator.
/// Of course, there glaring are exceptions to the rule.  Since some
/// data's location is transient in nature (Solution, State and Store
/// values, for example), the sources of these values must be provided
/// to the caller.  This is done via the getValue() function by
/// providing these containers.
///
/// Lastly, these operators must be defined to be runin parallel with
/// the caveat that the operators may not actually be constructable on
/// all processors.  I.E. that function's data may only exist on a
/// subset of the processors.  This means that during operator
/// construction, parallel consist dummy implementations must be
/// provided to processor not having the operator constructed directly.
///
/// To provide this functionality, the operators implement a ReduceOp
/// object that must behave the same way as the actually operator froma
/// parallel standpoint.  This is done by implementing the Op and
/// ReduceOp_ templates.  These templates take reduction and evaluation
/// policies that implement the data collection and the reduction
/// operations.  The ReduceNone and ReduceSum policy classes
/// implement the static function reduce() which implements the
/// reduction as a Noop or a summation, respectively.  Other reduction
/// operations can be defined, but are not currently needed.
///
/// When defining an operator, the operator class needs to implement the
/// get(op) static function which actually gets the value.
///
#ifndef Xyce_N_UTL_Op_h
#define Xyce_N_UTL_Op_h

#include <map>
#include <string>
#include <vector>
#include <iterator>

#include <Teuchos_SerialDenseMatrix.hpp>
  
#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>

#include <N_LAS_Vector.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Param.h>
#include <N_ANP_NoiseData.h>

namespace Xyce {
namespace Util {
namespace Op {

typedef size_t Identifier;      ///< Identifies the operator with a single value

//-----------------------------------------------------------------------------
// Function      : identifier_base
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 17 17:00:40 2014
//-----------------------------------------------------------------------------
///
/// Since operating systems like to move the executable around from to run,
/// an address difference must be used as this will continue to be the same
/// regardless of where the address map begins.
///
/// @return an Identifier to be used as a basis
///
Identifier identifier_base();

//-----------------------------------------------------------------------------
// Function      : identifier
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 06:51:28 2014
//-----------------------------------------------------------------------------
///
/// Returns a application unique identifier for the specified Operator class.
///
/// The identifier is actaully created by the linker by using the
/// class'ss id_() static function's address.  It is important that this
/// identifier be unique within the run and common to all processors
/// running the application in parallel.
///
/// @return Returns a application unique identifier for the specified Operator class.
///
///
template<class T>
inline Identifier identifier() {
  return T::id_() - identifier_base();
}

typedef std::map<std::string, Teuchos::SerialDenseMatrix<int, std::complex<double> > * > RFparamsData;

struct OpData 
{
  OpData()
  : currentIndex_(0),
    realSolutionVector_(0),
    imaginarySolutionVector_(0),
    stateVector_(0),
    realStoreVector_(0),
    imaginaryStoreVector_(0),
    realLeadCurrentVector_(0),
    imaginaryLeadCurrentVector_(0),
    realLeadCurrentDeltaVVector_(0),
    imaginaryLeadCurrentDeltaVVector_(0),
    objectiveVector_(0),
    dOdpDirectVector_(0),
    dOdpDirectScaledVector_(0),
    dOdpAdjointVector_(0),
    dOdpAdjointScaledVector_(0),
    onoise_(0.0),
    inoise_(0.0),
    noiseDataVec_(0),
    RFparams_(0)
  {}

  OpData(const OpData &op_data)
  : currentIndex_(op_data.currentIndex_),
    realSolutionVector_(op_data.realSolutionVector_),
    imaginarySolutionVector_(op_data.imaginarySolutionVector_),
    stateVector_(op_data.stateVector_),
    realStoreVector_(op_data.realStoreVector_),
    imaginaryStoreVector_(op_data.imaginaryStoreVector_),
    realLeadCurrentVector_(op_data.realLeadCurrentVector_),
    imaginaryLeadCurrentVector_(op_data.imaginaryLeadCurrentVector_),
    realLeadCurrentDeltaVVector_(op_data.realLeadCurrentDeltaVVector_),
    imaginaryLeadCurrentDeltaVVector_(op_data.imaginaryLeadCurrentDeltaVVector_),
    objectiveVector_(op_data.objectiveVector_),
    dOdpDirectVector_(op_data.dOdpDirectVector_),
    dOdpDirectScaledVector_(op_data.dOdpDirectScaledVector_),
    dOdpAdjointVector_(op_data.dOdpAdjointVector_),
    dOdpAdjointScaledVector_(op_data.dOdpAdjointScaledVector_),
    onoise_(op_data.onoise_),
    inoise_(op_data.inoise_),
    noiseDataVec_(op_data.noiseDataVec_),
    RFparams_(op_data.RFparams_)
  {}

  OpData(
    int                                 current_index,
    const Linear::Vector *              real_solution_vector,
    const Linear::Vector *              imaginary_solution_vector,
    const Linear::Vector *              state_vector,
    const Linear::Vector *              real_store_vector,
    const Linear::Vector *              imaginary_store_vector,
    const Linear::Vector *              real_LeadCurrent_vector = 0,
    const Linear::Vector *              imaginary_LeadCurrent_vector = 0,
    const Linear::Vector *              real_LeadCurrentDeltaV_vector = 0,
    const Linear::Vector *              imaginary_LeadCurrentDeltaV_vector = 0,
    const std::vector<double> *         sens_function = 0,
    const std::vector<double> *         sens_dOdPDirect = 0,
    const std::vector<double> *         sens_dOdPDirectScaled = 0,
    const std::vector<double> *         sens_dOdPAdjoint = 0,
    const std::vector<double> *         sens_dOdPAdjointScaled = 0,
    double                              onoise = 0.0,
    double                              inoise = 0.0,
    const std::vector<Xyce::Analysis::NoiseData*> * noiseDataVec = 0,
    const RFparamsData * RFparams = 0
    )
  : currentIndex_(current_index),
    realSolutionVector_(real_solution_vector),
    imaginarySolutionVector_(imaginary_solution_vector),
    stateVector_(state_vector),
    realStoreVector_(real_store_vector),
    imaginaryStoreVector_(imaginary_store_vector),
    realLeadCurrentVector_(real_LeadCurrent_vector),
    imaginaryLeadCurrentVector_(imaginary_LeadCurrent_vector),
    realLeadCurrentDeltaVVector_(real_LeadCurrentDeltaV_vector),
    imaginaryLeadCurrentDeltaVVector_(imaginary_LeadCurrentDeltaV_vector),
    objectiveVector_(sens_function),
    dOdpDirectVector_(sens_dOdPDirect),
    dOdpDirectScaledVector_(sens_dOdPDirectScaled),
    dOdpAdjointVector_(sens_dOdPAdjoint),
    dOdpAdjointScaledVector_(sens_dOdPAdjointScaled),
    onoise_(onoise),
    inoise_(inoise),
    noiseDataVec_(noiseDataVec),
    RFparams_(RFparams)
  {}

  int                           currentIndex_;
  const Linear::Vector *        realSolutionVector_;
  const Linear::Vector *        imaginarySolutionVector_;
  const Linear::Vector *        stateVector_;
  const Linear::Vector *        realStoreVector_;
  const Linear::Vector *        imaginaryStoreVector_;
  const Linear::Vector *        realLeadCurrentVector_;
  const Linear::Vector *        imaginaryLeadCurrentVector_;
  const Linear::Vector *        realLeadCurrentDeltaVVector_;
  const Linear::Vector *        imaginaryLeadCurrentDeltaVVector_;

  const std::vector<double> *   objectiveVector_;
  const std::vector<double> *   dOdpDirectVector_;
  const std::vector<double> *   dOdpDirectScaledVector_;
  const std::vector<double> *   dOdpAdjointVector_;
  const std::vector<double> *   dOdpAdjointScaledVector_;

  const double                  onoise_;
  const double                  inoise_;
  const std::vector<Xyce::Analysis::NoiseData*> * noiseDataVec_;

  const RFparamsData * RFparams_;
};

//-----------------------------------------------------------------------------
// Class         : Operator
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 06:52:21 2014
//-----------------------------------------------------------------------------
///
/// Operators are functors that calculate a single complex value.  This
/// class serves as the base class for those functors.  The operator()
/// function performs the operation by calling the derived classes
/// evaluate() function.  Note that this is a parallel operation, so the
/// Parallel::Machine argument to the operator() function serves and the
/// communicator.
///
class Operator
{
public:
  //-----------------------------------------------------------------------------
  // Function      : Operator
  // Purpose       : Constructor
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:54:51 2014
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the Operator.
  ///
  /// The name along with the arguments vector may be used by the deriving
  /// classes to locate data needed for the calculation.
  ///
  /// @invariant Sets the operator name.
  ///
  /// @param name         Name of the operator
  ///
  Operator(const std::string &name)
    : name_(name)
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Operator
  // Purpose       : Destructor
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:15 2014
  //-----------------------------------------------------------------------------
  ///
  /// Destroys the operator
  ///
  virtual ~Operator()
  {}

  //-----------------------------------------------------------------------------
  // Function      : id
  // Purpose       : Return a unique identifier for the operator
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:42 2014
  //-----------------------------------------------------------------------------
  ///
  /// Application unique and parallel common identifier for this operator.
  ///
  /// @return an application unique and parallel common identifier for this operator
  ///
  ///
  virtual Identifier id() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : operator()
  // Purpose       : Functor interface to the evaluate() function.
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:37 2014
  //-----------------------------------------------------------------------------
  ///
  /// Functor interface to the operator evaluate() function.
  ///
  /// @return result of the derived evaluate() function
  ///
  ///
  complex operator()(Parallel::Machine comm, const OpData &op_data) const
  {
    return evaluate(comm, op_data);
  }

  //-----------------------------------------------------------------------------
  // Function      : evaluate
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:48 2014
  //-----------------------------------------------------------------------------
  ///
  /// Performs the calculation.
  ///
  /// @param comm       Communicator
  ///
  /// @return complex result
  ///
  virtual complex evaluate(Parallel::Machine comm, const OpData &op_data) const = 0;

  virtual complex get(const OpData &op_data) const = 0;

  virtual complex reduce(Parallel::Machine comm, const complex &value) const = 0;

  virtual Parallel::ReduceInterface *reduce2(complex &value) const = 0;

  virtual complex eval(const complex &value) const = 0;

  //-----------------------------------------------------------------------------
  // Function      : getName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:57 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the name of the operator
  ///
  /// @return   Returns the name of the operator
  ///
  ///
  const std::string &getName() const
  {
    return name_;
  }

  //-----------------------------------------------------------------------------
  // Function      : addArg
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:33:42 2014
  //-----------------------------------------------------------------------------
  ///
  /// Add argument name to the argument list.
  ///
  /// @invariant Argument list is appended with the name
  ///
  /// @param arg        Argument name
  ///
  ///
  void addArg(const std::string &arg)
  {
    argList_.push_back(arg);
  }

  //-----------------------------------------------------------------------------
  // Function      : addArgs
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:36:12 2014
  //-----------------------------------------------------------------------------
  ///
  /// Adds several argument names to the argument list.
  ///
  /// @invariant Argument list is appeneded with the specified names
  ///
  /// @param begin      begin iterator of the names to be added
  /// @param end        end iterator of the names to be added
  ///
  ///
  template <class II>
  void addArgs(II begin, II end)
  {
    argList_.assign(begin, end);
  }

  //-----------------------------------------------------------------------------
  // Function      : &getArgs
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:37:50 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the argument names of the operator
  ///
  /// @return   the argument names of the operator
  ///
  const std::vector<std::string> &getArgs()
  {
    return argList_;
  }

protected:
  const std::string                     name_;
  std::vector<std::string>              argList_;
};


//-----------------------------------------------------------------------------
// Policy Class  : ReduceNone
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 08:53:06 2014
//-----------------------------------------------------------------------------
///
/// Implements the reduce() static function which simply returns that
/// value since no reduction is needed.
///
struct ReduceNone
{
  //-----------------------------------------------------------------------------
  // Function      : reduce
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:54:48 2014
  //-----------------------------------------------------------------------------
  ///
  /// ReduceNone policy does no reduction, simply return source result
  ///
  /// @param comm         communicator
  /// @param result       source result
  ///
  /// @return             source result
  ///
  static complex reduce(Parallel::Machine comm, complex result)
  {
    return result;
  }

  static Parallel::ReduceInterface *reduce2(complex &result);
};

//-----------------------------------------------------------------------------
// Policy Class  : ReduceSum
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 08:53:06 2014
//-----------------------------------------------------------------------------
///
/// Implements the reduce() static function which performs a MPI_SUM reduction.
///
struct ReduceSum
{
  //-----------------------------------------------------------------------------
  // Function      : reduce
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:54:48 2014
  //-----------------------------------------------------------------------------
  ///
  /// ReduceSum policy does a MPI_SUM all reduction.
  ///
  /// @param comm         communicator
  /// @param result       source result
  ///
  /// @return             source result
  ///
  static complex reduce(Parallel::Machine comm, complex result);

  static Parallel::ReduceInterface *reduce2(complex &result);
};

//-----------------------------------------------------------------------------
// Policy Class  : ReduceAverage
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 08:53:06 2014
//-----------------------------------------------------------------------------
///
/// Implements the reduce() static function which performs a MPI_SUM reduction.
///
struct ReduceAverage
{
  //-----------------------------------------------------------------------------
  // Function      : reduce
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:54:48 2014
  //-----------------------------------------------------------------------------
  ///
  /// ReduceAverage policy does a count of non-zero values and a MPI_SUM all reduction and then divides by the non-zero
  /// count.  A serious hack if zero is a valid value.
  ///
  /// @param comm         communicator
  /// @param result       source result
  ///
  /// @return             source result
  ///
  static complex reduce(Parallel::Machine comm, complex result);

  static Parallel::ReduceInterface *reduce2(complex &result);
};


//-----------------------------------------------------------------------------
// Policy Class  : EvalNoop
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 08:57:01 2014
//-----------------------------------------------------------------------------
///
/// Implements the eval() static function which simply returns the value.
///
struct EvalNoop
{
  //-----------------------------------------------------------------------------
  // Function      : eval
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:57:47 2014
  //-----------------------------------------------------------------------------
  ///
  /// EvalNoop policy returns the source result
  ///
  /// @param result     source result
  ///
  /// @return           source result
  ///
  ///
  static complex eval(complex result)
  {
    return result;
  }
};


//-----------------------------------------------------------------------------
// Class Template: ReduceOp_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 07:45:04 2014
//-----------------------------------------------------------------------------
///
/// Implements the dummy operator to ensure that parallel operators are
/// consistent.
///
/// Since an operator may not actual be created on all processors due to
/// lack of need or data, a dummy operator that is parallel consistent
/// with the processors that do create the operator must be provided.
/// This class provides the implementation for those objects.  The R
/// policy defines the reduction operator and the E policy provides the
/// additional calcuation after reduction, if needed.
///
/// The create() function creates the dummy operator
///
/// @param T         Name of the operator
/// @param R         Name of the operator
/// @param E         Name of the operator
///
template <class T, class R, class E = T>
class ReduceOp_ : public Operator
{
public:
  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:50:09 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates the dummy parallel consistent operator.
  ///
  /// @param name         Name of the operator
  ///
  /// @return             Pointer to the newly created operator
  ///
  ///
  static Operator *create(const std::string &name)
  {
    return new ReduceOp_(name);
  }

  //-----------------------------------------------------------------------------
  // Function      : ReduceOp_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:52:06 2014
  //-----------------------------------------------------------------------------
  ///
  /// Constructors the operator.
  ///
  /// @invariant Sets the operator name.
  ///
  /// @param name       Name of the operator 
  ///
  ReduceOp_(const std::string &name)
    : Operator(name)
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~ReduceOp_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:53:30 2014
  //-----------------------------------------------------------------------------
  ///
  /// Destroys the operator.
  ///
  virtual ~ReduceOp_()
  {}

  //-----------------------------------------------------------------------------
  // Function      : id
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:42 2014
  //-----------------------------------------------------------------------------
  ///
  /// Application unique and parallel common identifier for this operator.
  ///
  /// @return an application unique and parallel common identifier for this operator
  ///
  ///
  virtual Identifier id() const
  {
    return identifier<T>();
  }

protected:
  //-----------------------------------------------------------------------------
  // Function      : evaluate
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:19:04 2014
  //-----------------------------------------------------------------------------
  ///
  /// Evaluate the operator.
  ///
  /// The evaluation processor consist of getting the value, reducing it
  /// across the processors with the reduction policy and the applying
  /// evaluation policy.
  ///
  /// @param comm       Communicator
  ///
  /// @return           evaluated complex result
  ///
  ///
  virtual complex evaluate(Parallel::Machine comm, const OpData &op_data) const
  {
    return eval_(reduce_(comm, get_(op_data)));
  }

  virtual complex get(const OpData &op_data) const
  {
    return get_(op_data);
  }

  virtual complex reduce(Parallel::Machine comm, const complex &value) const
  {
    return reduce_(comm, value);
  }

  virtual Parallel::ReduceInterface *reduce2(complex &value) const
  {
    return reduce2_(value);
  }

  virtual complex eval(const complex &value) const
  {
    return eval_(value);
  }

  //-----------------------------------------------------------------------------
  // Function      : get_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:16:56 2014
  //-----------------------------------------------------------------------------
  ///
  /// Dummy implementation data collector.
  ///
  /// In general, zero is a good value, but for some reduction policies it is inappropriate
  ///
  /// @return complex(0.0, 0.0)
  ///
  ///
  complex get_(const OpData &op_data) const
  {
    return complex(0.0, 0.0);
  }

  //-----------------------------------------------------------------------------
  // Function      : reduce_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:18:28 2014
  //-----------------------------------------------------------------------------
  ///
  /// Call the reduction policy reduce() function
  ///
  /// @param comm       communicator
  /// @param value          source value
  ///
  /// @return           reduced value
  ///
  ///
  complex reduce_(Parallel::Machine comm, complex value) const
  {
    return R::reduce(comm, value);
  }
  Parallel::ReduceInterface *reduce2_(complex &value) const
  {
    return R::reduce2(value);
  }

  //-----------------------------------------------------------------------------
  // Function      : eval_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:21:52 2014
  //-----------------------------------------------------------------------------
  ///
  /// Call the evaluation policy eval() function.
  ///
  /// @param x          source value
  ///
  /// @return           evaluated value
  ///
  ///
  complex eval_(complex x) const
  {
    return E::eval(x);
  }

  //-----------------------------------------------------------------------------
  // Function      : set_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:16:56 2014
  //-----------------------------------------------------------------------------
  ///
  /// Dummy implementation data setter.
  ///
  ///
  void set_(const OpData &op_data, complex value) const
  {}
};

//-----------------------------------------------------------------------------
// Class Template: Op
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 07:45:04 2014
//-----------------------------------------------------------------------------
///
/// Implements the operator.
///
/// In general, the derived class, passed in as T, implements the get()
/// static function which take this as an argument.  
///
template<class T, class R, class E>
class Op : public Operator
{
public:
  //-----------------------------------------------------------------------------
  // Function      : id_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:35:35 2014
  //-----------------------------------------------------------------------------
  ///
  /// Defines the application unique and parallel common identifier.
  ///
  static Identifier id_()
  {
    static Identifier id = reinterpret_cast<Identifier>(&id_);

    return id;
  }

  typedef Op<T, R, E> Base;                     ///< Make construction of base class easy
  typedef ReduceOp_<T, R, E> ReduceOp;          ///< Dummy operator for non-creating processors

  //-----------------------------------------------------------------------------
  // Function      : Op
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:52:06 2014
  //-----------------------------------------------------------------------------
  ///
  /// Constructors the operator.
  ///
  /// @invariant Sets the operator name.
  ///
  /// @param name       Name of the operator 
  ///
  Op(const std::string &name)
    : Operator(name)
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Op_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 07:53:30 2014
  //-----------------------------------------------------------------------------
  ///
  /// Destroys the operator.
  ///
  virtual ~Op()
  {}

  //-----------------------------------------------------------------------------
  // Function      : id
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 06:55:42 2014
  //-----------------------------------------------------------------------------
  ///
  /// Application unique and parallel common identifier for this operator.
  ///
  /// @return an application unique and parallel common identifier for this operator
  ///
  ///
  virtual Identifier id() const
  {
    return identifier<T>();
  }

protected:
  //-----------------------------------------------------------------------------
  // Function      : evaluate
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:19:04 2014
  //-----------------------------------------------------------------------------
  ///
  /// Evaluate the operator.
  ///
  /// The evaluation processor consist of getting the value, reducing it
  /// across the processors with the reduction policy and the applying
  /// evaluation policy.
  ///
  /// @param comm       Communicator
  ///
  /// @return           evaluated complex result
  ///
  ///
  virtual complex evaluate(Parallel::Machine comm, const OpData &op_data) const
  {
    return eval_(reduce_(comm, get_(op_data)));
  }


  virtual complex get(const OpData &op_data) const
  {
    return get_(op_data);
  }

  virtual complex reduce(Parallel::Machine comm, const complex &value) const
  {
    return reduce_(comm, value);
  }

  virtual Parallel::ReduceInterface *reduce2(complex &value) const
  {
    return reduce2_(value);
  }

  virtual complex eval(const complex &value) const
  {
    return eval_(value);
  }

  //-----------------------------------------------------------------------------
  // Function      : get_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:16:56 2014
  //-----------------------------------------------------------------------------
  ///
  /// Data collector.
  ///
  /// Uses the base provided get() static function to get the value.
  ///
  /// @return value from the 
  ///
  ///
  complex get_(const OpData &op_data) const
  {
    return T::get(static_cast<const T &>(*this), op_data);
  }

  //-----------------------------------------------------------------------------
  // Function      : set_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:16:56 2014
  //-----------------------------------------------------------------------------
  ///
  /// Data collector.
  ///
  /// Uses the base provided get() static function to get the value.
  ///
  /// @return value from the 
  ///
  ///
  void set_(const OpData &op_data, complex value) const
  {
    T::set(static_cast<const T &>(*this), op_data, value);
  }

  //-----------------------------------------------------------------------------
  // Function      : reduce_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:18:28 2014
  //-----------------------------------------------------------------------------
  ///
  /// Call the reduction policy reduce() function
  ///
  /// @param comm       communicator
  /// @param x          source value
  ///
  /// @return           reduced value
  ///
  ///
  complex reduce_(Parallel::Machine comm, complex x) const
  {
    return R::reduce(comm, x);
  }

  Parallel::ReduceInterface *reduce2_(complex &x) const
  {
    return R::reduce2(x);
  }

  //-----------------------------------------------------------------------------
  // Function      : eval_
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:21:52 2014
  //-----------------------------------------------------------------------------
  ///
  /// Call the evaluation policy eval() function.
  ///
  /// @param x          source value
  ///
  /// @return           evaluated value
  ///
  ///
  complex eval_(complex x) const
  {
    return E::eval(x);
  }
};

//-----------------------------------------------------------------------------
// Class         : UndefinedOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 08:41:19 2014
//-----------------------------------------------------------------------------
///
/// UndefinedOp operator.
///
/// The UndefinedOp is created when the processor cannot provide an
/// operator due to lack of data.  During initialization, this will be
/// replaced with dummy operator, Op<T, R, E>::ReduceOp.
///
class UndefinedOp : public Op<UndefinedOp, ReduceNone, EvalNoop>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : UndefinedOp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 08:43:50 2014
  //-----------------------------------------------------------------------------
  ///
  /// 
  ///
  /// @invariant
  ///
  /// @param name 
  ///
  /// @return 
  ///
  ///
  UndefinedOp(const std::string &name = "Undefined")
    : Base(name)
  {}

  virtual ~UndefinedOp()
  {}

  static complex get(const UndefinedOp &op, const OpData &op_data)
  {
    return complex(0.0, 0.0);
  }
};

//-----------------------------------------------------------------------------
// Class         : ConstantOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : 
//-----------------------------------------------------------------------------
class ConstantOp : public Op<ConstantOp, ReduceNone, EvalNoop>
{

public:
  ConstantOp(const std::string &name, complex value)
    : Base(name),
      value_(value)
  {}

  virtual ~ConstantOp()
  {}

  static complex get(const ConstantOp &op, const OpData &op_data)
  {
    return op.value_;
  }

  const complex                 value_;
};

typedef std::vector<Operator *> OpList;

Operator *
makeOp(
  Parallel::Machine             comm,
  const BuilderManager &        op_builder_manager,
  ParamList::const_iterator &   param_it);

void
makeOps(
  Parallel::Machine                     comm,
  const BuilderManager &                op_builder_manager,
  const NetlistLocation &               netlist_location,
  ParamList::const_iterator             begin,
  ParamList::const_iterator             end,
  std::back_insert_iterator<OpList>     inserter);

complex getValue(
  Parallel::Machine             comm,
  const Util::Op::Operator &    op,
  const OpData &                op_data);

// void
// getValues(
//   Parallel::Machine             comm,
//   const Util::Op::OpList &      op_list,
//   const OpData &                op_data,
//   std::vector<double> &         result_list);

void
getValues(
  Parallel::Machine             comm,
  const Util::Op::OpList &      op_list,
  const OpData &                op_data,
  std::vector<complex> &        result_list);

} // namespace Op
} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Op_h
