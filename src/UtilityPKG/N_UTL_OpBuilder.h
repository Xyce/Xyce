//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_OpBuilder_h
#define Xyce_N_UTL_OpBuilder_h

#include <map>
#include <string>
#include <vector>

#include <N_UTL_fwd.h>
#include <N_UTL_Op.h>
#include <N_PDS_fwd.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Util {
namespace Op {

class BuilderManager;

//-----------------------------------------------------------------------------
// Class         : Builder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 09:04:11 2014
//-----------------------------------------------------------------------------
///
/// Builds an Op given the start of a ParamList.
///
/// The start iterator is advanced to the parameter just past the
/// consumed parameters need to assemble the Op.  If the parameter list
/// does not match the requirements of the Op, them the iterator is
/// unchanged and 0 is returned.
///
/// Interface for an operator builder.
///
class Builder
{
public:
  //-----------------------------------------------------------------------------
  // Function      : Builder
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:05:44 2014
  //-----------------------------------------------------------------------------
  ///
  /// Constructs a Builder
  ///
  Builder()
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Builder
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:06:05 2014
  //-----------------------------------------------------------------------------
  ///
  /// Destroys a Builder
  ///
  /// @invariant
  ///
  ///
  ///
  virtual ~Builder()
  {}

private:
  Builder(const Builder &);
  Builder &operator=(const Builder &);

public:
  //-----------------------------------------------------------------------------
  // Function      : createOp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:08:01 2014
  //-----------------------------------------------------------------------------
  ///
  /// Attempts to create an operator based on the parameters starting at it
  ///
  /// If an operator is not created, the return value will be 0 and the
  /// iterator wil not be advanced.  If the operator is created, the
  /// iterator will be advanced on element past the last utilized to
  /// create the operator.
  ///
  /// @return   pointer to new operator or 0 if no operator could be created
  ///
  Operator *createOp(ParamList::const_iterator &it) const;

  //-----------------------------------------------------------------------------
  // Function      : registerCreateFunctions
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:19:22 2014
  //-----------------------------------------------------------------------------
  ///
  /// Registers the create functions that this builder can create.
  ///
  /// Since dummy operators may need to be created for processor lacking
  /// data to identify or create the operator, create functions for
  /// these dummy operators must be registered so that they can be
  /// created later.
  ///
  /// @param op_builder_manager         builder manager to register functions to
  ///
  virtual void registerCreateFunctions(BuilderManager &op_builder_manager) const = 0;

private:
  //-----------------------------------------------------------------------------
  // Function      : Operator *makeOp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:22:11 2014
  //-----------------------------------------------------------------------------
  ///
  /// Building implementation that attempts to create an operator based
  /// on the parameters starting at it
  ///
  /// @return   pointer to new operator or 0 if no operator could be created
  ///
  virtual Operator *makeOp(ParamList::const_iterator &it) const = 0;
};

typedef Operator *(*CreateFunction)(const std::string &name);           ///< Function pointer for creating operators

//-----------------------------------------------------------------------------
// Class         : BuilderManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 10 09:28:08 2014
//-----------------------------------------------------------------------------
///
/// The BuilderManager maintains the builders and the dummy operator
/// creators.
///
class BuilderManager
{
public:
  typedef std::vector<const Builder *> BuilderVector;                   ///< Vector of operator builders
  typedef std::map<Identifier, CreateFunction> CreateMap;               ///< Map from operator identifier to dummy operator creators

  //-----------------------------------------------------------------------------
  // Function      : BuilderManager
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:33:10 2014
  //-----------------------------------------------------------------------------
  ///
  /// Constructs a builder manager.
  ///
  /// Registers the UndefinedOp operator.
  ///
  ///
  BuilderManager()
  {
    addCreateFunction<UndefinedOp>();
  }

  //-----------------------------------------------------------------------------
  // Function      : ~BuilderManager
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:34:00 2014
  //-----------------------------------------------------------------------------
  ///
  /// Destroys a builder manager.
  ///
  /// @invariant
  ///
  ///
  ///
  ~BuilderManager()
  {
    for (BuilderVector::iterator it = opBuilderVector_.begin(), end = opBuilderVector_.end(); it != end; ++it)
      delete *it;
  }

private:
  BuilderManager(const BuilderManager &);
  BuilderManager &operator=(const BuilderManager &);

public:
  //-----------------------------------------------------------------------------
  // Function      : addBuilder
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:34:41 2014
  //-----------------------------------------------------------------------------
  ///
  /// Adds the builders and requests that the builder registers all
  /// dummy operator creators that of operators that it may create.
  ///
  /// @invariant The builder is added and owned by the build manager
  ///
  /// @param op_builder         builder
  ///
  ///
  void addBuilder(const Builder *op_builder) {
    opBuilderVector_.push_back(op_builder);
    op_builder->registerCreateFunctions(*this);
  }

  //-----------------------------------------------------------------------------
  // Function      : addCreateFunction
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:37:03 2014
  //-----------------------------------------------------------------------------
  ///
  /// Adds a dummmy operator create function to the dummy operator
  /// creator map.
  ///
  /// @invariant The creator map has the new creator added.  Any
  /// duplicates are replaced, but with the same create function.
  ///
  template<class T>
  void addCreateFunction() {
    opCreateMap_[identifier<T>()] = T::ReduceOp::create;
  }

  //-----------------------------------------------------------------------------
  // Function      : createOp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:39:33 2014
  //-----------------------------------------------------------------------------
  ///
  /// Attempts to create an operator for the parameters starting at it.
  ///
  /// Iterates through all the registered operator creators, each
  /// attempting to create an operator based on the parameters starting
  /// at it.  If successful, the parameter iterator is advanced just
  /// past the last parameters used to create the operator.
  ///
  /// @param it         parameter iterator
  ///
  /// @return   pointer to new operator or 0 if no operator could be created
  ///
  ///
  Operator *createOp(ParamList::const_iterator &it) const;
  Operator *createOp(const std::string &name) const;

  //-----------------------------------------------------------------------------
  // Function      : findCreateFunction
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jul 10 09:43:56 2014
  //-----------------------------------------------------------------------------
  ///
  /// Returns the create function for the operator defined by the identifier.
  ///
  /// @param id         operator identifier
  ///
  /// @return           the dummy operator create function or 0 is not found
  ///
  CreateFunction findCreateFunction(Identifier id) const;

private:
  BuilderVector         opBuilderVector_;
  CreateMap             opCreateMap_;
};

} // namespace Op
} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_OpBuilder_h
