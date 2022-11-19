//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

//----------------------------------------------------------------------------
//
// Purpose        : Declare the circuit "context".
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 01/21/2003
//
//----------------------------------------------------------------------------

#ifndef N_IO_CIRCUITCONTEXT_H
#define N_IO_CIRCUITCONTEXT_H

// ---------- Standard Includes ----------
#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>

// trilinos includes
#include <Teuchos_RCP.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Param.h>

#include <N_IO_fwd.h>
#include <N_IO_FunctionBlock.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_UTL_Param.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_NetlistLocation.h>

#include <N_UTL_Pack.h>
#include <N_ERH_Message.h>

#include <N_PDS_fwd.h>

namespace Xyce {
namespace IO {

typedef Device::DeviceCountMap DeviceCountMap;

//----------------------------------------------------------------------------
// Class          : CircuitContext
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
class CircuitContext
{
  friend class Pack<CircuitContext>;
  
public:
  struct MutualInductance
  {
    friend class Pack<MutualInductance>;
    
    std::map<std::string,std::string> inductors;
    std::string coupling;
    std::string model;
    std::string firstInductor;

    // set of inductor names and associated terminals
    std::map< std::string, std::vector< std::string > > terminals;
    
    std::map< std::string, std::vector<Device::Param > > otherParams;

    // MI name
    std::string name;

    // sharedInductorTable and fastLookupTable key
    int sharedKey;

    MutualInductance() {}

    MutualInductance( DeviceBlock & device );
  }; // end of MutualIndutance struct

  // Constructors.
  CircuitContext(
    //Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group,
    Util::Op::BuilderManager &          op_builder_manager,
    const ParsingMgr &                  parsing_manager,
    std::list<CircuitContext*> &        context_list,
    CircuitContext *&                   current_circuit_context);

  // Destructor.
  virtual ~CircuitContext();

  // Begin a subcircuit to context.
  // Note: the current subcircuit context will remain active until
  // an endSubcircuitContext is issued.
  bool beginSubcircuitContext(
    const std::string & netlistFileName,
    TokenVector &       subcircuitLine);

  // End the current context, push it onto the previous contexts
  // list of contexts, and reset the currentContext_ pointer to the
  // previous context.
  bool endSubcircuitContext();

  // Add a model to the current context. Note, the circuit context
  // will assume responsibility for model pointers.
  void addModel(ParameterBlock * modelPtr);

  // Add a set of .PARAM parameters to the current context.
  void addParams(Util::ParamList::const_iterator paramIter, Util::ParamList::const_iterator paramEnd);

  // strings are parsed with quotes still attached.
  // this function removes the quotes
  void resolveQuote (Util::Param & parameter) const;
  
  // resolve tablefile( "filename ") into a table expression
  void resolveTableFileType (Util::Param & parameter) const;
  
  //resolve string("a string value") into a string 
  //RLS void resolveStringType (Util::Param & parameter) const;

  // Add a set of .GLOBAL_PARAM parameters to the current context.
  void addGlobalParams(Util::ParamList::const_iterator paramIter, Util::ParamList::const_iterator paramEnd);

  // Add a global node
  void addGlobalNode (std::string &gnode);

  // move params to globals, or vice-versa, as necessary.
  void  categorizeParams( std::list<Util::OptionBlock> &  optionsTable);

  // Add a .FUNC function to the current context.
  void addFunction(FunctionBlock const& function);

  // Setters and getters.
  void setName(std::string const& name);
  const std::string& getName() const;
  const std::string& getCurrentContextName() const;
  const std::string& getPrefix() const;
  unordered_map<std::string, std::string> *getNodeMapPtr() const;
  void setParentContextPtr( CircuitContext * const ptr );
  const CircuitContext * getParentContextPtr() const;
  CircuitContext * getCurrentContextPtr() const;

  // This method will look up "name" in the circuitContextTable
  // and return the CircuitContext pointer if one is found.
  // NOTE:  This method is for finding subcircuit DEFINITIONS and is
  // probably only useful for the top level CircuitContext object.
  CircuitContext * getSubcircuitContextPtr( const std::string& name ) const;

  // Get the node list for the current context.
  std::vector<std::string> const& getNodeList() const;

  // Get the node list for the current context.
  ModelMap const& getModelMap() const;

  // Increment the device count in the current context.
  void incrementDeviceCount(const std::string &device = "");
  void incrementLinearDeviceCount();

  // Add a subcircuit name and instance name to the lists. 
  // The instance list will be needed to get a total device count.
  void addInstance(std::string const& subcktName,
                   std::string const& instanceName,
                   std::string const& fileName, 
                   int const& lineNumber);

  // Return the subckt (model) and instance (name) list for this object, 
  // which indicates all the subcircuits used at this level.
  std::vector<std::string> const& getSubcktList() const;

  // Get the netlist locations for each instance.
  unordered_map<std::string,std::list<NetlistLocation> > const& getInstanceLocation() const;

  void setLocation( const NetlistLocation& loc ); 
  NetlistLocation getLocation() const; 

  // Resolve the parameters and functions in the curren context. Since this
  // operation is dependent on a given subcircuit instance, the resolution
  // for a given context may occur repeatedly.
  bool resolve(std::vector<Device::Param> const& subcircuitInstanceParams);

  // Set the current context to that corresponding to the given subcircuit
  // name. Save the previous context on the stack for later retrieval.
  // If there is no context corresponding to the subcircuit name, recursively
  // search parent contexts. Return true if found, false otherwise.
  bool setContext(std::string const& subcircuitName,
                  std::string const& subcircuitPrefix = "",
                  std::vector<std::string> const& instanceNodes = std::vector<std::string>(),
                  CircuitContext* previousContext = NULL) const;
  void setContext(CircuitContext* context) const;

  // Reset the context the context prior to the last invocation of
  // setContext.
  void restorePreviousContext() const;
  bool globalNode (const std::string &nodeName) const;

  // ERK. new version, with no exceptions strings (i.e. function arguments)
  bool resolveParameter(Util::Param& parameter) const;

  // ERK. new version, with no exceptions strings (i.e. function arguments)
  bool resolveGlobalParameter(Util::Param& parameter) const;

  // ERK. new function for new expression.
  bool resolveParameterThatIsAdotFunc(Util::Param& parameter, std::vector<std::string> funcArgs) const; 

  // Determine if expressionString has any unresolved strings and
  // resolve appropriately. Return true if all strings are resolved
  // otherwise return false.
  bool resolveStrings(Util::Expression & expression,
                      std::vector<std::string> exceptionStrings = std::vector<std::string>()) const;

  // Determine if expressionString has any unresolved functions and
  // resolve appropriately. Return true if all functions are resolved
  // otherwise return false.
  bool resolveFunctions(Util::Expression & expression) const;

  // Look for a parameter with tag parameterName in resolvedParams_.
  // Check current context and recursively check parent
  // contexts. Return the parameter if it is found, set the
  // parameter value to the empty string if it is not found.
  bool getResolvedParameter(Util::Param & parameter) const;

  // Look for a parameter with tag parameterName in resolvedGlobalParams_.
  // Check current context and recursively check parent
  // contexts. Return the parameter if it is found, set the
  // parameter value to the empty string if it is not found.
  bool getResolvedGlobalParameter(Util::Param & parameter) const;

  // Look for a function with tag functionName in resolvedFunctions_.
  // Check current context and recursively check parent
  // contexts. Return the function (as an Util::Param) if it is found,
  // set the Util::Param value to the empty string if it is not found.
  bool getResolvedFunction(Util::Param & parameter) const;

  void addMutualInductance( DeviceBlock & device )
  {
    currentContextPtr_->mutualInductances_.push_back( MutualInductance( device ) );
  }

  std::vector<MutualInductance> & getMutualInductances()
  {
    return currentContextPtr_->mutualInductances_;
  }

  std::vector< std::set< std::string > > & getSharedInductorTable()
  {
    return currentContextPtr_->sharedInductorTable_;
  }

  std::set< std::string > & getAllCoupledInductors()
  {
    return currentContextPtr_->allCoupledInductors_;
  }

  std::vector< std::vector< int > > & getAllIndexedMIs()
  {
    return currentContextPtr_->allIndexedMIs_;
  }

  int getNumMILines()
  {
    return currentContextPtr_->kLines_.size();
  }

  // Convert all MIs into tokenized device lines
  void bundleMIs();

  // Retrieve one tokenized device line
  TokenVector & getMILine( int i );

  bool haveMutualInductances()
  {
    return !( currentContextPtr_->mutualInductances_.empty() );
  }

  int totalMutualInductanceCount();

  // Search the models in the current context for the model of the
  // given name. If it is not found, recursively search each parent
  // context. Return a pointer to the parameterBlock for the model
  // if it is found, otherwise return NULL. Also, if the model is found,
  // construct the appropriate model prefix.
  bool findModel(std::string const& modelName,
                 ParameterBlock* & modelPtr,
                 std::string& modelPrefix) const;

  bool findModel(std::string const& modelName, ParameterBlock* & modelPtr) const;

  bool fullyResolveParam(Device::Param & param, double & value) const;

  bool findBinnedModel(std::string const& modelName,
                 ParameterBlock* & modelPtr,
                 std::string& modelPrefix,
                 const bool LWfound, const bool LNFINfound,
                 const double L, const double W, const double NFIN, std::string & binNumber ) const;

  bool findBinnedModel(std::string const& modelName, ParameterBlock* & modelPtr, 
      const bool LWfound, const bool LNFINfound,
      const double L, const double W, const double NFIN, std::string & binNumber ) const;

  // Check whether a subcircuit context is dependent on subcircuit parameters.
  // These are parameters on the .subckt line identified by "params:"
  // keyword. The result should be true if either the current subcircuit
  // context or any subcircuit context in the hierarchy containing the
  // current subcircuit context has subcircuit parameters.
  bool hasSubcircuitParams();

  // Calculate the total number of devices starting at current context
  // and including all subcircuit instances.  
  // NOTE:  This function also accumulates the used context list and device type count map.
  int getTotalDeviceCount();

  // Get the current number of devices at the current context.
  // NOTE:  The device count may not be correct unless deviceCountDone_ returns true.
  std::pair<int, bool> getDeviceCount() { return std::pair<int,bool>( deviceCount_, deviceCountDone_ ); }

  // Get the total number of linear devices starting at current context
  // and including all subcircuit instances.
  // NOTE:  This can be called after getTotalDeviceCount() is called, since the accumulation
  // happens in that method.
  int getTotalLinearDeviceCount() { return currentContextPtr_->linearDeviceCount_; }

  // A list of used context names, so that the rest can be pruned.
  // NOTE:  This is only determined during the getTotalDeviceCount method, otherwise the list is empty.
  std::vector<std::string>& getUsedContextList()
  {
    return currentContextPtr_->usedContextList_;
  }

  // A map of device types to the number of that type in this context.
  const DeviceCountMap &getDeviceCountMap()
  {
    return currentContextPtr_->localDeviceCountMap_;
  }

  const Util::ParamMap &getFunctions() const
  {
    return resolvedFunctions_;
  }

  const Util::UParamList &getParams() const
  {
    return resolvedParams_;
  }

  const Util::UParamList &getGlobals() const
  {
    return resolvedGlobalParams_;
  }

  bool getContextMultiplierSet() 
  {
    bool retval=false;
    if (currentContextPtr_) { retval = currentContextPtr_->getMultiplierSet(); }
    return retval;
  }

#if 0
  double getContextMultiplierValue() 
  {
    double retval=1.0;
    if (currentContextPtr_) { retval = currentContextPtr_->getMultiplierValue(); }
    return retval;
  }
#else
  Util::Param getContextMultiplierParam () 
  {
    Util::Param parameter("","");
    if (currentContextPtr_) 
    { 
      parameter = currentContextPtr_->getMultiplierParam (); 
    }
    else
    {  // will this ever get called?  Put in error trap as a test.
      parameter.setTag( std::string("M") ); // is this sufficient or should it be fully resolved?  check this.
      parameter.setVal(1.0);
      Report::DevelFatal()<< "Mistake in function getContextMultiplierParam" ;
    }
    return parameter;
  }
#endif

  void setMultiplierSet(bool tmp)  { multiplierSet_ = tmp; }
  bool getMultiplierSet() { return multiplierSet_; }

#if 0
  void setMultiplierValue(double val) { multiplierValue_ = val; }
  double getMultiplierValue() { return multiplierValue_; }
#else
  void setMultiplierParam (Util::Param & param) { multiplierParameter_ = param; }
  Util::Param getMultiplierParam () { return multiplierParameter_; }
#endif

  // Traverse the CircuitContext table and remove all subcircuit instances except for
  // the ones with the names provided.
  void pruneContexts( std::vector<std::string>& keepContexts );

  // Correct total number of devices after processing K-devices:
  void augmentTotalDeviceCount(int kLineCount, int coupledICount, int YDeviceCount);

private:
  CircuitContext();

  // Copy Constructor.
  CircuitContext(CircuitContext const& rhsCC);
  CircuitContext &operator=(const CircuitContext& rhsCC);
       
  // Reference to a Pointer to the current context;
  // ERK.  This was once static data, but is now non-static, and ultimately
  // owned by the IO_NetlistImportTool class.  For that reason, it is a
  // reference to a pointer, as it still needs to act static.
  CircuitContext*& currentContextPtr_;

  Util::Op::BuilderManager &            opBuilderManager_;

  const ParsingMgr &                    parsingMgr_;

  mutable CircuitContext * parentContextPtr_;

  // Stack of contexts to track context changes.
  std::list<CircuitContext*> & contextList_;

  // Context name, corresponds to subcircuit names.
  std::string name_;

  // Location of where this subcircuit is defined or used.
  NetlistLocation location_;

  bool deviceCountDone_;
  int deviceCount_, linearDeviceCount_;
  std::vector<std::string> subcktList_, instanceList_;
  unordered_map<std::string, std::list<NetlistLocation> > instanceLocation_;

  // This data can be dropped at serialization time.
  std::vector<std::string>        usedContextList_;
  DeviceCountMap                localDeviceCountMap_;

  std::vector<std::string> nodeList_;
  Util::ParamList subcircuitParameters_;

  unordered_map< std::string, CircuitContext* > circuitContextTable_;

  ModelMap models_;

  Util::UParamList unresolvedParams_;
  std::set<std::string> globalNodes_;
  Util::UParamList unresolvedGlobalParams_;
  std::vector<FunctionBlock> unresolvedFunctions_;

  std::vector<MutualInductance> mutualInductances_;

  // lookup tables used to create semiPackedMIs_
  std::vector< std::set< std::string > > sharedInductorTable_;
  std::set< std::string > allCoupledInductors_;
  std::vector< std::vector< int > > allIndexedMIs_;

  // tokenized MIs
  std::vector< TokenVector >kLines_;

  // Each of the following attributes is not set until pass 2, so they
  // do not need to be serialized.
  std::string subcircuitPrefix_;
  unordered_map<std::string, std::string> nodeMap_; // note: does not need to be serialized.
  bool resolved_;
  Util::UParamList resolvedParams_;
  Util::UParamList resolvedGlobalParams_;
  Util::ParamMap  resolvedFunctions_;

  Teuchos::RCP<Xyce::Util::baseExpressionGroup> expressionGroup_; ///< required for setting up expressions

  bool multiplierSet_;
  Util::Param multiplierParameter_;
};

//----------------------------------------------------------------------------
// Function       : packCircuitContext
// Purpose        :
// Special Notes  :
// Scope          : public, non-member function
// Creator        : 
// Creation Date  : 07/30/2015
//----------------------------------------------------------------------------
int packCircuitContext(CircuitContext* circuit_contexts, char* char_buffer, 
                       int char_buffer_size, Parallel::Communicator* pds_comm_ptr );

//----------------------------------------------------------------------------
// Function       : unpackCircuitContext
// Purpose        :
// Special Notes  :
// Scope          : public, non-member function
// Creator        : 
// Creation Date  : 07/30/2015
//----------------------------------------------------------------------------
bool unpackCircuitContext(CircuitContext* circuit_contexts, char* char_buffer, 
                          int bsize, Parallel::Communicator* pds_comm_ptr );

//----------------------------------------------------------------------------
// Function       : CircuitContext::setName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline void CircuitContext::setName(std::string const& name)
{
  name_ = name;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getName() const
{
  return name_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getPrefix
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getPrefix() const
{
  return currentContextPtr_->subcircuitPrefix_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getNodeMapPtr
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/05/2003
//----------------------------------------------------------------------------
inline unordered_map<std::string, std::string> * CircuitContext::getNodeMapPtr() const
{
  return &(currentContextPtr_->nodeMap_);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getCurrentContextName
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
inline const std::string& CircuitContext::getCurrentContextName() const
{
  return currentContextPtr_->name_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getNodeList
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
inline std::vector<std::string> const& CircuitContext::getNodeList() const
{
  return currentContextPtr_->nodeList_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getModelMap
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
inline ModelMap const& CircuitContext::getModelMap() const
{
  return currentContextPtr_->models_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::incrementDeviceCount
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
inline void CircuitContext::incrementDeviceCount(const std::string &device)
{
  currentContextPtr_->deviceCount_++;
  if (device != "")
  {
    currentContextPtr_->localDeviceCountMap_[device]++;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::incrementLinearDeviceCount
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 06/25/2014
//----------------------------------------------------------------------------
inline void CircuitContext::incrementLinearDeviceCount()
{
  currentContextPtr_->linearDeviceCount_++;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getSubcktList
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 09/10/2014
//----------------------------------------------------------------------------
inline std::vector<std::string> const& CircuitContext::getSubcktList() const
{
  return currentContextPtr_->subcktList_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getInstanceLocation
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 09/10/2014
//----------------------------------------------------------------------------
inline unordered_map<std::string, std::list<NetlistLocation> > const& CircuitContext::getInstanceLocation() const
{
  return currentContextPtr_->instanceLocation_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getLocation
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 09/10/2014
//----------------------------------------------------------------------------
inline NetlistLocation CircuitContext::getLocation() const
{
  return currentContextPtr_->location_;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setLocation
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 09/10/2014
//----------------------------------------------------------------------------
inline void CircuitContext::setLocation( const NetlistLocation& loc ) 
{
  location_ = loc;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline void CircuitContext::setParentContextPtr(
   CircuitContext * const ptr )
{
  parentContextPtr_ = ptr;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getParentContextPtr
// Purpose        : accessor
// Special Notes  :
// Scope          : public
// Creator        :
// Creation Date  :
//----------------------------------------------------------------------------
inline const CircuitContext * CircuitContext::getParentContextPtr()
  const
{
  return parentContextPtr_;
}

inline CircuitContext * CircuitContext::getCurrentContextPtr()
  const
{
  return currentContextPtr_;
}

inline CircuitContext * CircuitContext::getSubcircuitContextPtr( const std::string& name )
  const
{
  unordered_map< std::string, CircuitContext* >::const_iterator cctIter = circuitContextTable_.find( name );
  if ( cctIter != circuitContextTable_.end() )
  {
    return cctIter->second;
  }
  
  return 0;
}


} // namespace IO
} // namespace Xyce

#endif
