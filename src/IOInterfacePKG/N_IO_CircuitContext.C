//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Define the circuit "context".
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 01/21/2003
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <sstream>

#include <N_UTL_fwd.h>

#include <N_DEV_DeviceBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_Op.h>
#include <N_IO_ParameterBlock.h>
#include <N_PDS_Comm.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {

//----------------------------------------------------------------------------
// Function       : CircuitContext::CircuitContext
// Purpose        : Constructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
CircuitContext::CircuitContext(
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group,
  Util::Op::BuilderManager &    op_builder_manager,
  std::list<CircuitContext*> &  context_list,
  CircuitContext *&             current_context_pointer)
  : currentContextPtr_(current_context_pointer),
    opBuilderManager_(op_builder_manager),
    parentContextPtr_(NULL),
    contextList_(context_list),
    name_(""),
    deviceCountDone_(false),
    deviceCount_(0),
    linearDeviceCount_(0),
    subcircuitPrefix_(""),
    resolved_(false),
    resolvedParams_(),
    resolvedGlobalParams_(),
    expressionGroup_(group)
{
  if (currentContextPtr_ == NULL)
  {
    currentContextPtr_ = this;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::~CircuitContext
// Purpose        : Destructor
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/21/2003
//----------------------------------------------------------------------------
CircuitContext::~CircuitContext()
{
  // Delete each subcircuit context in this context.
  unordered_map< std::string, CircuitContext * >::iterator itsc = circuitContextTable_.begin();
  unordered_map< std::string, CircuitContext * >::const_iterator itsc_end = circuitContextTable_.end();

  for ( ; itsc != itsc_end; ++itsc )
  {
    delete itsc->second;
  }

  circuitContextTable_.clear();


  // We delete all our own stored model pointers.  Because we also
  // delete all our subcircuit contexts, *their* destructors take care of
  // deleting their model pointers.
  ModelMap::iterator modelIter = models_.begin();
  ModelMap::iterator modelIterEnd = models_.end();
  for (; modelIter != modelIterEnd; ++modelIter)
    delete modelIter->second;
  models_.clear();
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addInstance
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
void
CircuitContext::addInstance(std::string const& subcktName, std::string const& instanceName,
                            std::string const& fileName, int const& lineNumber)
{
  std::string subcktNameUpper(ExtendedString(subcktName).toUpper());
  std::string instanceNameUpper(ExtendedString(instanceName).toUpper());

  currentContextPtr_->subcktList_.push_back(subcktNameUpper);
  currentContextPtr_->instanceList_.push_back(instanceNameUpper);
  currentContextPtr_->instanceLocation_[subcktNameUpper].push_back( NetlistLocation(fileName, lineNumber) );
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::beginSubcircuitContext
// Purpose        : Add a subcircuit context to the current context.
// Special Notes  : This routine sets the current context to a newly
//                  created CircuitContext object. This context
//                  remains the current context until it is explicitly
//                  terminated with a call to endSubcircuitContext.
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
bool CircuitContext::beginSubcircuitContext(
  std::string const& netlistFileName,
  TokenVector & subcircuitLine)

{
  if (DEBUG_IO)
  {
    Xyce::dout() << "CircuitContext::beginSubcircuitContext" << std::endl
                 << "*******************************************" << std::endl
                 << " subcircuit file name: " << netlistFileName << std::endl
                 << "  size of line vector: " << subcircuitLine.size() << std::endl;
  }

  // remove parens
  TokenVector::iterator iterLine =
    subcircuitLine.begin();
  while( iterLine != subcircuitLine.end() )
  {
    if ( iterLine->string_ == "(" || iterLine->string_ == ")" )
    {
      subcircuitLine.erase( iterLine );
    }
    else
    {
      iterLine++;
    }
  }

  // Create a new circuit context for the subcircuit.
  CircuitContext* subcircuitContextPtr =
    new CircuitContext(expressionGroup_, opBuilderManager_, contextList_, currentContextPtr_);

  // Set the parent context, save the current context and reset it to the
  // newly created context.
  subcircuitContextPtr->parentContextPtr_ = currentContextPtr_;
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = subcircuitContextPtr;

  // Save location of this subcircuit definition
  NetlistLocation newSubcktLoc( netlistFileName, subcircuitLine[0].lineNumber_ );
  subcircuitContextPtr->setLocation( newSubcktLoc );

  // Extract the subcircuit data from subcircuitLine.
  int numFields = subcircuitLine.size();

  if (numFields < 2)
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "Subcircuit name required";
    return false;
  }

  if (numFields < 3)
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "At least one node required for subcircuit " << subcircuitLine[1].string_;
    return false;
  }

  if (DEBUG_IO)
  {
    Xyce::dout() << "  Subcircuit name: " << subcircuitLine[1].string_ << std::endl
                 << "  Number of fields on subckt line = " << numFields << std::endl;
  }

  // Extract and set the subcircuit name.
  ExtendedString field ( subcircuitLine[1].string_ );
  field.toUpper();
  subcircuitContextPtr->setName(field);

  // Extract the subcircuit external nodes.
  int i;
  for (i = 2; i < numFields; ++i)
  {
    ExtendedString fieldES (  subcircuitLine[i].string_ );
    fieldES.toUpper();
    // Exit the loop if there are netlist parameters
    if (fieldES == "PARAMS:")
    {
      break;
    }
    if (i < numFields-1 && subcircuitLine[i+1].string_ == "=")
    {
      --i;
      break;
    }
    if (fieldES == "0")
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_) 
        << "Ground node '0' appears in .SUBCKT line";
    }
    subcircuitContextPtr->nodeList_.push_back(fieldES);
  }

  // Extract the subcircuit parameters.
  Util::Param parameter("","");

  bool result = true;
  ++i; // Advance to start of parameters.
  while (i+2 < numFields)
  {
    ExtendedString fieldES = subcircuitLine[i].string_;
    fieldES.toUpper();
    if (!fieldES.possibleParam())
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[i].lineNumber_)
        << "Parameter name " << subcircuitLine[i].string_ <<  " contains illegal character(s)";
      result = false;
    }
    parameter.setTag( fieldES );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") )
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[i].lineNumber_)
        << "Parameter name " << parameter.uTag() << " not allowed in subcircuit parameter list for subcircuit "
        << getName();
      result = false;
    }

    if ( subcircuitLine[i+1].string_ != "=" )
    {
      Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
        << "Equal sign required between parameter and value in PARAM list for subcircuit " << getName();
      result = false;
    }

    else {
      i+=2; // Advance past "=" sign

      fieldES =  subcircuitLine[i].string_;
      fieldES.toUpper();
      parameter.setVal(std::string(fieldES));
      subcircuitContextPtr->subcircuitParameters_.push_back(parameter);
      if (DEBUG_IO)
      {
        Xyce::dout() << " Adding parameter " << parameter.uTag() << " with value " << fieldES << std::endl;
      }
    }
    ++i;
  }

  // check for truncated params: list
  if ( result && i < numFields )
  {
    Report::UserError0().at(netlistFileName, subcircuitLine[0].lineNumber_)
      << "Parameter list error in subcircuit " << getName();
    result = false;
  }

  if (DEBUG_IO)
  {
    Xyce::dout() << "End of CircuitContext::beginSubcircuitContext" << std::endl
                 << "*******************************************" << std::endl;
  }

  return result;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::endSubcircuitContext
// Purpose        : End the current context, push it onto the previous context's
//                  list of contexts, and reset the current context pointer to
//                  the previous context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 01/24/2003
//----------------------------------------------------------------------------
bool CircuitContext::endSubcircuitContext()
{
  if (DEBUG_IO)
  {
    Xyce::dout() << "CircuitContext::endSubcircuitContext" << std::endl;
  }

  // Get the previous context from the stack.
  if ( ! contextList_.empty() )
  {
    // Add the current context to the previous context's list of contexts.
    contextList_.front()->circuitContextTable_[ currentContextPtr_->name_ ] =
      currentContextPtr_;

    // Reset the current context to the previous context.
    currentContextPtr_ = contextList_.front();

    // remove from stack
    contextList_.pop_front();
    return true;
  }
  else
  {
    return false;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addModel
// Purpose        : Add a model to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addModel(ParameterBlock * modelPtr)
{
  ParameterBlock* tmpModel;
  if (findModel(modelPtr->getName(), tmpModel))
  {
    Report::UserWarning0 message;
    message << "Reading model named " << modelPtr->getName() << " in the ";

    if (getCurrentContextName() == "")
      message << "main circuit";
    else
      message << "subcircuit " << getCurrentContextName();

    message << " and found one or more models previously defined in this scope";
  }

  currentContextPtr_->models_[modelPtr->getName()] = modelPtr;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addParams
// Purpose        : Add a set of .PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addParams(
  Util::ParamList::const_iterator paramIter, 
  Util::ParamList::const_iterator paramEnd)
{
  Util::Param parameter;
  for ( ; paramIter != paramEnd; ++paramIter)
  {
    parameter = *paramIter;
    resolveQuote(parameter);
    resolveTableFileType(parameter);
    currentContextPtr_->unresolvedParams_.insert(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addGlobalParams
// Purpose        : Add a set of .GLOBAL_PARAM parameters to the current context.
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
void CircuitContext::addGlobalParams(
  Util::ParamList::const_iterator paramIter, 
  Util::ParamList::const_iterator paramEnd)
{
  Util::Param parameter;
  for ( ; paramIter != paramEnd; ++paramIter)
  {
    parameter = *paramIter;
    resolveQuote(parameter);
    resolveTableFileType(parameter);
    currentContextPtr_->unresolvedGlobalParams_.push_back(parameter);
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::categorizeParams
// Purpose        : move params to globals, or vice-versa, as necessary.  This
//                  should only ever be called on the top-level circuit context.
//
//                  This function was added to address issue #167, "Make
//                  globally-scoped .params available to the device package,
//                  so globally scoped normal params will be equivalent to
//                  .global_param"
//
// Special Notes  : As of this writing, this function only elevates (as
//                  needed) certain params from normal .param to .global_param.
//                  It does not currently lower .global_params down to .param,
//                  but this should be added.
//
//                  .param and .global_param are very similar, since the
//                  develoment of the new expression library, they can
//                  both change during calculation.  Earlier in Xyce development
//                  this was not true, and .param was treated as a constant.
//
//                  At this point, one of the biggest ways in which they are
//                  different is how they are treated by the device package.
//                  The device package tracks all parameters that depend
//                  on .global_params, which makes it possible to do "setParam"
//                  calls on them.  These calls are typically made to support
//                  things like STEP or SAMPLING.
//
//                  I experimented with converting all top-level .param
//                  be .global_params.  This "worked" but it caused efficiency
//                  problems in modern PDK netlists.  So, best solution is to keep
//                  as many parameters in the ".param" state as possible, and
//                  keep an absolute minimum of params in the ".global_param"
//                  state.
// Scope          :
// Creator        : Eric Keiter
// Creation Date  : 04/13/2021
//----------------------------------------------------------------------------
void  CircuitContext::categorizeParams( std::list<Util::OptionBlock> &  optionsTable)
{
  if (DEBUG_IO)
    std::cout << "CircuitContext::categorizeParams.  ERK:DEBUG option table contents:" <<std::endl;

  {
  std::list<Util::OptionBlock>::iterator  iter = optionsTable.begin();
  std::list<Util::OptionBlock>::iterator  end = optionsTable.end();
  bool sortingNeeded=false;
  for ( ; iter != end; ++iter)
  {
    if (DEBUG_IO)
    {
      Xyce::dout() << "Options name: ";
      Xyce::dout() << iter->getName() << std::endl;
    }

    if (
        (iter->getName() == "STEP") ||
        (iter->getName() == "DC") ||
        (iter->getName() == "DATA") ||
        (iter->getName() == "SAMPLING") ||
        (iter->getName() == "EMBEDDEDSAMPLING") ||
        (iter->getName() == "SENS") ||
        (iter->getName() == "PCE") ||
        (iter->getName() == "LOCA"))
        { sortingNeeded=true; break;}
  }
  if (!sortingNeeded) return;
  }

  {
  std::list<Util::OptionBlock>::iterator  iter = optionsTable.begin();
  std::list<Util::OptionBlock>::iterator  end = optionsTable.end();
  bool sortingNeeded=false;
  for ( ; iter != end; ++iter)
  {
    if ( (iter->getName() == "STEP") || (iter->getName() == "DC") || (iter->getName() == "DATA") )
    {
      Util::ParamList::iterator iterPar = iter->begin();
      Util::ParamList::iterator endPar  = iter->end ();
      for ( ; iterPar != endPar; ++iterPar)
      {
        if (DEBUG_IO)
        {
        Xyce::dout() << "STEP/DC param " << iterPar->tag() <<std::endl;
        }
        if (iterPar->tag() == "PARAM")
        {
          Util::Param parameter(iterPar->stringValue(), "");
          Util::UParamList::const_iterator urParamIter = unresolvedParams_.find( parameter );
          if ( urParamIter != unresolvedParams_.end() )
          {
            parameter = *urParamIter;
            unresolvedGlobalParams_.push_back(parameter);
            unresolvedParams_.erase(urParamIter);
          }
        }
      }
    }
    else if ( (iter->getName() == "SAMPLING")
           || (iter->getName() == "EMBEDDEDSAMPLING")
           || (iter->getName() == "SENS")
           || (iter->getName() == "PCE"))
    {
      Util::ParamList::iterator iterPar = iter->begin();
      Util::ParamList::iterator endPar  = iter->end ();
      for ( ; iterPar != endPar; ++iterPar)
      {
        if (DEBUG_IO)
        {
        Xyce::dout() << "SAMPLING/UQ param " << iterPar->tag() <<std::endl;
        }
        if ( std::string( iterPar->tag() ,0,5) == "PARAM")
        {
          Util::Param parameter(iterPar->stringValue(), "");
          Util::UParamList::const_iterator urParamIter = unresolvedParams_.find( parameter );
          if ( urParamIter != unresolvedParams_.end() )
          {
            parameter = *urParamIter;
            unresolvedGlobalParams_.push_back(parameter);
            unresolvedParams_.erase(urParamIter);
          }
        }
      }
    }
    else if ( (iter->getName() == "LOCA"))
    {
      Util::ParamList::iterator iterPar = iter->begin();
      Util::ParamList::iterator endPar  = iter->end ();
      for ( ; iterPar != endPar; ++iterPar)
      {
        if (DEBUG_IO)
        {
        Xyce::dout() << "LOCA param " << iterPar->tag() <<std::endl;
        }
        if ( std::string( iterPar->tag() ,0,8) == "CONPARAM")
        {
          Util::Param parameter(iterPar->stringValue(), "");
          Util::UParamList::const_iterator urParamIter = unresolvedParams_.find( parameter );
          if ( urParamIter != unresolvedParams_.end() )
          {
            parameter = *urParamIter;
            unresolvedGlobalParams_.push_back(parameter);
            unresolvedParams_.erase(urParamIter);
          }
        }
      }
    }
  }
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::addGlobalNode
// Purpose        : Add a global node from a .GLOBAL line
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/02/2009
//----------------------------------------------------------------------------
void CircuitContext::addGlobalNode( std::string & gnode)
{
  currentContextPtr_->globalNodes_.insert(gnode);
}


//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveQuote
// Purpose        : Remove quotes from string parameters
// Special Notes  : Strings are parsed with double quotes still at the beginning 
//                  and ending.  This function removes them
// Scope          :
// Creator        : 
// Creation Date  :
//----------------------------------------------------------------------------
void CircuitContext::resolveQuote(Util::Param & parameter) const
{
  if (parameter.isQuoted())
  {
    std::ifstream paramDataIn;
    std::string parameterData(parameter.stringValue().substr(1,parameter.stringValue().size()-2));
    std::string paramString = parameterData;
    parameter.setVal( paramString );
    return;
  }
}
 

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveTableFileType
// Purpose        : Resolve quoted parameters as soon as they are encountered.
// Special Notes  : This avoids every processor in a parallel run trying to
//                  access the same file.  Also exit handling does not work on
//                  error for parallel runs otherwise.
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  :
//----------------------------------------------------------------------------
void CircuitContext::resolveTableFileType(Util::Param & parameter) const
{
  if (parameter.isTableFileTypeQuoted())
  {
    // The parameter is time dependent with its time history defined by the set
    // of time-value pairs in the file given by the value of the parameter.
    // Open and read these values, using them to build a "Table" expression
    // for the value of the parameter.
    std::ifstream paramDataIn;
     
    // have to be careful because on the instance line 
    // tablefile("filename") stays tablefile("filename")
    // but on the model line 
    // tablefile("filename") becomes tablefile"filename"
    // because parens are striped off on the model lines.
    
    int paramLen = parameter.stringValue().size();
    const int tablefileLen = std::string("tablefile").size();
    // set offset to 1 assuming just a quote around file name.
    int offset=1;

    if( (parameter.stringValue()[tablefileLen] == '(') && (parameter.stringValue()[tablefileLen+1] == '"') &&
        (parameter.stringValue()[paramLen-2] == '"') && (parameter.stringValue()[paramLen-1] == ')') )
    {
      offset=2;
    }
    std::string parameterFile(parameter.stringValue().substr(tablefileLen+offset,paramLen-(tablefileLen+2*offset)));

    std::string tableFileString = "{tablefile(\"" + parameterFile + "\")}";
    parameter.setVal( Util::Expression(expressionGroup_,tableFileString) );

    return;
  }
}


//----------------------------------------------------------------------------
// Function       : CircuitContext::addFunction
// Purpose        : Add a .FUNC function to the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/03/2003
//----------------------------------------------------------------------------
void CircuitContext::addFunction(FunctionBlock const& function)
{
  currentContextPtr_->unresolvedFunctions_.push_back(function);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolve
// Purpose        : Resolve parameters in the current context.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/04/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolve( std::vector<Device::Param> const& subcircuitInstanceParams )
{
  Util::ParamList retryParams;
  Util::UParamList uretryParams;
  std::vector<FunctionBlock> retryFunctions;

  if (currentContextPtr_->subcircuitParameters_.empty() &&
      subcircuitInstanceParams.empty() &&
      currentContextPtr_->resolved_)
  {
    // If there are no subcircuit parameters and this context has already
    // been resolved, then no work to do.
    return true;
  }

  if (DEBUG_IO)
  {
    Xyce::dout() << " CircuitContext::resolve has something to do..." << std::endl;
  }

  currentContextPtr_->resolved_ = true;

  // Clear resolved params and funcs.
  currentContextPtr_->resolvedParams_.clear();
  currentContextPtr_->resolvedGlobalParams_.clear();
  currentContextPtr_->resolvedFunctions_.clear();

  Util::ParamList asYetUnresolvedSubcircuitParameters=currentContextPtr_->subcircuitParameters_;
  Util::UParamList asYetUnresolvedParameters=currentContextPtr_->unresolvedParams_;
  Util::ParamList asYetUnresolvedGlobalParameters=currentContextPtr_->unresolvedGlobalParams_;
  std::vector<FunctionBlock> asYetUnresolvedFunctions=currentContextPtr_->unresolvedFunctions_;

  bool resolvedSomethingThisLoop=true;
  bool somethingLeftToDo= (!(asYetUnresolvedSubcircuitParameters.empty()) ||
                           !(asYetUnresolvedParameters.empty()) ||
                           !(asYetUnresolvedGlobalParameters.empty()) ||
                           !(asYetUnresolvedFunctions.empty())  ||
                           !(subcircuitInstanceParams.empty()) );

  while (resolvedSomethingThisLoop && somethingLeftToDo)
  {
    resolvedSomethingThisLoop=false;
    somethingLeftToDo=false;
    retryParams.clear();
    uretryParams.clear();
    retryFunctions.clear();

    // Add subcircuitParameters_ to the set of resolved parameters.
    // currentContextPtr_->resolvedParams_.addParameters( currentContextPtr_->subcircuitParameters_ );
    Util::Param* paramPtr;
    Util::Param parameter;
    Util::UParamList::iterator ustart, uend, uparamIter;
    Util::ParamList::iterator paramIter;
    Util::ParamList::iterator start =
      asYetUnresolvedSubcircuitParameters.begin();
    Util::ParamList::iterator end =
      asYetUnresolvedSubcircuitParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      parameter = *paramIter;
      if (DEBUG_IO)
      {
        Xyce::dout() << " CircuitContext::resolve attempting to resolve" << parameter.uTag()<< std::endl;
      }

      if (parameter.getType() == Xyce::Util::STR && !parameter.isNumeric())
      {
        ExtendedString arg ( parameter.stringValue() );
        arg.toUpper();
        if (arg.possibleParam())
        {
          parameter.setVal(std::string("{" + arg + "}"));
        }
      }

      if (!resolveParameter(parameter))
      {
        if (DEBUG_IO)
        {
          Xyce::dout() << "Unable to resolve subcircuit param " << parameter.uTag() << std::endl;
        }

        retryParams.push_back(parameter);
        somethingLeftToDo=true;
      }
      else
      {
        if (DEBUG_IO)
        {
          Xyce::dout() << "resolveParameter returned true on parameter " << parameter.uTag() << " after resolution its type is " << parameter.getType() << "with value " ;
          switch (parameter.getType()) {
            case Xyce::Util::STR:
              Xyce::dout() << parameter.stringValue();
              break;
            case Xyce::Util::DBLE:
              Xyce::dout() << parameter.getImmutableValue<double>();
              break;
            case Xyce::Util::EXPR:
              Xyce::dout() << parameter.getValue<Util::Expression>().get_expression();
              break;
            default:
              Xyce::dout() << parameter.stringValue();
          }
          Xyce::dout() << std::endl;
        }

        currentContextPtr_->resolvedParams_.insert(parameter);
        resolvedSomethingThisLoop=true;
      }
      if (DEBUG_IO)
      {
        Xyce::dout() << " CircuitContext::resolve done attempting to resolve" << parameter.uTag()<< std::endl;
      }
    }
    asYetUnresolvedSubcircuitParameters=retryParams;
    retryParams.clear();

    // Reset the subcircuit parameter values with the instance parameter
    // values as needed.
    int i;
    int numInstanceParameters = subcircuitInstanceParams.size();
    for ( i = 0; i < numInstanceParameters; ++i )
    {
      if (DEBUG_IO)
      {
        Xyce::dout() << " CircuitContext::resolve resetting subcircuit instance parameter" << subcircuitInstanceParams[i].uTag()<< " with value " << std::endl;
        switch (subcircuitInstanceParams[i].getType()) {
          case Xyce::Util::STR:
            Xyce::dout() << subcircuitInstanceParams[i].stringValue();
            break;
          case Xyce::Util::DBLE:
            Xyce::dout() << subcircuitInstanceParams[i].getImmutableValue<double>();
            break;
          case Xyce::Util::EXPR:
          {
            Util::Expression foo(subcircuitInstanceParams[i].getValue<Util::Expression>());
            Xyce::dout() << "EXPR(" << foo.get_expression() << ")";
            break;
          }
          default:
            Xyce::dout() << subcircuitInstanceParams[i].stringValue();
        }
        Xyce::dout() << std::endl;
      }

      // Look for the parameter in resolvedParams_, issue
      // a warning and continue if not found. Otherwise, if found
      // set the value.
       
      std::pair<Util::UParamList::iterator, bool> rP_iter = currentContextPtr_->resolvedParams_.insert(static_cast<const Util::Param &>(subcircuitInstanceParams[i]));
      if ( rP_iter.second )
      {
        if (DEBUG_IO)
        {
          Xyce::dout() << " did not find in resolvedParams_, adding parameter " << subcircuitInstanceParams[i].uTag() << std::endl;
        }
      }
      else
      {
        if (DEBUG_IO)
        {
          Xyce::dout() << " found in resolvedParams_, setting parameter " << subcircuitInstanceParams[i].uTag() << std::endl;
        }
       
        Util::Param& par = const_cast<Util::Param&>(*rP_iter.first);     
        par.setVal(static_cast<const Util::Param &>(subcircuitInstanceParams[i]));
      }
    }

    // Resolve any .PARAM parameters in the current context and add to
    // the resolvedParams_.
    ustart = asYetUnresolvedParameters.begin();
    uend = asYetUnresolvedParameters.end();
    for (uparamIter = ustart; uparamIter != uend; ++uparamIter)
    {
      parameter = *uparamIter;

      if (!resolveParameter(parameter))
      {
        // save it for later, because it might use functions that haven't
        // been resolved yet.
        uretryParams.insert(parameter);
        somethingLeftToDo=true;
      }
      else
      {
        currentContextPtr_->resolvedParams_.insert(parameter);
        resolvedSomethingThisLoop=true;
      }
    }
    asYetUnresolvedParameters=uretryParams;
    uretryParams.clear();

    // Resolve any .GLOBAL_PARAM parameters in the current context and add to
    // the resolvedGlobalParams_.
    start = asYetUnresolvedGlobalParameters.begin();
    end = asYetUnresolvedGlobalParameters.end();
    for (paramIter = start; paramIter != end; ++paramIter)
    {
      parameter = *paramIter;
      if (DEBUG_IO)
      {
        Xyce::dout() << " CircuitContext::resolve Attempting to resolve global parameter " << parameter.uTag() <<std::endl;
      }

      if (!resolveParameter(parameter))
      {
        retryParams.push_back(parameter);
        somethingLeftToDo = true;
      }
      else
      {
        if (parameter.getType() == Xyce::Util::EXPR)
        {
          const std::vector<std::string> & nodes = parameter.getValue<Util::Expression>().getVoltageNodes();
          const std::vector<std::string> & instances = parameter.getValue<Util::Expression>().getDeviceCurrents();
          const std::vector<std::string> & variables = parameter.getValue<Util::Expression>().getVariables();
          const std::vector<std::string> & leads = parameter.getValue<Util::Expression>().getLeadCurrents();
          std::vector<std::string> specials;
          parameter.getValue<Util::Expression>().getSpecials(specials);

          if (!nodes.empty() || !instances.empty() || !leads.empty())
          {
            Report::UserError0 message;
            message << "The following are not allowed in global param expression: " << parameter.getValue<Util::Expression>().get_expression();
            
            // This should be caught earlier, but just in case it is checked here
            if (!nodes.empty())
            {
              message << std::endl << "node(s):";
              for (std::vector<std::string>::const_iterator s_i=nodes.begin() ; s_i!=nodes.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
            if (!instances.empty())
            {
              message << std::endl << "instance(s): ";
              for (std::vector<std::string>::const_iterator s_i=instances.begin() ; s_i!=instances.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
            if (!leads.empty())
            {
              message << std::endl << "lead(s): ";
              //for (std::vector<std::string>::iterator s_i=leads.begin() ; s_i!=leads.end() ; ++s_i)
              for (std::vector<std::string>::const_iterator s_i=leads.begin() ; s_i!=leads.end() ; ++s_i)
              {
                message << " " << *s_i;
              }
            }
          }

          if (!variables.empty())
          {
            // If variables are found, they must be previously defined global params
            //for (std::vector<std::string>::iterator s_i=variables.begin() ; s_i!=variables.end() ; ++s_i)
            for (std::vector<std::string>::const_iterator s_i=variables.begin() ; s_i!=variables.end() ; ++s_i)
            {
              if (findParameter(currentContextPtr_->resolvedGlobalParams_.begin(), currentContextPtr_->resolvedGlobalParams_.end(), *s_i) == NULL)
              {
                Report::UserError0() << "Unknown parameter (" << *s_i << ") found in global param expression: "
                                     << parameter.getValue<Util::Expression>().get_expression();
              }
            }
          }
          if (!specials.empty())
          {
            // This block is "conservative". TIME, TEMP and VT now work in global parameters that use
            // expressions.  If additional "specials" are added to Xyce then this check ensures
            // that they must also be enabled for use in expressions in global parameters.
            bool badSpecial = false;
            std::string message;
            for (std::vector<std::string>::iterator s_i=specials.begin() ; s_i!=specials.end() ; ++s_i)
            {
              if ( !((*s_i == "TIME") || (*s_i == "TEMP") || (*s_i == "VT") || (*s_i == "FREQ")) )
              {
                message += " " + *s_i;
                badSpecial = true;
              }
              if (badSpecial)
              {
                Report::UserError0() << "Unknown special var(s):" << message;
              }
            }
          }
        }
        currentContextPtr_->resolvedGlobalParams_.push_back(parameter);
        resolvedSomethingThisLoop=true;
      }
    }
    asYetUnresolvedGlobalParameters=retryParams;
    retryParams.clear();

    // Resolve functions in the current context.
    // ERK.  new expression code to follow.
    int numFunctions = asYetUnresolvedFunctions.size();
    for ( i = 0; i < numFunctions; ++i )
    {
      // Define the Util::Param for the function prototype (the function name
      // and arguments) and function body. Pass to the circuit for resolution.
      std::string functionName(asYetUnresolvedFunctions[i].functionName);
      std::string functionNameAndArgs(asYetUnresolvedFunctions[i].functionNameAndArgs);
      std::string functionBody(asYetUnresolvedFunctions[i].functionBody);

      std::vector<std::string> functionArgs =
        asYetUnresolvedFunctions[i].functionArgs;
      Util::Param functionParameter(functionNameAndArgs, functionBody);
      // this is only an error if we have no parameters left to resolve
      if (!resolveParameterThatIsAdotFunc(functionParameter, functionArgs))
      {
        retryFunctions.push_back(asYetUnresolvedFunctions[i]);
        somethingLeftToDo=true;
      }
      else
      {
        currentContextPtr_->resolvedFunctions_[functionName] = functionParameter;
        resolvedSomethingThisLoop=true;
      }
    }
    asYetUnresolvedFunctions=retryFunctions;
    retryFunctions.clear();
  }

  if (somethingLeftToDo)    // we failed to resolve everything, and the last loop did nothing.
  {
    for (Util::ParamList::iterator it = asYetUnresolvedSubcircuitParameters.begin(); it != asYetUnresolvedSubcircuitParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve .subckt parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (Util::UParamList::iterator it = asYetUnresolvedParameters.begin(); it != asYetUnresolvedParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (Util::ParamList::iterator it = asYetUnresolvedGlobalParameters.begin(); it != asYetUnresolvedGlobalParameters.end(); ++it)
    {
      Report::UserError0() << "Unable to resolve global parameter " << (*it).uTag() << " found in .PARAM statement";
    }

    for (std::vector<FunctionBlock>::iterator func_it = asYetUnresolvedFunctions.begin(); func_it != asYetUnresolvedFunctions.end(); ++func_it)
    {
      const std::string &functionNameAndArgs = (*func_it).functionNameAndArgs;
      const std::string &functionBody = (*func_it).functionBody;
      const std::vector<std::string> &functionArgs = (*func_it).functionArgs;

      Util::Param functionParameter(functionNameAndArgs, functionBody);
      if (!resolveParameterThatIsAdotFunc(functionParameter, (*func_it).functionArgs))
      {
        Report::UserError0() << functionNameAndArgs << " contains an undefined parameter or function.";
      }
      else
      {

        Util::Expression functionBodyExpression(expressionGroup_, functionParameter.stringValue());
        const std::vector<std::string> & strings = functionBodyExpression.getUnresolvedParams();
        for (std::vector<std::string>::const_iterator it = strings.begin(); it != strings.end(); ++it)
        {
          if (find(functionArgs.begin(), functionArgs.end(), *it) == functionArgs.end()
              && findParameter(currentContextPtr_->resolvedGlobalParams_.begin(), currentContextPtr_->resolvedGlobalParams_.end(), *it) == NULL)
          {
            Report::UserError0() << functionNameAndArgs << " contains unknown parameter " << (*it);
          }
        }
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::setContext
// Purpose        : Set the current context to that corresponding to the given
//                  subcircuit name. Save the previous context on the stack for
//                  later retrieval. If the given subcircuit name is not found,
//                  declare an error and abort.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
bool CircuitContext::setContext(
  const std::string &                   subcircuitName,
  const std::string &                   subcircuitPrefixIn,
  const std::vector<std::string> &      instanceNodes,
  CircuitContext *                      previousContext ) const
{
  bool success = false;

  unordered_map< std::string, CircuitContext* >::iterator ccIter =
    currentContextPtr_->circuitContextTable_.find( subcircuitName );

  if ( ccIter != currentContextPtr_->circuitContextTable_.end() )
  {
    // Save current context and reset the current context pointer.
    if (previousContext == NULL)
      contextList_.push_front(currentContextPtr_);
    else
      contextList_.push_front(previousContext);

    currentContextPtr_ = ccIter->second;

    // If subcircuitPrefix and instanceNodes were given, set the prefix
    // in the new context, and build the node map.
    currentContextPtr_->nodeMap_.clear();
    currentContextPtr_->subcircuitPrefix_ = subcircuitPrefixIn;

    if (!instanceNodes.empty())
    {
      std::vector<std::string>::const_iterator nodeIter = instanceNodes.begin();
      for (unsigned int i = 0; i < currentContextPtr_->nodeList_.size() && nodeIter != instanceNodes.end(); ++i, ++nodeIter)
      {
        if (currentContextPtr_->nodeMap_.find(currentContextPtr_->nodeList_[i]) !=
            currentContextPtr_->nodeMap_.end())
        {
          if (currentContextPtr_->nodeMap_[currentContextPtr_->nodeList_[i]] !=
              *nodeIter)
          {
            Report::UserError0()
              <<  "Duplicate nodes in .subckt " << subcircuitName << " point to different nodes in X line invocation";
            return false;
          }
        }
        else
        {
          if ((currentContextPtr_->nodeList_[i].size() >= 2 &&
               currentContextPtr_->nodeList_[i].substr(0,2) == "$G") ||
              globalNode(currentContextPtr_->nodeList_[i]))
          {
            if (currentContextPtr_->nodeList_[i] != *nodeIter)
            {
              Report::UserError0()
                << "Global node in subcircuit invocation must match same name in .subckt";
              return false;
            }
          }
          currentContextPtr_->nodeMap_[currentContextPtr_->nodeList_[i]] = *nodeIter;
        }
      }
    }

    return true;
  }

  // The context with the given name was not found in the current context,
  // recursively search parent contexts.
  if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    if (previousContext == NULL)
    {
      previousContext = currentContextPtr_;
    }

    currentContextPtr_ = currentContextPtr_->parentContextPtr_;
    success = setContext(subcircuitName, subcircuitPrefixIn, instanceNodes, previousContext);
  }

  // if (!success)
  // {
  //   Report::UserError0() << "Global node in subcircuit invocation must match same name in .subckt";
  // }

  return success;
}

void CircuitContext::setContext(CircuitContext* context) const
{
  contextList_.push_front(currentContextPtr_);
  currentContextPtr_ = context;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::restorePreviousContext
// Purpose        : Reset the context the context prior to the last invocation
//                  of setContext.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/07/2003
//----------------------------------------------------------------------------
void CircuitContext::restorePreviousContext() const
{
  // Restore the previous context from the list of contexts unless
  // the list is empty.
  if (!contextList_.empty())
  {
    currentContextPtr_ = contextList_.front();
    contextList_.pop_front();
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::globalNode
// Purpose        : 
// Special Notes  :
// Scope          :
// Creator        : Dave Shirley, PSSI
// Creation Date  : 10/20/2009
//----------------------------------------------------------------------------
bool CircuitContext::globalNode (const std::string &nodeName) const
{
  bool stat;

  if (!(currentContextPtr_->parentContextPtr_))
  {
    if (globalNodes_.find(nodeName) == globalNodes_.end())
      return false;
    else
      return true;
  }

  setContext(currentContextPtr_->parentContextPtr_);
  stat = globalNode ( nodeName );
  restorePreviousContext();

  return stat;
}

//----------------------------------------------------------------------------
// Function       : debugResolveParameterOutput1
// Purpose        : Helper debug output function so that 
//                  CircuitContext::resolveParameter is more readable
// Special Notes  :
// Scope          : public
// Creator        : Eric Keiter
// Creation Date  : 04/04/2021
//----------------------------------------------------------------------------
void debugResolveParameterOutput1(
      const std::vector<std::string> & nodes,
      const std::vector<std::string> & instances,
      const std::vector<std::string> & variables,
      const std::vector<std::string> & leads,
      const std::vector<std::string> & specials,
      bool isRandom
    )
{
  Xyce::dout()<<"CircuitContext::resolveParameter:  nodes, instances, leads, variables or specials not empty, or this has a random operator such as AGAUSS."<<std::endl;
  if (!nodes.empty()) { Xyce::dout()<<" Nodes: "<<std::endl;
    for (unsigned int ii=0; ii<nodes.size(); ++ii) {Xyce::dout()<<ii<<" : "<<nodes[ii]<<std::endl;}
  }
  if (!instances.empty()) { Xyce::dout()<<" Instances: "<<std::endl;
    for (unsigned int ii=0; ii<instances.size(); ++ii) {Xyce::dout()<<ii<<" : "<<instances[ii]<<std::endl;}
  }
  if (!leads.empty()) { Xyce::dout()<<" Leads: "<<std::endl;
    for (unsigned int ii=0; ii<leads.size(); ++ii) {Xyce::dout()<<ii<<" : "<<leads[ii]<<std::endl;}
  }
  if (!variables.empty()) { Xyce::dout()<<" Variables: "<<std::endl;
    for (unsigned int ii=0; ii<variables.size(); ++ii) {Xyce::dout()<<ii<<" : "<<variables[ii]<<std::endl;}
  }
  if (!specials.empty()) { Xyce::dout()<<" Specials: "<<std::endl;
    for (unsigned int ii=0; ii<specials.size(); ++ii) {Xyce::dout()<<ii<<" : "<<specials[ii]<<std::endl;}
  }

  if (isRandom) {Xyce::dout()<<" Depends on a random operator"<<std::endl;}
}

//----------------------------------------------------------------------------
// Function       : debugResolveParameterOutput2
// Purpose        : Helper debug output function so that 
//                  CircuitContext::resolveParameter is more readable
// Special Notes  :
// Scope          : public
// Creator        : Eric Keiter
// Creation Date  : 04/04/2021
//----------------------------------------------------------------------------
void debugResolveParameterOutput2 (Util::Param& parameter) 
{
    Xyce::dout() << " after setting the parameter " << parameter.uTag() 
      << ", its type is " << parameter.getType() << std::endl;
  switch (parameter.getType()) 
  {
    case Xyce::Util::STR:
      Xyce::dout()<<" "<<parameter.uTag()<<" is STR type; value =  "<< parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::DBLE:
      Xyce::dout()<<" "<<parameter.uTag()<<" is DBLE type; value =  "<< parameter.getImmutableValue<double>()<<std::endl;
      break;
    case Xyce::Util::EXPR:
      Xyce::dout()<<" "<<parameter.uTag()<<" is EXPR type; value =  "<<parameter.getValue<Util::Expression>().get_expression()<<std::endl;
      break;
    case Xyce::Util::BOOL:
      Xyce::dout()<<" "<<parameter.uTag()<<" is BOOL type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::STR_VEC:
      Xyce::dout()<<" "<<parameter.uTag()<<" is STR_VEC type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::INT_VEC:
      Xyce::dout()<<" "<<parameter.uTag()<<" is INT_VEC type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::DBLE_VEC:
      Xyce::dout()<<" "<<parameter.uTag()<<" is DBLE_VEC type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::DBLE_VEC_IND:
      Xyce::dout()<<" "<<parameter.uTag()<<" is DBLE_VEC_IND type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    case Xyce::Util::COMPOSITE:
      Xyce::dout()<<" "<<parameter.uTag()<<" is COMPOSITE type; value =  "<<parameter.stringValue()<<std::endl;
      break;
    default:
      Xyce::dout()<<" "<<parameter.uTag()<<" is default type (whatever that is): "<<parameter.stringValue()<<std::endl;
  }

  Xyce::dout()<<" and its value is ";
  switch (parameter.getType()) {
    case Xyce::Util::STR:
      Xyce::dout()<<parameter.stringValue();
      break;
    case Xyce::Util::DBLE:
      Xyce::dout()<<parameter.getImmutableValue<double>();
      break;
    case Xyce::Util::EXPR:
      Xyce::dout()<<parameter.getValue<Util::Expression>().get_expression();
      break;
    default:
      Xyce::dout()<<parameter.stringValue();
  }
  Xyce::dout()<<std::endl;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveParameter
// Purpose        : Simpler version, excluding exception strings.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters/Eric Keiter
// Creation Date  : 02/10/2003; 4/11/2020
//----------------------------------------------------------------------------
bool CircuitContext::resolveParameter(Util::Param& parameter) const
{
  if (hasExpressionTag(parameter) || parameter.hasExpressionValue() )
  {
    if (DEBUG_IO)
    {
      Xyce::dout() << "CircuitContext::resolveParameter parameter " << parameter.uTag()
                   << " has expression value " << std::endl;
    }

    // Extract the expression from the parameter value by stripping off
    // the enclosing braces.  Only strip if it's there!
    std::string expressionString;
    if (parameter.stringValue()[0] == '{')
    {
      expressionString = parameter.stringValue().substr(1, parameter.stringValue().size()-2);
    }
    else
    {
      expressionString = parameter.stringValue();
    }

    // Parse the expression:
    Util::Expression expression(expressionGroup_, expressionString);

    if (!expression.parsed()) { return false; }

    // Resolve the strings in the expression. Unresolved strings
    // may be parameters defined in a .PARAM statement or global
    // parameters defined in .GLOBAL_PARAM statement.
    bool stringsResolved = resolveStrings(expression);

    // Resolve functions in the expression.
    bool functionsResolved = resolveFunctions(expression);

#if 0
    parameter.setVal(expression); 
    if ( !(expression.getLeadCurrents().empty()) ) { return false; }

    if (DEBUG_IO) 
    {
       Xyce::dout() << "CircuitContext::resolveParameter: right before returns " << std::endl;
       debugResolveParameterOutput2(parameter); 
    }
    return stringsResolved && functionsResolved;
#else
    if ( !(expression.getLeadCurrents().empty()) )
    {
      parameter.setVal(expression);
      return false;
    }

    if (stringsResolved && functionsResolved)
    {
      // Check the expression for nodes or instance (probably a B-source
      // expression, in which case the parameter value should be
      // expressionString. Also check for "specials", the only special
      // allowed is "time" for time dependent parameters.
      const std::vector<std::string> & nodes = expression.getVoltageNodes();
      const std::vector<std::string> & instances = expression.getDeviceCurrents();
      const std::vector<std::string> & variables = expression.getVariables();
      const std::vector<std::string> & leads = expression.getLeadCurrents();
      std::vector<std::string> specials;
      expression.getSpecials(specials);
      bool isRandom = expression.isRandomDependent();

      if (!nodes.empty() || !instances.empty() || !leads.empty() ||
          !variables.empty() || !specials.empty() || isRandom)
      {
        if (DEBUG_IO) { debugResolveParameterOutput1( nodes, instances, variables, leads, specials, isRandom ); }

        parameter.setVal(expression);

        if (DEBUG_IO) 
        { 
           Xyce::dout() << "CircuitContext::resolveParameter: After all expression handling, get_expression returns "
               << expression.get_expression() << std::endl;
          debugResolveParameterOutput2(parameter); 
        }
      }
      else
      {
        // Reset the parameter value to the value of the expression.
        double value(0.0);
        expression.evaluateFunction ( value );
        parameter.setVal( value );
        // we have resolved the context so set it and the constant value to 
        // make later look ups easier.
        // parameter.addOp(Util::CONSTANT, new IO::ConstantOp(parameter.tag(), value));
        if (DEBUG_IO)
        {
          Xyce::dout() << " CircuitContext::resolveParameter --  Resetting parameter value from " << expressionString << " to " << value << " because it is resolved and not a function" << std::endl;
        }
      }
    }
    else
    {
      // Reset the parameter value to the value of the expression with
      // as much resolution as could be achieved.
      parameter.setVal(expression);
    }

    if (DEBUG_IO) 
    {
       Xyce::dout() << "CircuitContext::resolveParameter: right before returns " << std::endl;
       debugResolveParameterOutput2(parameter); 
    }

    return stringsResolved && functionsResolved;
#endif
  }

  // Handle quoted string parameters e.g. "some text" by removing quotes.
  resolveQuote(parameter);
  resolveTableFileType(parameter);
  
  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveParameterThatIsAdotFunc
// Purpose        : Parameters that are .funcs need special treatment
// Special Notes  : ERK.  Maybe make the return boolean more meaningful, 
//                  so downstream error checking is less necessary.
// Scope          : public
// Creator        : Eric Keiter
// Creation Date  : 04/11/2020
//----------------------------------------------------------------------------
bool CircuitContext::resolveParameterThatIsAdotFunc(Util::Param& parameter,
                                      std::vector<std::string> funcArgs) const
{
  if (hasExpressionTag(parameter) || parameter.hasExpressionValue() )
  {
    if (DEBUG_IO)
    {
      Xyce::dout() << "CircuitContext::resolveParameterThatIsAdotFunc parameter " << parameter.uTag()
                   << " has expression value " << std::endl;
    }

    // Extract the expression from the parameter value by stripping off
    // the enclosing braces.  Only strip if it's there!
    std::string expressionString;
    if (parameter.stringValue()[0] == '{')
    {
      expressionString = parameter.stringValue().substr(1, parameter.stringValue().size()-2);
    }
    else
    {
      expressionString = parameter.stringValue();
    }

    // Parse the expression:
    Util::Expression expression(expressionGroup_, expressionString,funcArgs);

    if (!expression.parsed()) { return false; }

    // Resolve the strings in the expression, which may be 
    // parameters defined by .PARAM or .GLOBAL_PARAM.  
    //
    // They may also be function arguments if the expression is the
    // body of a function defined by .FUNC.
    bool stringsResolved = resolveStrings(expression, funcArgs);

    // Resolve functions in the expression.
    bool functionsResolved = resolveFunctions(expression);

    parameter.setVal(expression); 
    if ( !(expression.getLeadCurrents().empty()) ) { return false; }

    if (DEBUG_IO)
    {
      Xyce::dout() << "CircuitContext::resolveParameterThatIsAdotFunc: right before returns " << std::endl;
       debugResolveParameterOutput2(parameter); 
    }
    return stringsResolved && functionsResolved;
  }
 
  // Handle quoted string parameters e.g. "some text" by removing quotes.
  resolveQuote(parameter);
  resolveTableFileType(parameter);
  
  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveStrings
// Purpose        : Determine if expression has any unresolved strings
//                  and resolve appropriately. Return true if all strings are
//                  resolved otherwise return false.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool CircuitContext::resolveStrings( Util::Expression & expression,
                                     std::vector<std::string> exceptionStrings) const
{
  // Strings in the expression must be previously resolved parameters
  // that appear in paramList else there is an error.
  bool unresolvedStrings = false;
  // Normally "strings" would be a reference.  But here it has to be a copy 
  // because it loops over them and accesses them in the function below.
  // If it is just a reference, then the vector keeps getting smaller each 
  // time a param/string is resolved, and it results in memory access errors.
  const std::vector<std::string> strings = expression.getUnresolvedParams(); 

  if ( !(strings.empty()) )
  {
    // If the expression is resolvable, each string in the current expression
    // must appear as a resolved parameter in netlistParameters. Get the value
    // if it appears there.
    ExtendedString parameterName("");
    int numStrings = strings.size();

    if (DEBUG_IO)
    {
      Xyce::dout() <<" CircuitContext::resolveStrings numStrings = " << numStrings << std::endl;
    }

    for (int i = 0; i < numStrings; ++i)
    {
      if (DEBUG_IO)
      {
        Xyce::dout() <<" CircuitContext::resolveStrings resolving " << strings[i] << std::endl;
      }

      // Skip the current string if it is in exception strings. This prevents
      // a function argument from being improperly resolved when there is
      // a parameter in a .param statement with the same name as the function
      // argument.
      if (!exceptionStrings.empty())
      {
        std::vector<std::string>::iterator stringIter = find(exceptionStrings.begin(),
                                                             exceptionStrings.end(),
                                                             strings[i]);
        if (stringIter != exceptionStrings.end()) { continue; }
      }

      // Look for the string in netlistParameters.
      parameterName = strings[i];
      parameterName.toUpper();

      Util::Param expressionParameter(parameterName, "");
      bool parameterFound = getResolvedParameter(expressionParameter);
      if (parameterFound)
      {
        if (DEBUG_IO)
        {
          Xyce::dout() <<" CircuitContext::resolveStrings string " << strings[i] << " is a resolved parameter " << expressionParameter.uTag() << " with type "
                       << expressionParameter.getType() << " and value ";
          switch (expressionParameter.getType()) {
            case Xyce::Util::STR:
              Xyce::dout() << expressionParameter.stringValue();
              break;
            case Xyce::Util::DBLE:
              Xyce::dout() << expressionParameter.getImmutableValue<double>();
              break;
            case Xyce::Util::EXPR:
              Xyce::dout() << "EXPR("<<expressionParameter.getValue<Util::Expression>().get_expression()<< ")";
              break;
            default:
              Xyce::dout() << expressionParameter.stringValue();
          }
          Xyce::dout() << std::endl;
        }

        if ( expressionParameter.getType() == Xyce::Util::STR ||
             expressionParameter.getType() == Xyce::Util::DBLE )
        {
          enumParamType paramType=DOT_PARAM;
          if (!expression.make_constant(strings[i], expressionParameter.getImmutableValue<double>(),paramType))
          {
            Report::UserWarning0() << "Problem converting parameter " << parameterName << " to its value.";
          }
        }
        else if (expressionParameter.getType() == Xyce::Util::EXPR)
        {
          std::string expressionString=expression.get_expression();

          // ERK.  Add an error test for nodes that cannot be attached below.  Something like:
          //
          //  Report::UserWarning0() << "Problem inserting expression " << expressionParameter.getValue<Util::Expression>().get_expression()
          //                         << " as substitute for " << parameterName << " in expression " << expressionString;
          const std::vector<std::string> & variables = expressionParameter.getValue<Util::Expression>().getVariables ();

          enumParamType paramType=DOT_PARAM;
          if (variables.empty()) paramType=DOT_PARAM;
          else paramType=SUBCKT_ARG_PARAM;
          expression.attachParameterNode(strings[i], expressionParameter.getValue<Util::Expression>(),paramType); 
        }
      }
      else
      {
        parameterFound = getResolvedGlobalParameter(expressionParameter);
        if (parameterFound)
        {
          if (DEBUG_IO)
          {
            Xyce::dout() <<" CircuitContext::resolveStrings string " << strings[i] << " is a resolved global parameter " << expressionParameter.uTag() << " with type "
                         << expressionParameter.getType() << " and value ";
            switch (expressionParameter.getType()) {
              case Xyce::Util::STR:
                Xyce::dout() << expressionParameter.stringValue();
                break;
              case Xyce::Util::DBLE:
                Xyce::dout() << expressionParameter.getImmutableValue<double>();
                break;
              case Xyce::Util::EXPR:
                Xyce::dout() << "EXPR("<<expressionParameter.getValue<Util::Expression>().get_expression()<< ")";
                break;
              default:
                Xyce::dout() << expressionParameter.stringValue();
            }
            Xyce::dout() << std::endl;
          }

          // ERK right thing to do, but won't work until set_vars/order_vars, etc are removed, 
          // and a better group is set up.
          if (expressionParameter.getType() == Xyce::Util::EXPR)
          {
            Util::Expression & expToBeAttached = expressionParameter.getValue<Util::Expression>();
            expression.attachParameterNode(strings[i], expToBeAttached);
          }
          else
          {
            if (!expression.make_var(strings[i])) 
            {
              Report::UserWarning0() << "Problem converting parameter " << parameterName <<" to its value";
            }
          }
        }
        else
        {
          if (Util::isBool(strings[i]))
          {
            bool stat = false;
            enumParamType paramType=DOT_PARAM;
            if (Util::Bval(strings[i]))
              stat = expression.make_constant(strings[i], static_cast<double>(1),paramType);
            else
              stat = expression.make_constant(strings[i], static_cast<double>(0),paramType);
            if (!stat)
            {
              Report::UserWarning0() << "Problem converting parameter " << parameterName << " to its value";
            }
          }
          else
          {
            unresolvedStrings = true;
          }
        }
      }
    }
  }
  return !unresolvedStrings;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::resolveFunctions
// Purpose        : Determine if expression has any unresolved functions
//                  and resolve appropriately. Return true if all functions
//                  are resolved otherwise return false.
// Special Notes  : 
// Scope          :
// Creator        : Lon Waters, Eric Keiter
// Creation Date  : 02/11/2003, 2020
//----------------------------------------------------------------------------
bool CircuitContext::resolveFunctions(Util::Expression & expression) const
{
  bool unresolvedFunctions = false;
  std::vector<std::string> funcNames;
  expression.getFuncNames(funcNames);
  if ( funcNames.size() > 0 )
  {
    for (int ii = 0; ii < funcNames.size(); ++ii)
    {
      // Look for the function in resolvedFunctions_.
      Util::Param functionParameter(funcNames[ii], "");
      bool functionfound = getResolvedFunction(functionParameter);
      if (functionfound)
      { // pull out the RHS expression and attach
        if(functionParameter.getType() == Xyce::Util::EXPR)
        {
          Util::Expression & expToBeAttached = functionParameter.getValue<Util::Expression>();
          expression.attachFunctionNode(funcNames[ii], expToBeAttached);
        }
        else 
        { 
          Report::DevelFatal()<< "functionParameter " <<  funcNames[ii] << " is not EXPR type!!!"; 
        }
      }
      else
      {
        unresolvedFunctions = true;
      }
    }
  }
  return !unresolvedFunctions;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getResolvedParameter
// Purpose        : Look for a parameter (set in the netlist by .param) 
//                  in the current context's set of resolved parameters. 
//                  Check the current context and recurively check parent 
//                  contexts.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/11/2003
//----------------------------------------------------------------------------
bool CircuitContext::getResolvedParameter(Util::Param & parameter) const
{
  bool success = false;

  Util::UParamList::const_iterator rP_iter = currentContextPtr_->resolvedParams_.find( parameter );
 
  if ( rP_iter != currentContextPtr_->resolvedParams_.end() )
  {
    // Found a parameter with given name, set the value
    // of parameter and return.
    parameter.setVal(*rP_iter);
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedParameter(parameter);
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getResolvedGlobalParameter
// Purpose        : Look for a parameter (set in the netlist by .global_param)
//                  in the current context's set of resolved parameters. 
//                  Check the current context and recurively check parent 
//                  contexts.
// Special Notes  :
// Scope          : public
// Creator        : Dave Shirley, PSSI
// Creation Date  : 11/17/2005
//----------------------------------------------------------------------------
bool CircuitContext::getResolvedGlobalParameter(Util::Param & parameter) const
{
  bool success = false;

  const Util::Param* parameterPtr = findParameter(currentContextPtr_->resolvedGlobalParams_.begin(), currentContextPtr_->resolvedGlobalParams_.end(), parameter.tag());
  if (parameterPtr != NULL)
  {
    // Found a parameter with given name, set the value
    // of parameter and return.
    parameter.setVal(*parameterPtr);
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedGlobalParameter(parameter);
    restorePreviousContext();
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::getResolvedFunction
// Purpose       : Look for a function (set in the netlist by .func) 
//                 in resolvedFunctions_.  Check current context and recursively 
//                 check each parent context.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/27/2001
//-----------------------------------------------------------------------------
bool CircuitContext::getResolvedFunction(Util::Param & parameter) const
{
  bool success = false;

  Util::ParamMap::const_iterator it = currentContextPtr_->resolvedFunctions_.find(parameter.uTag());

  if (it != currentContextPtr_->resolvedFunctions_.end())
  {
    parameter = it->second;
    success = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = getResolvedFunction(parameter);
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::findModel
// Purpose        : Search the models in the current context for the model of
//                  the given name. If it is not found, recursively
//                  search each parent context. Return a pointer to
//                  the parameterBlock for the model if it is found,
//                  otherwise return NULL. Also, if the model is found,
//                  construct the appropriate model prefix.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/12/2003
//----------------------------------------------------------------------------
bool CircuitContext::findModel(
  const std::string &           modelName,
  ParameterBlock* &             modelPtr,
  std::string &                 modelPrefix) const
{
  bool success = false;

  ModelMap::const_iterator modelIter = currentContextPtr_->models_.find(modelName);
  if (modelIter != currentContextPtr_->models_.end())
  {
    modelPtr = modelIter->second;
    if (modelPtr->hasExpressionValuedParams())
    {
      modelPrefix = currentContextPtr_->getPrefix();
    }
    else
    {
      std::string prefix = currentContextPtr_->getCurrentContextName();
      if (prefix == "")
        modelPrefix = "";
      else
        modelPrefix = prefix + ":";
    }

    return true;
  }

  // The model was not found in the current circuit's model list,
  // recursively search the parent circuit's model list.
  if ( currentContextPtr_->parentContextPtr_ != NULL )
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = findModel( modelName, modelPtr, modelPrefix );
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::findModel
// Purpose        : Overloaded version of findModel for cases when the model
//                  prefix is not needed.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
bool CircuitContext::findModel(
  const std::string &           modelName,
  ParameterBlock* &        modelPtr) const
{
  bool success;
  std::string temp;
  success = findModel(modelName, modelPtr, temp);

  return success;
}


//----------------------------------------------------------------------------
// Function       : FindPrefix
// Purpose        : Find the first model in the ModelMap whose name starts 
//                  with the specified prefix.  The ModelMap is an ordered map,
//                  so other model names which have this prefix will immediately 
//                  follow it, in order in the container.
// Special Notes  : 
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//----------------------------------------------------------------------------
ModelMap::const_iterator FindBinningName(const ModelMap& modelMap, const std::string & binningModelName) 
{
  std::string tmpName(binningModelName);
  Util::toUpper(tmpName);
  ModelMap::const_iterator i = modelMap.lower_bound(tmpName);
  if (i != modelMap.end()) 
  {
    const std::string & key = i->first;
    if (key.compare(0, tmpName.size(), tmpName) == 0) // Really a prefix?
    {
      return i;
    }
  }
  return modelMap.end();
}

//----------------------------------------------------------------------------
// Function       : is_equal
// Purpose        :
// Special Notes  : Added to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//----------------------------------------------------------------------------
bool is_equal(double result, double expectedResult)
{
  return fabs(result - expectedResult) < 1e-15;
}

//----------------------------------------------------------------------------
// Function       : InBinRange
// Purpose        :
// Special Notes  : Added to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//----------------------------------------------------------------------------
bool InBinRange(double value, double min, double max)
{
  // the standard binning rule is: min <= value < max 
  return is_equal(value, min) || (min < value && value < max);
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::fullyResolveParam
// Purpose        : Goes further than the "resolveParameter" function.
// Special Notes  : Added to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/16/2018
//----------------------------------------------------------------------------
bool CircuitContext::fullyResolveParam(Device::Param & param, double & value) const
{
  bool successfullyResolved = false;

  if (Util::isValue(param.stringValue()))
  {
    value = Util::Value(param.stringValue());
    successfullyResolved=true;
  }
  else
  {
    if (param.hasExpressionValue() || param.isQuoted() ||
        param.isTableFileTypeQuoted() || param.isStringTypeQuoted())
    {
      if(resolveParameter(param))
      {
        // ERK. This is to make this (binning) work with L,W params being set via 
        // global_params.
        // 
        // if "param" is L or W, and it is an expression which depends on a global, then 
        // the "getImmutableValue" function will produce an error.  But, if the the globals
        // are NOT used to change values (such as via .STEP), then there is nothing wrong
        // with this.  So, I have decided this use case should issue a warning,
        // rather than have a fatal error.
        if (param.hasExpressionValue())
        {
          Util::Expression &expression = const_cast<Util::Expression &>(param.getValue<Util::Expression>());
          double dVal;
          expression.evaluateFunction (dVal);
          value = dVal;
        }
        else
        { // if not an expression, do things the old-fashioned way
          value = param.getImmutableValue<double>();
        }
        successfullyResolved=true;
      }
    }
  }

  return successfullyResolved;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::findBinnedModel
// Purpose        : 
// Special Notes  : Special case to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//----------------------------------------------------------------------------
bool CircuitContext::findBinnedModel(
  const std::string &           modelName,
  ParameterBlock* &             modelPtr,
  std::string &                 modelPrefix,
  const bool LWfound,
  const bool LNFINfound,
  const double L,
  const double W,
  const double NFIN, std::string & binNumber) const
{
  bool success = false;

  ModelMap::const_iterator modelIter = FindBinningName(currentContextPtr_->models_, modelName);
  ModelMap::const_iterator end = currentContextPtr_->models_.end();

  std::string tmpname=modelName; Util::toUpper(tmpname);
  for(;modelIter!=end;modelIter++)
  {
    const std::string& key = modelIter->first;
    bool done=false;
    std::string tmpName(modelName); Util::toUpper(tmpName);
    if (key.compare(0, tmpName.size(), tmpName) == 0) // check if prefix still matches, and then check if the correct bin
    {
      Device::Param parLmax(std::string("LMAX"), "");
      Device::Param parLmin(std::string("LMIN"), "");

      ParameterBlock & tmpModelParamBlock = *(modelIter->second);
      Device::Param * Lmaxptr = tmpModelParamBlock.findParameter(parLmax);
      Device::Param * Lminptr = tmpModelParamBlock.findParameter(parLmin);

      // for the case of L & W geometric parameters 
      if (LWfound)
      {
        Device::Param parWmax(std::string("WMAX"), "");
        Device::Param parWmin(std::string("WMIN"), "");

        Device::Param * Wmaxptr = tmpModelParamBlock.findParameter(parWmax);
        Device::Param * Wminptr = tmpModelParamBlock.findParameter(parWmin);

        if (Lmaxptr != NULL && Lminptr != NULL && Wmaxptr != NULL && Wminptr != NULL) 
        { 
          double Lmax = 1.0, Lmin = 0.0, Wmax = 1.0, Wmin = 0.0;
          bool resolvedLmax=false, resolvedLmin=false, resolvedWmax=false, resolvedWmin=false;

          resolvedLmax = fullyResolveParam(*Lmaxptr,Lmax);
          resolvedLmin = fullyResolveParam(*Lminptr,Lmin);
          resolvedWmax = fullyResolveParam(*Wmaxptr,Wmax);
          resolvedWmin = fullyResolveParam(*Wminptr,Wmin);

          if (resolvedLmax && resolvedLmin && resolvedWmax && resolvedWmin)
          {
            if (InBinRange(L, Lmin, Lmax) && InBinRange(W, Wmin, Wmax)) 
            {
              binNumber = key.substr( (tmpName.size()+1), (key.size()-1) );
              done=true;
            }
          }
        }
      }


      // for the case of L & NFIN geometric parameters 
      if (LNFINfound)
      {
        Device::Param parNFINmax(std::string("NFINMAX"), "");
        Device::Param parNFINmin(std::string("NFINMIN"), "");

        Device::Param * NFINmaxptr = tmpModelParamBlock.findParameter(parNFINmax);
        Device::Param * NFINminptr = tmpModelParamBlock.findParameter(parNFINmin);

        if (Lmaxptr != NULL && Lminptr != NULL && NFINmaxptr != NULL && NFINminptr != NULL) 
        { 
          double Lmax = 1.0, Lmin = 0.0, NFINmax = 1.0, NFINmin = 0.0;
          bool resolvedLmax=false, resolvedLmin=false, resolvedNFINmax=false, resolvedNFINmin=false;

          resolvedLmax = fullyResolveParam(*Lmaxptr,Lmax);
          resolvedLmin = fullyResolveParam(*Lminptr,Lmin);
          resolvedNFINmax = fullyResolveParam(*NFINmaxptr,NFINmax);
          resolvedNFINmin = fullyResolveParam(*NFINminptr,NFINmin);

          if (resolvedLmax && resolvedLmin && resolvedNFINmax && resolvedNFINmin)
          {
            if (InBinRange(L, Lmin, Lmax) && InBinRange(NFIN, NFINmin, NFINmax)) 
            {
              binNumber = key.substr( (tmpName.size()+1), (key.size()-1) );
              done=true;
            }
          }
        }
      }

    }
    else // if prefix no longer matches, we are done, and have probably failed.
    {
      done=true;
    }

    if (done) break;
  }

  if (modelIter != currentContextPtr_->models_.end())
  {
    modelPtr = modelIter->second;
    if (modelPtr->hasExpressionValuedParams())
    {
      modelPrefix = currentContextPtr_->getPrefix();
    }
    else
    {
      std::string prefix = currentContextPtr_->getCurrentContextName();
      if (prefix == "")
        modelPrefix = "";
      else
        modelPrefix = prefix + ":";
    }

    return true;
  }

  // The model was not found in the current circuit's model list,
  // recursively search the parent circuit's model list.
  if (currentContextPtr_->parentContextPtr_ != NULL )
  {
    setContext(currentContextPtr_->parentContextPtr_);
    success = findBinnedModel( modelName, modelPtr, modelPrefix, LWfound, LNFINfound, L, W, NFIN, binNumber);
    restorePreviousContext();
  }

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::findBinnedModel
// Purpose        : 
// Special Notes  : Special case for model binning
// Scope          : public
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//----------------------------------------------------------------------------
bool CircuitContext::findBinnedModel(
  const std::string &           modelName,
  ParameterBlock* &        modelPtr,
  const bool LWfound, const bool LNFINfound,
  const double L, const double W, const double NFIN, std::string & binNumber) const
{
  bool success;
  std::string temp;
  success = findBinnedModel(modelName, modelPtr, temp, LWfound, LNFINfound, L, W, NFIN, binNumber);

  return success;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::hasSubcircuitParams
// Purpose        : Check whether a subcircuit context is dependent on
//                  subcircuit parameters. These are parameters on the
//                  .subckt line identified by "params:" keyword. The result
//                  should be true if either the current subcircuit context
//                  or any subcircuit context in the hierarchy containing the
//                  current subcircuit context has subcircuit parameters.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/15/2003
//----------------------------------------------------------------------------
bool CircuitContext::hasSubcircuitParams()
{
  bool foundSubcircuitParams = false;

  if (!currentContextPtr_->subcircuitParameters_.empty())
  {
    foundSubcircuitParams = true;
  }
  else if (currentContextPtr_->parentContextPtr_ != NULL)
  {
    setContext(currentContextPtr_->parentContextPtr_);
    foundSubcircuitParams = hasSubcircuitParams();
    restorePreviousContext();
  }

  return foundSubcircuitParams;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::getTotalDeviceCount
// Purpose        : Calculate the total number of devices starting at current
//                  context and including all subcircuit instances.
// Special Notes  : This function also accumulates the used context list and 
//                  device type count map.
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/13/2003
//----------------------------------------------------------------------------
int CircuitContext::getTotalDeviceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;
  int linearCount = currentContextPtr_->linearDeviceCount_;

  // Determine the device count associated with each subcircuit instance.
  if (!currentContextPtr_->deviceCountDone_)
  {
    DeviceCountMap tmpDCMap = currentContextPtr_->localDeviceCountMap_;

    std::vector<std::string>::iterator instanceIter;
    std::vector<std::string>::iterator start = currentContextPtr_->subcktList_.begin();
    std::vector<std::string>::iterator end = currentContextPtr_->subcktList_.end();

    for (instanceIter = start; instanceIter != end; ++instanceIter)
    {
      bool result = setContext(*instanceIter);
      if (result)
      {
        // MUST perform getTotalDeviceCount before getTotalLinearDeviceCount.
        count += getTotalDeviceCount();
        linearCount += currentContextPtr_->getTotalLinearDeviceCount();

        // Add in devices to map for counting.
        DeviceCountMap::iterator lDCMb = currentContextPtr_->localDeviceCountMap_.begin();
        DeviceCountMap::iterator lDCMe = currentContextPtr_->localDeviceCountMap_.end();
        DeviceCountMap::iterator lDCMit = lDCMb;
        for ( ; lDCMit != lDCMe; ++lDCMit) 
        {
          if ( tmpDCMap[lDCMit->first] )
          {
            tmpDCMap[lDCMit->first] += lDCMit->second;
          }
          else
          {
            tmpDCMap[lDCMit->first] = lDCMit->second;
          }
        }
      }
      
      restorePreviousContext();

      if (result)
      {
        // Just push back the used subcircuit, multiple entries will be filtered out later.
        (currentContextPtr_->usedContextList_).push_back( *instanceIter );
      }  
    }

    // Propagate this list up to the parent, eventually the top circuit context will have
    // all the names of the used subcircuits.
    if (currentContextPtr_->parentContextPtr_ != NULL)
    {
      std::vector<std::string>::iterator ucIter;
      std::vector<std::string>::iterator uc_start = (currentContextPtr_->usedContextList_).begin();
      std::vector<std::string>::iterator uc_end = (currentContextPtr_->usedContextList_).end();

      for (ucIter = uc_start; ucIter!=uc_end; ++ucIter)
      {
        // Just push back the used subcircuit, multiple entries will be filtered out later.
        ((currentContextPtr_->parentContextPtr_)->usedContextList_).push_back( *ucIter );
      }
    }

    currentContextPtr_->deviceCount_ = count;
    currentContextPtr_->linearDeviceCount_ = linearCount;
    currentContextPtr_->localDeviceCountMap_ = tmpDCMap;
    currentContextPtr_->deviceCountDone_ = true;
  }
    
  return count;
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::pruneContexts
// Purpose        : Remove any context that is not specified in the passed in 
//                : list.  These were obtained through the computation of the
//                : total number of devices.
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 06/20/2014
//----------------------------------------------------------------------------
void CircuitContext::pruneContexts( std::vector<std::string>& keepContexts )
{
  unordered_map< std::string, CircuitContext * >::iterator itsc2, itsc = circuitContextTable_.begin();

  // Delete each context object in the table that is not on the provided list.
  while ( itsc != circuitContextTable_.end() )
  {
    itsc->second->pruneContexts( keepContexts );

    // Safely remove any contexts that were not used.
    if ( !std::binary_search( keepContexts.begin(), keepContexts.end(), itsc->first ) )
    {
      // Delete the pointer to the object and then remove it from the map.
      // Make sure that itsc is a valid iterator after the entry is removed from the map.
      itsc2 = itsc;
      itsc++;
      delete itsc2->second;
      circuitContextTable_.erase( itsc2 );
    }
    else
    {
      itsc++;
    }
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::totalMutualInductanceCount
// Purpose        : Calculate the total number of MIs starting at current
//                  context and including all subcircuit instances.
// Special Notes  :
// Scope          : public
// Creator        : Rob Hoekstra
// Creation Date  : 08/27/2004
//----------------------------------------------------------------------------
int CircuitContext::totalMutualInductanceCount()
{
  // Get device count for the current context.
  int count = currentContextPtr_->mutualInductances_.size();

  // Determine the device count associated with each subcircuit instance.
  std::vector<std::string>::iterator instanceIter;
  std::vector<std::string>::iterator start = currentContextPtr_->subcktList_.begin();
  std::vector<std::string>::iterator end = currentContextPtr_->subcktList_.end();

  for (instanceIter = start; instanceIter != end; ++instanceIter)
  {
    bool result = setContext(*instanceIter);
    if (result)
    {
      count += totalMutualInductanceCount();
    }

    restorePreviousContext();
  }

  return count;
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::MutualInductance
// Purpose       : Collect data from MI line "Ktrans L1 L2 " etc.
// Special Notes : This function just collects the data from a line starting
//                 with K for a linear or non-linear mutual inductor.
//                 it does not try to find the component inductors.
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
CircuitContext::MutualInductance::MutualInductance( DeviceBlock & device )
{
  int numParameters = device.getNumberOfInstanceParameters();
  Device::Param parameter;
  bool first = true;
  for ( int i = 0; i < numParameters; ++i )
  {
    parameter = device.getInstanceParameter(i);

    if ( parameter.tag() != "COUPLING" )
    {
      std::string inductorName (parameter.uTag());
      if( first )
      {
        firstInductor = inductorName;
        first = false;
      }

      inductors[inductorName] = "0.0";
    }
    else
      coupling = parameter.stringValue();
  }

  model = device.getModelName();

  name = device.getInstanceName().getEncodedName();
  sharedKey = 0;

}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::bundleMIs
// Purpose       : Convert all MIs into tokenized device lines
// Special Notes : Collects the component inductors of a mutual indcutor
//                 into a full set of data (i.e. type of mutual indcutor
//                 and all the component inductors.)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void CircuitContext::bundleMIs()
{
  std::string type;
  unsigned int i, j, mTableSize, mTableRowSize;
  StringToken field;
  std::vector< std::string > tmpInductorList, tmpCouplingList, tmpOtherParamsList;
  TokenVector tmpLine;

  // retrieve number of K lines in this (sub)circuit
  mTableSize = currentContextPtr_->allIndexedMIs_.size();
  
  for( i = 0; i < mTableSize; ++i )
  {
    // reset lists, tables and indices
    tmpInductorList.clear();
    tmpCouplingList.clear();
    tmpOtherParamsList.clear();

    // reset tmp vars
    std::string tmpName("");
    std::string tmpModel("");
    tmpLine.clear();

    std::vector< std::set< std::string > > & tmpTable =
      currentContextPtr_->getSharedInductorTable();

    mTableRowSize = currentContextPtr_->allIndexedMIs_[i].size();

    for( j = 0; j < mTableRowSize; ++j )
    {
      // select current mutual inductance
      MutualInductance & mutind =
        currentContextPtr_->mutualInductances_[
          currentContextPtr_->allIndexedMIs_[i][j]];

      // collect name segments
      if ( tmpName == "" )
        tmpName = mutind.name;
      else
        tmpName += "_" + mutind.name;

      // flag nonlinear coupling
      if( mutind.model != "" )
      {
        tmpModel = mutind.model;
        type = "N";
      }
      // else linear coupling
      else
      {
        type = "L";
      }

      // assemble inductor list
      std::map< std::string, std::string >::iterator mIter = mutind.inductors.begin();
      for( ; mIter != mutind.inductors.end(); ++mIter )
      {
        // add inductor data if not already present
        if( ( type == "N" ) ||
            ( tmpTable[mutind.sharedKey].find( (*mIter).first ) !=
              tmpTable[mutind.sharedKey].end() ) )
        {
          // retrieve the inductor name
          tmpInductorList.push_back( (*mIter).first );

          //Do some error-checking to make sure that the inductor has been
          //defined elsewhere!

          std::map<std::string,std::vector<std::string> >::iterator inducIter;
          inducIter = mutind.terminals.find((*mIter).first);
          if (inducIter == mutind.terminals.end())
          {
            Report::UserError0() << "Undefined inductor " << (*mIter).first << " in mutual inductor " << mutind.name << " definition.";
          }
          else
          {
            // retrieve the inductor nodes
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][0] );
            tmpInductorList.push_back( mutind.terminals[(*mIter).first][1] );
          }

          // retrieve the inductance value
          tmpInductorList.push_back((*mIter).second);
          
          // add in initial condition value 
          std::vector<Device::Param >::iterator mParamIter =  mutind.otherParams[mIter->first].begin();
          for( ; mParamIter != mutind.otherParams[mIter->first].end(); ++mParamIter)
          {
            // this is really backwards.  The device line has already been parsed into Param objects
            // but we need to put it back into a string to pack up a mutual inductor -- ick!
            // only do this for params that were actually given
            //if( mParamIter->given())
            //{
              //std::stringstream paramstr;
              //paramstr << mParamIter->tag() << "=" << mParamIter->stringValue();
              tmpInductorList.push_back( mParamIter->stringValue() );
            //}
        
          }

          // remove inductor from list to prevent duplicates
          if( type == "L" )
          {
            tmpTable[mutind.sharedKey].erase( (*mIter).first );
          }

          /*
          // manually check for nonlin/lin inductor overlap & exit
          else
          {
          // FIXME
          }
          */
        }
        
        // append to coupling list; names are a bit redundant, but easy to track
        tmpCouplingList.push_back( (*mIter).first );
      }

      // append coupling value to string
      std::stringstream cnvtr;
      cnvtr << mutind.coupling;
      tmpCouplingList.push_back( cnvtr.str() );
    }

    // fully contruct line from components and type
    if( tmpName != "" )
    {
      // append type
      field.string_ = "YMI" + type;
      tmpLine.push_back( field );

      // append name
      field.string_ = tmpName;
      tmpLine.push_back( field );

      // append list of (inductor, terminals, inductance)+
      for( j = 0; j < tmpInductorList.size(); ++j )
      {
        field.string_ = tmpInductorList[j];
        tmpLine.push_back( field );
      }

      // append coupling list
      for( j = 0; j < tmpCouplingList.size(); ++j )
      {
        field.string_ = tmpCouplingList[j];
        tmpLine.push_back( field );
      }
      
      // append model name
      if( type == "N" )
      {
        field.string_ = tmpModel;
        tmpLine.push_back( field );
      }

      // store in current context
      currentContextPtr_->kLines_.push_back( tmpLine );
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::getMILine
// Purpose       : Retrieve one tokenized device line
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
TokenVector &
CircuitContext::getMILine( int i )
{
  if( ( i < 0 ) || ( i > (int)currentContextPtr_->kLines_.size() ) )
  {
    // bounds checking error exit
    Report::UserError() << "Request exceeds number of mutual inductances in this subcircuit";
  }

  return currentContextPtr_->kLines_[i];
}

//----------------------------------------------------------------------------
// Function       : CircuitContext::augmentTotalDeviceCount
// Purpose        : Augment total device count after we process k-lines
//
// Special Notes  :
// Scope          : public
// Creator        : Keith Santarelli
// Creation Date  : 09/22/08
//----------------------------------------------------------------------------
void CircuitContext::augmentTotalDeviceCount(int kLineCount,
                                             int coupledICount,
                                             int YDeviceCount)
{
  // Get device count for the current context.
  int count = currentContextPtr_->deviceCount_;
  count += YDeviceCount - kLineCount - coupledICount;

  int linearCount = currentContextPtr_->linearDeviceCount_;
  linearCount -= coupledICount;

  // Remove the inductors from the count
  if (coupledICount)
  {
    currentContextPtr_->localDeviceCountMap_["L"] -= coupledICount;
  }

  if (count < 0)
  {
    Report::DevelFatal() << "Augmented number of devices is less than 0.";
  }
  else
  {
    currentContextPtr_->deviceCount_ = count;
    currentContextPtr_->linearDeviceCount_ = linearCount;
  }
}

} // namespace IO

//-----------------------------------------------------------------------------
// Function      : CircuitContext::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
int Pack<IO::CircuitContext>::packedByteCount(
  const IO::CircuitContext &    circuit_context)
{
  int byteCount = 0;
  int size, i;

  // count name
  byteCount += sizeof( int );
  byteCount += circuit_context.name_.length();

  // count device count
  byteCount += sizeof( int );

  // count number of devices
  byteCount += sizeof( int );

  // count number of characters for device types
  Device::DeviceCountMap::const_iterator it_dcM = circuit_context.localDeviceCountMap_.begin();
  Device::DeviceCountMap::const_iterator it_dceM = circuit_context.localDeviceCountMap_.end();
  for( ; it_dcM != it_dceM; ++it_dcM )
  {
    byteCount += sizeof( int );  
    byteCount += it_dcM->first.length();
    byteCount += sizeof( int );
  }

  // count models
  IO::ModelMap::const_iterator it_spbM = circuit_context.models_.begin();
  IO::ModelMap::const_iterator it_speM = circuit_context.models_.end();
  byteCount += sizeof( int );
  for( ; it_spbM != it_speM; ++it_spbM )
  {
    // ---- count key
    byteCount += sizeof( int );
    byteCount += it_spbM->first.length();

    // ---- count data
    byteCount += Xyce::packedByteCount(*it_spbM->second);
  }

  // count unresolved functions
  size = circuit_context.unresolvedFunctions_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += Xyce::packedByteCount(circuit_context.unresolvedFunctions_[i]);
  }

  // count instance list
  std::vector< std::string >::const_iterator it_stbL = circuit_context.subcktList_.begin();
  std::vector< std::string >::const_iterator it_steL = circuit_context.subcktList_.end();
  byteCount += sizeof( int );
  for( ; it_stbL != it_steL; ++it_stbL )
  {
    byteCount += sizeof( int );
    byteCount += it_stbL->length();
  }

  // count node list
  size = circuit_context.nodeList_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += sizeof( int );
    byteCount += circuit_context.nodeList_[ i ].length();
  }

  // count subcircuit parameters
  size = circuit_context.subcircuitParameters_.size();
  byteCount += sizeof( int );
  for (Util::ParamList::const_iterator it = circuit_context.subcircuitParameters_.begin(), end = circuit_context.subcircuitParameters_.end(); it != end; ++it)
  {
    byteCount += Pack<Util::Param>::packedByteCount(*it);
  }

  // count unresolved params
  size = circuit_context.unresolvedParams_.size();
  byteCount += sizeof( int );
  for (Util::UParamList::const_iterator it = circuit_context.unresolvedParams_.begin(), end = circuit_context.unresolvedParams_.end(); it != end; ++it)
  {
    byteCount += Pack<Util::Param>::packedByteCount(*it);
  }

  // pack global node names
  std::set<std::string>::const_iterator globalNodes_i, globalNodes_end;
  globalNodes_i = circuit_context.globalNodes_.begin();
  globalNodes_end = circuit_context.globalNodes_.end();
  byteCount += sizeof( int );
  for ( ; globalNodes_i != globalNodes_end ; ++globalNodes_i )
  {
    byteCount += sizeof( int );
    byteCount += globalNodes_i->length();
  }

  // count global params
  size = circuit_context.unresolvedGlobalParams_.size();
  byteCount += sizeof( int );
  for (Util::ParamList::const_iterator it = circuit_context.unresolvedGlobalParams_.begin(), end = circuit_context.unresolvedGlobalParams_.end(); it != end; ++it)
  {
    byteCount += Pack<Util::Param>::packedByteCount(*it);
  }

  //count mutual inductances
  size = circuit_context.mutualInductances_.size();
  byteCount += sizeof( int );
  for( i = 0; i < size; ++i )
  {
    byteCount += Xyce::packedByteCount(circuit_context.mutualInductances_[i]);
  }

  // pack list of coupled inductors
  std::set<std::string>::const_iterator coupledL_i, coupledL_end;
  coupledL_i = circuit_context.allCoupledInductors_.begin();
  coupledL_end = circuit_context.allCoupledInductors_.end();
  byteCount += sizeof( int );
  for ( ; coupledL_i != coupledL_end; ++coupledL_i )
  {
    byteCount += sizeof( int );
    byteCount += coupledL_i->length();
  }

  // count subcircuit contexts
  size = circuit_context.circuitContextTable_.size();
  byteCount += sizeof( int );
  unordered_map< std::string, IO::CircuitContext * >::const_iterator itsc;
  unordered_map< std::string, IO::CircuitContext * >::const_iterator itsc_end = circuit_context.circuitContextTable_.end();
  for ( itsc = circuit_context.circuitContextTable_.begin(); itsc != itsc_end; ++itsc )
  {
    byteCount += sizeof( int );
    byteCount += itsc->first.length();
    byteCount += Xyce::packedByteCount(*itsc->second);
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::pack
// Purpose       : Packs circuit context into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::CircuitContext>::pack(
  const IO::CircuitContext &    circuit_context,
  char *                        buf,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm )
{
  int size, length, i;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+Pack<IO::CircuitContext>::packedByteCount( circuit_context );
#endif

  // pack name
  length = circuit_context.name_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( circuit_context.name_.c_str(), length, buf, bsize, pos );

  // pack device count
  comm->pack( &circuit_context.deviceCount_, 1, buf, bsize, pos );

  // count number of devices
  size = circuit_context.localDeviceCountMap_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  Device::DeviceCountMap::const_iterator it_dcM = circuit_context.localDeviceCountMap_.begin();
  Device::DeviceCountMap::const_iterator it_dceM = circuit_context.localDeviceCountMap_.end();
  for( ; it_dcM != it_dceM; ++it_dcM ) 
  {
    length = it_dcM->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_dcM->first.c_str(), length, buf, bsize, pos );
    comm->pack( &it_dcM->second, 1, buf, bsize, pos );
  } 

  // pack models_
  IO::ModelMap::const_iterator it_spbM = circuit_context.models_.begin();
  IO::ModelMap::const_iterator it_speM = circuit_context.models_.end();
  size = circuit_context.models_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( ; it_spbM != it_speM; ++it_spbM )
  {
    // ---- pack key
    length = it_spbM->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_spbM->first.c_str(), length, buf, bsize, pos );

    // ---- pack data
    Xyce::pack(*it_spbM->second, buf, bsize, pos, comm );
  }

  // pack unresolved functions
  size = circuit_context.unresolvedFunctions_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    Xyce::pack(circuit_context.unresolvedFunctions_[i], buf, bsize, pos, comm );
  }

  // pack instance list
  std::vector< std::string >::const_iterator it_stbL = circuit_context.subcktList_.begin();
  std::vector< std::string >::const_iterator it_steL = circuit_context.subcktList_.end();
  size = circuit_context.subcktList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( ; it_stbL != it_steL; ++it_stbL )
  {
    length = it_stbL->length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( it_stbL->c_str(), length, buf, bsize, pos );
  }

  // pack node list
  size = circuit_context.nodeList_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    length = circuit_context.nodeList_[ i ].length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( circuit_context.nodeList_[ i ].c_str(), length, buf, bsize, pos );
  }

  // pack subcircuit parameters
  size = circuit_context.subcircuitParameters_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for (Util::ParamList::const_iterator it = circuit_context.subcircuitParameters_.begin(), end = circuit_context.subcircuitParameters_.end(); it != end; ++it)
  {
    Pack<Util::Param>::pack(*it, buf, bsize, pos, comm );
  }

  // pack unresolved params
  size = circuit_context.unresolvedParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for (Util::UParamList::const_iterator it = circuit_context.unresolvedParams_.begin(), end = circuit_context.unresolvedParams_.end(); it != end; ++it)
  {
    Pack<Util::Param>::pack(*it, buf, bsize, pos, comm );
  }

  // pack global node names
  size = circuit_context.globalNodes_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  std::set<std::string>::const_iterator globalNodes_i, globalNodes_end;
  globalNodes_i = circuit_context.globalNodes_.begin();
  globalNodes_end = circuit_context.globalNodes_.end();
  for ( ; globalNodes_i != globalNodes_end ; ++globalNodes_i )
  {
    length = globalNodes_i->length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( globalNodes_i->c_str(), length, buf, bsize, pos );
  }

  // pack global params
  size = circuit_context.unresolvedGlobalParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for (Util::ParamList::const_iterator it = circuit_context.unresolvedGlobalParams_.begin(), end = circuit_context.unresolvedGlobalParams_.end(); it != end; ++it)
  {
    Pack<Util::Param>::pack(*it, buf, bsize, pos, comm );
  }

  // pack mutual inductances
  size = circuit_context.mutualInductances_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( i = 0; i < size; ++i )
  {
    Xyce::pack(circuit_context.mutualInductances_[i], buf, bsize, pos, comm );
  }

  // pack coupled inductor names
  size = circuit_context.allCoupledInductors_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  std::set<std::string>::const_iterator coupledL_i, coupledL_end;
  coupledL_i = circuit_context.allCoupledInductors_.begin();
  coupledL_end = circuit_context.allCoupledInductors_.end();
  for ( ; coupledL_i != coupledL_end; ++coupledL_i)
  {
    length = coupledL_i->length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( coupledL_i->c_str(), length, buf, bsize, pos );
  }
  
  // pack circuitContextTable_
  size = circuit_context.circuitContextTable_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  unordered_map< std::string, IO::CircuitContext * >::const_iterator itsc = circuit_context.circuitContextTable_.begin();
  unordered_map< std::string, IO::CircuitContext * >::const_iterator itsc_end = circuit_context.circuitContextTable_.end();
  for ( ; itsc != itsc_end; ++itsc )
  {
    length = itsc->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( itsc->first.c_str(), length, buf, bsize, pos );
    Xyce::pack(*itsc->second, buf, bsize, pos, comm );
  }

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    Report::DevelFatal() << "CircuitContext::pack - predicted pos (" << predictedPos 
                         << ") does not match actual pos (" << pos << ")";
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : CircuitContext::unpack
// Purpose       : Unpacks circuit context from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::CircuitContext>::unpack(
  IO::CircuitContext &     circuit_context,
  char *                   pB,
  int                      bsize,
  int &                    pos,
  Parallel::Communicator * comm )
{
  int size = 0;
  int length = 0;

  // unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  circuit_context.name_ = std::string( ( pB + pos ), length );
  pos += length;

  // unpack device count
  comm->unpack( pB, bsize, pos, &circuit_context.deviceCount_, 1 );

  // unpack device types
  comm->unpack( pB, bsize, pos, &size, 1 );

  for (int i = 0; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string aString(std::string( ( pB + pos ), length ));
    pos += length;
    comm->unpack( pB, bsize, pos, &length, 1 );

    circuit_context.localDeviceCountMap_[ aString ] = length;
  } 

  // unpack models_
  comm->unpack( pB, bsize, pos, &size, 1 );

  for (int i = 0; i < size; ++i )
  {

    // ---- unpack key
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string aString(std::string( ( pB + pos ), length ));
    pos += length;

    // ---- unpack data
    circuit_context.models_[ aString ] = new IO::ParameterBlock();
    Xyce::unpack(*circuit_context.models_[aString], pB, bsize, pos, comm );
  }

  // unpack unresolved functions (vector of FunctionBlock)
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i )
  {
    IO::FunctionBlock aFuncBlk;
    Xyce::unpack(aFuncBlk, pB, bsize, pos, comm );
    circuit_context.unresolvedFunctions_.push_back( aFuncBlk );
  }

  // unpack instance list
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    circuit_context.subcktList_.push_back( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack node list
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    circuit_context.nodeList_.push_back( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack subcircuit parameters
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    Util::Param aParam;
    Pack<Util::Param>::unpack(aParam, pB, bsize, pos, comm );
    circuit_context.subcircuitParameters_.push_back( aParam );
  }

  // unpack unresolved params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    Util::Param aParam;
    Pack<Util::Param>::unpack(aParam, pB, bsize, pos, comm );
    circuit_context.unresolvedParams_.insert( aParam );
  }

  // unpack global node names
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    circuit_context.globalNodes_.insert( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack unresolved params vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    Util::Param aParam;
    Pack<Util::Param>::unpack(aParam, pB, bsize, pos, comm );
    circuit_context.unresolvedGlobalParams_.push_back( aParam );
  }

  // unpack mutual inductances vector
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    IO::CircuitContext::MutualInductance MI;
    Xyce::unpack(MI, pB, bsize, pos, comm );
    circuit_context.mutualInductances_.push_back( MI );
  }

  // unpack coupled inductors
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0;  i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    circuit_context.allCoupledInductors_.insert( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack circuitContextTable_
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int i = 0; i < size; ++i)
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string tmp( ( pB + pos ), length );
    pos += length;
    std::pair< unordered_map< std::string, IO::CircuitContext *>::iterator, bool > p =
      circuit_context.circuitContextTable_.insert(
        std::pair< std::string, IO::CircuitContext *>(
          tmp,
          new IO::CircuitContext(circuit_context.expressionGroup_, circuit_context.opBuilderManager_, circuit_context.contextList_, circuit_context.currentContextPtr_ ) ) );

    // set the parent context of my children to me
    p.first->second->setParentContextPtr( &circuit_context );
    Xyce::unpack(*p.first->second, pB, bsize, pos, comm );
  }

  circuit_context.currentContextPtr_ = &circuit_context;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 8/27/04
//-----------------------------------------------------------------------------
template<>
int
Pack<IO::CircuitContext::MutualInductance>::packedByteCount(
  const IO::CircuitContext::MutualInductance &  mutual_inductance)
{
  int byteCount = 0;

  // coupling value
  byteCount += sizeof(int);
  byteCount += mutual_inductance.coupling.length();

  // model name
  byteCount += sizeof(int);
  byteCount += mutual_inductance.model.length();

  // first inductor name
  byteCount += sizeof(int);
  byteCount += mutual_inductance.firstInductor.length();

  // inductor info
  byteCount += sizeof(int);
  std::map<std::string,std::string>::const_iterator iterI = mutual_inductance.inductors.begin();
  std::map<std::string,std::string>::const_iterator  endI = mutual_inductance.inductors.end();
  for( ; iterI != endI; ++iterI )
  {
    byteCount += sizeof(int);
    byteCount += iterI->first.length();
    byteCount += sizeof(int);
    byteCount += iterI->second.length();
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::pack
// Purpose       : Packs MI into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::CircuitContext::MutualInductance>::pack(
  const IO::CircuitContext::MutualInductance &  mutual_inductance,
  char *                                        buf,
  int                                           bsize,
  int &                                         pos,
  Parallel::Communicator *                      comm ) 
{
  int size, length;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+Pack<IO::CircuitContext::MutualInductance>::packedByteCount( mutual_inductance );
#endif

  // coupling value
  length = mutual_inductance.coupling.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( mutual_inductance.coupling.c_str(), length, buf, bsize, pos );

  // model name
  length = mutual_inductance.model.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( mutual_inductance.model.c_str(), length, buf, bsize, pos );

  // first inductor name
  length = mutual_inductance.firstInductor.length();
  comm->pack( &length, 1, buf, bsize, pos );
  if( length ) comm->pack( mutual_inductance.firstInductor.c_str(), length, buf, bsize, pos );

  // pack inductors
  size = mutual_inductance.inductors.size();
  comm->pack( &size, 1, buf, bsize, pos );
  std::map<std::string,std::string>::const_iterator iterI = mutual_inductance.inductors.begin();
  std::map<std::string,std::string>::const_iterator  endI = mutual_inductance.inductors.end();
  for( ; iterI != endI; ++iterI )
  {
    length = iterI->first.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( iterI->first.c_str(), length, buf, bsize, pos );
    length = iterI->second.length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( iterI->second.c_str(), length, buf, bsize, pos );
  }
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    Report::DevelFatal() << "CircuitContext::MutualInductance::pack - predicted pos (" << predictedPos 
                         << ") does not match actual pos (" << pos << ")";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : CircuitContext::MutualInductance::unpack
// Purpose       : Unpacks MI from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 08/27/04
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::CircuitContext::MutualInductance>::unpack(
  IO::CircuitContext::MutualInductance &        mutual_inductance,
  char *                                        pB,
  int                                           bsize,
  int &                                         pos,
  Parallel::Communicator *                      comm )
{
  int size, length, i;

  // coupling value
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    mutual_inductance.coupling = std::string( ( pB + pos ), length );
    pos += length;
  }

  // model name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    mutual_inductance.model = std::string( ( pB + pos ), length );
    pos += length;
  }

  // first inductor name
  comm->unpack( pB, bsize, pos, &length, 1 );
  if( length )
  {
    mutual_inductance.firstInductor = std::string( ( pB + pos ), length );
    pos += length;
  }

  // unpack inductors
  mutual_inductance.inductors.clear();
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( i = 0; i < size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string name(std::string( ( pB + pos ), length ));
    pos += length;
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string value(std::string( ( pB + pos ), length ));
    pos += length;
    mutual_inductance.inductors[name] = value;
  }
}

} // namespace Xyce

//-----------------------------------------------------------------------------
// Function      : packCircuitContext
// Purpose       : send circuit context to all procs
// Special Notes :
// Scope         : non-member function
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Xyce::IO::packCircuitContext(Xyce::IO::CircuitContext* circuit_contexts, char* char_buffer, 
                                 int char_buffer_size, Xyce::Parallel::Communicator* pds_comm_ptr )
{
  int bsize = 0;

  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;

    // pack circuit contexts
    Xyce::pack(*circuit_contexts, char_buffer, char_buffer_size, pos, pds_comm_ptr );

    bsize=pos;
  }

  return bsize;
}


//-----------------------------------------------------------------------------
// Function      : unpackCircuitContext
// Purpose       : receive circuit context from proc 0
// Special Notes :
// Scope         : non-member function
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Xyce::IO::unpackCircuitContext(Xyce::IO::CircuitContext* circuit_contexts, char* char_buffer, 
                                    int bsize, Xyce::Parallel::Communicator* pds_comm_ptr )
{
  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;

    circuit_contexts->setParentContextPtr( NULL );

    Xyce::unpack(*circuit_contexts, char_buffer, bsize, pos, pds_comm_ptr);

    circuit_contexts->resolve( std::vector<Device::Param>() );
  }

  return true;
}

