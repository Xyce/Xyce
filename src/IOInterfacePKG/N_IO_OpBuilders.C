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
// Purpose        : Build various operators, like V(), I(), N(), P(), W(),
//                  DNO() and DNI()
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  : 08/04/14
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_Op.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>
#include <N_UTL_DeviceNameConverters.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {

namespace {

//-------------------------------------------------------------------------- 
// Function      : parameterNameAndArgs 
// Purpose       : Return a parameter's names and arguments
// Special Notes :  
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
void parameterNameAndArgs(
  std::string &                         name,
  std::vector<std::string> &            args,
  Util::ParamList::const_iterator &     it)
{
  const std::string &param_tag = (*it).tag();

  // don't enter this if statement if param_tag == "SENS", or is a Y-device
  if ( (*it).getType() == Util::INT && (param_tag[0] == 'V' || param_tag[0] == 'I' || param_tag[0] == 'N' ||
                                        param_tag[0] == 'P' || param_tag[0] == 'W' || param_tag[0] == 'D' ||
			                param_tag[0] == 'S' || param_tag[0] == 'Y' || param_tag[0] == 'Z') )
  {
    std::ostringstream oss;
    oss << param_tag << "(";
    int arg_count = (*it).getImmutableValue<int>();
    for (int i = 0; i < arg_count; ++i)
    {
      ++it;
      if (i != 0)
        oss << ",";
      oss << (*it).tag();
      args.push_back((*it).tag());
    }
    oss << ")";
    name = oss.str();
  }
}

//-------------------------------------------------------------------------- 
// Function      : findNodeIndex 
// Purpose       : Determine whether a requested node name is "valid".
//                 That includes node names in the solution vector, the
//                 aliasNodeMap and also the Ground node (0).
//                 The return values are:
//                   a) -2 if the node is not found on this processor.
//                   b) -1 if it is a Ground node.
//                   c) the node index (0 ...N), otherwise.
// Special Notes :  
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
int findNodeIndex(
  const std::string &       name,
  const NodeNameMap &       node_map,
  const AliasNodeMap &      alias_map)
{
  // The return value will be -2 if the specified node name is not
  // found AND the specified node name is not Ground (0). A value 
  // of -2 may be returned even for a valid node name, in parallel,
  // because of how node_map is used in parallel.
  int nodeIndex = -2;

  // Handle Ground (0).  If .PREPROCESS REPLACEGROUND was used,
  // then this function assumes that all of the GND nodes have
  // been replaced with 0 by this point in Xyce startup.
  if ( name == "0")
  { 
    // A node index of -1 will be used in the various get() functions
    // to denote the Ground node.
    nodeIndex = -1;
  } 
  else
  {
    NodeNameMap::const_iterator node_it = node_map.find(name);
    if (node_it == node_map.end()) 
    {
      // If the specified node name is not found in the node_map then
      // look for it in the AliasNodeMap. An example where this can 
      // happen is for a subcircuit interface node.
      AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
      if (alias_node_it != alias_map.end())
      {
        // (*alias_node_it).second will be the "real name" of the alias node.
        node_it = node_map.find((*alias_node_it).second);
      }
    }

    // get the node index if it exists on this processor.
    if (node_it != node_map.end()) 
    {
      nodeIndex = (*node_it).second;
    } 
  }

  return nodeIndex;
}

} // namespace <unnamed>

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitTemperatureOpBuilder 
// Purpose       : This creates an OutputMgrTemperatureOp
// Special Notes :  
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitTemperatureOpBuilder : public Util::Op::Builder
{
  CircuitTemperatureOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitTemperatureOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrTemperatureOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "TEMP") 
    {
      new_op  = new OutputMgrTemperatureOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitTimeOpBuilder 
// Purpose       : This creates an OutputMgrTimeOp
// Special Notes :  
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitTimeOpBuilder : public Util::Op::Builder
{
  CircuitTimeOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitTimeOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrTimeOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "TIME") 
    {
      new_op  = new OutputMgrTimeOp(param_tag, outputManager_, 1.0);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitOutputNoiseOpBuilder 
// Purpose       : This creates an OutputMgrOutputNoiseOp
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitOutputNoiseOpBuilder : public Util::Op::Builder
{
  CircuitOutputNoiseOpBuilder(const OutputMgr & output_manager,
                              const Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)
  {}

  virtual ~CircuitOutputNoiseOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrOutputNoiseOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "ONOISE")
    {
      if (!analysisManager_.getNoiseFlag())
      {
        Report::UserError0() << "ONOISE operator only supported for .NOISE analyses";
        return new_op;
      }

      new_op  = new OutputMgrOutputNoiseOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  const Analysis::AnalysisManager &    analysisManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitInputNoiseOpBuilder 
// Purpose       : This creates an OutputMgrInputNoiseOp 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitInputNoiseOpBuilder : public Util::Op::Builder
{
  CircuitInputNoiseOpBuilder(const OutputMgr & output_manager,
                             const Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)

  {}

  virtual ~CircuitInputNoiseOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrInputNoiseOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "INOISE")
    {
      if (!analysisManager_.getNoiseFlag())
      {
        Report::UserError0() << "INOISE operator only supported for .NOISE analyses";
        return new_op;
      }

      new_op  = new OutputMgrInputNoiseOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  const Analysis::AnalysisManager &    analysisManager_;
};

//--------------------------------------------------------------------------
// Structure     : Util::Op::Builder::CircuitNoiseContOpBuilder
// Purpose       : This creates an either an OutputMgrOutputNoiseContOp
//               : or an OutputMgrInputNoiseContOp
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 11/20/17
//--------------------------------------------------------------------------
struct CircuitNoiseContOpBuilder : public Util::Op::Builder
{
  CircuitNoiseContOpBuilder(const OutputMgr & output_manager,
                            const Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)
  {}

  virtual ~CircuitNoiseContOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrOutputNoiseContOp>();
    builder_manager.addCreateFunction<OutputMgrInputNoiseContOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'D' && args.size() > 0)
    {
      if (!analysisManager_.getNoiseFlag())
      {
        Report::UserError0() << "DNO and DNI operators only supported for .NOISE analyses";
        return new_op;
      }

      // The DNO() and DNI() operators come in two forms.  For example, DNO(Q1) and DN(Q1,rc).
      // So, we need to find the index of the device (in the noiseDataVec_ of the NOISE object),
      // and then the index(es) of the requested noise type for that device, if the device has
      // multiple types sources.
      int devIndex = -1;
      std::vector<int> typeIndex;
      NodeNameMap::const_iterator nd_it =  outputManager_.getNoiseDeviceNameMap().find(args[0]+"_ND");
      if (nd_it != outputManager_.getNoiseDeviceNameMap().end())
           devIndex = (*nd_it).second;

      // Now find the index(es) of the noise type, for the specified device, if one was requested.
      if (args.size() == 2)
      {
	NodeNameMap::const_iterator nt_it = outputManager_.getNoiseTypeNameMap().begin();
        int suffix = 0;

        while ( !(nt_it == outputManager_.getNoiseTypeNameMap().end()) )
	{
          // For ADMS devices, that may have duplicate entries for a given noise type, the entries
          // were "suffixed" with _0, _1, _2, ... .  If there were no duplicate entries, then just
          // the _0 suffix was used.  So, the DNO and DNI operators use a vector-of-ints for
          // typeIndex to handle this.
          std::ostringstream s;
          s << suffix;

	  nt_it = outputManager_.getNoiseTypeNameMap().find("noise_" + args[0] + "_" + args[1] + "_" + s.str());
          if (nt_it != outputManager_.getNoiseTypeNameMap().end())
	  {
            typeIndex.push_back((*nt_it).second);
            ++suffix;
          }
        }
      }

      // Only make the op if the relevant indices have been found.
      // If not then this op will be flagged as in error, since new_op is still a null pointer
      if ( (devIndex != -1) && ( ( args.size() == 1 ) || ((args.size() == 2) && !typeIndex.empty()) ) )
      {
        if (param_tag == "DNI")
        {
          new_op  = new OutputMgrInputNoiseContOp(name, devIndex, typeIndex, outputManager_);
        }
        else if (param_tag == "DNO")
        {
          new_op  = new OutputMgrOutputNoiseContOp(name, devIndex, typeIndex, outputManager_);
        }
      }
    }

    if (new_op)
          new_op->addArg(args[0]);

    return new_op;
  }

private:
  const OutputMgr &                    outputManager_;
  const Analysis::AnalysisManager &    analysisManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitFrequencyOpBuilder 
// Purpose       : This creates an OutputMgrFrequencyOp 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitFrequencyOpBuilder : public Util::Op::Builder
{
  CircuitFrequencyOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~CircuitFrequencyOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrFrequencyOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if ( param_tag == "FREQ" ){
      new_op  = new OutputMgrFrequencyOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::StepSweepOpBuilder 
// Purpose       : This creates an OutputMgrStepSweepOp 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct StepSweepOpBuilder : public Util::Op::Builder
{
  StepSweepOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~StepSweepOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrStepSweepOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    for (size_t i = 0; i < outputManager_.getStepSweepVector().size(); ++i)
    {
      if (param_tag == outputManager_.getStepSweepVector()[i].name)
      {
        new_op = new OutputMgrStepSweepOp(param_tag, outputManager_, i);
        break;
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::DCSweepOpBuilder 
// Purpose       : This creates an OutputMgrDCSweepOp 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct DCSweepOpBuilder : public Util::Op::Builder
{
  DCSweepOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~DCSweepOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrDCSweepOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    for (size_t i = 0; i < outputManager_.getDCSweepVector().size(); ++i)
    {
      if (param_tag == outputManager_.getDCSweepVector()[i].name)
      {
        new_op = new OutputMgrDCSweepOp(param_tag, outputManager_, i);
        break;
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::DCSweepCurrentValueOpBuilder 
// Purpose       : This creates an OutputMgrDCSweepCurrentValueOp 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct DCSweepCurrentValueOpBuilder : public Util::Op::Builder
{
  DCSweepCurrentValueOpBuilder(const OutputMgr & output_manager)
    : outputManager_(output_manager)
  {}

  virtual ~DCSweepCurrentValueOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<OutputMgrDCSweepCurrentValueOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "sweep") {
      new_op  = new OutputMgrDCSweepCurrentValueOp(param_tag, outputManager_);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
};

//--------------------------------------------------------------------------
// Structure     : Util::Op::Builder::StepNumOpBuilder
// Purpose       : This creates a StepNumOp
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 08/19/14
//--------------------------------------------------------------------------
struct StepNumOpBuilder : public Util::Op::Builder
{
  StepNumOpBuilder(const OutputMgr &output_manager)
    :outputManager_(output_manager)
  {}

  virtual ~StepNumOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<StepNumOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "STEPNUM") {
      new_op  = new StepNumOp(param_tag, outputManager_);
      new_op->addArg(param_string);
    }

    return new_op;
  }

private:
  const OutputMgr &       outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::InternalVariableOpBuilder 
// Purpose       : Builds various "solution operators".  See the list
//                 below in registerCreateFunctions().  It includes N(), 
//                 StateOp() and StoreOp().  The code below also lists
//                 NR(), NI(), NM(), NP() and NDB(). However, those operators
//                 won't pass netlist parsing and aren't fully implemented
//                 yet. See SON Bug 880 for more details. 
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct InternalVariableOpBuilder : public Util::Op::Builder
{
  InternalVariableOpBuilder(const OutputMgr & output_manager, Parallel::Machine comm)
    : outputManager_(output_manager),
      comm_(comm)
  {}

  virtual ~InternalVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseDegOp>();
    builder_manager.addCreateFunction<SolutionPhaseRadOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<StateOp>();
    builder_manager.addCreateFunction<StoreOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'N' && args.size() == 1)
    {
      // An index of -2 means the node was not found on this processor.  -1 is the "index"
      // for the Ground node.  An index >= 0 means that the node was found on this processor.
      int index = findNodeIndex(args[0], outputManager_.getSolutionNodeMap(), outputManager_.getAliasNodeMap());
      
      // Only make the op if the node was found on some processor.  
      int maxIndex=index;
      Parallel::AllReduce(comm_, MPI_MAX, &maxIndex, 1); 
      if (maxIndex > -2)
      {
        if (param_tag == "N" )
        {
          new_op = new SolutionOp(name, index);
        }
        else if (param_tag == "NR" )
        {
          new_op = new SolutionRealOp(name, index);
        }
        else if (param_tag == "NI" )
        {
          new_op = new SolutionImaginaryOp(name, index);
        }
        else if (param_tag == "NM" )
        {
          new_op = new SolutionMagnitudeOp(name, index);
        }
        else if (param_tag == "NP" )
        {
          if (outputManager_.getPhaseOutputUsesRadians())
	  {
            new_op = new SolutionPhaseRadOp(name, index);
          }
          else
          {
            new_op = new SolutionPhaseDegOp(name, index);
          }
        }
        else if (param_tag == "NDB" )
        {
          new_op = new SolutionDecibelsOp(name, index);
        }
      }
      else
      {
        int index = findNodeIndex(args[0], outputManager_.getStateNodeMap(), outputManager_.getAliasNodeMap());
        
        // Only make the op if the node was found on some processor.
        int maxIndex=index;
        Parallel::AllReduce(comm_, MPI_MAX, &maxIndex, 1);
        if (maxIndex > -2)
        {
          new_op = new StateOp(name, index);
        }
        else
        {
          NodeNameMap::const_iterator it = outputManager_.getStoreNodeMap().find(args[0]);
          if (it != outputManager_.getStoreNodeMap().end())
          {
            new_op = new StoreOp(name, (*it).second);
          }
          else
          {
            new_op = new Util::Op::UndefinedOp(param_tag);
          }
        }
      }

      if (new_op)
        new_op->addArg(args[0]);
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  Parallel::Machine     comm_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::VoltageVariableOpBuilder 
// Purpose       : Builds various "solution operators" for voltage variables.
//                 See the list below in registerCreateFunctions().  This
//                 includes the operators for V(), VR(), VI(), VM(), VP()
//                 and VDB(). It includes, for example, V(a) and also the
//                 corresponding "voltage difference" form of V(a,b).
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct VoltageVariableOpBuilder : public Util::Op::Builder
{
  VoltageVariableOpBuilder(const OutputMgr & output_manager, Parallel::Machine comm)
    : outputManager_(output_manager),
      comm_(comm)
  {}

  virtual ~VoltageVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionRealOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseDegOp>();
    builder_manager.addCreateFunction<SolutionPhaseRadOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<VoltageDifferenceOp>();
    builder_manager.addCreateFunction<VoltageDifferenceRealOp>();
    builder_manager.addCreateFunction<VoltageDifferenceImaginaryOp>();
    builder_manager.addCreateFunction<VoltageDifferenceMagnitudeOp>();
    builder_manager.addCreateFunction<VoltageDifferencePhaseDegOp>();
    builder_manager.addCreateFunction<VoltageDifferencePhaseRadOp>();
    builder_manager.addCreateFunction<VoltageDifferenceDecibelsOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'V' && args.size() > 0)
    {

      // Solution variable
      if (args.size() == 1)
      {
        // An index of -2 means the node was not found on this processor.  -1 is the "index"
        // for the Ground node.  An index >= 0 means that the node was found on this processor.
        int index = findNodeIndex(args[0], outputManager_.getSolutionNodeMap(), outputManager_.getAliasNodeMap());

        // Only make the op if the node was found on some processor.
        int maxIndex=index;
        Parallel::AllReduce(comm_, MPI_MAX, &maxIndex, 1);
        if ( maxIndex > -2 )
	{
          if (param_tag == "V" )
          {
            new_op = new SolutionOp(name, index);
          }
          else if (param_tag == "VR" )
          {
            new_op = new SolutionRealOp(name, index);
          }
          else if (param_tag == "VI" )
          {
            new_op = new SolutionImaginaryOp(name, index);
          }
          else if (param_tag == "VM" )
          {
            new_op = new SolutionMagnitudeOp(name, index);
          }
          else if (param_tag == "VP" )
          {
	    if (outputManager_.getPhaseOutputUsesRadians())
	     {
               new_op = new SolutionPhaseRadOp(name, index);
             }
             else
             {
               new_op = new SolutionPhaseDegOp(name, index);
             }
          }
          else if (param_tag == "VDB" )
          {
            new_op = new SolutionDecibelsOp(name, index);
          }

          if (new_op)
            new_op->addArg(args[0]);
        }
      }

      // Voltage Difference
      else if (args.size() == 2)
      {
        // An index of -2 means the node was not found on this processor.  -1 is the "index"
        // for the Ground node.  An index >= 0 means that the node was found on this processor.
        int index1 = findNodeIndex(args[0], outputManager_.getSolutionNodeMap(), outputManager_.getAliasNodeMap());
        int index2 = findNodeIndex(args[1], outputManager_.getSolutionNodeMap(), outputManager_.getAliasNodeMap());       
        
        // Only make the op if both nodes were found on some processor, where the nodes may have
        // been found on different processors in parallel.
        int maxIndex1=index1;
        int maxIndex2=index2;
        Parallel::AllReduce(comm_, MPI_MAX, &maxIndex1, 1);
        Parallel::AllReduce(comm_, MPI_MAX, &maxIndex2, 1);
        if ( (maxIndex1 > -2) && (maxIndex2 > -2) )
        {
          if (param_tag == "V" )
          {
            new_op = new VoltageDifferenceOp(name, index1, index2);
          }
          else if (param_tag == "VR" )
          {
            new_op = new VoltageDifferenceRealOp(name, index1, index2);
          }
          else if (param_tag == "VI" )
          {
            new_op = new VoltageDifferenceImaginaryOp(name, index1, index2);
          }
          else if (param_tag == "VM" )
          {
            new_op = new VoltageDifferenceMagnitudeOp(name, index1, index2);
          }
          else if (param_tag == "VP" )
          {
            if (outputManager_.getPhaseOutputUsesRadians())
	    {
              new_op = new VoltageDifferencePhaseRadOp(name, index1, index2);
            }
            else
            {
              new_op = new VoltageDifferencePhaseDegOp(name, index1, index2);
            }
          }
          else if (param_tag == "VDB" )
          {
            new_op = new VoltageDifferenceDecibelsOp(name, index1, index2);
          }

          if (new_op)
            new_op->addArgs(args.begin(), args.end());
        }
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  Parallel::Machine     comm_;
};

//-------------------------------------------------------------------------- 
// Allowed set of current operators, besides I()
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
namespace {
static const char * const func_names[] = {"II", "IR", "IP", "IM", "IDB"};
}

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CurrentVariableOpBuilder 
// Purpose       : Builds for various "solution operators" and "store
//                 operators" for lead current variables.  See the list
//                 below in registerCreateFunctions().  This includes the 
//                 operators for I(), IR(), II(), IM(), IP() and IDB(). 
//                 The BranchData forms are used to get values for lead 
//                 currents that are not part of the solution vector.
//                 Note also that special handling is needed for YPDE, C 
//                 and L devices.
// Special Notes :
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CurrentVariableOpBuilder : public Util::Op::Builder
{

  CurrentVariableOpBuilder(const OutputMgr & output_manager,
                           const Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)
  {}

  virtual ~CurrentVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SolutionOp>();
    builder_manager.addCreateFunction<SolutionRealOp>();
    builder_manager.addCreateFunction<SolutionImaginaryOp>();
    builder_manager.addCreateFunction<SolutionMagnitudeOp>();
    builder_manager.addCreateFunction<SolutionPhaseDegOp>();
    builder_manager.addCreateFunction<SolutionPhaseRadOp>();
    builder_manager.addCreateFunction<SolutionDecibelsOp>();
    builder_manager.addCreateFunction<StoreOp>();
    builder_manager.addCreateFunction<StoreRealOp>();
    builder_manager.addCreateFunction<StoreImaginaryOp>();
    builder_manager.addCreateFunction<StoreMagnitudeOp>();
    builder_manager.addCreateFunction<StorePhaseDegOp>();
    builder_manager.addCreateFunction<StorePhaseRadOp>();
    builder_manager.addCreateFunction<StoreDecibelsOp>();
    builder_manager.addCreateFunction<BranchDataCurrentOp>();
    builder_manager.addCreateFunction<BranchDataCurrentRealOp >();
    builder_manager.addCreateFunction<BranchDataCurrentImaginaryOp>();
    builder_manager.addCreateFunction<BranchDataCurrentMagnitudeOp>();
    builder_manager.addCreateFunction<BranchDataCurrentPhaseDegOp>();
    builder_manager.addCreateFunction<BranchDataCurrentPhaseRadOp>();
    builder_manager.addCreateFunction<BranchDataCurrentDecibelsOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if (param_tag[0] == 'I' && !args.empty())
    {
      // Node name could be circuit_context:DeviceTypeDeviceName while internally it should be DeviceType:circuit_context:DeviceName.
      std::string modifiedName = Util::xyceDeviceNameToSpiceName(args[0]);

      // could be a device lead current "DEV_I" or a branch current.
      // so we don't have to duplicate solution vars(branch currents) in the
      // store vector, look for each type.
      bool param_func = std::find(func_names, func_names + sizeof(func_names)/sizeof(func_names[0]), param_tag) != func_names + sizeof(func_names)/sizeof(func_names[0]);
      std::string store_name = modifiedName + ":DEV_" + (param_func ? "I" : param_tag);  // if it is in the state/store vec.
      std::string leadCurrent_name = modifiedName + ":BRANCH_D";
      if (!param_func && (param_tag.length() > 1) )
      {
        leadCurrent_name += param_tag[1];
      } 

      // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
      // we try to find devices based on store_name in the following blocks of code,
      // so do this modification now.
      std::string::size_type space = store_name.find_first_of(" ");
      if (space != std::string::npos)
      {
        if (space == 4 && store_name.substr(0, 4) == "YPDE")
        {
          store_name.replace(4, 1, ":");
        }
      }
      
      NodeNameMap::const_iterator it;
      // Search solution vector first, specifically excluding the C Device.  For the C device, 
      // we will search in the lead current vector.  If the IC parameter is specified for a
      // capacitor, then a voltage source is used to enforce that IC value at the DCOP.  However,
      // the contribution from that V-source has been included in the lead-current vector.  
      if ( !new_op && !startswith_nocase(modifiedName, "C") ) 
      {
        std::string solution_name = modifiedName + "_BRANCH";         // if it is in the solution vec.
        const NodeNameMap &x = outputManager_.getSolutionNodeMap();
        it = x.find(solution_name);
        if (it != outputManager_.getSolutionNodeMap().end())
        {
          int index = (*it).second;
          if (param_tag == "I" )
          {
            new_op = new SolutionOp(name, index);
          }
          else if (param_tag == "IR" )
          {
            new_op = new SolutionRealOp(name, index);
          }
          else if (param_tag == "II" )
          {
            new_op = new SolutionImaginaryOp(name, index);
          }
          else if (param_tag == "IM" )
          {
            new_op = new SolutionMagnitudeOp(name, index);
          }
          else if (param_tag == "IP" )
          {
            if (outputManager_.getPhaseOutputUsesRadians())
	    {
              new_op = new SolutionPhaseRadOp(name, index);
            }
            else
            {
              new_op = new SolutionPhaseDegOp(name, index);
            }
          }
          else if (param_tag == "IDB" )
          {
            new_op = new SolutionDecibelsOp(name, index);
          }
        }
      }

      // Search solution
      if (!new_op && startswith_nocase(modifiedName, "L")) // Mutual inductor special name
      {
        const NodeNameMap &x = outputManager_.getSolutionNodeMap();
        it = x.find(modifiedName);
        if (it != outputManager_.getSolutionNodeMap().end())
        {
          int index = (*it).second;
          if (param_tag == "I" )
          {
            new_op = new SolutionOp(name, index);
          }
          else if (param_tag == "IR" )
          {
            new_op = new SolutionRealOp(name, index);
          }
          else if (param_tag == "II" )
          {
            new_op = new SolutionImaginaryOp(name, index);
          }
          else if (param_tag == "IM" )
          {
            new_op = new SolutionMagnitudeOp(name, index);
          }
          else if (param_tag == "IP" )
          {
            if (outputManager_.getPhaseOutputUsesRadians())
	    {
              new_op = new SolutionPhaseRadOp(name, index);
            }
            else
            {
              new_op = new SolutionPhaseDegOp(name, index);
            }
          }
          else if (param_tag == "IDB" )
          {
            new_op = new SolutionDecibelsOp(name, index);
          }
        }
      }

      // The remaining places to look for lead currents are NOT valid for either .AC
      // or .NOISE analyses.  So, return new_op=0 at this point for those two cases.
      // Note: the code can't also return a UserError0() message at this point.  Otherwise,
      // I(V1) and I(L1) will stop working in parallel for both the .AC and .NOISE cases
      // since a valid new_op is only non-NULL on one processor in parallel.
      if (!new_op && (analysisManager_.getNoiseFlag() || analysisManager_.getACFlag()) )
      {
        return new_op;
      }

      // These cases are valid for .DC and .TRAN
      if( !new_op) 
      {
        it = outputManager_.getBranchVarsNodeMap().find(leadCurrent_name );
        // Search lead current vector 
        if (it != outputManager_.getBranchVarsNodeMap().end())
        {
          int index = (*it).second;
          //new_op = new BranchDataCurrentOp(name, index);
            // need to support the following for other lead current usage 
          if (param_tag == "IR" )
          {
            new_op = new BranchDataCurrentRealOp(name, index);
          }
          else if (param_tag == "II" )
          {
            new_op = new BranchDataCurrentImaginaryOp(name, index);
          }
          else if (param_tag == "IM" )
          {
            new_op = new BranchDataCurrentMagnitudeOp(name, index);
          }
          else if (param_tag == "IP" )
          {
            if (outputManager_.getPhaseOutputUsesRadians())
	    {
              new_op = new BranchDataCurrentPhaseRadOp(name, index);
            }
            else
            {
              new_op = new BranchDataCurrentPhaseDegOp(name, index);
            }
          }
          else if (param_tag == "IDB" )
          {
            new_op = new BranchDataCurrentDecibelsOp(name, index);
          }
          else // IC, IE, IB
          {
            new_op = new BranchDataCurrentOp(name, index);
          }
        }
      } 
 
      if( !new_op )
      {
        // Search store
        NodeNameMap::const_iterator it = outputManager_.getStoreNodeMap().find(store_name);
        if (it != outputManager_.getStoreNodeMap().end())
        {
          int index = (*it).second;
          if (param_tag == "IR" )
          {
            new_op = new StoreRealOp(name, index);
          }
          else if (param_tag == "II" )
          {
            new_op = new StoreImaginaryOp(name, index);
          }
          else if (param_tag == "IM" )
          {
            new_op = new StoreMagnitudeOp(name, index);
          }
          else if (param_tag == "IP" )
          {
            if (outputManager_.getPhaseOutputUsesRadians())
	    {
              new_op = new StorePhaseRadOp(name, index);
            }
            else
            {
              new_op = new StorePhaseDegOp(name, index);
            }
          }
          else if (param_tag == "IDB" )
          {
            new_op = new StoreDecibelsOp(name, index);
          }
          else // IC, IE, IB
          {
            new_op = new StoreOp(name, index);
          }
        }
      }
 
      if (new_op)
      {
        new_op->addArg(args[0]);
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  const Analysis::AnalysisManager &    analysisManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::PowerVariableOpBuilder 
// Purpose       : Builds various operators" for power calculations.  See 
//                 the list below in registerCreateFunctions(). This includes  
//                 includes the operators for P() and W().  Note that
//                 different operators are needed for each of the various 
//                 multi-terminal devices like J, Q, M, T and Z.
// Special Notes :
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct PowerVariableOpBuilder : public Util::Op::Builder
{

  PowerVariableOpBuilder(const OutputMgr & output_manager,
                         const Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)
  {}

  virtual ~PowerVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<BranchDataCurrentOp>();
    builder_manager.addCreateFunction<BranchDataCurrentRealOp >();
    builder_manager.addCreateFunction<BranchDataCurrentImaginaryOp>();
    builder_manager.addCreateFunction<BranchDataCurrentMagnitudeOp>();
    builder_manager.addCreateFunction<BranchDataCurrentPhaseDegOp>();
    builder_manager.addCreateFunction<BranchDataCurrentPhaseRadOp>();
    builder_manager.addCreateFunction<BranchDataCurrentDecibelsOp>();
    builder_manager.addCreateFunction<BranchDataPosNegPowerOp>();
    builder_manager.addCreateFunction<BranchDataBJTPowerOp>();
    builder_manager.addCreateFunction<BranchDataJFETPowerOp>();
    builder_manager.addCreateFunction<BranchDataMESFETPowerOp>();
    builder_manager.addCreateFunction<BranchDataMOSFETPowerOp>();
    builder_manager.addCreateFunction<BranchDataTRAPowerOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if ((param_tag == "P") || (param_tag == "W"))
    {
      // P() and W() are not supported for either .AC or .NOISE analyses for any device type.
      if ( !new_op && (analysisManager_.getNoiseFlag() || analysisManager_.getACFlag()) )
      {
        Report::UserError0() << "P() and W() are not supported for .AC and .NOISE analyses";
        return new_op;
      }

      // This check may no longer be needed, now that P() and W() work in a .MEASURE statement.
      // Leaving it in for now.
      if (args.size() == 0)
      {
        Report::UserFatal0() << "Error building operator: " << name << " missing variable name";
      }
      // Node name could be circuit_context:DeviceTypeDeviceName while internally it should be DeviceType:circuit_context:DeviceName.
      std::string modifiedName = Util::xyceDeviceNameToSpiceName(args[0]);
            
      // could be a device lead current "DEV_I" or a branch current.
      // so we don't have to duplicate solution vars(branch currents) in the
      // store vector, look for each type.
      bool param_func = std::find(func_names, func_names + sizeof(func_names)/sizeof(func_names[0]), "I") != func_names + sizeof(func_names)/sizeof(func_names[0]);
      std::string store_name = modifiedName + ":BRANCH_D" ; // + (param_func ? "I" : "I");  // if it is in the state/store vec.
      // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
      // we try to find devices based on store_name in the following blocks of code,
      // so do this modification now.
      std::string::size_type space = store_name.find_first_of(" ");
      if (space != std::string::npos)
      {
        if (space == 4 && store_name.substr(0, 4) == "YPDE")
        {
          store_name.replace(4, 1, ":");
        }
      }

      // Search store
      if( modifiedName[0] == 'Q' )
      {
        // BJT type.  
        // this should be encapsulated into the BJT class as some point
        NodeNameMap::const_iterator itB = outputManager_.getBranchVarsNodeMap().find(store_name+"B");
        NodeNameMap::const_iterator itE = outputManager_.getBranchVarsNodeMap().find(store_name+"E");
        NodeNameMap::const_iterator itC = outputManager_.getBranchVarsNodeMap().find(store_name+"C");
        NodeNameMap::const_iterator itS = outputManager_.getBranchVarsNodeMap().find(store_name+"S");
        if ( (itB != outputManager_.getBranchVarsNodeMap().end()) &&
             (itE != outputManager_.getBranchVarsNodeMap().end()) &&
             (itC != outputManager_.getBranchVarsNodeMap().end()))
		{
		  int indexB = (*itB).second;
		  int indexC = (*itC).second;
		  int indexE = (*itE).second;
		  int indexS = -1;
                  if (itS != outputManager_.getBranchVarsNodeMap().end())
                    indexS = (*itS).second;
		  
		  new_op = new BranchDataBJTPowerOp(name, indexB, indexC, indexE, indexS);
		}
      
      }
      else if( modifiedName[0] == 'J' )
      {
        // JFET type.  
        // this should be encapsulated into the JFET class as some point
        NodeNameMap::const_iterator itD = outputManager_.getBranchVarsNodeMap().find(store_name+"D");
        NodeNameMap::const_iterator itG = outputManager_.getBranchVarsNodeMap().find(store_name+"G");
        NodeNameMap::const_iterator itS = outputManager_.getBranchVarsNodeMap().find(store_name+"S");
        if ( (itD != outputManager_.getBranchVarsNodeMap().end()) &&
             (itG != outputManager_.getBranchVarsNodeMap().end()) &&
             (itS != outputManager_.getBranchVarsNodeMap().end()))
		{
		  int indexD = (*itD).second;
		  int indexG = (*itG).second;
		  int indexS = (*itS).second;
 		  new_op = new BranchDataJFETPowerOp(name, indexD, indexG, indexS);
		}
      
      }
      else if( modifiedName[0] == 'M' )
      {
        // MOSFET type.  
        // this should be encapsulated into the MOSFET class as some point
        NodeNameMap::const_iterator itD = outputManager_.getBranchVarsNodeMap().find(store_name+"D");
        NodeNameMap::const_iterator itG = outputManager_.getBranchVarsNodeMap().find(store_name+"G");
        NodeNameMap::const_iterator itS = outputManager_.getBranchVarsNodeMap().find(store_name+"S");
        NodeNameMap::const_iterator itB = outputManager_.getBranchVarsNodeMap().find(store_name+"B");
        if ( (itD != outputManager_.getBranchVarsNodeMap().end()) &&
             (itG != outputManager_.getBranchVarsNodeMap().end()) &&
             (itS != outputManager_.getBranchVarsNodeMap().end()))
	{
	  int indexD = (*itD).second;
	  int indexG = (*itG).second;
	  int indexS = (*itS).second;
	  int indexB = -1;
          if (itB != outputManager_.getBranchVarsNodeMap().end())
	  {
            indexB = (*itB).second;
	  }	  
	  new_op = new BranchDataMOSFETPowerOp(name, indexD, indexG, indexS, indexB);
        }
      }
      else if( modifiedName[0] == 'T')
      {
        // T Device (TRA or Lossless Transmission Line)
        NodeNameMap::const_iterator it1 = outputManager_.getBranchVarsNodeMap().find(store_name+"1");
        NodeNameMap::const_iterator it2 = outputManager_.getBranchVarsNodeMap().find(store_name+"2");
        if ( (it1 != outputManager_.getBranchVarsNodeMap().end()) &&
             (it2 != outputManager_.getBranchVarsNodeMap().end())) 
	{
	  int index1 = (*it1).second;
	  int index2 = (*it2).second;
 	  new_op = new BranchDataTRAPowerOp(name, index1, index2);
	}
      }
      else if( modifiedName[0] == 'Z' )
      {
        // MESFET type.  
        // this should be encapsulated into the MESFET class as some point
        NodeNameMap::const_iterator itD = outputManager_.getBranchVarsNodeMap().find(store_name+"D");
        NodeNameMap::const_iterator itG = outputManager_.getBranchVarsNodeMap().find(store_name+"G");
        NodeNameMap::const_iterator itS = outputManager_.getBranchVarsNodeMap().find(store_name+"S");
        if ( (itD != outputManager_.getBranchVarsNodeMap().end()) &&
             (itG != outputManager_.getBranchVarsNodeMap().end()) &&
             (itS != outputManager_.getBranchVarsNodeMap().end()))
	{
	  int indexD = (*itD).second;
	  int indexG = (*itG).second;
	  int indexS = (*itS).second;
 	  new_op = new BranchDataMESFETPowerOp(name, indexD, indexG, indexS);
	}
      }
      else
      {
        // default for devices for which power = I*V
	NodeNameMap::const_iterator it = outputManager_.getBranchVarsNodeMap().find(store_name);
	if (it != outputManager_.getBranchVarsNodeMap().end())
	{
	  int index = (*it).second;
          new_op = new BranchDataPosNegPowerOp(name, index);
        }
      }

      if (new_op)
      {
        new_op->addArg(args[0]);
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  const Analysis::AnalysisManager &   analysisManager_;
};

//--------------------------------------------------------------------------
// Structure     : Util::Op::Builder::RFparamsVariableOpBuilder
// Purpose       : Builds various operators for RF parameter output.
//                 See the list below in registerCreateFunctions().  This
//                 includes the operators for S(), SR(), SI(), SM(), SP()
//                 and SDB(). It also includes the corresponding Y and Z
//                 parameter output.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 07/01/19
//--------------------------------------------------------------------------
struct RFparamsVariableOpBuilder : public Util::Op::Builder
{
  RFparamsVariableOpBuilder(const OutputMgr & output_manager,
                            Analysis::AnalysisManager & analysis_manager)
    : outputManager_(output_manager),
      analysisManager_(analysis_manager)
  {}

  virtual ~RFparamsVariableOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<RFparamsOp>();
    builder_manager.addCreateFunction<RFparamsRealOp>();
    builder_manager.addCreateFunction<RFparamsImaginaryOp>();
    builder_manager.addCreateFunction<RFparamsMagnitudeOp>();
    builder_manager.addCreateFunction<RFparamsPhaseDegOp>();
    builder_manager.addCreateFunction<RFparamsPhaseRadOp>();
    builder_manager.addCreateFunction<RFparamsDecibelsOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    std::vector<std::string> args;
    std::string name;
    parameterNameAndArgs(name, args, it);

    if ( ((param_tag[0] == 'S') || (param_tag[0] == 'Y') || (param_tag[0] == 'Z')) && args.size() == 2)
    {
      if (!analysisManager_.getACLinFlag())
      {
        Report::UserError0() << "S(), Y() and Z() operators only supported for .LIN analyses";
        return new_op;
      }

      int rowIdx= atoi(args[0].c_str());
      int colIdx= atoi(args[1].c_str());
      if ( rowIdx < 1  || colIdx < 1 )
      {
        Report::UserError0() << "Indices for S(), Y() and Z() operators must be > 0";
        return new_op;
      }

      const std::string type = name.substr(0,1);
      if ( (param_tag == "S") || (param_tag == "Y") || (param_tag == "Z"))
      {
        new_op = new RFparamsOp(name, type, rowIdx, colIdx);
      }
      else if ( (param_tag == "SR") || (param_tag == "YR") || (param_tag == "ZR") )
      {
        new_op = new RFparamsRealOp(name, type, rowIdx, colIdx);
      }
      else if ( (param_tag == "SI") || (param_tag == "YI") || (param_tag == "ZI") )
      {
        new_op = new RFparamsImaginaryOp(name, type, rowIdx, colIdx);
      }
      else if ( (param_tag == "SM") || (param_tag == "YM") || (param_tag == "ZM") )
      {
        new_op = new RFparamsMagnitudeOp(name, type, rowIdx, colIdx);
      }
      else if ( (param_tag == "SP" ) || (param_tag == "YP") || (param_tag == "ZP") )
      {
        if (outputManager_.getPhaseOutputUsesRadians())
	{
          new_op = new RFparamsPhaseRadOp(name, type, rowIdx, colIdx);
        }
        else
	{
          new_op = new RFparamsPhaseDegOp(name, type, rowIdx, colIdx);
        }
      }
      else if ( (param_tag == "SDB" ) || (param_tag == "YDB" ) || (param_tag == "ZDB") )
      {
        new_op = new RFparamsDecibelsOp(name, type, rowIdx, colIdx);
      }

      if (new_op)
      {
        new_op->addArgs(args.begin(), args.end());
        analysisManager_.setRFParamsRequested(type);
      }
    }

    return new_op;
  }

private:
  const OutputMgr &     outputManager_;
  Analysis::AnalysisManager &    analysisManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::ExpressionOpBuilder 
// Purpose       : Makes an ExpressionOp or ConstantOp.
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct ExpressionOpBuilder : public Util::Op::Builder
{
  ExpressionOpBuilder(Parallel::Machine comm, const OutputMgr & output_manager)
    : comm_(comm),
      outputManager_(output_manager)
  {}

  virtual ~ExpressionOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<ExpressionOp>();
    builder_manager.addCreateFunction<Util::Op::ConstantOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    int param_type = (*it).getType();

    if (Util::hasExpressionTag(*it))
    {
      new_op = new ExpressionOp(param_tag, param_tag, comm_, outputManager_);
    }

    return new_op;
  }

private:
  const Parallel::Machine       comm_;
  const OutputMgr &             outputManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::MeasurementOpBuilder 
// Purpose       : Builds a MeasurementOp.
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct MeasurementOpBuilder : public Util::Op::Builder
{
  MeasurementOpBuilder(const Measure::Manager & measure_manager)
    : measureManager_(measure_manager)
  {}

  virtual ~MeasurementOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<MeasureOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    const Measure::Base *measure = measureManager_.find(param_tag);
    if (measure)
    {
      new_op = new MeasureOp(param_tag, *measure);
    }

    return new_op;
  }

private:
  const Measure::Manager &             measureManager_;
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::CircuitIndexOpBuilder 
// Purpose       : Builds a CircuitIndexOp.
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct CircuitIndexOpBuilder : public Util::Op::Builder
{
  CircuitIndexOpBuilder()
  {}

  virtual ~CircuitIndexOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<CurrentIndexOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();

    if (param_tag == "INDEX") 
    {
      new_op  = new CurrentIndexOp(param_tag);
    }

    return new_op;
  }
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::SensitivityOpBuilder 
// Purpose       : Builds various "sensitivity operators".  See the
//                 list below in registerCreateFunctions(). This includes 
//                 the operators for both direct and adjoint sensitivities.
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct SensitivityOpBuilder : public Util::Op::Builder
{
  SensitivityOpBuilder()
  {}

  virtual ~SensitivityOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SensitivityObjFunctionOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectScaledOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointScaledOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    // Certain variables were added to this list with the temporary tag "SENS",
    // in the OutputMgr::registerSens function.  They also contain 
    // OP information, which can be used here to construct the proper name
    // for each derivative.
    if (param_tag == "SENS")
    {
      std::string function;
      std::string parameter;
      ptrdiff_t type;
      int index;

      Util::Marshal min(param_string);

      min >> function >> parameter >> type >> index;

      std::string name = "d" + function + "/d(" + parameter + ")";
      if (type == Util::Op::identifier<SensitivityObjFunctionOp>())
      {
        new_op = new SensitivityObjFunctionOp(function, index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpDirectOp>())
      {
        new_op = new SensitivitydOdpDirectOp(name + "_Dir", index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpDirectScaledOp>())
      {
        new_op = new SensitivitydOdpDirectScaledOp(name + "_Dir_scaled", index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpAdjointOp>())
      {
        new_op = new SensitivitydOdpAdjointOp(name + "_Adj", index);
      }
      if (type == Util::Op::identifier<SensitivitydOdpAdjointScaledOp>())
      {
        new_op = new SensitivitydOdpAdjointScaledOp(name + "_Adj_scaled", index);
      }
    }

    return new_op;
  }
};

//-------------------------------------------------------------------------- 
// Structure     : Util::Op::Builder::TransientAdjointOpBuilder 
// Purpose       : Builder for various "sensitivity operators".  See
//                 the list below in registerCreateFunctions().
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
struct TransientAdjointOpBuilder : public Util::Op::Builder
{
  TransientAdjointOpBuilder()
  {}

  virtual ~TransientAdjointOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<SensitivityObjFunctionOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectOp>();
    builder_manager.addCreateFunction<SensitivitydOdpDirectScaledOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointOp>();
    builder_manager.addCreateFunction<SensitivitydOdpAdjointScaledOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const 
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    // Certain variables were added to this list with the temporary tag "SENS",
    // in the OutputMgr::registerSens function.  They also contain 
    // OP information, which can be used here to construct the proper name
    // for each derivative.
    if (param_tag == "SENS")
    {
      std::string function;
      std::string parameter;
      ptrdiff_t type;
      int index;

      Util::Marshal min(param_string);

      min >> function >> parameter >> type >> index;

      std::string name = "d" + function + "/d(" + parameter + ")";
      if (type == Util::Op::identifier<SensitivityObjFunctionOp>())
      {
        new_op = new SensitivityObjFunctionOp(function, index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpDirectOp>())
      {
        new_op = new SensitivitydOdpDirectOp(name + "_Dir", index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpDirectScaledOp>())
      {
        new_op = new SensitivitydOdpDirectScaledOp(name + "_Dir_scaled", index);
      }
      else if (type == Util::Op::identifier<SensitivitydOdpAdjointOp>())
      {
        new_op = new SensitivitydOdpAdjointOp(name + "_Adj", index);
      }
      if (type == Util::Op::identifier<SensitivitydOdpAdjointScaledOp>())
      {
        new_op = new SensitivitydOdpAdjointScaledOp(name + "_Adj_scaled", index);
      }
    }

    return new_op;
  }
};

//-------------------------------------------------------------------------- 
// Function      : registerOpBuilders 
// Purpose       : register various builders with the op_builder_manager
// Special Notes : Error handling for some Ops depends on the analysis mode.
//                 For those Op types, a const reference to the analysis_manager
//                 is passed in.
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
void registerOpBuilders(Util::Op::BuilderManager &op_builder_manager, Parallel::Machine comm, OutputMgr &output_manager,
                        Analysis::AnalysisManager &analysis_manager)
{
  op_builder_manager.addBuilder(new CircuitTemperatureOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitTimeOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitNoiseContOpBuilder(output_manager,analysis_manager));
  op_builder_manager.addBuilder(new CircuitOutputNoiseOpBuilder(output_manager,analysis_manager));
  op_builder_manager.addBuilder(new CircuitInputNoiseOpBuilder(output_manager,analysis_manager));
  op_builder_manager.addBuilder(new CircuitFrequencyOpBuilder(output_manager));
  op_builder_manager.addBuilder(new CircuitIndexOpBuilder());
  op_builder_manager.addBuilder(new SensitivityOpBuilder());
  op_builder_manager.addBuilder(new TransientAdjointOpBuilder());
  op_builder_manager.addBuilder(new ExpressionOpBuilder(comm, output_manager));
//  op_builder_manager.addBuilder(new StepSweepOpBuilder(output_manager));
  op_builder_manager.addBuilder(new DCSweepOpBuilder(output_manager));
  op_builder_manager.addBuilder(new DCSweepCurrentValueOpBuilder(output_manager));
  op_builder_manager.addBuilder(new StepNumOpBuilder(output_manager));
  op_builder_manager.addBuilder(new InternalVariableOpBuilder(output_manager,comm));
  op_builder_manager.addBuilder(new VoltageVariableOpBuilder(output_manager,comm));
  op_builder_manager.addBuilder(new CurrentVariableOpBuilder(output_manager,analysis_manager));
  op_builder_manager.addBuilder(new PowerVariableOpBuilder(output_manager,analysis_manager));
  op_builder_manager.addBuilder(new RFparamsVariableOpBuilder(output_manager,analysis_manager));
}

//-------------------------------------------------------------------------- 
// Function      : registerOpBuilders 
// Purpose       : register MeasurementOpBuilder with the op_builder_manager
// Special Notes : 
// Creator       : Dave Baur
// Creation Date : 08/04/14 
//--------------------------------------------------------------------------
void registerOpBuilders(Util::Op::BuilderManager &op_builder_manager, Parallel::Machine comm, Measure::Manager &measure_manager)
{
  op_builder_manager.addBuilder(new MeasurementOpBuilder(measure_manager));
}

} // namespace IO
} // namespace Xyce
