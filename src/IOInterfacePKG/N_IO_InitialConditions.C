//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iterator>

#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LAS_Vector.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::InitialConditionsManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
InitialConditionsManager::InitialConditionsManager(
  const std::string &         netlist_filename)
  : netlistFilename_(netlist_filename),
    outputOnceAlreadyFlag_(false)
{}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::registerSave
// Purpose       : registers set of variables to set for .SAVE.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/11/07
//-----------------------------------------------------------------------------
bool InitialConditionsManager::registerSave(const Util::OptionBlock & option_block)
{
  if (DEBUG_IO)
    Xyce::dout() << " SAVE Registered::" << std::endl;

  icData_.saveFlag_ = true;

  ExtendedString sval("");

  for (Util::ParamList::const_iterator it = option_block.begin(), 
      end = option_block.end(); it != end; ++it)
  {
    if (DEBUG_IO)
      Xyce::dout() << "(*it).tag = " << (*it).tag() << std::endl;

    if ((*it).tag() == "TYPE")
    {
      sval = (*it).stringValue();
      sval.toUpper();

      if (sval == "NODESET" || sval == ".NODESET")
      {
        icData_.saveFileType_ = ".NODESET";
      }
      else if (sval == "IC" || sval == ".IC")
      {
        icData_.saveFileType_ = ".IC";
      }
      else
      {
        Report::UserWarning0() 
          << "Unrecognized type specified on .SAVE command.  Defaulting to .NODESET";
      }
    }
    else if ((*it).tag() == "FILE")
    {
      icData_.saveOutputFile_ = (*it).stringValue();
    }
    else if ((*it).tag() == "TIME")
    {
      icData_.saveTime_ = (*it).getImmutableValue<double>();
    }
    else if ((*it).tag() == "LEVEL")
    {
      sval = (*it).stringValue();
      sval.toUpper();

      if (sval == "ALL")
      {
        // do nothing
      }
      else if (sval == "NONE")
      {
        // none means don't output anything.  Pretend the .SAVE line isn't in the netlist.
        icData_.saveFlag_ = false;
        icData_.saveFileLevel_ = "NONE";
      }
      else if (sval == "TOP")
      {
        Report::UserWarning0() << "LEVEL=TOP in .SAVE line not supported.  Defaulting to ALL.";
      }
      else
      {
        Report::UserWarning0() << "Unrecognized LEVEL " << sval << " specified in .SAVE command.  Defaulting to ALL.";
      }
    }
    else
    {
      Report::UserWarning0() << "Parameter " << (*it).tag() << " not recognized in .SAVE command";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::registerIC
// Purpose       : registers set of variables to set for .IC or .DCVOLT
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool InitialConditionsManager::registerIC(const Util::OptionBlock & OB)
{
  icData_.ICflag_ = true;

  // save the ICblock_.  This will be needed later.
  // We allow multiple .IC blocks in the netlist, so there needs to be an
  // STL vector of these.
  ICblockVec_.push_back(OB);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::registerNodeSet
// Purpose       : registers set of variables to set for .NODESET
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool InitialConditionsManager::registerNodeSet(const Util::OptionBlock & OB)
{
 icData_.nodesetflag_ = true;

  // save the nodesetblock_.  This will be needed later.
  // We allow multiple .nodeset blocks in the netlist, so there needs to be an
  // STL vector of these.
  nodesetblockVec_.push_back(OB);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::getICData
// Purpose       : Provide OP data for LOCA
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/04/06
//-----------------------------------------------------------------------------
InitialConditionsData::NodeNamePairMap & InitialConditionsManager::getICData(
  int & op_found, int & icType)
{
  op_found = icData_.op_found_;
  icType = icData_.icType_;

  return icData_.opData_;
}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::setupIC_or_NODESET
// Purpose       :
// Special Notes : Assumes that all_nodes and alias_nodes are set up correctly.
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool setupIC_or_NODESET(
  Parallel::Machine                             comm,
  const NodeNameMap &                           all_nodes,
  const AliasNodeMap &                          alias_nodes,
  Linear::Vector &                              nextSolnVec,
  Linear::Vector &                              flagSolnVec,
  int                                           icType,
  std::vector<Util::OptionBlock> &              initBlockVec,
  InitialConditionsData::NodeNamePairMap &      op_data,
  int &                                         op_found_,
  int &                                         total_soln_)
{
  int lid = 0;
  bool success = false;

  int totalNodes = all_nodes.size();

  nextSolnVec.putScalar(0.0);
  flagSolnVec.putScalar(-1.0);

  for (NodeNameMap::const_iterator it = all_nodes.begin(),
      nodes_end = all_nodes.end(); it != nodes_end ; ++it)
  {
    flagSolnVec[(*it).second] = 0;
  }

  std::set<std::string> notFoundInCkt;
  std::set<std::string> matched;
  for (std::vector<Util::OptionBlock>::const_iterator it = initBlockVec.begin(), 
      end = initBlockVec.end(); it != end; ++it)
  {
    // param loop
    for (Util::ParamList::const_iterator it_param = (*it).begin(), 
        endPar = (*it).end(); it_param != endPar; ++it_param)
    {
      // Check if this .IC or .NODESET block is to be filled out by subcircuits.
      if ( (*it_param).tag() == "SUBCKT" )
        break;

      // the first tag is always "V".  At this point I variables are not supported,
      // and there is an error trap for them in the extractICData function.
      ++it_param;
      std::string node = (*it_param).tag();
      ++it_param;

      double value = 0.0;
      ExtendedString strValueNodeTag((*it_param).tag());
      strValueNodeTag.toUpper();
      if (strValueNodeTag=="VALUE")
      {
        value = (*it_param).getImmutableValue<double>();
      }
      else
      {
        Report::UserFatal0() << "Problems processing initial condition for node "
                             << node << " to " << strValueNodeTag << " value";
      }

      int found = false; // will be set to true if the node is found in all_nodes
      int aliasNodeFound = false; // will be set to true if the node is found in alias_nodes

      // the real node name (e.g., 2) of a subcircuit interface alias node (e.g., X1:A)
      // found in the alias_nodes map.
      std::string realNodeName;

      NodeNameMap::const_iterator iterCI = all_nodes.find(node);
      if (iterCI != all_nodes.end())
      {
        lid = iterCI->second;
        nextSolnVec[lid] = value;
        flagSolnVec[lid] = 1;
        found = true;
      }
      else
      {
        // See SRN Bug 1962 for more details but, for subcircuit interface nodes, 
        // alias_nodes would have key-value pairs like <"X1:A","2">.  So, we
        // want the "real node name" (say 2) for X1:A.  Based on the "real 
        // node name", we can then look up its lid on this processor. 
        AliasNodeMap::const_iterator iterAN = alias_nodes.find(node);
        if (iterAN != alias_nodes.end())
	{
          // to be conservative, check that the "real node name" actually exists
          // in the all_nodes map.
          NodeNameMap::const_iterator iterNode = all_nodes.find(iterAN->second);
          if (iterNode != all_nodes.end()) 
          {
            realNodeName = iterNode->first;
            lid = iterNode->second;
	    nextSolnVec[lid] = value;
            flagSolnVec[lid] = 1;
            aliasNodeFound = true;
	  }
        }
      }

      if (DEBUG_IO) 
      {
        if (found || aliasNodeFound)
        {
          Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
          Xyce::dout()
            << icType << " found, and set: V(" << node << "):  solution["<<lid<<"] = " << value << std::endl;
        }
      }

      Parallel::AllReduce(comm, MPI_SUM, &found, 1);
      Parallel::AllReduce(comm, MPI_SUM, &aliasNodeFound, 1);

      if (found || aliasNodeFound)
      {
        success = true;

        if ( found && (matched.find(node) == matched.end()) )
        {
          matched.insert(node);
        }
        else if ( aliasNodeFound && (matched.find(realNodeName) == matched.end()) )
        {
          matched.insert(realNodeName);
        }
      }
      else
      {
        notFoundInCkt.insert(node);
      }
    }
  }

  nextSolnVec.importOverlap();
  flagSolnVec.importOverlap();

  op_found_ = matched.size();
  total_soln_ = totalNodes;

  std::set<std::string>::iterator matched_i = matched.begin();
  std::set<std::string>::iterator matched_end = matched.end();
  for (; matched_i != matched_end ; ++matched_i)
  {
    // only add this to the opData map if it is local to the processor.
    // The allNodes object only contains local stuff.  Otherwise there
    // isn't much point.
    NodeNameMap::const_iterator iterCI = all_nodes.find(*matched_i);
    if (iterCI != all_nodes.end())
    {
      op_data[*matched_i].first = iterCI->second;
      op_data[*matched_i].second = nextSolnVec[op_data[*matched_i].first];
    }
  }

  if (DEBUG_IO)
  {
    // Identify nodes in the circuit which were not set.
    std::set<std::string> local;
    for (NodeNameMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end ; ++it)
    {
      if (matched.find((*it).first) == matched.end())
        local.insert((*it).first);
    }

    Util::Marshal mout;
    mout << local;

    std::vector<std::string> dest;
    Parallel::GatherV(comm, 0, mout.str(), dest);

    std::set<std::string> notSpecified;
    if (Parallel::rank(comm) == 0) 
    {
      for (int p = 0; p < Parallel::size(comm); ++p) 
      {
        Util::Marshal min(dest[p]);

        std::vector<std::string> x;
        min >> x;
        std::copy(x.begin(), x.end(), std::inserter(notSpecified, notSpecified.begin()));
      }
    }
  }

  if (Parallel::rank(comm) == 0)
  {
    if (DEBUG_IO)
    {
      if (totalNodes > matched.size())
      {
        Xyce::dout() << icType << ":  Initialized " << matched.size() 
          << " of a possible " << totalNodes << " nodes" << std::endl;
      }
      else
      {
        Xyce::dout() << icType << ":  All " << totalNodes << " nodes initialized" 
          << std::endl;
      }
    }

    if (notFoundInCkt.size() > 0)
    {
      Xyce::lout() << "Netlist warning: Initial conditions specified at nodes not present in circuit." << std::endl
                   << "May be error in .IC or .NODESET line. Ignoring nodes:" << std::endl;
      for (std::set<std::string>::iterator nm_i = notFoundInCkt.begin(), 
                end = notFoundInCkt.end(); nm_i != end; ++nm_i)
      {
        lout() << *nm_i << std::endl;
      }
    }

    // if (DEBUG_IO)
    // {
    //   // Note: for a typical .IC line, this list of unspecified nodes will be
    //   // a very long list.  So, don't output except in debug mode.
    //   if (notSpecified.size() > 0)
    //   {
    //     dout() << icType << ":  Nodes present in circuit, but not specified on ." 
    //            << icType << " line(ignoring):" << std::endl;

    //     for (std::set<std::string>::iterator nm_i = notSpecified.begin(), 
    //               end = notSpecified.end(); nm_i != end; ++nm_i)
    //     {
    //       dout() << *nm_i << std::endl;
    //     }
    //   }
    // }
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : InitialConditionsManager::outputIC_or_NODESET
// Purpose       : Outputs either an *.ic file, either in .ic or .nodeset format.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
void outputIC_or_NODESET(
  Parallel::Machine             comm,
  std::ofstream &               os,
  const std::string &           saveFileType_,
  InitialConditionsData::NodeNamePairMap               all_nodes,
  const Linear::Vector &          solution)
{
  // Collect the local string for all the nodes owned by each processor.
  std::ostringstream localString;

  for (InitialConditionsData::NodeNamePairMap::iterator it = all_nodes.begin(), 
      name_end = all_nodes.end(); it != name_end ; ++it)
  {
    ExtendedString tmpName((*it).first);
        tmpName.toUpper();
        std::string::size_type pos = tmpName.rfind("BRANCH");
    if (pos == std::string::npos)
    {
      int index =(*it).second.first;
      (*it).second.second = solution[index];
      localString << saveFileType_ << " V(" << (*it).first << ") = " << (*it).second.second << std::endl;
    }
  }

  // Gather the strings for each processor.
  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, localString.str(), dest);

  // Print the string, one at a time.
  if (Parallel::rank(comm) == 0) 
  {
    InitialConditionsData::NodeNamePairMap global;
    for (int p = 0; p < Parallel::size(comm); ++p) 
    {
      os << dest[p];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setupInitialConditions
// Purpose       : This function is an umbrella function, under which the
//                 various types of initial conditions are potentially set up.
//                 These include, but aren't limited to .IC and .NODESET.
//
//                 In addition to setting these things up, it does some
//                 nominal checking to make sure that more than one initial
//                 condition hasn't been specified.
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool InitialConditionsManager::setupInitialConditions(
  Parallel::Machine             comm,
  const NodeNameMap &           allNodes_,
  const AliasNodeMap &          alias_nodes,
  Linear::Vector &                solnVec,
  Linear::Vector &                flagVec)
{
    bool dotIC = false;
    bool NODESET = false;
    icData_.icType_ = InitialConditionsData::IC_TYPE_UNDEFINED;

    if (icData_.ICflag_ && icData_.nodesetflag_)
    {
      Report::UserFatal0()
        << "Cannot set both .IC and .NODESET simultaneously.";
    }
    else if (icData_.ICflag_)
    {
      // check for .IC
      if (icData_.ICflag_)
      {
        dotIC = setupIC_or_NODESET(comm, allNodes_, alias_nodes, solnVec, flagVec, 
            InitialConditionsData::IC_TYPE_IC, ICblockVec_, 
            icData_.opData_, icData_.op_found_, icData_.total_soln_);
      }
      if (dotIC)
      {
        icData_.icType_ = InitialConditionsData::IC_TYPE_IC;
      }
    }
    else if (icData_.nodesetflag_)
    {
      // check for .NODESET
      if (icData_.nodesetflag_)
      {
        NODESET = setupIC_or_NODESET(comm, allNodes_, alias_nodes, solnVec, flagVec, 
            InitialConditionsData::IC_TYPE_NODESET, nodesetblockVec_, 
            icData_.opData_, icData_.op_found_, icData_.total_soln_);
      }
      if (NODESET)
      {
        icData_.icType_ = InitialConditionsData::IC_TYPE_NODESET;
      }
    }

    return dotIC || NODESET;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputDCOP
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void InitialConditionsManager::outputDCOP(
  Parallel::Machine             comm,
  const NodeNameMap &           all_nodes,
  const Linear::Vector &        solution)
{
  if (outputOnceAlreadyFlag_)
  {
    return;
  }

  InitialConditionsData::NodeNamePairMap        solution_node_names;

  for (NodeNameMap::const_iterator it = all_nodes.begin(), end = all_nodes.end(); it != end; ++it)
    solution_node_names[(*it).first] = std::make_pair((*it).second, 0.0);

  if (icData_.saveFlag_)
  {
    std::ofstream saveOutputStream;

    if (Parallel::rank(comm) == 0)
    {
      std::string outputFileName;
      if (icData_.saveOutputFile_  == "")
      {
        outputFileName = netlistFilename_ + ".ic";
      }
      else
      {
        outputFileName = icData_.saveOutputFile_;
      }

      saveOutputStream.open(outputFileName.c_str());

      if (!saveOutputStream.good())
      {
        Report::UserFatal() << "Cannot create Save File "  << outputFileName;
      }
    }
    Report::safeBarrier(comm);

    // check for .IC
    outputIC_or_NODESET(comm, saveOutputStream, icData_.saveFileType_, solution_node_names, solution);
    if (Parallel::rank(comm) == 0)
    {
      saveOutputStream.close();
    }
  }

  outputOnceAlreadyFlag_ = true;
}

namespace {

//-----------------------------------------------------------------------------
// Function      : extractICData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool extractICData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("IC", 
      Util::OptionBlock::ALLOW_EXPRESSIONS, 
      netlist_filename, 
      parsed_line[0].lineNumber_);

  const std::string& subcktName = circuit_block.getName();

  bool addBlock = true;
  int numFields = parsed_line.size();

  if ( numFields < 2 )
  {
    Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Ignored .IC and/or .DCVOLT, no arguments provided.";
    addBlock = false;
  }

  if (DEBUG_IO)
  {
    for (int ieric=0;ieric<parsed_line.size();++ieric)
    {
      Xyce::dout() << "parsed_line["<<ieric<<"] = " << parsed_line[ieric].string_ << std::endl;
    }
  }

  std::ostringstream msg;
  Util::Param parameter;
  ExtendedString field("");
  int p_err=0;
  int parameterStartPos  = 1;
  int position = parameterStartPos;

  // .IC can have two formats:
  //   (1)   .ic  V(a)=1.0   V(2)=2.0
  //   (2)   .ic  a 1.0   2 2.0
  // Check for case 1 by looking for = in the line.
  bool formatOne(false);

  while ( position < numFields )
  {
    if (parsed_line[position].string_ == "=")
    {
      formatOne = true;
    }
    ++position;
  }

  position = parameterStartPos;

  if (formatOne)
  {
    while ( position < numFields )
    {
      if (position+5 < numFields && parsed_line[position+1].string_ == "(")
      {
        // current variables.  ERK: Note, originally I thought to have it
        // be possible to set both voltages and currents.  However, it
        // appears (I think) that other codes just set voltages.
        const std::string & posString = parsed_line[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg << "Unsupported current specification.";
            p_err = position;
          }
        }
        // voltage variables:
        else if( posString == "V" || posString == "v")
        {
          if( parsed_line[position+3].string_ == ")" )
          {
            if (subcktName != "")
            {
              parameter.setTag("SUBCKT");
              parameter.setVal(subcktName);
              option_block.addParam( parameter );
            }

            parameter.setTag("V");
            parameter.setVal( 1.0 );
            option_block.addParam( parameter );

            // node name:
            field = parsed_line[position+2].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            option_block.addParam( parameter );

            // Check that this is not a global node in a subcircuit.
            if ( field.substr(0,2) == "$G" && subcktName != "" )
            {
              Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
                << "Ignored .IC and/or .DCVOLT on global node " << field << " inside subcircuit " 
                << subcktName << ", move statement to global scope.";
              addBlock = false;
            }

            // value:
            field = parsed_line[position+5].string_;
            field.toUpper();
            parameter.setTag( "VALUE" );
            parameter.setVal( field );
            option_block.addParam( parameter );

            position += 6;
          }
          else if( parsed_line[position+6].string_ == "," )
          {
            msg << "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg << "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg << "Unrecognized .IC and/or .DCVOLT specification";
        p_err = position;
      }

      if (!msg.str().empty())
      {
        msg << " in .IC and/or .DCVOLT near ";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg << parsed_line[p_err].string_ + " ";
          ++p_err;
        }
        Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
          << msg.str();
      }
    }
  }
  else
  {
    while ( position < numFields )
    {
      if (subcktName != "")
      {
        parameter.setTag("SUBCKT");
        parameter.setVal(subcktName);
        option_block.addParam( parameter );
      }

      parameter.setTag("V");
      parameter.setVal( 1.0 );
      option_block.addParam( parameter );

      // node name:
      field = parsed_line[position].string_;
      field.toUpper();
      std::cout << "field = " << field << std::endl;      
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      // Check that this is not a global node in a subcircuit.
      if ( field.substr(0,2) == "$G" && subcktName != "" )
      {
        Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Ignored .IC and/or .DCVOLT on global node " << field << " inside subcircuit " 
          << subcktName << ", move statement to global scope.";
        addBlock = false; 
      }

      // value:
      field = parsed_line[position+1].string_;
      field.toUpper();
      parameter.setTag( std::string("VALUE") );
      parameter.setVal( field );
      option_block.addParam( parameter );

      position += 2;
    }
  }

  if (addBlock)
  {
    circuit_block.addOptions(option_block);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : extractNodeSetData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/07/2007
//-----------------------------------------------------------------------------
bool
extractNodeSetData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("NODESET", 
      Util::OptionBlock::ALLOW_EXPRESSIONS, 
      netlist_filename, 
      parsed_line[0].lineNumber_);

  const std::string& subcktName = circuit_block.getName();

  bool addBlock = true;
  int numFields = parsed_line.size();

  if ( numFields < 2 )
  {
    Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Ignored .NODESET, no arguments provided.";
    addBlock = false;
  }

  if (DEBUG_IO)
  {
    for (int ieric=0;ieric<parsed_line.size();++ieric)
    {
      Xyce::dout() << "parsed_line["<<ieric<<"] = " << parsed_line[ieric].string_ << std::endl;
    }
  }

  std::ostringstream msg;
  Util::Param parameter;
  ExtendedString field("");
  int p_err=0;
  int parameterStartPos  = 1;
  int position = parameterStartPos;

  // .NODESET can have two formats:
  //   (1)   .nodeset  V(a)=1.0   V(2)=2.0
  //   (2)   .nodeset  a 1.0   2 2.0
  // Check for case 1 by looking for = in the line.
  bool formatOne(false);

  while ( position < numFields )
  {
    if (parsed_line[position].string_ == "=")
    {
      formatOne = true;
    }
    ++position;
  }


  position = parameterStartPos;

  if (formatOne)
  {
    while ( position < numFields )
    {
      if (position+5 < numFields && parsed_line[position+1].string_ == "(")
      {
        // current variables:
        const std::string & posString = parsed_line[position].string_;

        if (posString == "I" || posString == "i"
         || (posString.size() == 2 &&
             (posString[0] == 'I' || posString[0] == 'i')))
        {
          {
            msg << "Unsupported current specification.";
            p_err = position;
          }
        }
        // voltage variables:
        else if( posString == "V" || posString == "v")
        {
          if( parsed_line[position+3].string_ == ")" )
          {
            if (subcktName != "")
            {
              parameter.setTag("SUBCKT");
              parameter.setVal(subcktName);
              option_block.addParam( parameter );
            }

            parameter.setTag("V");
            parameter.setVal( 1.0 );
            option_block.addParam( parameter );

            // node name:
            field = parsed_line[position+2].string_;
            field.toUpper();
            parameter.setTag( field );
            parameter.setVal( 0.0 );
            option_block.addParam( parameter );

            // Check that this is not a global node in a subcircuit.
            if ( field.substr(0,2) == "$G" && subcktName != "" )
            {
              Report::UserWarning().at(netlist_filename, parsed_line[0].lineNumber_)
                << "Ignored .NODESET on global node " << field << " inside subcircuit " 
                << subcktName << ", move statement to global scope.";
              addBlock = false;
            }

            // value:
            field = parsed_line[position+5].string_;
            field.toUpper();
            parameter.setTag( "VALUE" );
            parameter.setVal( field );
            option_block.addParam( parameter );
            
            position += 6;
          }
          else if( parsed_line[position+6].string_ == "," )
          {
            msg << "Voltage differences not supported.";
            p_err = position;
          }
          else
          {
            msg << "Unrecognized parenthetical specification.";
            p_err = position;
          }
        }
        else
        {
          msg << "Unrecognized parenthetical specification";
          p_err = position;
        }
      }
      else // This is either an expression, or a STEP parameter.
      {
        msg << "Unrecognized .NODESET specification";
        p_err = position;
      }

      if (!msg.str().empty())
      {
        msg << " in .NODESET near ";
        position = p_err+4;
        p_err -= 2;
        if (p_err<0)
          p_err = 0;
        if (position >= numFields)
          position = numFields;
        while (p_err < position)
        {
          msg << parsed_line[p_err].string_ + " ";
          ++p_err;
        }
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << msg.str();
      }
    }
  }
  else
  {
    while ( position < numFields )
    {
      if (subcktName != "")
      {
        parameter.setTag("SUBCKT");
        parameter.setVal(subcktName);
        option_block.addParam( parameter );
      }

      parameter.setTag("V");
      parameter.setVal( 1.0 );
      option_block.addParam( parameter );

      // node name:
      field = parsed_line[position].string_;
      field.toUpper();
      parameter.setTag( field );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      // Check that this is not a global node in a subcircuit.
      if ( field.substr(0,2) == "$G" && subcktName != "" )
      {
        Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Ignored .NODESET on global node " << field << " inside subcircuit " 
          << subcktName << ", move statement to global scope.";
        addBlock = false;
      }

      // value:
      field = parsed_line[position+1].string_;
      field.toUpper();
      parameter.setTag( std::string("VALUE") );
      parameter.setVal( field );
      option_block.addParam( parameter );
      
      position += 2;
    }
  }

  if (addBlock)
  {
    circuit_block.addOptions(option_block);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : extractSaveData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool
extractSaveData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("SAVE", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();
  int linePosition = 0;

  Util::Param parameter("", "");
  ExtendedString ES("");

  ++linePosition;
  while (linePosition < numFields)
  {
    ES = parsed_line[linePosition].string_ ;
    ES.toUpper ();

    if (ES == "TYPE" || ES == "FILE" || ES == "TIME" || ES == "LEVEL" )
    {
      parameter.setTag( ES );
      // error out on syntax such as:
      //    TYPE
      // which is missing the =<val>
      if (++linePosition >= numFields) 
      {
        Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Unrecognized field in .SAVE line";
        return false;
      }
      // error out on syntax such as:
      //    TYPE=
      // which is missing the <val>
      if (parsed_line[linePosition].string_ == "=")
      {
        if (++linePosition >= numFields) 
	{
          Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
            << "Unrecognized field in .SAVE line";
          return false;
        }
      }
      parameter.setVal( parsed_line[linePosition].string_ );
      option_block.addParam( parameter );
      ++linePosition;
    }
    else
    {
      Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
        << "Unrecognized field in .SAVE line";
      return false;
    }
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(InitialConditionsManager & initial_conditions_manager, PkgOptionsMgr &options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("OP_IO");

  parameters.insert(Util::ParamMap::value_type("INPUT", Util::Param("INPUT", "")));
  parameters.insert(Util::ParamMap::value_type("OUTPUT", Util::Param("OUTPUT", "")));

  options_manager.addCommandParser(".DCVOLT", extractICData);
  options_manager.addCommandParser(".IC", extractICData);
  options_manager.addCommandParser(".NODESET", extractNodeSetData);
  options_manager.addCommandParser(".SAVE", extractSaveData);

  options_manager.addCommandProcessor("IC", 
      IO::createRegistrationOptions(initial_conditions_manager, &InitialConditionsManager::registerIC));
  options_manager.addCommandProcessor("NODESET", 
      IO::createRegistrationOptions(initial_conditions_manager, &InitialConditionsManager::registerNodeSet));
  options_manager.addCommandProcessor("SAVE", 
      IO::createRegistrationOptions(initial_conditions_manager, &InitialConditionsManager::registerSave));

  return true;
}

} // namespace IO
} // namespace Xyce

