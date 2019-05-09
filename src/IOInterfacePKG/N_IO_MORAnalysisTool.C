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
// Purpose        : Defines the MORAnalysisTool class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 09/10/2014
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_CIR_Xyce.h>
#include <N_IO_fwd.h>
#include <N_IO_MORAnalysisTool.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_PDS_Comm.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CmdParse.h>
#include <N_IO_ParsingHelpers.h>
#include <N_IO_ParameterBlock.h>

#include <N_UTL_Param.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_NoCase.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::reduceLinearSubcircuits()
// Purpose       : Find the linear subcircuits in the CircuitContext and reduce
//                 them according to the MOR options.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
bool MORAnalysisTool::reduceLinearSubcircuits()
{
  bool success = true;

  // Get the top level instance list, if empty there are no subcircuits to reduce.
  const std::vector< std::string >& topLevelInstanceList = circuitContext_.getSubcktList();
  if (topLevelInstanceList.empty())
  {
    return false;
  }

  // Get the .MOR line, if there are any arguments other than the analysis statement,
  // we will produce a warning and eliminate the MOR statement from the desired analysis.
  bool isSaveRedSys = false, isAutoSize = false;
  std::vector< NetlistLocation > optionsLocation;
  std::list<Util::OptionBlock>::iterator op_analysis_it = optionsTable_.begin();
  for ( ; op_analysis_it != optionsTable_.end(); op_analysis_it++ )
  {
    if (op_analysis_it->getName() == "MOR")
    {
      // Check the arguments in this line.
      if (op_analysis_it->size() > 0)
      {
        Report::UserWarning() << 
          "Specifying port list for MOR reduction with an analysis type is not supported, .MOR line is being ignored."
          << std::endl;
        return false;
      }

      // Collect where the MOR line is in the input files.
      optionsLocation.push_back( op_analysis_it->getNetlistLocation() );
    }
    
    // While were going through the option block list, grab any user-desired subcircuits.
    if (op_analysis_it->getName() == "MOR_OPTS")
    {
      for (Util::ParamList::const_iterator it = (*op_analysis_it).begin(), end = (*op_analysis_it).end(); it != end; ++it )
      {
        if ( std::string((*it).tag(),0,7) == "SUBCKTS" )
        {
          std::string subcktName = (*it).stringValue();
          Util::toUpper( subcktName );
          desiredSubckts_.push_back( subcktName );
        }
        else if ( std::string((*it).tag(),0,12) == "ANALYSISONLY" 
                  && (*it).stringValue()=="1" )
        {
          analysisOnly_ = true;
        }
        else
        {
          // Add this option to the MOR_OPTS string for later use.
          // NOTE:  We are only allowing autosizing with autodetection, so any "SIZE" will be ignored.
          //        We are also making sure the reduced system is written to file for later use.
          if ( (*it).tag() != "SIZE" )
          {
            morOptionsLine_ += (*it).tag();
            morOptionsLine_ += "=";
            morOptionsLine_ += (*it).stringValue();
            morOptionsLine_ += " ";
            if ((*it).tag()=="SAVEREDSYS" && (*it).stringValue()=="1")
            {
              isSaveRedSys = true;
            }
            if ((*it).tag()=="AUTOSIZE" && (*it).stringValue()=="1")
            {
              isAutoSize = true;
            }
          }
        }
      }

      // Collect where the MOR_OPTS line is in the input files.
      optionsLocation.push_back( op_analysis_it->getNetlistLocation() );
    }
  }

  // It is necessary to save the reduced systems, so add the option if it doesn't exist.
  if (!isSaveRedSys)
  {
    morOptionsLine_ += " SAVEREDSYS=1";
  }
  if (!isAutoSize)
  {
    morOptionsLine_ += " AUTOSIZE=1";
  }

  // Store the location of the MOR options so they can be removed from the final netlist.
  usedSubcktLocations_["OPTIONSLOCATION"] = optionsLocation;

  // Go through the hierarchy in the CircuitContext and collect all reducibleSubcircuits
  // (and their children for extraction)
  std::string emptyName;
  collectReducibleSubcircuits( emptyName ); 

  if ( desiredSubckts_.size() > 0 )
  {
    std::map< std::string, std::map< std::string, NetlistLocation > > newRedSubcktNames;
    for ( unsigned i=0; i<desiredSubckts_.size(); i++ )
    {
      if (DEBUG_MOR)
      {
        Xyce::dout() << "User requested subcircuit " << desiredSubckts_[i] << " be reduced." << std::endl;
      }
      std::map< std::string, std::map< std::string, NetlistLocation > >::iterator 
        rsnIter = redSubcktNames_.find( desiredSubckts_[i] );
      if (rsnIter == redSubcktNames_.end())
      {
        Report::UserWarning() << "MOR cannot be applied to user-specified subcircuit " 
                              << desiredSubckts_[i] << std::endl;
      }
      else
      {
        newRedSubcktNames[ rsnIter->first ] = rsnIter->second;
      }
    }

    // Now copy over pruned version of reducible subcircuit list, if everything is pruned then we failed.
    if ( newRedSubcktNames.empty() )
    {
      return false;
    }
    else
    {
      redSubcktNames_ = newRedSubcktNames;
    }
  }
 
  // Print all reducible subcircuits and their children.
  if (DEBUG_MOR)
  {
    Xyce::dout() << "MOR analysis found " << redSubcktNames_.size() << " individual subcircuits to reduce!" << std::endl;
    std::map< std::string, std::map< std::string, NetlistLocation > >::iterator rsnIter = redSubcktNames_.begin();
    for ( ; rsnIter != redSubcktNames_.end(); rsnIter++ )
    {
      Xyce::dout() << "Reducible subckt " << rsnIter->first  << " : ";
      std::map< std::string, NetlistLocation >::iterator innerIter = rsnIter->second.begin();
      for ( ; innerIter != rsnIter->second.end(); innerIter++ )
        Xyce::dout() << innerIter->first << " ";
      Xyce::dout() << std::endl;
    }
  }

  // If user only requests the analysis be performed, but nothing reduced, exit here.
  if ( analysisOnly_ )
  {
    return false;
  }
 
  // Generate netlist files for each of the reducible circuits. 
  generateMORFiles();

  // Run model order reduction via Xyce on these files.
  runFiles( morFilenames_ );

  // Generate reintegrated netlist.
  std::string newNetlistName;
  generateReintegratedFiles( newNetlistName );

  // Run netlist with regenerated files.
  std::vector< std::string > topNetlist( 1, newNetlistName );
  runFiles( topNetlist );

  return success;
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::generateMORFiles
// Purpose       : Output reducible subcircuits to file for MOR analysis.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::generateMORFiles()
{
  std::string topLevelNetlist = circuitContext_.getLocation().getFilename();
  std::map< std::string, std::map< std::string, NetlistLocation > >::iterator rsnIter = redSubcktNames_.begin();
  for ( ; rsnIter != redSubcktNames_.end(); rsnIter++ )
  {
     std::string newMORFilename = topLevelNetlist + "_" + rsnIter->first + ".cir";
     morFilenames_.push_back( newMORFilename );
     if (DEBUG_MOR)
     {
       Xyce::dout() << "New filename: " << newMORFilename << std::endl;
     }

     std::ofstream morFile(newMORFilename.c_str());

     // Some error checking in case we can't open the file.
     if(morFile.fail())
     {
       Report::UserError0() << "MOR analysis autodetection cannot open the output file " << newMORFilename;
     }

     //Create title line (not the same as title line of original netlist file!)
     std::string header("XYCE-generated MOR Netlist file for subcircuit " + rsnIter->first);
     morFile << header << std::endl << std::endl;

/*
     // Write out model lines for reducible subcircuit and children to morFile.
     // NOTE:  This needs to be done differently because the parsedLine is no longer
     //        provided by the model blocks.
     std::list< ModelMap > modelLines = redSubcktModels_[ rsnIter->first ];
     std::list< ModelMap >::iterator mlIter = modelLines.begin();
     for ( ; mlIter != modelLines.end(); mlIter++ )
     {
       ModelMap::iterator iter2 = mlIter->begin();
       for ( ; iter2 != mlIter->end(); iter2++ )
       {
         bool isParen = false;
         bool isEqual = false;
         int lineSize = (iter2->second)->parsedLine.size();
         morFile << (iter2->second)->parsedLine[0].string_;
         std::cout << (iter2->second)->parsedLine[0].string_;
         for ( int i = 1; i < lineSize; ++i )
         {
           if (!isParen)
           {
             morFile << " " << (iter2->second)->parsedLine[i].string_;
             std::cout << " " << (iter2->second)->parsedLine[i].string_;
           }
           else
           {
             morFile << (iter2->second)->parsedLine[i].string_;
             std::cout<< (iter2->second)->parsedLine[i].string_;
             if (isEqual)
               morFile << " ";
           }
 
           if ((iter2->second)->parsedLine[i].string_ == "(")
             isParen = true;
           if ((iter2->second)->parsedLine[i].string_ == "=")
           {
             isEqual = true;
           }
           else
           {
             isEqual = false;
           }
         }
         morFile << std::endl << std::endl;
         std::cout << std::endl << std::endl;
       }
     }
*/
     // Write out reducible subcircuit and children to morFile.
     std::map< std::string, NetlistLocation >::iterator subcktIter = rsnIter->second.begin();
     for ( ; subcktIter != rsnIter->second.end(); subcktIter++ )
     {
        // Get ifstream for each original file.
        FileSSFPair& ssfPtr = ssfMap_[subcktIter->second.getFilename()];
        if (ssfPtr.first == 0)
          Report::UserError0() << "Could not find the input filestream for file " << subcktIter->second.getFilename();

        // Rewind file back to beginning.
        (ssfPtr.second)->setLocation(0);
        (ssfPtr.second)->setLineNumber(1);

        // Proceed through file to required line number (not efficient, but ok).
        std::string chompLine;
        chompLine.reserve(256);

        // Chomp all the lines before the subcircuit definition.
        for (int i=0; i<subcktIter->second.getLineNumber()-1; i++)
        {
          IO::readLine( *(ssfPtr.first), chompLine );
        }

        // Write all requested subcircuit lines to morFile.
        bool keepGoing = true;
        while (keepGoing)
        {
          IO::readLine( *(ssfPtr.first), chompLine );
          morFile << chompLine << std::endl;
      
          // Check if this is the end of the subcircuit 
          Util::toUpper( chompLine ); 
          if ( chompLine.find( ".ENDS" ) != std::string::npos )
          {
            morFile << std::endl; 
            keepGoing = false;
          }
        }
     }

     // Add voltage sources for each port in the port list.
     std::map< std::string, std::vector< std::string > >::iterator
       nodeListIter = redSubcktNodes_.find( rsnIter->first );
     std::string portList;

     std::vector< std::string >::iterator nodeIter = nodeListIter->second.begin();
     for ( ; nodeIter != nodeListIter->second.end(); nodeIter++ )
     {
       // Create a voltage source line: 'V<node_name> <node_name> 0 0'
       morFile << "V" << *nodeIter << " " << *nodeIter << " 0 0 " << std::endl;

       // Accumulate port names;
       portList += " ";
       portList += *nodeIter;
     }
     morFile << std::endl;

     // Add subcircuit instance for reducible subcircuit.
     morFile << "X1" << portList << " " << rsnIter->first << std::endl << std::endl;  

     // Add .MOR line with top subcircuit port list.
     morFile << ".MOR" << portList << std::endl << std::endl;
          
     // Add MOR options line.
     morFile << ".OPTIONS MOR_OPTS " << morOptionsLine_ << std::endl;

     // Must put and end of file line in for parser.
     morFile << ".END" << std::endl;

     // Close the file, we are done writing.
     morFile.close();
  }
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::collectReducibleSubcircuits
// Purpose       : Find the linear subcircuit names and locations in the input files.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::collectReducibleSubcircuits( const std::string& name )
{
  // Search for a reducible linear subcircuit if the name is an empty string.
  // Otherwise, we have found a reducible linear subcircuit and this is a subcircuit
  // used by that subcircuit.  So, the netlist location for this subcircuit needs to
  // be added the list so that the top-level linear subcircuit can be resolved.
  if ( name.empty() )
  {
    // Get current model map, just in case we need it for the reducible subcircuit.
    ModelMap currentModelMap = circuitContext_.getModelMap();

    // Get the instance list from the current CircuitContext.
    std::vector< std::string > subcktList = circuitContext_.getSubcktList();
    std::vector< std::string >::const_iterator ilIter = subcktList.begin();
    for ( ; ilIter != subcktList.end(); ilIter++ )
    {
      // Set the context for this subcircuit.
      circuitContext_.setContext( *ilIter );

      if (DEBUG_MOR)
      {
        Xyce::dout() << "For subcircuit " << *ilIter << " the total device count is : " << circuitContext_.getTotalDeviceCount() << std::endl;
        Xyce::dout() << "For subcircuit " << *ilIter << " the total linear device count is : " << circuitContext_.getTotalLinearDeviceCount() << std::endl;
      }
      if (circuitContext_.getTotalDeviceCount()==circuitContext_.getTotalLinearDeviceCount())
      {
        // If only it were this easy, but we have to look for sources.
        DeviceCountMap devCountMap = circuitContext_.getDeviceCountMap();
        // Because there are WAY too many "linear" devices, perform a find_first_not_of for RLCs.
     
        bool isOnlyRLC = true;
        DeviceCountMap::const_iterator dcmIter = devCountMap.begin();
        for ( ; dcmIter != devCountMap.end(); dcmIter++ )
        {
          if (dcmIter->first != "R" && dcmIter->first != "L" && dcmIter->first != "C")
          {
            isOnlyRLC = false;
            break;
          }  
        }

        // We are now certain that this is a reducible subckt, so recurse and get all of
        // the netlist file data from the children, if there are any.
        if (isOnlyRLC)
        {
          if (DEBUG_MOR)
          {
            Xyce::dout() << "We have found a reducible subcircuit " << *ilIter << std::endl;
          }
          // Create new entry in the map for this linear subckt and enter its location.
          if (redSubcktNames_.find( *ilIter ) == redSubcktNames_.end())
          {
            (redSubcktNames_[ *ilIter ])[ *ilIter ] = circuitContext_.getLocation();
            redSubcktNodes_[ *ilIter ] = circuitContext_.getNodeList();
            if (!currentModelMap.empty())
              redSubcktModels_[ *ilIter ].push_back( currentModelMap );

            // Now recursively call this method to collect all the children's names.
            collectReducibleSubcircuits( *ilIter );
          }
        }
      }
      else
      {
         // Continue recursion to find other linear subcircuits
         std::string emptyName;
         collectReducibleSubcircuits( emptyName );
      }
 
      // Restore context before next loop. 
      circuitContext_.restorePreviousContext();
    }
  }
  else
  {
     // A linear subcircuit has been found, just place the netlist location for any
     // supporting subckts in the redSubcktNames_ map.
     if (DEBUG_MOR)
     {
       Xyce::dout() << "Collecting children of subcircuit : " << name << std::endl;
     } 
     // Find the map associated with the largest linear subckt (passed in as "name").
     std::map< std::string, std::map< std::string, NetlistLocation > >::iterator
       rsIter = redSubcktNames_.find( name );

     // Get the instance list from the current CircuitContext.
     std::vector< std::string > subcktList = circuitContext_.getSubcktList();
     std::vector< std::string >::const_iterator ilIter = subcktList.begin();

     for ( ; ilIter != subcktList.end(); ilIter++ )
     { 
        // Set the context for this subcircuit.
        circuitContext_.setContext( *ilIter );

        // Get current model map, just in case we need it for the reducible subcircuit.
        ModelMap currentModelMap = circuitContext_.getModelMap();
      
        if (DEBUG_MOR)
        {
          Xyce::dout() << "Looking at child subcircuit " << *ilIter << std::endl; 
        }
        if ( rsIter->second.find( *ilIter ) == (rsIter->second).end() ) 
        {
           if (DEBUG_MOR)
           {
             Xyce::dout() << "Adding " << *ilIter << std::endl;
           }
           (rsIter->second)[ *ilIter ] = circuitContext_.getLocation();
           if (!currentModelMap.empty())
             redSubcktModels_[ name ].push_back( currentModelMap );
 
           // Now recursively call this method to collect all the children's locations.
           collectReducibleSubcircuits( name );
        }
      
        // Restore context before next loop. 
        circuitContext_.restorePreviousContext();
     }
  } 
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::runFiles
// Purpose       : Run a new simulator on list of input files.
// Special Notes : This is used to perform MOR analysis on reducible subcircuits.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::runFiles( const std::vector< std::string >& files )
{
  CmdParse tmpCmdLine( commandLine_ );

  // Loop over all generated netlists and run them through Xyce.
  std::vector< std::string >::const_iterator fIter = files.begin();
  for ( ; fIter != files.end(); fIter++ )
  {
    // Copy the new MOR file into the character string.
    tmpCmdLine.setNetlist( *fIter );
    
    // Create a new simulator
    N_CIR_Xyce newSimulator;

    // Run new simulator on MOR file. 
    newSimulator.run(tmpCmdLine.argc(), tmpCmdLine.argv());
  }
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::generateReintegratedFiles
// Purpose       : Output reintegrated circuit files.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::generateReintegratedFiles( std::string& newNetlistName )
{
  // Create a vector of all reduced subcircuits.
  std::vector< std::string > reducedNames;
  std::map< std::string, std::map< std::string, NetlistLocation > >::iterator
    rsnIter = redSubcktNames_.begin();
  for ( ; rsnIter != redSubcktNames_.end(); rsnIter++ )
  {
    reducedNames.push_back( rsnIter->first );
  }

  // Find which used subcircuits will be replaced with YROM devices and their locations.
  collectReducedSubcircuits( reducedNames );  

  // Find list of include files and locations that also need updating.
  std::vector< NetlistLocation > includeLocations;
  collectIncludeFileLocations( includeLocations );
  usedSubcktLocations_["INCLUDEFILES"] = includeLocations;

  // Now create a map of old filenames to new filenames.
  std::map< std::string, std::string > newFilenames;
  std::map< std::string, std::vector< NetlistLocation > >::iterator usnIter = usedSubcktLocations_.begin();
  for ( ; usnIter != usedSubcktLocations_.end(); usnIter++ ) 
  {
    if (DEBUG_MOR)
    {
      Xyce::dout() << "Reduced subcircuit: " << usnIter->first << std::endl;
    }
    for (int i=0; i<usnIter->second.size(); i++)
    {
      if ( newFilenames.find( (usnIter->second)[i].getFilename() ) == newFilenames.end() )
      {
        newFilenames[ (usnIter->second)[i].getFilename() ] = (usnIter->second)[i].getFilename() + "_ROM";
      } 
      if (DEBUG_MOR)
      {
        Xyce::dout() << "File : " << (usnIter->second)[i].getFilename() << ", Line : " << (usnIter->second)[i].getLineNumber() << std::endl;
      }
    }
  }

  // Now generate the reintegrated files.
  std::map< std::string, std::string >::iterator nfIter = newFilenames.begin();
  for ( ; nfIter != newFilenames.end(); nfIter++ )
  {
     // Collect a list of each line that will change in this file.
     std::vector< std::pair< int, std::string > > newLines;
     
     std::map< std::string, std::vector< NetlistLocation > >::iterator usnIter = usedSubcktLocations_.begin();
     for ( ; usnIter != usedSubcktLocations_.end(); usnIter++ ) 
     {
       for (int i=0; i<usnIter->second.size(); i++)
       {
         if ( nfIter->first == (usnIter->second)[i].getFilename() )
         { 
           newLines.push_back( std::pair<int, std::string>( (usnIter->second)[i].getLineNumber(), usnIter->first) );
         }
       }
     }
     
     // Sort this list based on line number.   
     std::sort(newLines.begin(), newLines.end(), sort_by_line());

     if (DEBUG_MOR)
     {
       Xyce::dout() << "For original file : " << nfIter->first << std::endl;
       for (int j=0; j<newLines.size(); j++)
       {
         Xyce::dout() << "Changing line number: " << newLines[j].first << " which is a line of type " << newLines[j].second << std::endl;
       }
     }

     // Open each new file, one at a time and rewrite them.
     std::ofstream romFile( (nfIter->second).c_str() );

     // Some error checking in case we can't open the file.
     if(romFile.fail())
     {
       Report::UserError0() << "ROM reintegration cannot open the output file " << nfIter->second;
     }

     // Get ifstream for each original file.
     FileSSFPair& ssfPtr = ssfMap_[nfIter->first];
     if (ssfPtr.first == 0)
       Report::UserError0() << "Could not find the input filestream for file " << nfIter->first;

     // Rewind file back to beginning.
     (ssfPtr.second)->setLocation(0);
     (ssfPtr.second)->setLineNumber(1);

     // Copy over unaffected lines, process only affected ones.
     std::string chompLine;
     chompLine.reserve(256);

        // Write all requested subcircuit lines to morFile.
     int currentLine = 1;
     TokenVector separatedLine;
     for (int i=0; i<(int)newLines.size(); i++)
     {
       if (DEBUG_MOR)
       {
         Xyce::dout() << "Modifying line : " << newLines[i].first << " of type " << newLines[i].second << std::endl;
       }
       // Copy over all lines until we get to the line that needs modification. 
       for (int j=currentLine; j<newLines[i].first; j++, currentLine++)
       {
         IO::readLine( *(ssfPtr.first), chompLine );
         romFile << chompLine << std::endl;
         if (DEBUG_MOR)
         {
           Xyce::dout() << "Copying over : " << chompLine << std::endl;
         }
       }
 
       // Read line being replaced.
       IO::readLine( *(ssfPtr.first), chompLine );
       if (DEBUG_MOR)
       {
         Xyce::dout() << "Replacing line: " << chompLine << std::endl;
       }

       // Replace this line and increment the counters.
       // (I think there is an issue here if the line is continued with a "+".) 
       if ( newLines[i].second == "INCLUDEFILES" )
       {
         // This will work except for whitespace and comments on the same line.
         romFile << chompLine << "_ROM" << std::endl;
       }
       else if ( newLines[i].second == "OPTIONSLOCATION" )
       {
         // Do nothing here to remove the options lines in the reintegrated file.
       }
       else
       { 
         // Have the SpiceSeparatedFieldTool split the current line.
         IO::splitLine( chompLine, separatedLine );
         
         // Now generate the new line that uses the reduced order model from the tokenized line.
         // NOTE:  Line structure is "YROM <device name> <port list> BASE_FILENAME=<file name>"
         std::string newLine("YROM ");

         // Add name for reduced order model, reuse subckt name (as it should have been unique).
         newLine += "ROM";
         newLine += separatedLine[0].string_.substr(1);
         for (unsigned int k = 1; k < separatedLine.size(); ++k)
         {
           // Copy everything except the linear subcircuit that has been reduced.
           Util::toUpper( separatedLine[k].string_ );
           if ( separatedLine[k].string_ != newLines[i].second )
           {
             newLine += separatedLine[k].string_;
           }
           else
           {
             break;
           }
         }
       
         // Now add filename for subcircuit that was reduced 
         std::string romFilename = "BASE_FILENAME=" 
            + circuitContext_.getLocation().getFilename() + "_" + newLines[i].second + ".cir";
         newLine += romFilename;

         // Write this line to the file.  
         romFile << newLine << std::endl;
       }
       currentLine++;
     }
     
     // Copy over rest of file.  
     while ( !((ssfPtr.first)->eof()) )
     {  
       IO::readLine(*(ssfPtr.first), chompLine );
       romFile << chompLine << std::endl;
     }

     // Close the file. 
     romFile.close();
  }

  // Now return the new top level netlist name.
  // If the file was not regenerated then this is the same as the top level netlist.
  newNetlistName = circuitContext_.getLocation().getFilename();
  std::map< std::string, std::string >::iterator topIter = newFilenames.find( newNetlistName );
  if ( topIter != newFilenames.end() )
  {
    newNetlistName = topIter->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::collectReducedSubcircuits
// Purpose       : Collect the name and location of all reduced subcircuits.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::collectReducedSubcircuits( const std::vector< std::string >& reducedNames )
{
  // Get the instance list and their locations from the current CircuitContext.
  std::vector< std::string > subcktList = circuitContext_.getSubcktList();

  // Make sure the list of instances is unique.
  std::sort( subcktList.begin(), subcktList.end() );
  std::vector<int>::iterator it;
  subcktList.erase(std::unique(subcktList.begin(), subcktList.end()), subcktList.end());

  // Get the instance location list from the current CircuitContext.
  unordered_map< std::string, std::list<NetlistLocation> > 
    instanceLocation = circuitContext_.getInstanceLocation();

  std::vector< std::string >::const_iterator ilIter = subcktList.begin();
  for (int i=0 ; ilIter != subcktList.end(); i++, ilIter++ )
  {
    // Set the context for this subcircuit.
    circuitContext_.setContext( *ilIter );

    std::vector< std::string >::const_iterator
      rnIter = std::find( reducedNames.begin(), reducedNames.end(), *ilIter );
    if ( rnIter != reducedNames.end() )
    {
      std::list< NetlistLocation > iLoc = instanceLocation[ *rnIter ];
      std::list< NetlistLocation >::iterator iLocIter = iLoc.begin();
      for ( ; iLocIter != iLoc.end(); iLocIter++ )
      {
        usedSubcktLocations_[ *rnIter ].push_back( *iLocIter );
      }
    } 
    else
    {
      // Recurse, we didn't find any here
      collectReducedSubcircuits( reducedNames );
    }

    // Restore context before next loop. 
    circuitContext_.restorePreviousContext();
  }
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::collectIncludeFileLocations
// Purpose       : Collect the name and location of all affected include files.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::collectIncludeFileLocations( std::vector< NetlistLocation >& includeFiles )
{
  std::map< std::string, std::vector< NetlistLocation > >::iterator usnIter = usedSubcktLocations_.begin();
  for ( ; usnIter != usedSubcktLocations_.end(); usnIter++ ) 
  {
    for (int i=0; i<usnIter->second.size(); i++)
    {
      if (DEBUG_MOR)
      {
        Xyce::dout() << "File that includes " << (usnIter->second)[i].getFilename() << " is : " << iflMap_[ (usnIter->second)[i].getFilename() ].location.getFilename() << std::endl;
      }
      bool addIncFiles = true;
      std::string newIncludeFile = (usnIter->second)[i].getFilename();
      while ( addIncFiles )
      {
        std::string currIncFile = iflMap_[ newIncludeFile ].location.getFilename(); 
        if ( currIncFile != "" )
        {
          bool isUnique = true;
          for (int j=0; j<includeFiles.size(); j++)
          {
            if ( currIncFile == includeFiles[j].getFilename() )
            { 
              isUnique = false;
              break;
            } 
          }
          if (isUnique)
          {
            includeFiles.push_back( iflMap_[ newIncludeFile ].location );
            newIncludeFile = iflMap_[ newIncludeFile ].location.getFilename();
          }
          else
          {
            addIncFiles = false;
          }
        }
        else
        {
          addIncFiles = false;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MORAnalysisTool::removeMOROptions
// Purpose       : Remove any .MOR and .OPTIONS MOR_OPTS lines from the OptionBlock
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void MORAnalysisTool::removeMOROptions()
{
  // Remove .MOR and .OPTIONS MOR_OPTS from option table before continuing simulation
  // MOR was determined to not be viable.
  std::list<Util::OptionBlock>::iterator op_analysis_it = optionsTable_.begin();
  for ( ; op_analysis_it != optionsTable_.end(); )
  {
     if (op_analysis_it->getName() == "MOR" || op_analysis_it->getName() == "MOR_OPTS")
     {
       // Remove this iterator from the list and the next one is returned.
       op_analysis_it = optionsTable_.erase( op_analysis_it );
     }
     else
     {
       op_analysis_it++;
     }
  }
}


} // namespace IO
} // nfamespace Xyce
