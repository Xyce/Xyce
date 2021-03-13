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

//-----------------------------------------------------------------------------
// Purpose       : Class for re-reading Xyce file output, of simulation results,
//                 that is in FORMAT=<STD|CSV|NOINDEX>
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/06/2012
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include<N_IO_OutputPrn.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputPrn::OutputPrn() 
// Purpose       : Constructor for PRN and CSV file IO
// Special Notes : This function also works for .CSV files, since they have the 
//                 same format as comma-delimited .PRN files, except that the 
//                 .PRN files have an additional "footer" line of something like
//                 "End of Xyce(TM) Simulation" as their last line.  That 
//                 difference has no effect on this re-measure code. 
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
OutputPrn::OutputPrn() 
{
}


//-----------------------------------------------------------------------------
// Function      : OutputPrn::~OutputPrn() 
// Purpose       : Destructor for PRN and CSV file IO
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
OutputPrn::~OutputPrn() 
{
}


// Note The following two functions should be made more general 
// to extract strings given some set of deliminators chars 
// and the same for doubles.  Then the could be moved to the base class
// and re-used for other text based formats.

//-----------------------------------------------------------------------------
// Function      : OutputPrn::getOutputVarNames() 
// Purpose       : Get header line of PRN or CSV file for var names.
// Special Notes : return false if we cannot read from file.
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
bool OutputPrn::getOutputVarNames( std::vector< std::string > & varNames )
{
  bool retVal = true;
  std::stringstream extractedVarName;
  bool doneWithReadLine = false;
  bool withinWord=false;
  int numBraces=0;
  int numParens=0;
  while( !doneWithReadLine )
  {
    char characterRead=0;
    istreamPtr_->get( characterRead );
    if( characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
    }
    
    // numBraces is used to determine whether the parsing is inside of 
    // an expresion in the first line of the .PRN or .CSV file
    if( characterRead == '{' )
      ++numBraces;
    if( characterRead == '}' )
      --numBraces;

    // numParens is then used to determine whether the parsing is inside of
    // a V(a,b) syntax that is outside of an expression.
    if( characterRead == '(' )
      ++numParens;
    if( characterRead == ')' )
      --numParens;
    
    // This IF statement can be read as:
    //    a) Are we in an expression, or
    //    b) Are we inside of a V(a,b) syntax that is outside of an expression, or
    //    c) is the current character NOT one of the allowed delimiters for a
    //       .PRN or .CSV file?
    //
    // If any of this is true then the IF clause is TRUE, and the current characterRead 
    // will be part of the same varNames entry.
    //
    // In Xyce 6.5, or earlier, this IF statement only checked for whitespace as
    // delimters. For Xyce 6.6, that was changed to include comma as a delimiter
    // to account for remeasure of .CSV files, and also .PRN files made with the 
    //     Xyce -delim COMMA <netlist> syntax.
    if( (numBraces >0) || (numParens>0) ||
        (characterRead != ' ' && characterRead != '\t' && 
         characterRead != '\r' && characterRead != '\n' &&
         characterRead != ',' ) )
    {
      withinWord = true;
      // need to remap spaces with "{ }" to "_"
      if(characterRead != ' ' && characterRead != '\t' && 
         characterRead != '\r' && characterRead != '\n' )
        extractedVarName.put(characterRead);
      else
        extractedVarName.put('_');
    }
    else
    {
      if( withinWord )
      {
        // just found first white space (or comma) after word, when we
        // were not in an expression, or a V(a,b) syntax outside of an
        // expression.  So, transfer extractedVarName to varNames array 
        std::string name;
        extractedVarName >> name;
        varNames.push_back( name );
        extractedVarName.clear();
        withinWord=false;
      }
    }
  }
  
  if( varNames.size() == 0 )
  {
    // nothing read so return false 
    retVal = false;
  }  

  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : OutputPrn::getOutputNextVarValuesParallel() 
// Purpose       : Get a line of data from .PRN or .CSV file.  
// Special Notes : Return false if we hit the end of the file.
// Creator       : Rich Schiek
// Creation Date : 4/16/13
//-----------------------------------------------------------------------------
int OutputPrn::getOutputNextVarValuesParallel( Linear::Vector * varValues )
{
  int retVal = 1;
  // endMarkerText1 is the text written at the end of a .PRN file for a
  // .TRAN simulation if .STEP is not used.  endMarkerText2 is the text 
  // when .STEP is used with .TRAN for a .PRN file.  For a .CSV file,
  // we read until end-of-file.
  const std::string endMarkerText1("End of Xyce(TM) Simulation");
  const std::string endMarkerText2("End of Xyce(TM) Parameter Sweep");
  bool doneWithReadLine = false;
  
  // read an entire line to check if this is the end of the file 
  // or one of the endMarkerText strings

  std::string aLine;
  std::getline( *istreamPtr_, aLine );
  if( (istreamPtr_->eof()) || (aLine.compare(endMarkerText1)==0) ||
        (aLine.compare(endMarkerText2)==0)    )
  {
    // hit end of file. return false 
    doneWithReadLine = true;
    retVal = 0;
    return retVal;
  }
 
  std::stringstream theLineStream( aLine ); 
  std::stringstream extractedVarValue;
  bool withinWord=false;
  // validNumberChars also include +-Ee. to account for the numbers
  // being in scientific notation 
  const std::string validNumberChars="+-Ee.0123456789";
  // assume we start filling the array at index 0.
  int varIndex = 0;
  while( !doneWithReadLine )
  {
    char characterRead=0;
    theLineStream.get( characterRead );
    if(  theLineStream.eof() || characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
    }

    // Columns in a Xyce .PRN file are typically delimited with whitespace.
    // However, they can be delimited with comma, if the Xyce -delim COMMA <netlist> 
    // syntax was used to make the .PRN file.  They are also delimited with commas
    // for .CSV files.
    if( characterRead != ' ' && characterRead != '\t' && 
        characterRead != '\r' && characterRead != '\n' &&
        characterRead != ',')
    {
      // need to ensure that characters read are valid numbers 
      // 0123456789
      if( withinWord || (validNumberChars.find( characterRead ) != std::string::npos ) )
      {
        withinWord = true;
        extractedVarValue.put(characterRead);
      }
    }
    else
    {
      if( withinWord ) 
      {
        // just found first white space after word so 
        // transfer extractedVarName to varValues array 
        double value;
        extractedVarValue >> value;
        (*varValues)[varIndex] = value ;
        varIndex++;
        extractedVarValue.clear();
        withinWord=false;
      }
    }
    
    if( doneWithReadLine && withinWord )
    {
      // this is to capture the last element on the line
      // where there isn't any more whitespace after the last char.
      double value;
      extractedVarValue >> value;
      (*varValues)[varIndex] = value ;
      varIndex++;
      extractedVarValue.clear();
      withinWord=false;
    }
     
  }
 
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : OutputPrn::getOutputNextVarValuesSerial() 
// Purpose       : Get a line of data from .PRN or .CSV file. This version is used 
//               : within the ERROR measure.   
// Special Notes : Return false if we hit the end of the file.
// Creator       : Pete Sholander
// Creation Date : 8/22/16
//-----------------------------------------------------------------------------
int OutputPrn::getOutputNextVarValuesSerial( std::vector<double> * varValues )
{
  int retVal = 1;
  // endMarkerText1 is the text written at the end of a .PRN file for a
  // .TRAN simulation if .STEP is not used.  endMarkerText2 is the text 
  // when .STEP is used with .TRAN for a .PRN file.  For a .CSV file,
  // we read until end-of-file.
  const std::string endMarkerText1("End of Xyce(TM) Simulation");
  const std::string endMarkerText2("End of Xyce(TM) Parameter Sweep");
  bool doneWithReadLine = false;
  
  // read an entire line to check if this is the end of the file 
  // or one of the endMarkerText strings

  std::string aLine;
  std::getline( *istreamPtr_, aLine );
  if( (istreamPtr_->eof()) || (aLine.compare(endMarkerText1)==0) ||
        (aLine.compare(endMarkerText2)==0)    )
  {
    // hit end of file. return false 
    doneWithReadLine = true;
    retVal = 0;
    return retVal;
  }
 
  std::stringstream theLineStream( aLine ); 
  std::stringstream extractedVarValue;
  bool withinWord=false;
  // validNumberChars also include +-Ee. to account for the numbers
  // being in scientific notation 
  const std::string validNumberChars="+-Ee.0123456789";
  // assume we start filling the array at index 0.
  int varIndex = 0;
  while( !doneWithReadLine )
  {
    char characterRead=0;
    theLineStream.get( characterRead );
    if(  theLineStream.eof() || characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
    }

    // Columns in a Xyce .PRN file are typically delimited with whitespace.
    // However, they can be delimited with comma, if the Xyce -delim COMMA <netlist> 
    // syntax was used to make the .PRN file.  They are also delimited with commas
    // for .CSV files.
    if( characterRead != ' ' && characterRead != '\t' && 
        characterRead != '\r' && characterRead != '\n' &&
        characterRead != ',')
    {
      // need to ensure that characters read are valid numbers 
      // 0123456789
      if( withinWord || (validNumberChars.find( characterRead ) != std::string::npos ) )
      {
        withinWord = true;
        extractedVarValue.put(characterRead);
      }
    }
    else
    {
      if( withinWord ) 
      {
        // just found first white space after word so 
        // transfer extractedVarName to varValues array 
        double value;
        extractedVarValue >> value;
        (*varValues).push_back(value) ;
        varIndex++;
        extractedVarValue.clear();
        withinWord=false;
      }
    }
    
    if( doneWithReadLine && withinWord )
    {
      // this is to capture the last element on the line
      // where there isn't any more whitespace after the last char.
      double value;
      extractedVarValue >> value;
      (*varValues).push_back(value) ;
      varIndex++;
      extractedVarValue.clear();
      withinWord=false;
    }
     
  }
 
  return retVal;
}

} // namespace IO
} // namespace Xyce
