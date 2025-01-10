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

//-----------------------------------------------------------------------------
// Purpose       : Class for re-reading Xyce file output, of simulation results,
//                 that is in FORMAT=PROBE.
// Special Notes :
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/10/2015
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include<N_IO_OutputCsd.h>
#include <N_IO_ParsingHelpers.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputCsd::OutputCsd() 
// Purpose       : Constructor for CSD file IO
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 12/10/15
//-----------------------------------------------------------------------------
OutputCsd::OutputCsd() 
{
}


//-----------------------------------------------------------------------------
// Function      : OutputCsd::~OutputCsd() 
// Purpose       : Destructor for CSD file IO
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 12/10/15
//-----------------------------------------------------------------------------
OutputCsd::~OutputCsd() 
{
}


// Note The following two functions should be made more general 
// to extract strings given some set of deliminators chars 
// and the same for doubles.  Then the could be moved to the base class
// and re-used for other text based formats.

//-----------------------------------------------------------------------------
// Function      : OutputCsd::getOutputVarNames() 
// Purpose       : Get header line of CSD file for var names.
// Special Notes : return false if we cannot read from file.
// Creator       : Pete Sholander
// Creation Date : 12/10/15
//-----------------------------------------------------------------------------
bool OutputCsd::getOutputVarNames( std::vector< std::string > & varNames )
{
  bool retVal = true;
  std::stringstream extractedVarName;
  bool doneWithReadLine = false;
  bool withinQuotes=false;
  int numBraces=0;

  // end of data is denoted by a line with only #; on it
  const std::string endDataBlockText("#;");
  // list of variable names is proceeded by a line with only
  // #N on it.
  bool foundPoundN = false;

  // parse forward to find the line with the variable names.  
  // It is preceded by a line that contains only #N.  Note
  // also that the variables are enclosed in single quotes
  // An example is:
  // #N
  // 'V(1)' 'V(2)' 'V(3)' 'V(4)' 'V(5)' 
  std::string aLine;
  while (!foundPoundN && std::getline( *istreamPtr_, aLine ))
  {
    if( aLine.compare(endDataBlockText) ==0 )
    {
      // hit end-of-data-block string (which is #;) without
      // finding the #N line that starts the header section of
      // each data block.  So, return false 
      retVal = 0;
      return retVal;
    }
    else if ( (aLine.length() > 1) && (aLine.compare(0,2,"#N")==0) )
    {
      foundPoundN = true;
    }
  }

  // reached end of .CSD file without finding either #N or #; lines
  if (!foundPoundN)
  {
    retVal = 0;
    return retVal;
  }

  // parse the column header names from the line following #N line
  // add TIME first, since the #N line only lists the variables
  varNames.push_back( "TIME" );
  while( !doneWithReadLine )
  {
    char characterRead=0;
    istreamPtr_->get( characterRead );
    if( characterRead == '\n' || characterRead == '\r' )
    {
      doneWithReadLine=true;
      break;
    }
    
    if( characterRead == '\'' && withinQuotes == false )
    {
      withinQuotes=true;
      istreamPtr_->get( characterRead );
      if( characterRead == '\'') 
      {
        // this is the error-case of back-to-back ''
        Report::UserError0() << "Header column name missing on #N line in re-measured .CSD file";
      }
    }
    
    if ( withinQuotes )  
    {
      if (characterRead == '\'')
      {
        if (numBraces == 1)
        {
          // this is the error-case of an unmatched {
          Report::UserError0() << "Unmatched { in #N line in in re-measured .CSD file";
        }
        else
        {
          withinQuotes = false;
          std::string name;
          extractedVarName >> name;
          varNames.push_back( name );
          extractedVarName.clear();
        }
      }
      else 
      {
        if( characterRead == '{' )
        {
          ++numBraces;
        }
        else if( characterRead == '}' )
        {
        --numBraces;
        }
      
        if (numBraces >0)  
        {
          // need to remap spaces within "{ }" to "_"
          if(characterRead != ' ' && characterRead != '\t' && 
             characterRead != '\r' && characterRead != '\n' )
             extractedVarName.put(characterRead);
          else
            extractedVarName.put('_');
        }
        else
        {
          extractedVarName.put(characterRead);
        }
      } 
    }
  }

  if( varNames.size() == 0 )
  {
    // nothing read, so return false 
    Report::UserError0() << "No varible names found on #N line in re-measured .CSD file";
    retVal = false;
  }  
  else if (numBraces == 1)
  {
    // this is the error-case of an unmatched {
    Report::UserError0() << "Unmatched { in #N line in re-measured .CSD file";
    retVal=false;
  }
  else if (withinQuotes == true)
  {
    // this is the error-case of an unmatched '
    Report::UserError0() << "Unmatched ' in #N line in re-measured .CSD file";
    retVal=false;
  }

  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : OutputCsd::getOutputNextVarValuesParallel() 
// Purpose       : Get a line of data from CSD file.  
// Special Notes : Return false if we hit the end of the file.
// Creator       : Pete Sholander
// Creation Date : 12/10/15
//-----------------------------------------------------------------------------
int OutputCsd::getOutputNextVarValuesParallel( Linear::Vector * varValues )
{
  int retVal = 1;
  const std::string endDataBlockText("#;");
  std::string word;

  std::vector<std::string> parsedLine;
  int numDataPoints;
  int varIndex = 0;
  
  // read an entire line, and then check if this is the end of the file 
  // or the endMarkerText. readLine() comes from N_IO_ParsingHelpers.C
  std::string aLine;
  readLine(*istreamPtr_, aLine );
  if( istreamPtr_->eof() )
  {
    // hit end of file. return false 
    retVal = 0;
    return retVal;
  }
  else if ( aLine.compare(endDataBlockText)==0 ) 
  {
    // Each data block ends with a #; line.
    // Function returns 0 if it hit end-of-file.
    // It returns 1 if it successfully skipped over a 
    // header block at the start of a block of .STEP data,
    // found a #C line 
    if ( !handleEndOfCsdFileDataBlock(parsedLine) )
    {
      // hit end of file. return false 
      retVal = 0;
      return retVal;
    }
  }
  else
  {
    // turn the input lines into tokens
    splitCsdFileLine(aLine, parsedLine );
  }

  // first line, of each block of time-value data, should start with 
  // a #C line, which has the time value and the number of data points.  
  // It is followed by a list of values (0.00000000e+00:1) that may 
  // span multiple lines.  The :1 after a data value says this is the 
  // value for the first data column.  An example is:
  // #C 0.00000000e+00 5
  // 0.00000000e+00:1   -5.00000000e-01:2   5.00000000e-01:3   0.00000000e+00:4
  // 0.00000000e+00:5   
  // #C 1.00000000e-05 5
  // 6.27773973e-02:1   -4.37569034e-01:2   4.37569034e-01:3   1.25554795e-01:4
  // -1.25554795e-01:5 

  if ( (parsedLine.size() == 3) && (parsedLine[0] == "#C") )
  {
    // time value is in token 1.  Number of Data Points is in token 2.
    (*varValues)[varIndex] = atof(parsedLine[1].c_str());
    numDataPoints = atoi(parsedLine[2].c_str());
    varIndex++;
  }
  else
  {
    Report::DevelFatal() << "Error reading #C line in remeasured CSD file";
    retVal = 0;
    return retVal;
  }
 
  std::size_t colonIdx;
  const std::string validNumberChars="+-Ee.0123456789";
  while (varIndex < numDataPoints+1)
  {
    
    // get next line
    readLine(*istreamPtr_, aLine );
    if( istreamPtr_->eof() )
    {
      // hit end of file. return false 
      retVal = 0;
      return retVal;
    }
    else if ( aLine.compare(endDataBlockText)==0 ) 
    {
      // Each data block ends with a #; line.
      // Function returns 0 if it hit end-of-file.
      // It returns 1 if it successfully skipped over a 
      // header block at the start of a block of .STEP data 
      if ( !handleEndOfCsdFileDataBlock(parsedLine) )
      {
        // hit end of file. return false 
        retVal = 0;
        return retVal;
      }
    }
    else
    {
      splitCsdFileLine(aLine, parsedLine );
    }

    for (int idx=0;idx<parsedLine.size();idx++)
    {
      // need to remove the characters after the colon
      colonIdx = parsedLine[idx].find(":");
      word = parsedLine[idx].substr(0,colonIdx);
      // check that the remaining string only contains legal characters for
      // a double.
      if (word.find_first_not_of(validNumberChars) != std::string::npos)
      {
        Report::DevelFatal() << "Error reading data line in remeasured CSD file";
        retVal = 0;
        return retVal;
      }
      (*varValues)[varIndex] = atof(word.c_str());
      varIndex++;
    }
  }
 
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : OutputCsd::getNextVarValuesSerial() 
// Purpose       : Get a line of data from CSD file.  This version is used 
//               : within the ERROR measure.  
// Special Notes : Return false if we hit the end of the file.
// Creator       : Pete Sholander
// Creation Date : 8/22/16
//-----------------------------------------------------------------------------
int OutputCsd::getOutputNextVarValuesSerial( std::vector<double> * varValues )
{
  int retVal = 1;
  const std::string endDataBlockText("#;");
  std::string word;

  std::vector<std::string> parsedLine;
  int numDataPoints;
  int varIndex = 0;
  
  // read an entire line, and then check if this is the end of the file 
  // or the endMarkerText. readLine() comes from N_IO_ParsingHelpers.C
  std::string aLine;
  readLine(*istreamPtr_, aLine );
  if( istreamPtr_->eof() )
  {
    // hit end of file. return false 
    retVal = 0;
    return retVal;
  }
  else if ( aLine.compare(endDataBlockText)==0 ) 
  {
    // Each data block ends with a #; line.
    // Function returns 0 if it hit end-of-file.
    // It returns 1 if it successfully skipped over a 
    // header block at the start of a block of .STEP data,
    // found a #C line 
    if ( !handleEndOfCsdFileDataBlock(parsedLine) )
    {
      // hit end of file. return false 
      retVal = 0;
      return retVal;
    }
  }
  else
  {
    // turn the input lines into tokens
    splitCsdFileLine(aLine, parsedLine );
  }

  // first line, of each block of time-value data, should start with 
  // a #C line, which has the time value and the number of data points.  
  // It is followed by a list of values (0.00000000e+00:1) that may 
  // span multiple lines.  The :1 after a data value says this is the 
  // value for the first data column.  An example is:
  // #C 0.00000000e+00 5
  // 0.00000000e+00:1   -5.00000000e-01:2   5.00000000e-01:3   0.00000000e+00:4
  // 0.00000000e+00:5   
  // #C 1.00000000e-05 5
  // 6.27773973e-02:1   -4.37569034e-01:2   4.37569034e-01:3   1.25554795e-01:4
  // -1.25554795e-01:5 

  if ( (parsedLine.size() == 3) && (parsedLine[0] == "#C") )
  {
    // time value is in token 1.  Number of Data Points is in token 2.
    (*varValues).push_back( atof(parsedLine[1].c_str()) );
    numDataPoints = atoi(parsedLine[2].c_str());
    varIndex++;
  }
  else
  {
    Report::DevelFatal() << "Error reading #C line in remeasured CSD file";
    retVal = 0;
    return retVal;
  }
 
  std::size_t colonIdx;
  const std::string validNumberChars="+-Ee.0123456789";
  while (varIndex < numDataPoints+1)
  {
    
    // get next line
    readLine(*istreamPtr_, aLine );
    if( istreamPtr_->eof() )
    {
      // hit end of file. return false 
      retVal = 0;
      return retVal;
    }
    else if ( aLine.compare(endDataBlockText)==0 ) 
    {
      // Each data block ends with a #; line.
      // Function returns 0 if it hit end-of-file.
      // It returns 1 if it successfully skipped over a 
      // header block at the start of a block of .STEP data 
      if ( !handleEndOfCsdFileDataBlock(parsedLine) )
      {
        // hit end of file. return false 
        retVal = 0;
        return retVal;
      }
    }
    else
    {
      splitCsdFileLine(aLine, parsedLine );
    }

    for (int idx=0;idx<parsedLine.size();idx++)
    {
      // need to remove the characters after the colon
      colonIdx = parsedLine[idx].find(":");
      word = parsedLine[idx].substr(0,colonIdx);
      // check that the remaining string only contains legal characters for
      // a double.
      if (word.find_first_not_of(validNumberChars) != std::string::npos)
      {
        Report::DevelFatal() << "Error reading data line in remeasured CSD file";
        retVal = 0;
        return retVal;
      }
      (*varValues).push_back( atof(word.c_str()) );
      varIndex++;
    }
  }
 
  return retVal;
}

//----------------------------------------------------------------------------
// Function       : splitCsdFileLine
// Purpose        : Split a line in a .CSD file into a tokenized vector of strings.
// Special Notes  :
// Scope          : public
// Creator        : Pete Sholander
// Creation Date  : 12/14/2015
//----------------------------------------------------------------------------
void OutputCsd::splitCsdFileLine( const std::string& charLine, std::vector<std::string>& line )
{
  int lineLength = charLine.length();
  int currPtr = 0;
  char c;
  const std::string nonid(" \t\n\r\0");
  line.clear();

  while ( currPtr < lineLength )
  {
    std::string string_;
    string_.reserve(16);

    c=charLine[currPtr];
    if (nonid.find(c) == nonid.npos) //not a whitespace character
    {
      string_ += c;
      currPtr++;
      while ( currPtr < lineLength )
      {
        c=charLine[currPtr];
        if (nonid.find(c) == nonid.npos)
        {
          string_ += c;
          c = charLine[currPtr++ ];
        }
        else
        {
          currPtr++;
          break;
        }
      }
      
      if (string_.length() > 0)
      {
        line.push_back(string_);
      }
    }
    else
    {
      currPtr++;
    }   
  }
}

//----------------------------------------------------------------------------
// Function       : handleEndOfCsdFileDataBlock
// Purpose        : Exits at end-of-file, or skip over the header block at 
//                  the start of each block of .STEP data in a .CSD file.
// Special Notes  :
// Scope          : public
// Creator        : Pete Sholander
// Creation Date  : 12/18/2015
//----------------------------------------------------------------------------
bool OutputCsd::handleEndOfCsdFileDataBlock(std::vector<std::string>& parsedLine)
{
  std::string aLine;
  bool foundPoundC = false;
  bool retVal = 1;

  // need to skip past the #; line, and read the next line.
  readLine(*istreamPtr_, aLine );

  if( istreamPtr_->eof() )
  {
    // hit end of file. return false 
    retVal = 0;
    return retVal;
  }
  else 
  {
    // look for the #H line that starts each data block.  
    // Exit if it's not the next line.
    if (aLine.compare(0,2,"#H")!=0)
    {
      Report::DevelFatal() << "Did not find #H at start of STEP data in remeasured CSD file";
      retVal = 0;
      return retVal;
    }
  }

  while ( !foundPoundC && std::getline( *istreamPtr_, aLine ) )
  {
    if ( (aLine.length() > 1) && (aLine.compare(0,2,"#C")==0) )
    {
      foundPoundC = true;
      splitCsdFileLine(aLine, parsedLine );
      retVal = 1;
      return retVal;
    }
  }

  // reached end of .CSD file without finding another #C line
  if (!foundPoundC)
  {
    retVal = 0;
    return retVal;
  }

  return retVal;
}

} // namespace IO
} // namespace Xyce
