//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Filename      : N_IO_SpiceSeparatedFieldTool.H
//
// Purpose       : This file defines the N_IO_SpiceSeparatedFieldTool
//                 class. (for use in parsing the netlist)
//
// Special Notes :
//
// Creator       : Alan Lundin, SNL
//
// Creation Date : 08/31/99
//
//-------------------------------------------------------------------------

#ifndef SpiceSeparatedFieldTool_H
#define SpiceSeparatedFieldTool_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>
#include <iosfwd>

#include <N_UTL_Pack.h>

namespace Xyce {
namespace IO {

struct StringToken
{
  friend class Pack<StringToken>;

  StringToken(): lineNumber_(0), string_("")
  {}

  ~StringToken()
  {}

  size_t lineNumber_;
  std::string string_;
};

class SpiceSeparatedFieldTool
{
public:

  // Constructor
  SpiceSeparatedFieldTool(std::ifstream & input, std::string const & fileName, 
                          const std::vector< std::pair< std::string, std::string > > & externalParams );

  // Destructor
  ~SpiceSeparatedFieldTool();

  int getLine(std::vector<StringToken> & line, bool replgndvar=false,
              const std::vector<std::string> findFirstEntry = {});
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line".
  //  KRS, 10/11/07:  this new version of getLine has been added specifically
  // for the preprocess phase so as to ensure that no changes to the netlist
  // file are made during preprocessing (i.e., if '.PREPROCESS REPLACEGND 
  // true' is in the netlist file halfway through the netlist, we don't want
  // occurences of 'GND' to be replaced by '0' halfway through the netlist.
  // Technically, doing so shouldn't cause errors, but this is more of a 
  // precaution for unforseen circumstances.

  int getLineWithComments(std::vector<StringToken> & line);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O str
  // O- A string containing the next line of input.
  // Read a line and split it into fields that are stored in the
  // the private attribute "Fields". The line is stored in the
  // private member attribute "Line".
  // KRS, 12/05/07:  this is the original version of getLine, except that it
  // does NOT skip comments and blank lines (used for netlist copying).
  
  int peekAtNextLine(char& nextChar);
  // R int
  // R- Returns 1 if end-of-file has not been reached, 0 otherwise.
  // O char
  // O- A character containing the next character of input.
  // Read the first character of the next line that isn't a comments or blank.
  // This only peeks at the character, if the line needs to be processed, call
  // getLine, otherwise call skipToEndOfLine(). 

  void changeCursorLineNumber(int token);
  // Increase or decrease the line number of cursor by the given amount.
  // Increase if token is positive; decrease if token is negative.
  // If the result cursor line number is less than 1, set it to 1.

  // accessors
  void setLineNumber( int loc ) 
  { 
    ( loc <= 0 ) ? cursorLineNum_ = 1 : cursorLineNum_ = loc;
  }
  int getLineNumber() { return cursorLineNum_; }

  // Return the current character position in file.
  std::streampos getFilePosition() const;

  // Set the location in the ifstream at which the next input operation
  // will begin.
  bool setLocation(std::streampos const& startLocation);

  const std::string & getFileName() const {return fileName_;};

  void substituteExternalParams(std::vector<StringToken>& line);

  void skipToEndOfLine();

private:
  std::ifstream & in_;

  std::string fileName_;
  size_t cursorLineNum_; //The physical line number of the cursor
  std::vector< std::pair< std::string, std::string > > externalParams_;

  // We'll keep a (C-style) vector of pointers to sets of allowed
  // vector-composite names.  These will also be stored in a pair of maps
  // from model or instance names to allowed sets
  static const int numCompositeSets_=8;
  unordered_set<std::string> *allowedVectorCompositeSets_[numCompositeSets_];
  unordered_map<std::string,unordered_set<std::string>* > modelsWithVC_;
  unordered_map<std::string,unordered_set<std::string>* > yDevicesWithVC_;
  
  // R bool
  // R-
  bool NextChar_(char & c);
  void skipCommentsAndBlankLines_();
  void initializeVCMaps_();
  
  bool withinQuote_;
};

} // namespace IO
} // namespace Xyce

#endif // SpiceSeparatedFieldTool_H 
