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
// Purpose        : Define an N_IO_FunctionBlock instantiations of which are
//                  associated with netlist .FUNC lines.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 12/26/2001
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <algorithm>

// ----------   Xyce Includes   ----------
#include <N_IO_CircuitBlock.h>
#include <N_UTL_ExtendedString.h>
#include <N_IO_FunctionBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Expression.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::FunctionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
FunctionBlock::FunctionBlock(
    std::string const& fileName,
    TokenVector const& parsedInputLine)
 : netlistLocation_(fileName, parsedInputLine[0].lineNumber_)
{
  int len;

  len = parsedInputLine[parsedInputLine.size()-1].string_.size();
  if (parsedInputLine[parsedInputLine.size()-1].string_.substr(0,1) != "{" ||
      parsedInputLine[parsedInputLine.size()-1].string_.substr(len-1,1) != "}") {
    Report::UserFatal0().at(netlistLocation_)
      << "In .func line for function: " << parsedInputLine[1].string_ << ", expression must be enclosed by curly braces";
  }

  // Now extract the data from this parsed line.
  extractData(parsedInputLine);
}

//----------------------------------------------------------------------------
// Function       : FunctionBlock::FunctionBlock
// Purpose        : copy constructor
// Special Notes  : 
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/26/2003
//----------------------------------------------------------------------------
FunctionBlock::FunctionBlock( FunctionBlock const& rhsFB )
  : functionName(rhsFB.functionName),
    functionNameAndArgs(rhsFB.functionNameAndArgs),
    functionArgs(rhsFB.functionArgs),
    functionBody(rhsFB.functionBody),
    netlistLocation_(rhsFB.netlistLocation_)
{
}

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::print
// Purpose       : Output the details of a function block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
void FunctionBlock::print()
{
  Xyce::dout() << std::endl
               << "Function Information" << std::endl
               << "--------------------" << std::endl
               << "  name: " << functionName << std::endl
               << "  name and args: " << functionNameAndArgs << std::endl
               << "  body: " << functionBody << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::extractData
// Purpose       : Extract function data from parsed line formed from a
//                 netlist .func statement.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/26/2001
//-----------------------------------------------------------------------------
bool FunctionBlock::extractData(TokenVector const& parsedInputLine)
{
  ExtendedString ES1("");

  // Set the function name.
  ES1 = parsedInputLine[1].string_;
  ES1.toUpper();
  functionName = ES1;

  int arg_start = 2; // Start position of the function's argument list.
  int iend = parsedInputLine.size();
  int arg_end = iend - 2; // End position of the argument list.

  // The argument list must be enclosed by parentheses.
  if ( (parsedInputLine[arg_start].string_ != "(") || 
       (parsedInputLine[arg_end].string_ != ")") )
  {
    Report::UserFatal0().at(netlistLocation_.getFilename(), parsedInputLine[arg_start].lineNumber_)
      << ".FUNC argument list must be enclosed by parentheses in function " << functionName;
  }

  // Collect up the function name and arguments (with enclosing parentheses),
  // since the expression class needs functions in this form.
  ES1 = functionName + "(";

  // Get the list of arguments.
  for ( int i = arg_start+1; i <= arg_end-1; ++i )
  {
    ES1 += parsedInputLine[i].string_;

    if (parsedInputLine[i].string_ != ",") 
      functionArgs.push_back(ExtendedString(parsedInputLine[i].string_).toUpper());
  }

  // Add the closing parenthese.
  ES1 += ")";

  ES1.toUpper();
  functionNameAndArgs = ES1;

  // Get the function body.
  ES1 =  parsedInputLine[iend - 1].string_;
  ES1.toUpper();
  functionBody = ES1;

  return true; // Only get here on success.
}


} // namespace IO

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
int
Pack<IO::FunctionBlock>::packedByteCount(
  const IO::FunctionBlock &     function_block)
{
  int byteCount = 0;
  int size, j;

  // count functionName
  byteCount += sizeof( int );
  byteCount += function_block.functionName.length();

  // count functionNameAndArgs
  byteCount += sizeof( int );
  byteCount += function_block.functionNameAndArgs.length();

  // count functionArgs
  size = function_block.functionArgs.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j)
  {
    byteCount += sizeof( int );
    byteCount += function_block.functionArgs[j].length();
  }

  // count functionBody
  byteCount += sizeof( int );
  byteCount += function_block.functionBody.length();  

  //----- count netlistLocation_
  byteCount += 2*sizeof(int);

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::pack
// Purpose       : Packs function block into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::FunctionBlock>::pack(
  const IO::FunctionBlock &     function_block,
  char *                        buf,
  int                           bsize,
  int &                         pos,
  N_PDS_Comm *                  comm )
{
  int size, length;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+Pack<IO::FunctionBlock>::packedByteCount( function_block );
#endif

  // pack functionName;
  length = function_block.functionName.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( function_block.functionName.c_str(), length, buf, bsize, pos );

  // pack functionNameAndArgs;
  length = function_block.functionNameAndArgs.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( function_block.functionNameAndArgs.c_str(), length, buf, bsize, pos );

  // pack functionArgs
  size = function_block.functionArgs.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( int k = 0; k < size; ++k)
  {
    length = function_block.functionArgs[ k ].length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( function_block.functionArgs[ k ].c_str(), length, buf, bsize, pos );
  }

  // pack functionBody;
  length = function_block.functionBody.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( function_block.functionBody.c_str(), length, buf, bsize, pos );

  //----- pack netlistLocation_
  int file_number = function_block.netlistLocation_.getFileNumber();
  comm->pack(&file_number, 1, buf, bsize, pos );
  int line_number = function_block.netlistLocation_.getLineNumber();
  comm->pack(&line_number, 1, buf, bsize, pos );

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    Report::DevelFatal() << "FunctionBlock::pack - predicted pos (" << predictedPos
                         << ") does not match actual pos (" << pos << ")";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : FunctionBlock::unpack
// Purpose       : Unpacks function block from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::FunctionBlock>::unpack(
  IO::FunctionBlock &   function_block,
  char *                pB,
  int                   bsize,
  int &                 pos,
  N_PDS_Comm *          comm)
{
  int length, size, j;

  // unpack function name
  comm->unpack( pB, bsize, pos, &length, 1 );
  function_block.functionName = std::string( ( pB + pos ), length );
  pos += length;

  // unpack function name and args
  comm->unpack( pB, bsize, pos, &length, 1 );
  function_block.functionNameAndArgs = std::string( ( pB + pos ), length );
  pos += length;

  // unpack function args
  comm->unpack( pB, bsize, pos, &size, 1 );
  for( j = 0; j < size; ++j )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    function_block.functionArgs.push_back( std::string( ( pB + pos ), length ) );
    pos += length;
  }

  // unpack function body
  comm->unpack( pB, bsize, pos, &length, 1 );
  function_block.functionBody = std::string( ( pB + pos ), length );
  pos += length;

  //----- unpack netlistLocation_
  int file_number = 0;
  comm->unpack( pB, bsize, pos, &file_number, 1 );
  function_block.netlistLocation_.setFileNumber(file_number);
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  function_block.netlistLocation_.setLineNumber(line_number);
}

} // namespace Xyce
