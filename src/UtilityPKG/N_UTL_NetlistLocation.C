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
// Purpose       :
//
// Special Notes :
//
// Creator       : David G. Baur
//
// Creation Date : 1/24/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <algorithm>

#include <N_UTL_fwd.h>
#include <N_UTL_Pack.h>
#include <N_UTL_NetlistLocation.h>

#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>


namespace Xyce {

NetlistLocation::NetlistLocation()
: fileNumber_(0),
  lineNumber_(0)
{
  std::string empty;
  fileNumber_ = Filename::getFileNumber(empty);
}

NetlistLocation::NetlistLocation(const std::string &filename, const int line_number)
: fileNumber_(Filename::getFileNumber( filename )),
  lineNumber_(line_number)
{}
  
NetlistLocation & NetlistLocation::setFilename(const std::string &filename)
{ 
  fileNumber_ = Filename::getFileNumber( filename );
  return *this;
}
  
const std::string & NetlistLocation::getFilename() const
{ 
  return Filename::getFilename( *this );
}

bool operator<(const NetlistLocation &s0, const NetlistLocation &s1)
{
  return s0.getFilename() < s1.getFilename() || (!(s1.getFilename() < s0.getFilename()) && s0.getLineNumber() < s1.getLineNumber());
}

std::ostream &operator<<(std::ostream &os, const NetlistLocation &x)
{
  os << "file " << x.getFilename() << " at or near line " << x.getLineNumber();
  return os;
}

Util::Marshal &operator<<(Util::Marshal &mout, const NetlistLocation &netlist_location) 
{
  return mout << netlist_location.getFilename() << netlist_location.getLineNumber();
}

Util::Marshal &operator>>(Util::Marshal &min, NetlistLocation &netlist_location) {
  std::string filename;
  int line_number;
  min >> filename >> line_number;
  netlist_location.setFilename(filename);
  netlist_location.setLineNumber(line_number);
  return min;
}

// This is a simple struct to hold the map between the filename and the filenumber.
struct FileData
{ 
  Filename::FilenameVector    filenameVector_;             ///<Filename registration
  
  ~FileData() {}
};

//-----------------------------------------------------------------------------
// Function      : getFileData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Feb  5 11:49:42 2014
//-----------------------------------------------------------------------------
/// returns the configuration data singleton.
///
/// @return reference to the configuration data singleton.
///
FileData &getFileData()
{ 
  static FileData data_;
  
  return data_;
}

//-----------------------------------------------------------------------------
// Function      : Filename::getFilenameVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : Fri Mar 14 13:23:31 2014
//-----------------------------------------------------------------------------
const Filename::FilenameVector &
Filename::getFilenameVector()
{
  return getFileData().filenameVector_;
}

//-----------------------------------------------------------------------------
// Function      : Filename::getFilename
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : Fri Mar 14 13:23:52 2014 
//-----------------------------------------------------------------------------
const std::string& 
Filename::getFilename(const NetlistLocation& loc)
{
  return getFileData().filenameVector_[loc.getFileNumber()];
}

//-----------------------------------------------------------------------------
// Function      : Filename::getFilenumber
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : Fri Mar 14 13:23:52 2014 
//-----------------------------------------------------------------------------
int 
Filename::getFileNumber(const std::string& filename)
{
  FilenameVector::iterator it;
  FilenameVector& fileVector = getFileData().filenameVector_;  
  it = std::find(fileVector.begin(), fileVector.end(), filename);
  if ( it != fileVector.end() )
  {
    return std::distance( fileVector.begin(), it );
  }
  else
  {
    fileVector.push_back( filename ); 
    return (fileVector.size()-1);
  }
}

//-----------------------------------------------------------------------------
// Function      : FilenameVector::packedByteCount
// Purpose       : Counts bytes needed to pack filename vector.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 02/08/2019
//-----------------------------------------------------------------------------
template<>
int
Pack<Filename::FilenameVector>::packedByteCount(
  const Filename::FilenameVector & filename_obj)
{
  int byteCount = 0;

  //----- count name
  byteCount += sizeof(int);

  //----- count filename size
  Filename::FilenameVector::const_iterator it_fV = filename_obj.begin();
  Filename::FilenameVector::const_iterator it_fVE = filename_obj.end();
  for ( ; it_fV != it_fVE; ++it_fV )
  {
    byteCount += sizeof(int);
    byteCount += (*it_fV).length();
  } 
 
  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : FilenameVector::pack
// Purpose       : Packs FilenameVector into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 02/08/2019
//-----------------------------------------------------------------------------
template<>
void
Pack<Filename::FilenameVector>::pack(const Filename::FilenameVector &filename_obj,
                                     char * buf, int bsize, int & pos, Parallel::Communicator * comm)
{ 
  if (DEBUG_IO)
    Xyce::dout() << "Packing FilenameVector "<< std::endl;
  
  int size, length;
 
  //----- pack length of filename vector
  size = filename_obj.size();
  comm->pack( &size, 1, buf, bsize, pos );
 
  //----- pack filenames
  Filename::FilenameVector::const_iterator it_fV = filename_obj.begin();
  Filename::FilenameVector::const_iterator it_fVE = filename_obj.end();
  for ( ; it_fV != it_fVE; ++it_fV )
  {
    length = (*it_fV).length();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( (*it_fV).c_str(), length, buf, bsize, pos );
  }
 
  if (DEBUG_IO)
    Xyce::dout() << "Packed " << pos << " of " << bsize << " bytes for FilenameVector " << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FilenameVector::unpack
// Purpose       : Unpacks FilenameVector into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 02/08/2019
//-----------------------------------------------------------------------------
template<>
void
Pack<Filename::FilenameVector>::unpack(Filename::FilenameVector &filename_obj,
                                       char * pB, int bsize,int & pos, Parallel::Communicator * comm)
{
  int size = 0;
  int length = 0;
  int i;

  //----- unpack length of filename vector
  comm->unpack( pB, bsize, pos, &size, 1 );

  if (DEBUG_IO)
    Xyce::dout() << "Unpacking FilenameVector of length " << size << std::endl;

  //----- unpack filenames
  for ( int i=0; i<size; ++i )
  {
    comm->unpack( pB, bsize, pos, &length, 1 );
    std::string newFilename( (pB+pos), length );
    pos += length;
    
    // Now check this filename against current filename vector.
    // Calling getFileNumber() inserts newFilename, if it doesn't already exist, and returns the index for it.
    int newIndex = Filename::getFileNumber( newFilename );
    if (newIndex != i)
      std::cout << "Processor " << comm->procID() << " has index " << newIndex << " instead of " << i << " for filename " << newFilename << std::endl;
  }

  if (DEBUG_IO)
    Xyce::dout() << "Unpacked " << pos << " of " << bsize << " bytes for FilenameVector " << std::endl;
}

} // namespace Xyce

