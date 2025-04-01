//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_Outputter.h>
#include <N_IO_OutputMgr.h>

#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : getTecplotTimeDateStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:20:41 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
///
/// @return 
///
///
std::string getTecplotTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME= \" %I:%M:%S %p %b %d, %Y \" ", localtime( &now));

  return std::string( timeDate);
}

//-----------------------------------------------------------------------------
// Function      : tecplotTimeHeader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:20:35 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
/// @param os 
/// @param print_title 
/// @param title 
/// @param op_list 
/// @param output_manager 
///
///
void tecplotTimeHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);


  if (print_title)
  {
    std::string localTitle(title);
    std::size_t loc=0;
    while( (loc = localTitle.find_first_of( '"',loc ) ) != std::string::npos )
    {
      localTitle.insert(loc, 1, '\\');
      loc+=2;
    }
    os << "TITLE = \"" << localTitle << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
      os << "\" " << (*it)->getName() << "\" " << std::endl;

    // output some AUXDATA
    os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;

    if (!output_manager.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << output_manager.getCircuitTemp() << " \"" << std::endl;
    }

  } // print header calls=0

  os << "ZONE F=POINT ";

  if ( output_manager.getStepSweepVector().empty())
  {
    //os << "T=\"Xyce data - " <<  output_manager.getNetlistFilename() << " \" ";
    os << "T=\"" <<  output_manager.getNetlistFilename() << " \" ";
  }
  else
  {
    os << "T= \" ";
    int maxParams=10;
    if (maxParams > output_manager.getStepSweepVector().size()) { maxParams = output_manager.getStepSweepVector().size(); }
    for (int ip=0;ip<maxParams;ip++)
    {
      std::vector<Analysis::SweepParam>::const_iterator it = output_manager.getStepSweepVector().begin() + ip;
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }
  os << std::endl;

  // put in the various sweep parameters as auxdata, as long as there are not too many of them:
  if (!output_manager.getStepSweepVector().empty())
  {
    int maxParams=10;
    if (maxParams > output_manager.getStepSweepVector().size()) { maxParams = output_manager.getStepSweepVector().size(); }
    for (int ip=0;ip<maxParams;ip++)
    {
      std::vector<Analysis::SweepParam>::const_iterator it = output_manager.getStepSweepVector().begin() + ip;
      // convert any ":", "%" or "!" in the name to a "_", so as not to confuse tecplot.
      std::string name(it->name);
      std::replace(name.begin(), name.end(), '%', '_');
      std::replace(name.begin(), name.end(), ':', '_');
      std::replace(name.begin(), name.end(), '!', '_');
      os << "AUXDATA " << name << " = " << "\" " << it->currentVal << "\" ";
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : tecplotFreqHeader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:21:00 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
/// @param os 
/// @param print_title 
/// @param title 
/// @param op_list 
/// @param output_manager 
///
///
void tecplotFreqHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);

  if (print_title)
  {
    std::string localTitle(title);
    std::size_t loc=0;
    while( (loc = localTitle.find_first_of( '"',loc ) ) != std::string::npos )
    {
      localTitle.insert(loc, 1, '\\');
      loc+=2;
    }
    os << " TITLE = \" Xyce Frequency Domain data, " << localTitle << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
    {
      os << "\" ";
      os << (*it)->getName() ;
      os << "\" " << std::endl;
    }
    os << "DATASETAUXDATA ";
    os << getTecplotTimeDateStamp();
    os << std::endl;

    if (!output_manager.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << output_manager.getCircuitTemp() << " \"" << std::endl;
    }
  }

  // output some AUXDATA
  os << "ZONE F=POINT  ";

  if ( output_manager.getStepSweepVector().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    int maxParams=10;
    if (maxParams > output_manager.getStepSweepVector().size()) { maxParams = output_manager.getStepSweepVector().size(); }
    for (int ip=0;ip<maxParams;ip++)
    {
      std::vector<Analysis::SweepParam>::const_iterator it = output_manager.getStepSweepVector().begin() + ip;
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata, as long as there are not too many of them:
  if (!output_manager.getStepSweepVector().empty())
  {
    int maxParams=10;
    if (maxParams > output_manager.getStepSweepVector().size()) { maxParams = output_manager.getStepSweepVector().size(); }
    for (int ip=0;ip<maxParams;ip++)
    {
      std::vector<Analysis::SweepParam>::const_iterator it = output_manager.getStepSweepVector().begin() + ip;
      // convert any ":", "%" or "!" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(it->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      replace(tmpName.begin(), tmpName.end(), '!', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << it->currentVal << "\" ";
    }
    os << std::endl;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
