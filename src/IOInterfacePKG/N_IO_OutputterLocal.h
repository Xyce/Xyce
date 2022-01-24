///-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterLocal_h
#define Xyce_N_IO_OutputterLocal_h

#include <N_UTL_Math.h>
#include <list>
#include <string>
#include <vector>
#include <iterator>

#include <N_IO_Outputter.h>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Demangle.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {

namespace Outputter {

// ERK 4/11/2016: the following 3 classes served to group together outputters 
// for same print types.  They were only there to provide implementations the 
// useless purely virtual functions in the Xyce::IO::Outputter base class.   
// For example, in the "TimeInterface" class, the "outputNoise" function was
// implemented here, and did nothing.
//
// I recently made these functions so they were no longer purely virtual, as
// most of them should never be called, and I was tired of all the typing.
// But, I chose to leave these classes in place for now, in case they are useful
// later, in a later IO package re-org.

class TimeInterface : public Interface
{

};

class FrequencyInterface : public Interface
{

};

class HBInterface : public Interface
{

};


 void fixupColumns(Parallel::Machine comm,
                   const Util::Op::BuilderManager &op_builder_manager,
                   PrintParameters &print_parameters,
                   Util::Op::OpList &op_list);

 void fixupColumnsFromStrVec(Parallel::Machine comm,
                             PrintParameters &print_parameters,
                             std::vector<std::string> & colNames);

 void createOps(Parallel::Machine                comm,
                const Util::Op::BuilderManager & op_builder_manager,
                bool                             expandComplexTypes,
                const double                     time_scale_factor,
                const NetlistLocation &          netlist_location,
                Util::ParamList::const_iterator  begin,
                Util::ParamList::const_iterator  end,
                std::back_insert_iterator<Util::Op::OpList> inserter);
  
std::string outputFilename(const std::string &filename,
                           const std::string &default_extension,
                           const std::string &suffix,
                           const std::string &net_list_filename,
                           const std::string &overrideRawFilename,
                           const bool &formatSupportsOverrideRaw,
                           const std::string &dashoFilename,
                           const bool &fallbackPrintLine);

std::ostream &printHeader(std::ostream &os,
                          const PrintParameters &print_parameters);
std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list,
                          const std::string &delimiter);
std::ostream &printHeader(std::ostream &os, const Table::Column &column);
std::ostream &printValue(std::ostream &os, const Table::Column &column,
                         const std::string &delimiter, const int column_index,
                         double value);

//-----------------------------------------------------------------------------
// Function      : filter
// Purpose       : Applies a filter to a double value, returning 0 if the
//                 absolute value of the data is less than the filter.
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 11/25/2013
//-----------------------------------------------------------------------------
inline double filter(double value, double filter)
{
  return std::fabs(value) < filter ? 0.0 : value;
}

} // namespace Outputter

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterLocal_h
