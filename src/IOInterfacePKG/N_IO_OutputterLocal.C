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

#include <string>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_Outputter.h>
#include <N_IO_OutputterLocal.h>
#include <N_IO_Op.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------

namespace {

//-----------------------------------------------------------------------------
// Function      : hasExtension
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 17 14:49:26 2014
//-----------------------------------------------------------------------------
///
/// Return true if path ends with extension
///
/// @param path         path to check if ends with extension
/// @param extension    extension searched for
///
/// @return true if path ends with extension
///
///
bool hasExtension(const std::string &path, const std::string &extension)
{
  return path.size() >= extension.size() &&
    path.compare(path.size() - extension.size(), extension.size(), extension) == 0;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : outputFilename
// Purpose       : return a string for the name of the output file.
// Special Notes : The "precedence" for constructing the output "base filename 
//                 is to use the name given by the -r command line option, if
//                 specified, and the print format supports raw file output.
//                 The second choice is to use the name given by the -o
//                 command line option, if specified. (See SON Bug 911 for more
//                 details.  The third choice is the name, if any, given by FILE=  
//                 on the .print line.  Otherwise, use the netlist name as 
//                 the base filename and add suffixes or extensions as 
//                 determined by the analysis/print type.
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 6/14/2013
//-----------------------------------------------------------------------------
std::string
outputFilename(
  const std::string &           filename,
  const std::string &           default_extension,
  const std::string &           suffix,
  const std::string &           net_list_filename,
  const std::string &           overrideRawFilename,
  const bool &                  formatSupportsOverrideRaw,
  const std::string &           dashoFilename,
  const bool &                  fallbackPrintLine)
{
  std::string output_filename;

  // implement the precedence specified above in Special Notes
  if ( !overrideRawFilename.empty() && formatSupportsOverrideRaw )
  {
    output_filename = overrideRawFilename;
  }
  else if (!dashoFilename.empty())
  {
    output_filename = dashoFilename;
  }
  else
  {
    if (filename.empty())
    {
      output_filename = net_list_filename;
    }
    else
    {
      // Use filename (from FILE= on .PRINT line) if this is NOT a fallback print 
      // line.  If it is a fallback then add the default extension to filename.
      output_filename = !fallbackPrintLine ? filename : filename + default_extension;
    }
  }

  std::string::size_type pos = output_filename.find('.');
  std::string extension;

  // If filename ends with the default extension, rip it off and set the extension to the default extension.
  if (hasExtension(output_filename, default_extension))
  {
    output_filename = output_filename.substr(0, output_filename.size() - default_extension.size());
    extension = default_extension;
  }
  // Otherwise, if the filename is the netlist file, or it ends in '.cir' or it doesn't 
  // contain any extension, set the extension to the default extension.
  else if ((output_filename == net_list_filename) || hasExtension(output_filename, ".cir"))
  {
    if (default_extension.size() > 0)
      extension = default_extension;
    else
      extension=".xyce_unknown_output";
  }
  else if (pos != std::string::npos) {
    extension = output_filename.substr(pos);
    output_filename = output_filename.substr(0, pos);
  }

  // error condition if FILE= qualifier for a .PRINT line that doesn't support 
  // raw output is being sent to the file specified with the -r command line option
  if (!overrideRawFilename.empty() && !formatSupportsOverrideRaw &&
      ((output_filename + suffix + extension) == overrideRawFilename) )
  {
    Report::UserFatal() << "Conflict between FILE= on .PRINT line and -r command line option";
  }

  // Stick the suffix between the filename and extension
  return output_filename + suffix + extension;
}

//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : Given print parameters and a stream, print the header
// Special Notes : top level function
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const PrintParameters &print_parameters)
{
  return printHeader(os, print_parameters.table_.columnList_, print_parameters.delimiter_);
}

//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : Given stream, column list, and delimiter, print header
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list, const std::string &delimiter)
{
  for (Table::ColumnList::const_iterator it = column_list.begin(); it != column_list.end(); ++it)
  {
    if (it != column_list.begin())
      os << (delimiter.empty() ? " " : delimiter);

    printHeader(os, (*it));
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : print a single column of the header on the given stream.
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const Table::Column &column)
{
  std::string name = column.name_;
  if (name == "INDEX")
    name = "Index";

  size_t left_padding = 0;
  size_t right_padding = 0;

  // if (column.width_ < name.size())
  //   column.width_ = name.size();

  if (column.width_ > name.size())
  {
    switch (column.justification_)
    {
      case Table::JUSTIFICATION_LEFT:
        right_padding = column.width_ - left_padding - name.size();
        break;
      case Table::JUSTIFICATION_CENTER:
        left_padding = (column.width_ - column.name_.size())/2;
        right_padding = column.width_ - left_padding - name.size();
        break;
      case Table::JUSTIFICATION_RIGHT:
        left_padding = column.width_ - column.name_.size();
        break;
      case Table::JUSTIFICATION_NONE:
        // this empty case is here just to shut up warnings from clang about unhandled enum cases
        break;
    }
  }

  os << std::setw(left_padding) << "" << std::setw(0) << name << std::setw(right_padding) << "";

  return os;
}

//-----------------------------------------------------------------------------
// Function      : createOps
// Purpose       : given an output manager and a (print) parameter list, and
//                 a back_inserter iterator for an OpList, generate all the
//                 "Ops" needed to obtain the values in the print parameter
//                 list, and add to the end of the OpList associated
//                 with the iterator.
// Special Notes : If we're in frequency domain, solution var access
//                 such as V(A) gets expanded into two ops, one for real part,
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
createOps(
  Parallel::Machine                             comm,
  const Util::Op::BuilderManager &              op_builder_manager,
  bool                                          expandComplexTypes,
  const double                                  time_scale_factor,
  const NetlistLocation &                       netlist_location,
  Util::ParamList::const_iterator               begin,
  Util::ParamList::const_iterator               end,
  std::back_insert_iterator<Util::Op::OpList>   inserter)
{
  Util::Op::OpList tempOpList;

  makeOps(comm, op_builder_manager, netlist_location, begin, end, std::back_inserter(tempOpList));

  for (Util::Op::OpList::const_iterator it = tempOpList.begin(); it != tempOpList.end(); ++it)
  {
    if (time_scale_factor != 1.0 && (*it)->id() == Util::Op::identifier<OutputMgrTimeOp>())
    {
      *inserter++ = new OutputMgrTimeOp((*it)->getName(), (static_cast<const OutputMgrTimeOp *>(*it))->outputMgr_, time_scale_factor);
      delete *it;
    }
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<SolutionOp>())
    {
      std::string solutionName = (*it)->getName();
      int index = -1;
      const SolutionOp *op = dynamic_cast<const SolutionOp *>(*it);
      if (op)
        index = op->index_;

      delete *it;
      *inserter++ = new SolutionRealOp("Re(" + solutionName + ")", index);
      *inserter++ = new SolutionImaginaryOp("Im(" + solutionName + ")", index);
    }
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<StoreOp>())
    {
      std::string storeName = (*it)->getName();
      int index = -1;
      const StoreOp *op = dynamic_cast<const StoreOp *>(*it);
      if (op)
        index = op->index_;

      delete *it;
      *inserter++ = new StoreRealOp("Re(" + storeName + ")", index);
      *inserter++ = new StoreImaginaryOp("Im(" + storeName + ")", index);
    }
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<BranchDataCurrentOp>())
    {
      std::string branchDataCurrentName = (*it)->getName();
      int index = -1;
      const BranchDataCurrentOp *op = dynamic_cast<const BranchDataCurrentOp *>(*it);
      if (op)
        index = op->index_;

      delete *it;
      *inserter++ = new BranchDataCurrentRealOp("Re(" + branchDataCurrentName + ")", index);
      *inserter++ = new BranchDataCurrentImaginaryOp("Im(" + branchDataCurrentName + ")", index);
    }
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<VoltageDifferenceOp>())
    {
      std::string solutionName = (*it)->getName();
      int index1 = -1;
      int index2 = -1;
      const VoltageDifferenceOp *op =
        dynamic_cast<const VoltageDifferenceOp *>(*it);
      if (op)
      {
        index1 = op->index1_;
        index2 = op->index2_;
      }

      delete *it;
      *inserter++ = new VoltageDifferenceRealOp("Re(" + solutionName + ")", index1, index2);
      *inserter++ = new VoltageDifferenceImaginaryOp("Im(" + solutionName + ")", index1, index2);
    }
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<RFparamsOp>())
    {
      std::string RFparamsName = (*it)->getName();
      std::string type;
      int index1 = -1;
      int index2 = -1;
      const RFparamsOp *op =
        dynamic_cast<const RFparamsOp *>(*it);
      if (op)
      {
        type = op->type_;
        index1 = op->index1_;
        index2 = op->index2_;
      }

      delete *it;
      *inserter++ = new RFparamsRealOp("Re(" + RFparamsName + ")", type, index1, index2);
      *inserter++ = new RFparamsImaginaryOp("Im(" + RFparamsName + ")", type, index1, index2);
    }
#if 0
    // Not sure if this is good idea.    This was added as a candidate fix for
    // bug 1032 on the SON, entitled "The V(a,b) operator may not print out 
    // correctly when used in expressions for .AC analyses".  But having added
    // it and played with it, I now think it is not a good choice.
    //
    // The problem with the below implementation is that it applies the Re and Im 
    // prefixes unconditionally to an expression-based output.  Some expressions 
    // will clearly have both real and imaginary components and some never will.
    //
    // For example, one would never need Im(ONOISE), as ONOISE is always real-valued.
    //
    // The new expression library gives you the power to request these things. It would
    // allow .print  AC {Re(V(A))} {Img(V(A))}
    //
    else if (expandComplexTypes && (*it)->id() == Util::Op::identifier<ExpressionOp>())
    {
      const ExpressionOp *op = dynamic_cast<const ExpressionOp *>(*it);
      if (op)
      {
        *inserter++ = new ExpressionRealOp( *op );
        *inserter++ = new ExpressionImaginaryOp( *op );
        delete *it;
      }
    }
#endif
    else
      *inserter++ = *it;
  }
}

//-----------------------------------------------------------------------------
// Function      : fixupColumns
// Purpose       : given an output manager and a (print) parameter list, and
//                 an OpList, fill the list with the ops needed via createOps,
//                 and add additional output columns for things like
//                 index, time, frequency.  Set formatting options for these
//                 additional columns, set delimiter.
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void fixupColumns(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager, PrintParameters &print_parameters, Util::Op::OpList &op_list)
{
  createOps(comm, op_builder_manager, print_parameters.expandComplexTypes_, print_parameters.outputTimeScaleFactor_, print_parameters.netlistLocation_, print_parameters.variableList_.begin(), print_parameters.variableList_.end(), std::back_inserter(op_list));

  Table::Justification justification = print_parameters.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

  for (Util::Op::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
  {
    if ((*it)->id() == Util::Op::identifier<StepNumOp>())
    {
      print_parameters.table_.addColumn("STEPNUM", std::ios_base::fixed, 8, 0, Table::JUSTIFICATION_LEFT);
    }
    else if ((*it)->id() == Util::Op::identifier<CurrentIndexOp>())
    {
      print_parameters.table_.addColumn("INDEX", std::ios_base::fixed, 5, 0, Table::JUSTIFICATION_LEFT);
    }
    else
    {
      print_parameters.table_.addColumn((*it)->getName(), print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : fixupColumnsFromStrVec
// Purpose       : Inputs are a (print) parameter list, and a string vector of
//                 column names rather than an opList.  This function then
//                 sets the formatting options for those named columns.  It
//                 also sets the delimiter.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/27/19
//-----------------------------------------------------------------------------
void fixupColumnsFromStrVec(Parallel::Machine comm, PrintParameters &print_parameters, std::vector<std::string> & colNames)
{
  Table::Justification justification = print_parameters.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;
  for (std::vector<std::string>::const_iterator it = colNames.begin() ; it != colNames.end(); ++it)
  {
    print_parameters.table_.addColumn(*it, print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
  }
}

//-----------------------------------------------------------------------------
// Function      : printValue
// Purpose       :
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printValue(std::ostream &os, const Table::Column &column, const std::string &delimiter, const int column_index, double value)
{
  if (delimiter.empty())
  {
    if (column_index != 0)
      os <<  " ";
    os << std::resetiosflags(std::ios_base::floatfield) << std::setiosflags(column.format_)
       << std::resetiosflags(std::ios_base::adjustfield)
       << std::setiosflags(column.justification_ == Table::JUSTIFICATION_LEFT ? std::ios_base::left : std::ios_base::right)
       << std::setprecision(column.precision_) << std::setw(column.width_)
       << value;
  }
  else
  {
    if (column_index != 0)
      os << delimiter;
    os << std::resetiosflags(std::ios_base::floatfield) << std::setiosflags(column.format_)
      << std::setw(0) << std::setprecision(column.precision_) << value;
  }

  return os;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
