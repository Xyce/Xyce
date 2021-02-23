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
// Purpose       : Base class for re-reading Xyce file output of simulation results
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 12/5/2012
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_IO_OutputFileBase.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Vector.h>
#include <N_UTL_DeviceNameConverters.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::OutputFileBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
OutputFileBase::OutputFileBase() :
  ostreamPtr_(NULL),
  outputFileBaseName_(""),
  fileSuffix_(""),
  simulationSuffix_(""),
  fileFormatName_("FormatUndefined"),
  appendOutputFlag_(false)    
{};


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::~OutputFileBase
// Purpose       : Destructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
OutputFileBase::~OutputFileBase()
{
  if ( (ostreamPtr_ != &std::cout) && (ostreamPtr_ != NULL) ) 
  {
    Report::DevelFatal().in("OutputFileBase::~OutputFileBase()") << "Non-null ostreamPtr_ from " << fileFormatName_ << " derived class.";
  }
  
  if( !(istreamPtr_ != & std::cin) && (istreamPtr_ != 0)  )
  { 
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
}



//-----------------------------------------------------------------------------
// Function      : OutputFileBase::openFile
// Purpose       : open output file for writing.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------
void OutputFileBase::openFile( std::string basename, std::string simulationSuffix )
{
  simulationSuffix_ = simulationSuffix;
  outputFileBaseName_ = basename;
  fullFileName_ = outputFileBaseName_ + simulationSuffix_ + fileSuffix_;
  
  if( ostreamPtr_ == 0 )
  {
    // only try to open a new ostream if the pointer is NULL
    if (fullFileName_ == "CONSOLE")
    {
      ostreamPtr_ = &std::cout;
    }
    else
    {
      if( appendOutputFlag_ == true ) 
      {
        ostreamPtr_ = new std::ofstream(fullFileName_.c_str(), std::ios_base::app );
      }
      else
      {
        ostreamPtr_ = new std::ofstream(fullFileName_.c_str());
      }
    }
    if( !ostreamPtr_ )
    {
      std::string msg = "Could not open file, \""  + fullFileName_ + "\" for output.";
      Report::UserFatal() << msg;
    }
  }
  else
  {
    Report::DevelFatal().in("void OutputFileBase::openFile( std::string basename, std::string simulationSuffix )")
      << "ostreamPtr_ was not NULL.";
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::closeOutput
// Purpose       : close the output file.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 12/05/12
//-----------------------------------------------------------------------------   
void OutputFileBase::closeOutput()
{
  // Try just closing the output file.  Report any errors.
  if ( ostreamPtr_ != &std::cout && ostreamPtr_ ) 
  {
    delete ostreamPtr_;
    ostreamPtr_ = 0;
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::openFileForRead
// Purpose       : opens an existing output file for reading.
// Special Notes : returns false if open fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::openFileForRead( std::string filename )
{
  bool retVal=true;
  istreamPtr_ = new std::ifstream(filename.c_str());
  retVal = istreamPtr_->is_open();
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::getOutputVarNames
// Purpose       : Read the header for far names
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::getOutputVarNames(std::vector<std::string> & varNames)
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::getOutputNextVarValuesParallel
// Purpose       : Get one row of simulation data
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
int OutputFileBase::getOutputNextVarValuesParallel(Linear::Vector * varValues)
{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::getOutputNextVarValuesSerial
// Purpose       : Get one row of simulation data
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical Systems Modeling
// Creation Date : 08/22/16
//----------------------------------------------------------------------------- 
int OutputFileBase::getOutputNextVarValuesSerial(std::vector<double> * varValues)
{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : OutputFileBase::closeFileForRead
// Purpose       : close the old output file.
// Special Notes : returns false if close fails
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/15/13
//----------------------------------------------------------------------------- 
bool OutputFileBase::closeFileForRead()
{
  bool retVal=true;
  if( istreamPtr_ != 0 )
  {
    istreamPtr_->close();
    retVal = !(istreamPtr_->is_open());
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
  if( !(istreamPtr_ != & std::cin) && (istreamPtr_ != 0)  )
  { 
    delete istreamPtr_;
    istreamPtr_ = 0;
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : OutputFileBase::convertOutputNamesToSolVarNames
// Purpose       : Converts names as they appear in an output header to
//                 variable names as they would appear in Xyce's solution 
//                 vector. E.g.: v(a) to "a", I(a) to "a_dev" and Ix(a) to "a_devx"
// Special Notes : Voltage difference v(a,b) and expressions cannot be reduced 
//                 to solution var names, so they are left as is.  Also, in Xyce, 
//                 currents through voltage sources and inductors are solution variables. 
//                 Currents through other devices may be either "lead currents" or store
//                 variables. This routine, which is specific to re-measure, ignores 
//                 those distinctions, and treats all I() operators as lead currents.
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/17/13
//----------------------------------------------------------------------------- 
void OutputFileBase::convertOutputNamesToSolVarNames( std::vector< std::string > & varNames,
                                                      std::vector<int> & solSymIdx, 
                                                      std::vector<int> & brVarSymIdx )
{
  int numNames = varNames.size();
  for( int i=0; i<numNames; i++ )
  {
    // Now I need to make the var names from the header match what would normally 
    // be found in the solutionNodeMap or the BranchVarsNodeMap.  
    // Essentially: V(name) becomes name,  I(name) must be transformed into the 
    // "lead current name" expected by the I() operator class.  

    // Expressions on the print line such ax {V(x)-V(y)} can't be reduced further 
    // than the expression they are.  This is a complication we will have to deal 
    // with at some point.
    
    if( varNames[i].length() > 1 )
    {
      std::string::size_type begBraceLoc = varNames[i].find_first_of('{');  
      std::string::size_type endBraceLoc = varNames[i].find_first_of('}');
      if( begBraceLoc != std::string::npos &&  endBraceLoc != std::string::npos )
      {
        // an expression so leave it as it is.
      }
      else 
      {
        // flag used to differentiate between V() and I() operators when the indices
        // are stored in either the solSymIdx or brVarSymIdx vectors
        bool solSym =true; 

        std::string::size_type begloc = varNames[i].find_first_of('(');  
        std::string::size_type endloc = varNames[i].find_first_of(')');
        if( begloc != std::string::npos &&  endloc != std::string::npos )
        {
		  // This takes care of v(a) and i(a) but not v(a,b) as that
	          // just becomes a,b -- can't do more with that because of
                  // how voltage difference operators are made.
		  std::string nodeName;
		  nodeName.assign(varNames[i], begloc+1, endloc-begloc-1);
		  if (varNames[i][0] == 'I' || varNames[i][0] == 'i' )
		  {
                    solSym = false;
                    // This is a bit of a hack, but treat all I() as lead currents
                    // for the purposes of re-measure.  However, we then need to use 
                    // the "Spice Name" in order to handle devices in subcircuits.
                    // We then treat it as a lead current by appending :BRANCH_D to its
                    // Spice name.  So, examples are:
                    //
                    //  Operator        nodeName
                    //    I(V4)           V4:BRANCH_D
                    //    I(X1:V4)        V:X1:4:BRANCH_D
                    //    IC(Q2)          Q2:BRANCH_DC
                    //    IC(X1:Q2)       Q:X1:2:BRANCH_DC
                    nodeName = Util::xyceDeviceNameToSpiceName(nodeName);
		    nodeName.append(":BRANCH_D");

		    // If begloc != 1 then there is an additional lead suffix 
		    // on the I(a) as in Ib(a).  An example would be lead currents
                    // for the BJT device, as shown above. 
		    if( begloc != 1 )
		    {
		     nodeName.append(varNames[i], 1, 1);
		    }
		  }
		  varNames[i].assign( nodeName );
        }
        else
        {
          // couldn't find '(' and ')' characters as part of V( x ) or I( y )
          // so just leave the name as it is.
        }

        // put the index into the correct vector, depending on whether it is being
        // treated like a solution symbol (from a V() operator) or a branch symbol
        // (from an I() operator). 
        if (solSym)
          solSymIdx.push_back(i);
        else
          brVarSymIdx.push_back(i);
      }
    }
    else
    {
      // unexpectedly short name, just leave it as is.
    }
  }

}

} // namespace IO
} // namespace Xyce
