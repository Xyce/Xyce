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
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Peter Sholander, SNL, Electrical Models & Simulation
//
// Creation Date  : 06/01/18
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_IBIS.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_IO_DeviceBlock.h>
#include <N_IO_ParsingHelpers.h>
#include <sstream>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_Expression.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_FeatureTest.h>

#include <N_DEV_ExpressionGroupWrapper.h>

namespace Xyce {
namespace Device {

namespace IBIS {

void Traits::loadInstanceParameters(ParametricData<IBIS::Instance> &p)
{
  // Set up configuration constants:
  // Set up double precision variables:
  p.addPar("FILE", "", &IBIS::Instance::fileName_)
    .setGivenMember(&IBIS::Instance::fileName_given)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("File Name");

  p.addPar("MODEL", "", &IBIS::Instance::modelName_)
    .setGivenMember(&IBIS::Instance::modelName_given)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("Model Name");

  p.addPar ("NODELIST",std::vector<std::string>(),&IBIS::Instance::nodeList_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("");

  p.addPar("GNDCLAMPTBL", 0.0, &IBIS::Instance::gndClampI_)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("VI table for GND Clamp");

  p.addPar("PWRCLAMPTBL", 0.0, &IBIS::Instance::powerClampI_)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("VI table for Power Clamp");

  p.addPar("PULLDOWNTBL", 0.0, &IBIS::Instance::pulldownI_)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("VI table for Pulldown");

  p.addPar("PULLUPTBL", 0.0, &IBIS::Instance::pullupI_)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("VI table for Pullup");

  p.addPar("I", 0.0, &IBIS::Instance::I)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_AMP)
    .setDescription("Current for current source");

  p.addPar("TEMP", 0.0, &IBIS::Instance::temp)
    .setUnit(U_DEGC)
    .setCategory(CAT_NONE)
    .setDescription("Device temperature");
}

void Traits::loadModelParameters(ParametricData<IBIS::Model> &p)
{}

#define Xyce_NONPOINTER_MATRIX_LOAD 1

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Peter Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/01/18
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       BMiter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(BMiter),
    Exp_ptr_pc(0),
    Exp_ptr_gc(0),
    Exp_ptr_pu(0),
    Exp_ptr_pd(0),
    expNumVars(0),
    expNumVars_pc(0),
    expNumVars_gc(0),
    expNumVars_pu(0),
    expNumVars_pd(0),
    expVal(0),
    expVal_pc(0),
    expVal_gc(0),
    expVal_pu(0),
    expVal_pd(0),
    IB(IB),
    fileName_(""),
    modelName_(""),
    ibisCommentChar_('|'),
    ibisVer_("-1"),
    ibisComponent_(""),
    pkgRLCVec_(),
    Model_type(IBIS_MODEL_INVALID),
    Polarity(IBIS_POLARITY_INVALID),
    Vinl(0),
    Vinl_given(false),
    Vinh(0),
    Vinh_given(false),
    Vmeas(0),
    Vmeas_given(false),
    Rref(0),
    Rref_given(false),
    Vref(0),
    Vref_given(false),
    C_comp(0),
    C_comp_given(false),
    GND_Clamp_Reference(0),
    GND_Clamp_Reference_given(false),
    ramp_dV_dt_r_num(0),
    ramp_dV_dt_r_num_given(false),
    ramp_dV_dt_r_den(0),
    ramp_dV_dt_r_den_given(false),
    ramp_dV_dt_f_num(0),
    ramp_dV_dt_f_num_given(false),
    ramp_dV_dt_f_den(0),
    ramp_dV_dt_f_den_given(false),
    bufferModel_(),
    gndClampI_(0),
    powerClampI_(0),
    pulldownI_(0),
    pullupI_(0),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_branch_data(0),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    APosEquPosNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0)
{
  numIntVars   = 1;
  numExtVars   = IB.numExtVars;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();

  // read .ibs file
  readIbsFile();

  // hack since only two buffer models are supported now
  if ( !( (bufferModel_.Model_type == IBIS_INPUT) || 
          (bufferModel_.Model_type == IBIS_OUTPUT) ))
  {
    UserError(*this) << "Only supported IBIS buffer types are: Input and Output" << std::endl; 
  }

  // Set any non-constant parameter defaults:
  if ( given("GNDCLAMPTBL") && given("PWRCLAMPTBL") ) 
  {
    // current sources don't have current as part of the solution vector
    // so store them in the store vector
    setNumStoreVars(1);
  }
  else
  {
    UserError(*this) << "IBIS buffer model mising Ground Clamp or Power Clamp table";
  }

  // Set any non-constant parameter defaults:
  if ( (!given("PULLUPTBL") || !given("PULLDOWNTBL")) && (bufferModel_.Model_type == IBIS_OUTPUT) ) 
  {
    UserError(*this) << "IBIS output buffer model mising Pullup or Pulldown table";
  }

  numIntVars = 0;

  std::vector<Depend>::const_iterator d;
  std::vector<Depend>::const_iterator begin = getDependentParams().begin();
  std::vector<Depend>::const_iterator end = getDependentParams().end();

  for  (d = begin ; d != end ; ++d)
  {
    if (d->name == "GNDCLAMPTBL")
    {
      expNumVars = d->numVars;
      expNumVars_gc = d->numVars;
      Exp_ptr_gc = d->expr;

      if (expNumVars < 1)
      {
        UserError(*this) << "Error making Ground Clamp Table";
      }
    }
    else if (d->name == "PWRCLAMPTBL")
    {
      expNumVars_pc = d->numVars;
      Exp_ptr_pc = d->expr;
     
      if (d->numVars < 1)
      {
        UserError(*this) << "Error making Power Clamp Table";
      }
    }
    else if (d->name == "PULLUPTBL")
    {
      expNumVars_pu = d->numVars;
      Exp_ptr_pu = d->expr;
     
      if (d->numVars < 1)
      {
        UserError(*this) << "Error making Pullup Table";
      }
    }
    else if (d->name == "PULLDOWNTBL")
    {
      expNumVars_pd = d->numVars;
      Exp_ptr_pd = d->expr;
     
      if (d->numVars < 1)
      {
        UserError(*this) << "Error making Pulldown Table";
      }
    }
  }

  if( jacStamp.empty() )
  {
    jacStamp.resize( 2 );
    jacStamp[0].resize(expNumVars);
    jacStamp[1].resize(expNumVars);
    for( int i = 0; i < expNumVars; ++i )
    {
      jacStamp[0][i] = i+2;
      jacStamp[1][i] = i+2;
    }
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/30/03
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Instance::updateTemperature(const double & temp_tmp)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Peter Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/01/18
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Peter Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/01/18
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  IBISInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.
  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
  }

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The IBIS device uses exactly one "branch data vector" element, where
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current if this is a current source
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 01/17/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}


//-----------------------------------------------------------------------------
// Function      : Instance::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & Instance::getDepSolnVars()
{
  return DeviceInstance::getDepSolnVars();
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs(
  const std::vector< std::vector<int> > & jacLIDVec)
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  
  APosEquExpVarOffsets.resize( expNumVars );
  ANegEquExpVarOffsets.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    APosEquExpVarOffsets[i] = jacLIDVec[0][i];
    ANegEquExpVarOffsets[i] = jacLIDVec[1][i];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquExpVarPtrs.resize( expNumVars );
  fNegEquExpVarPtrs.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    fPosEquExpVarPtrs[i]  = &(dFdx[li_Pos][ APosEquExpVarOffsets[i] ]);
    fNegEquExpVarPtrs[i]  = &(dFdx[li_Neg][ ANegEquExpVarOffsets[i] ]);
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  // Evaluate Expression with corrected time derivative values
  if (expNumVars_gc != 0)
  {
    Exp_ptr_gc->evaluate( expVal_gc, expVarDerivs);
  }

  // Test derivatives, if too big, zero out
  for (int i = 0; i < expNumVars_gc; ++i)
  {
    double maxMag = 1.0e+10;
    if (expVarDerivs[i] > maxMag || expVarDerivs[i] < -maxMag)
    {
      static Report::MessageCode id;

      Report::UserWarning(id) << "In device " << getName() << ": Expression derivative for variable number " << i << " |" << expVarDerivs[i] << "| exceeds " << maxMag << ", value reduced";

      expVarDerivs[i] = (expVarDerivs[i] > 0) ? maxMag : -maxMag;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 ibis instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double source(0.0), v_pos(0.0), v_neg(0.0), i_bra(0.0);
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  source = expVal_gc;

  fVec[li_Pos] += source;
  fVec[li_Neg] += -source;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    leadF[li_branch_data] = source;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  double coef = 1.0;

  if( expNumVars )
  {
    for( int i = 0; i < expNumVars; ++i )
    {
      dFdx[li_Pos][APosEquExpVarOffsets[i]] += expVarDerivs[i];
      dFdx[li_Neg][ANegEquExpVarOffsets[i]] -= expVarDerivs[i];
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/02
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  
}

//-------------------------------------------------------------------------- 
// Purpose      : Allowed set of characters for use as the IBIS comment 
//                character
//Special Notes : The sequences \\ and \" in the string valid_ibis_comment_chars
//                put the characters backslash and double-quote into that string.
// Creator       : Pete Sholander
// Creation Date : 08/07/18 
//--------------------------------------------------------------------------
static const std::string valid_ibis_comment_chars("!\"#$%&'()*,:;<>?@\\^`{|}~");


//-----------------------------------------------------------------------------
// Function      : Instance::readIbsFile
// Purpose       : This is the initial version of a parser to read in a .ibs
//                 file.  At some point, it should likely be re-done in
//                 something like Flex and Bison.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/05/18
//-----------------------------------------------------------------------------
bool Instance::readIbsFile()
{
  // where we are in the file, and whether parsing had an error
  int lineNum = 0;
  bool psuccess = true;

  // a fixed length of line for some read-through.
  //const size_t maxLineLength=122;
  //char aLine[maxLineLength];

  bool skipPastModelLine = false;
  bool tableEnd = false;

  // try to open the data file
  std::ifstream inputFile;

  // Error out if the user-specified file does not exist, cannot
  // be opened, or is a directory name rather than a file name.  
  // See SON Bug 785 and SRN Bug 2100 for more details.
  if ( !(Util::checkIfValidFile(fileName_)) )
  {
    Report::UserError() << ".ibs file \"" << fileName_ << "\" for device " 
                        << getName() << " could not be found.";
    return false;
  }

  // open the .ibs file
  inputFile.open( fileName_.c_str(),std::ifstream::in);
  if( !inputFile.good() ) 
  {
    Report::UserError() << "Warning .ibs file \"" << fileName_ << "\" for device " 
                        << getName() << " could not be opened.";
    return false;
  }
 
  // parse the file
  std::string aLine;
  IO::TokenVector parsedLine;
  readIBISFileLine(inputFile,aLine,lineNum);
  while( (!inputFile.eof()) || (aLine.substr(0,5) == "[END]") )
  {
    if( (aLine[0] != ibisCommentChar_) )
    {
      // these are "control lines"
      if (aLine[0] == '[')
      {
        splitIBISFileLine(aLine,parsedLine);
        if (aLine.substr(0,10) == "[IBIS ver]")
        {
          if (parsedLine.size() < 3)
	  { 
            Report::UserError() << "Invalid [IBIS ver] line in file \"" << fileName_ 
                                << "\" for device " << getName() << " at line " << lineNum;
            return false;
          }
          else
	  {
            ibisVer_ = parsedLine[2].string_;
          }
        }
        else if (aLine.substr(0,14) == "[Comment Char]" )
        {
          // The new comment character (say semi-colon) will be defined by a line like:
          //    [Comment Char] ;_char
          // where _char is post-pended to the new comment character
          if ( (parsedLine.size() < 3) || (parsedLine[2].string_.substr(1,5) != "_char") )
	  {
            Report::UserError() << "Invalid [Comment Char] line in file \"" << fileName_ 
                                << "\" for device " << getName() << " at line " << lineNum;
            return false;
          }
          
          // Make sure that the new comment characters is an allowed one, per
          // the IBIS standard
          if (valid_ibis_comment_chars.find(parsedLine[2].string_[0]) != std::string::npos )
	  {
	    ibisCommentChar_ = parsedLine[2].string_[0];
          }
          else
	  {
            Report::UserError() << "Invalid comment character on [Comment Char] line in file \"" 
                                << fileName_ << "\" for device " << getName() << " at line " << lineNum;
            return false;
          }         
        }
        else if (aLine.substr(0,11) == "[Component]")
        {
          if ( parsedLine.size() < 2 )
	  {
            Report::UserError() << "Invalid [Component] line in file \"" << fileName_ 
                                << "\" for device " << getName() << " at line " << lineNum;
            return false;
          }
          else
	  {
            ibisComponent_ = parsedLine[1].string_;
          }

          // Parse the typ, min and max values of R_pkg, L_pkg and C_pkg.
          // The typ value are mandatory.  The min and max values may be NA.
	  parsedLine.clear();
          readIBISFileLine(inputFile,aLine,lineNum);
          IBIS::pkgRLC typVec, minVec, maxVec;
          typVec.type_ = IBIS_TYP;
          minVec.type_ = IBIS_MIN;
          maxVec.type_ = IBIS_MAX;
          bool typRFound(false), minRFound(false), maxRFound(false);
          bool typLFound(false), minLFound(false), maxLFound(false);
          bool typCFound(false), minCFound(false), maxCFound(false);
          while( !inputFile.eof() && aLine.substr(0,5) != "[Pin]")
	  {
            if( (aLine[0] != ibisCommentChar_) )
	    {
              if ( aLine.substr(0,5) == "R_pkg")
              {
                splitIBISFileLine(aLine,parsedLine);
                if ( parsedLine.size() < 4 )
	        {
                  Report::UserError() << "Invalid R_pkg line in file \"" << fileName_ 
                                      << "\" for device " << getName() << " at line " << lineNum;
                  return false;
                }
                else
	        {
                  if (parsedLine[1].string_ == "NA")
                  {
                    Report::UserError() << "Invalid R_pkg typ value in file \"" << fileName_ 
                                        << "\" for device " << getName() << " at line " << lineNum;
                    return false;
                  }
                  else
	          {
                    bool ps = ibisStrToVal(parsedLine[1].string_, typVec.R_pkg, typRFound, lineNum); 
                    if (!ps) return false;
                  }

                  if (parsedLine[2].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, minVec.R_pkg, minRFound, lineNum);
                    if (!ps) return false; 
                  }
                  if (parsedLine[3].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[3].string_, maxVec.R_pkg, maxRFound, lineNum);
                    if (!ps) return false; 
                  }
                }
                parsedLine.clear();
              }
              else if ( aLine.substr(0,5) == "L_pkg")
	      {
                splitIBISFileLine(aLine,parsedLine);
                if ( parsedLine.size() < 4 )
	        {
                  Report::UserError() << "Invalid L_pkg line in file \"" << fileName_ 
                                      << "\" for device " << getName() << " at line " << lineNum;
                  return false;
                }
                else
	        {
                  if (parsedLine[1].string_ == "NA")
                  {
                    Report::UserError() << "Invalid L_pkg typ value in file \"" << fileName_ 
                                        << "\" for device " << getName() << " at line " << lineNum;
                    return false;
                  }
                  else
	          {
                    bool ps = ibisStrToVal(parsedLine[1].string_, typVec.L_pkg, typLFound, lineNum);
                    if (!ps) return false; 
                  }

                  if (parsedLine[2].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, minVec.L_pkg, minLFound, lineNum);
                    if (!ps) return false; 
                  }
                  if (parsedLine[3].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[3].string_, maxVec.L_pkg, maxLFound, lineNum);
                    if (!ps) return false; 
                  }
                }
                parsedLine.clear();
              }
              else if ( aLine.substr(0,5) == "C_pkg")
	      {
                splitIBISFileLine(aLine,parsedLine);
                if ( parsedLine.size() < 4 )
	        {
                  Report::UserError() << "Invalid C_pkg line in file \"" << fileName_ 
                                      << "\" for device " << getName() << " at line " << lineNum;
                  return false;
                }
                else
	        {
                  if (parsedLine[1].string_ == "NA")
                  {
                    Report::UserError() << "Invalid C_pkg typ value in file \"" << fileName_ 
                                        << "\" for device " << getName() << " at line " << lineNum;
                    return false;
                  }
                  else
                  {
                    bool ps = ibisStrToVal(parsedLine[1].string_, typVec.C_pkg, typCFound, lineNum);
                    if (!ps) return false; 
                  }

                  if (parsedLine[2].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, minVec.C_pkg, minCFound, lineNum);
                    if (!ps) return false; 
                  }
                  if (parsedLine[3].string_ != "NA")
                  {
                    bool ps = ibisStrToVal(parsedLine[3].string_, maxVec.C_pkg, maxCFound, lineNum);
                    if (!ps) return false; 
                  }
                } 
                parsedLine.clear();

                // typ values must be found
                if (typRFound && typLFound && typCFound)
                {
                  pkgRLCVec_.push_back(typVec);
                }
                else
		{
                  Report::UserError() << "Invalid values for typical RLC_pkg values in file \"" 
                                      << fileName_ << "\" for device " << getName() << " at line " << lineNum;
                  return false;
                }

                // min and max values are optional
                if (minRFound && minLFound && minCFound)
		{
                  pkgRLCVec_.push_back(minVec);
                }
                if (maxRFound && maxLFound && maxCFound)
		{
                  pkgRLCVec_.push_back(maxVec);
                }

                break; // done with processing for package R,L and C
              }
            }
            readIBISFileLine(inputFile,aLine,lineNum);
          }

          parsedLine.clear();
          readIBISFileLine(inputFile,aLine,lineNum);
          // have we seen the [Pin] line yet?
          bool pinLineFound = false; 
          // the pin information should be followed by the first [Model] section
          while( !inputFile.eof() && aLine.substr(0,7) != "[Model]")
	  {
            if( aLine.substr(0,5) == "[Pin]" )
	    {
              pinLineFound = true;
            }
	    else if( (aLine[0] != ibisCommentChar_) && (pinLineFound == true) )
	    {
              splitIBISFileLine(aLine,parsedLine);
              if ( parsedLine.size() < 3 )
	      {
                Report::UserError() << "Invalid [Pin] line in file \"" << fileName_ 
                                    << "\" for device " << getName() << " at line " << lineNum;
                return false;
              }
              else
	      {
                Pin p;
                p.pinNum = atoi(parsedLine[0].string_.c_str());
                p.signal_name = parsedLine[1].string_;
                p.model_name = parsedLine[2].string_;
                if (parsedLine.size() >= 6)
	        {
                  bool psR = ibisStrToVal(parsedLine[3].string_, p.R_pin, p.R_pin_given, lineNum);
                  bool psC = ibisStrToVal(parsedLine[4].string_, p.L_pin, p.L_pin_given, lineNum);
                  bool psL = ibisStrToVal(parsedLine[5].string_, p.C_pin, p.C_pin_given, lineNum);
                  if (!psR || !psC || !psL) return false; 
                }
                pinVec_.push_back(p);
              }
            }
            parsedLine.clear();
            readIBISFileLine(inputFile,aLine,lineNum);
          }
        }
        else if (aLine.substr(0,7) == "[Model]")
	{
          splitIBISFileLine(aLine,parsedLine);
          if ( parsedLine.size() < 2 )
	  {
            Report::UserError() << "Invalid [Model] line in file \"" << fileName_ 
                                << "\" for device " << getName() << " at line " << lineNum;
            return false;
          }
          else if (parsedLine[1].string_ == modelName_)
	  {
	    parsedLine.clear();
            readIBISFileLine(inputFile,aLine,lineNum);
             
            while( !inputFile.eof() && aLine.substr(0,7) != "[Model]")
	    {
              if( (aLine[0] != ibisCommentChar_) )
	      {
                if( aLine.substr(0,10) == "Model_type" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 2 )
	          {
                    Report::UserError() << "Invalid [Model] line in file \"" << fileName_ 
                                        << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
                  {
                    bufferModel_.Model_type = setIBISModelType(parsedLine[1].string_);  
                  }
                }
                else if ( aLine.substr(0,8) == "Polarity" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 2 )
	          {
                    Report::UserError() << "Invalid Polarity line in file \"" << fileName_ 
                                        << "\" for device " << getName() << " at line " << lineNum;
                    return false;
	          }
                  else
		  {
                    bufferModel_.Polarity = setIBISModelPolarity(parsedLine[1].string_);  
                  }
                }
                else if ( aLine.substr(0,4) == "Vinl" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 3 )
	          {
                    Report::UserError() << "Invalid Vinl line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
		    bool ps = ibisStrToVal(parsedLine[2].string_, bufferModel_.Vinl, bufferModel_.Vinl_given, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,4) == "Vinh" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 3 )
	          {
                    Report::UserError() << "Invalid Vinl line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, bufferModel_.Vinh, bufferModel_.Vinh_given, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,5) == "Vmeas" )
		{
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 3 )
	          {
                    Report::UserError() << "Invalid Vmeas line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, bufferModel_.Vmeas, bufferModel_.Vmeas_given, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,4) == "Rref" )
		{
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 3 )
	          {
                    Report::UserError() << "Invalid Rref line in file \"" << fileName_ 
				  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, bufferModel_.Rref, bufferModel_.Rref_given, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,4) == "Vref" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 3 )
	          {
                    Report::UserError() << "Invalid Vref line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = ibisStrToVal(parsedLine[2].string_, bufferModel_.Vref, bufferModel_.Vref_given, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,6) == "C_comp" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 4 )
	          {
                    Report::UserError() << "Invalid C_comp line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = makeTmmVec(parsedLine, bufferModel_.C_comp, 1, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,15) == "[Voltage Range]" )
		{
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 5 )
	          {
                    Report::UserError() << "Invalid [Voltage Range] line in file \"" << fileName_ 
			  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = makeTmmVec(parsedLine, bufferModel_.Voltage_Range, 2, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,19) == "[Temperature Range]" )
	        {
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 5 )
	          {
                    Report::UserError() << "Invalid [Temperature Range] line in file \"" << fileName_ 
				  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = makeTmmVec(parsedLine, bufferModel_.Temperature_Range, 2, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,21) == "[GND Clamp Reference]" )
		{
                  splitIBISFileLine(aLine,parsedLine);
                  if ( parsedLine.size() < 6 )
	          {
                    Report::UserError() << "Invalid [Voltage Range] line in file \"" << fileName_ 
					  << "\" for device " << getName() << " at line " << lineNum;
                    return false;
		  }
                  else
		  {
                    bool ps = makeTmmVec(parsedLine, bufferModel_.GND_Clamp_Reference, 3, lineNum);
                    if (!ps) return false;
                  }
                }
                else if ( aLine.substr(0,11) == "[GND Clamp]" )
		{
                  std::string expString = "TABLE {V(" + nodeList_[2] + "," + 
                                                        nodeList_[3] + ")} = ";

                  // skip to next line in file and then read table entries which end with the
                  // next command line that starts with [
                  parsedLine.clear();
                  readIBISFileLine(inputFile,aLine,lineNum);             
                  while ( aLine[0] != '[' ) 
		  {
                    while ( aLine[0] == ibisCommentChar_ )
		    {
                      readIBISFileLine(inputFile,aLine,lineNum);
                    }
		    if( aLine[0] != '[' )
		    {
                      splitIBISFileLine(aLine,parsedLine);
                      double voltage, current;
                      ExtendedString posVal = parsedLine[0].string_;

                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        voltage = posVal.IBISValue();
                      }
                      posVal = parsedLine[1].string_;     
                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        current = posVal.IBISValue();
                      }

		      expString += "( " + parsedLine[0].string_ + " , " + parsedLine[1].string_ + " )";
                      parsedLine.clear();
                      readIBISFileLine(inputFile,aLine,lineNum);
                    }
                  }
                  Param tbl_param( "", "" );
                  tbl_param.setTag( "GNDCLAMPTBL" );
                  tbl_param.setVal( Util::Expression(getSolverState().getGroupWrapper()->expressionGroup_,expString) );
                  tbl_param.setGiven( true );
		  IB.params.push_back(tbl_param);

                  tableEnd = true;
                }
                else if ( aLine.substr(0,13) == "[POWER Clamp]" )
		{
                  std::string expString = "TABLE {V(" + nodeList_[2] + "," + 
                                                        nodeList_[3] + ")} = ";
                  
                  // skip to next line in file and then read table entries which end with the
                  // next command line
                  parsedLine.clear();
                  readIBISFileLine(inputFile,aLine,lineNum);                 
                  while ( aLine[0] != '[' ) 
		  {
                    while ( aLine[0] == ibisCommentChar_)
		    {
                      readIBISFileLine(inputFile,aLine,lineNum);
		    }
                    if( aLine[0] != '[' )
		    {
                      splitIBISFileLine(aLine,parsedLine);
                      double voltage, current;
                      ExtendedString posVal = parsedLine[0].string_;

                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        voltage = posVal.IBISValue();
                      }
                      posVal = parsedLine[1].string_;     
                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        current = posVal.IBISValue();
                      }

                      expString += "( " + parsedLine[0].string_ + " , " + parsedLine[1].string_ + " )";
                      parsedLine.clear();
                      readIBISFileLine(inputFile,aLine,lineNum);
                    }
                  }
                  Param tbl_param( "", "" );
                  tbl_param.setTag( "PWRCLAMPTBL" );
                  tbl_param.setVal( Util::Expression(getSolverState().getGroupWrapper()->expressionGroup_,expString) );
                  tbl_param.setGiven( true );
		  IB.params.push_back(tbl_param);

                  tableEnd = true;
                }
	        else if ( aLine.substr(0,10) == "[Pulldown]" )
	        {
                  std::string expString = "TABLE {V(" + nodeList_[2] + "," + 
                                                        nodeList_[3] + ")} = ";

                  // skip to next line in file and then read table entries which end with the
                  // next command line
                  parsedLine.clear();
                  readIBISFileLine(inputFile,aLine,lineNum);                 
                  while ( aLine[0] != '[' ) 
		  {
                    while ( aLine[0] == ibisCommentChar_)
		    {
                      readIBISFileLine(inputFile,aLine,lineNum);
		    }
                    if( aLine[0] != '[' )
		    {
                      splitIBISFileLine(aLine,parsedLine);
                      double voltage, current;
                      ExtendedString posVal = parsedLine[0].string_;

                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        voltage = posVal.IBISValue();
                      }
                      posVal = parsedLine[1].string_;     
                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        current = posVal.IBISValue();
                      }

                      expString += "( " + parsedLine[0].string_ + " , " + parsedLine[1].string_ + " )";
                      parsedLine.clear();
                      readIBISFileLine(inputFile,aLine,lineNum);
                    }
                  }
                  Param tbl_param( "", "" );
                  tbl_param.setTag( "PULLDOWNTBL" );
                  tbl_param.setVal( Util::Expression(getSolverState().getGroupWrapper()->expressionGroup_,expString) );
                  tbl_param.setGiven( true );
                  IB.params.push_back(tbl_param);

                  tableEnd = true;
	        }
                else if ( aLine.substr(0,8) == "[Pullup]" )
	        {
                  std::string expString = "TABLE {V(" + nodeList_[2] + "," + 
                                                        nodeList_[3] + ")} = ";

                  // skip to next line in file and then read table entries which end with the
                  // next command line
                  parsedLine.clear();
                  readIBISFileLine(inputFile,aLine,lineNum);                 
                  while ( aLine[0] != '[' ) 
		  {
                    while ( aLine[0] == ibisCommentChar_)
		    {
                      readIBISFileLine(inputFile,aLine,lineNum);
		    }
                    if( aLine[0] != '[' )
                    {
                      splitIBISFileLine(aLine,parsedLine);
                      double voltage, current;
                      ExtendedString posVal = parsedLine[0].string_;

                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        voltage = posVal.IBISValue();
                      }
                      posVal = parsedLine[1].string_;     
                      if (!posVal.isValue())
                      {
                        Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                                            << "\" for device " << getName() << " at line " << lineNum;
                        return false;
                      }
                      else
                      { 
                        current = posVal.IBISValue();
                      }

                      expString += "( " + parsedLine[0].string_ + " , " + parsedLine[1].string_ + " )";
                      parsedLine.clear();
                      readIBISFileLine(inputFile,aLine,lineNum);
                    }
                  }
                  Param tbl_param( "", "" );
                  tbl_param.setTag( "PULLUPTBL" );
                  tbl_param.setVal( Util::Expression(getSolverState().getGroupWrapper()->expressionGroup_,expString) );
                  tbl_param.setGiven( true );
		  IB.params.push_back(tbl_param);

                  tableEnd = true;
	        }    
              } 
              
              parsedLine.clear();
              if (!tableEnd)
	      {
                readIBISFileLine(inputFile,aLine,lineNum);
              }
              else
	      {
                tableEnd=false;
              }
            }
	  }
	  else
	  {
	    skipPastModelLine = true;
	  }
        }  
      }
    }
    parsedLine.clear();
    if ( (aLine.substr(0,7) != "[Model]") || (skipPastModelLine == true) ) 
    {
      readIBISFileLine(inputFile,aLine,lineNum);
    }
  }

  inputFile.close();

  // update any instance parameters that were set in the .ibs file
  setParams (IB.params);

  // do some error checking
  if (pkgRLCVec_.size() == 0)
  {
    Report::UserError() << "No valid set of R_pkg, C_pkg and L_pkg values found in file \"" 
                        << fileName_ << "\" for device " << getName();
    return false;
  }

  if (pinVec_.size() == 0)
  {
    Report::UserError() << "No valid Pin lines found in file \"" << fileName_ 
                        << "\" for device " << getName();
    return false;
  }

  if ( bufferModel_.Model_type == IBIS_MODEL_INVALID )
  {
    Report::UserError() << "No valid model of type " <<  modelName_ <<  " found in \"" 
                        << fileName_  << "\" for device " << getName();
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::ibisStrToVal
// Purpose       : Turn a string in a .ibs file into a value of type double.
//                 It also sets the flag valFound to true if that conversion
//                 works.  (Note: the valFound flag is used to set the isGiven
//                 flag for a given private variable.)  If the conversion 
//                 doesn't work then the function returns false.
// Special Notes : This function is needed because the IBIS spec has a different 
//                 set of string identifiers for the units, and those identifiers 
//                 are case sensitive.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
bool Instance::ibisStrToVal(const std::string& valStr, double& val, bool& valFound, 
                            const int lineNum)
{
  bool success = false;  // assume failure
  valFound = false;

  ExtendedString posVal = valStr;
  // At some point, isValue() should likely be specialized to only recognize 
  // the string identifiers that are legal in IBIS rather than legal in Xyce.
  if (!posVal.isValue())
  {
    Report::UserError() << "string can not be converted to a value in file \"" << fileName_ 
                        << "\" for device " << getName() << " at line " << lineNum;
    success = false;
  }
  else
  {
    val = posVal.IBISValue();
    valFound = true;
    success = true;
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : Instance::makeTmmVec
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/27/18
//-----------------------------------------------------------------------------
bool Instance::makeTmmVec(const IO::TokenVector& pl, std::vector<tmmParam>& tpv, 
                          const int offset, const int lineNum)
{
  bool success = false;  // assume failure

  if (pl[offset].string_ == "NA")
  {
    Report::UserError() << "Invalid typ value in file \"" << fileName_ 
                        << "\" for device " << getName() << " at line " << lineNum;
    success = false;
  }
  else
  {
    tmmParam tp;
    bool ps = ibisStrToVal(pl[offset].string_, tp.val_, tp.given_, lineNum);
    if (ps)
    {
      tp.type_= IBIS_TYP;
      tpv.push_back(tp);
      success = true;
    }
    else 
    {
      return false; 
    }
  }

  if (pl[offset+1].string_ != "NA")
  {
    tmmParam tp;
    bool ps = ibisStrToVal(pl[offset+1].string_, tp.val_, tp.given_, lineNum);
    if (ps)
    {
      tp.type_= IBIS_MIN;
      tpv.push_back(tp);
      success = true;
    }
    else 
    {
      return false; 
    }
  }

  if (pl[offset+2].string_ != "NA")
  {
    tmmParam tp;
    bool ps = ibisStrToVal(pl[offset+2].string_, tp.val_, tp.given_, lineNum);
    if (ps)
    {
      tp.type_= IBIS_MAX;
      tpv.push_back(tp);
      success = true;
    }
    else 
    {
      return false; 
    }
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIBISModelType
// Purpose       :
// Special Notes : It might make sense to replace these if-else blocks with a 
//                 map at some point.  For each Xyce release, I might modify
//                 this if-else statement to only those IBISModelType values
//                 that are supported in that release.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
IBIS::IBISModelType Instance::setIBISModelType(std::string modelTypeStr) 
{
  IBISModelType mType = IBIS_MODEL_INVALID;

  if (modelTypeStr == "Input")
  {
    mType = IBIS_INPUT;
  }
  else if (modelTypeStr == "Output")
  {
    mType = IBIS_OUTPUT;
  }
  else if (modelTypeStr == "I/O")
  {
    mType = IBIS_IO;
  }
  else if (modelTypeStr == "3-state")
  {
    mType = IBIS_3STATE;
  }
  else if (modelTypeStr == "Open_drain")
  {
    mType = IBIS_OPEN_DRAIN;
  }
  else if (modelTypeStr == "I/0_open_drain")
  {
    mType = IBIS_IO_OPEN_DRAIN;
  }
  else if (modelTypeStr == "Open_sink")
  {
    mType = IBIS_OPEN_SINK;
  }
  else if (modelTypeStr == "I/O_open_sink")
  {
    mType = IBIS_IO_OPEN_SINK;
  }
   else if (modelTypeStr == "open_source")
  {
    mType = IBIS_OPEN_SOURCE;
  }
  else if (modelTypeStr == "I/O_open_source")
  {
    mType = IBIS_IO_OPEN_SOURCE;
  }
  else if (modelTypeStr == "Input_ECL")
  {
    mType = IBIS_INPUT_ECL;
  }
  else if (modelTypeStr == "Output_ECL")
  {
    mType = IBIS_OUTPUT_ECL;
  }
  else if (modelTypeStr == "IO_ECL")
  {
    mType = IBIS_IO_ECL;
  }
  else if (modelTypeStr == "3-state_ECL")
  {
    mType = IBIS_3STATE_ECL;
  }
  else if (modelTypeStr == "Terminator")
  {
    mType = IBIS_TERMINATOR;
  }
   else if (modelTypeStr == "Series")
  {
    mType = IBIS_IO_ECL;
  }
  else if (modelTypeStr == "Series_switch")
  {
    mType = IBIS_SERIES;
  }
  else if (modelTypeStr == "Terminator")
  {
    mType = IBIS_SERIES_SWITCH;
  }
  else if (modelTypeStr == "Input_diff")
  {
    mType = IBIS_INPUT_DIFF;
  }
  else if (modelTypeStr == "Output_diff")
  {
    mType = IBIS_OUTPUT_DIFF;
  }
  else if (modelTypeStr == "I/O_diff")
  {
    mType = IBIS_IO_DIFF;
  }
  else if (modelTypeStr == "3-state_diff")
  {
    mType = IBIS_3STATE_DIFF;
  }
  else
  {
    Report::UserError() << "Invalid model type " << modelTypeStr << " in \"" << fileName_ 
                         << "\" for device " << getName();
    mType = IBIS_MODEL_INVALID;
  }

  return mType;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIBISModelType
// Purpose       :
// Special Notes : It might make sense to replace these if-else blocks with a 
//                 map at some point.  For each Xyce release, I should truncate
//                 this if-else statement to only those IBISModelPolarity
//                 values that are supported in that release.  There's also
//                 issue of checking "compatibility" between the IBISModelType
//                 and IBISModelPolarity values.  I'm not sure that all possible
//                 combinations are allowed by the IBIS standard.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
IBIS::IBISModelPolarity Instance::setIBISModelPolarity(std::string polarityTypeStr) 
{
  IBISModelPolarity pType = IBIS_POLARITY_INVALID;

  if (polarityTypeStr == "Inverting")
  {
    pType = IBIS_POLARITY_INVERTING;
  }
  else if (polarityTypeStr == "Non-Inverting")
  {
    pType = IBIS_POLARITY_NONINVERTING;
  }
  else
  {
    Report::UserError() << "Invalid polarity " << polarityTypeStr << " in \"" << fileName_ 
                         << "\" for device " << getName();
    pType = IBIS_POLARITY_INVALID;
  }

  return pType;
}

//-----------------------------------------------------------------------------
// Function      : Instance::splitIBISFileLine
// Purpose       : Handle in-line comments, and then split a line from a .ibs
//                 file into a TokenVector
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/08/18
//-----------------------------------------------------------------------------
void Instance::splitIBISFileLine(const std::string& aLine, IO::TokenVector & parsedLine)
{
  // trim off any in-line comments
  std::string trimmedLine(aLine);
  size_t pos = trimmedLine.find(ibisCommentChar_);
  if (pos != std::string::npos) trimmedLine.erase(trimmedLine.begin()+pos,trimmedLine.end());

  // now split the line into "tokens", and also trim out extra white space from those
  // tokens.
  IO::splitLineNoWS(trimmedLine,parsedLine);  

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::readIBISFileLine
// Purpose       : Read in a line from a .ibs file, and also increment the 
//                 lineNum counter
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/08/18
//-----------------------------------------------------------------------------
void Instance::readIBISFileLine(std::istream & inputFile, std::string& aLine, int& lineNum)
{
  IO::readLine(inputFile,aLine);    
  lineNum++;

  return;
}
 
// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool Model::processParams()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
bool Model::processInstanceParams()
{
  return true;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

// IBIS Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
    return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi =  *(*it);

    // Evaluate Expression with corrected time derivative values
    if (bi.expNumVars_gc != 0)
    {
      bi.Exp_ptr_gc->evaluate( bi.expVal_gc, bi.expVarDerivs);
    }

    // Test derivatives, if too big, zero out
    for (int k = 0; k < bi.expNumVars_gc; ++k)
    {
      double maxMag = 1.0e+10;
      if (bi.expVarDerivs[k] > maxMag || bi.expVarDerivs[k] < -maxMag)
      {
        static Report::MessageCode id;

        Report::UserWarning(id) << "In device " << bi.getName() << ": Expression derivative for variable number " << k << " |" << bi.expVarDerivs[k] << "| exceeds " << maxMag << ", value reduced";
        bi.expVarDerivs[k] = (bi.expVarDerivs[k] > 0) ? maxMag : -maxMag;
      }
    }
  }

  return true;
  }

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);
    double v_pos(0.0), v_neg(0.0), i_bra(0.0);
    double source = bi.expVal_gc;

    fVec[bi.li_Pos] += source;
    fVec[bi.li_Neg] += -source;
    if( bi.loadLeadCurrent )
    {
      leadF[bi.li_branch_data] = source;
      junctionV[bi.li_branch_data] = solVec[bi.li_Pos] - solVec[bi.li_Neg];
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & bi = *(*it);
    double coef = 1.0;

    if( bi.expNumVars )
    {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      for( int k = 0; k < bi.expNumVars; ++k )
      {
        *bi.fPosEquExpVarPtrs[k] += bi.expVarDerivs[k];
        *bi.fNegEquExpVarPtrs[k] -= bi.expVarDerivs[k];
      }
#else
      for( int j = 0; j < bi.expNumVars; ++j )
      {
        dFdx[bi.li_Pos][bi.APosEquExpVarOffsets[j]] += bi.expVarDerivs[j];
        dFdx[bi.li_Neg][bi.ANegEquExpVarOffsets[j]] -= bi.expVarDerivs[j];
      }
#endif
    }
  }
  return true;
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("IBIS", 1);
}

} // namespace IBIS
} // namespace Device
} // namespace Xyce
