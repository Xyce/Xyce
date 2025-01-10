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
// Purpose        : This class defines a user-defined membrane mechanism.
//
// Special Notes  :
//
// Creator        : Christy Warrender, Cognitive Modeling
//
// Creation Date  : 12/14/2010
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include<iostream>
#include<algorithm>

// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_DEV_MembraneUserDefined.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExpressionGroupWrapper.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::MembraneUserDefined
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
MembraneUserDefined::MembraneUserDefined(const SolverState & ss1, double cMem, double gMem, double vRest,
    	std::vector<std::string> & currentEqus, std::vector<std::string> & indepVars, std::vector<std::string> & fEqs,
    	std::vector<std::string> & qEqs, std::vector<std::string> & extraFunctions, std::vector<std::string> & extraParameters) :
    MembraneModel(ss1),
    cMem_(cMem),
    gMem_(gMem),
    vRest_(vRest),
    currentEqus_(currentEqus),
    indepVars_(indepVars),
    fEqs_(fEqs),
    qEqs_(qEqs),
    extraFunctions_(extraFunctions),
    extraParameters_(extraParameters)
{
  // make expressions from the passed in strings
  convertStringsToExpression( currentEqus_, currentEqusExpRCP_ );
  convertStringsToExpression( indepVars_, indepVarsExpRCP_ );
  convertStringsToExpression( fEqs_, fEqsExpRCP_ );
  convertStringsToExpression( qEqs_, qEqsExpRCP_ );
  convertStringsToExpression( extraFunctions_, extraFunctionsExpRCP_ );
  convertStringsToExpression( extraParameters_, extraParametersExpRCP_ );

  // try to reduce equations to just the unknowns (indepVars and V)
  consolidateExpressions();

  numIndependentVars_ = indepVarsExpRCP_.size() + 1;  // the +1 is for the voltage at each node which
                                                      // is always there even if this is a passive cable
  systemJacOffset_.resize(numIndependentVars_);       // this will hold a map to for the local jacobian
                                                      // offsets on each row.
}

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::setJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
void MembraneUserDefined::setJacStamp( int numExtVars, int segmentNumber, int vOffset,
                                             std::vector< std::vector< int > > & segmentJacStamp )
{
  // caller sets up size of row of jac stamp to numIndependentVars_ + extra's needed for
  // its modeling.  So for a cable based device this is Vpre, Vseg, (other membrane vars), Vnext.
  // in general this is numIndependentVars_ + 2 (Vpre and Vnext).  In this routine we fill in
  // just what is needed for the membrane model

  int offset = numExtVars + numIndependentVars_*segmentNumber;

  // note, any external vars come before the internal vars.  Dependance on Vprev and Vnext (if
  // they are really there) are handled by the caller as in a cable equation model.

  // general Jacobian strcuture:
  //             Vpre    V     user var 1   user var 2   user var 3   Vnext
  // kcl          yes   yes   	   ?             ?            ?   	   yes
  // user var 1          ?         yes           ?            ?
  // user var 2          ?         ?             yes          ?
  // user var 3          ?         ?             ?            yes
  //
  // need to parse user-defined equations to determine dependence of each internal variable on the
  // others
  std::map<std::string,int>::iterator currVarOffset = indepVarOffset_.begin();
  std::map<std::string,int>::iterator endVarOffset = indepVarOffset_.end();

  while( currVarOffset != endVarOffset )
  {
    std::set<std::string> depVarNames;
    depVarNames.clear();
    if( currVarOffset->first == "V" )
    {
      // voltage is different as it depends on a set of current equations
      for( int i=0; i<currentEqusVarNames_.size(); i++ )
      {
        depVarNames.insert( currentEqusVarNames_[i].begin(), currentEqusVarNames_[i].end() );
      }
    }
    else
    {
      // for other variables we look at the F & Q equation
      depVarNames.insert( fEqsEqusVarNames_[ currVarOffset->second - 1 ].begin(), fEqsEqusVarNames_[ currVarOffset->second - 1 ].end() );
      depVarNames.insert( qEqsEqusVarNames_[ currVarOffset->second - 1 ].begin(), qEqsEqusVarNames_[ currVarOffset->second - 1 ].end() );
    }

    // now we have a set of the variable names, but the order is unknown.  Use the offsetToIndepVar_ map to make
    // a new map in the right order - meaning the order of nonzero terms in the jacobian - and in sequence.
    // cew - taking next bit out.  systemJacOffsets_ is shared by all segments, so we can't use it to
    // keep track of the differences in which variables V depends on - instead, add vOffset to the values
    // stored here when they're used, both in setting up jacStamp below, and in loadDAEdFdx (it's
    // not used in loadDAEdQdx - we assume no contributions to V from user-defined equations).
    // Also, if the currVarOffset->first (the name of the variable we're processing now) =="V" then
    // we have to add vOffset as this allows for "Vprevious" on the same row of the jacobian
    // In other words, we're assuming that the internal segment variables don't depend on anything
    // outside the membrane model, but we know that V does, and that may affect internal variables'
    // indices in the segmentJacStamp
    //int voltageOffset=0;
    //if( currVarOffset->first == "V" )
    //{
      //voltageOffset=vOffset;
    //}
    std::set<std::string>::iterator depVarNamesEnd = depVarNames.end();
    int numEnteries = offsetToIndepVar_.size();
    std::map< std::string, int > tmpIndVarMap;
    int numFound = 0;
    for( int i=0; i<numEnteries; i++)
    {
      if( depVarNames.find( offsetToIndepVar_[i] ) != depVarNamesEnd )
      {
        //tmpIndVarMap[ offsetToIndepVar_[i] ] = voltageOffset + numFound++;
        tmpIndVarMap[ offsetToIndepVar_[i] ] = numFound++;
        Xyce::dout() << "Assigned tmpIndVarMap[ " << offsetToIndepVar_[i] << " ]= " << (numFound-1) << std::endl;
      }
    }

    systemJacOffset_[currVarOffset->second ]= tmpIndVarMap;
    Xyce::dout() << "tmpIndVarMap assigned to variable " << currVarOffset->first << std::endl;

    // get the number of things on this row
    int numvars = systemJacOffset_[currVarOffset->second].size();
    // Xyce::dout() << "systemJacOffset_[ " << currVarOffset->second << "].size() = " << numvars << std::endl;

    // the owning device allocates this row for us in the case of V
    // so don't resize it in that case
    if( currVarOffset->first != "V" )
    {
      segmentJacStamp[offset + currVarOffset->second ].resize( numvars );
    }

    // now loop over the set of independent vars to set jacStamp dependance
    std::set<std::string>::iterator depVarNamesCurr = depVarNames.begin();
    while( depVarNamesCurr != depVarNamesEnd )
    {
      if( currVarOffset->first == "V" )
      {
        segmentJacStamp[offset + currVarOffset->second ][ vOffset+systemJacOffset_[currVarOffset->second ][*depVarNamesCurr ]] = offset + indepVarOffset_[*depVarNamesCurr];
      }
      else
      {
        segmentJacStamp[offset + currVarOffset->second ][ systemJacOffset_[currVarOffset->second ][*depVarNamesCurr ]] = offset + indepVarOffset_[*depVarNamesCurr];
      }
      depVarNamesCurr++;
    }
    currVarOffset++;
  }

  // print out jacstamp for debugging

  Xyce::dout() << "MembraneUserDefined::setJacStamp() jacStamp = " << std::endl;
  for( int i = 0; i<segmentJacStamp.size(); ++i )
  {
    Xyce::dout() << "jacStamp[ " << i << " ] = { ";
    for( int j=0; j<segmentJacStamp[i].size(); ++j)
    {
       Xyce::dout() << segmentJacStamp[i][j];
      if( j != ( segmentJacStamp[i].size() -1 ) )
      {
        Xyce::dout() << ", ";
      }
    }
    Xyce::dout() << " }" << std::endl;
  }
  Xyce::dout() << Xyce::section_divider << std::endl;

  Xyce::dout() << "MembraneUserDefined::setJacStamp() systemJacOffset_ = " << std::endl;
  for( int i = 0; i<systemJacOffset_.size(); ++i )
  {
    Xyce::dout() << "systemJacOffset[ " << i << " ] = { ";
    std::map<std::string,int>::iterator currVarOffset = systemJacOffset_[i].begin();
    std::map<std::string,int>::iterator endVarOffset = systemJacOffset_[i].end();
    while( currVarOffset != endVarOffset )
    {
        Xyce::dout() << currVarOffset->first << ", " << currVarOffset->second << std::endl;
        currVarOffset++;
    }
    Xyce::dout() << " }" << std::endl;
  }
  Xyce::dout() << Xyce::section_divider << std::endl;



}

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
void MembraneUserDefined::loadDAEQVector( int segmentNumber,
                                                     std::vector< int > & lidIndexVector,
                                                     Linear::Vector * solnVecPtr,
                                                     Linear::Vector * daeQVecPtr,
                                                     double segArea)
{
  int index = segmentNumber * numIndependentVars_;

  // contribution of membrane capacitance to segment voltage
  (*daeQVecPtr)[lidIndexVector[index]] += cMem_ * segArea * (*solnVecPtr)[lidIndexVector[index]];
  int numCurrentExp = currentEqusExpRCP_.size();
  Xyce::dout() << "loadDAEQVector:  entry for index " << index << " : " << cMem_ * segArea * (*solnVecPtr)[lidIndexVector[index]] << std::endl;

  // add Q terms for user-defined vars
  int numExp = qEqsExpRCP_.size();
  for( int i=0; i<numExp; i++)
  {
    // get the number of vars in this expression and update their values
    int numvars = qEqsEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      qEqsEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ qEqsEqusVarNames_[i][j] ] ] ];
    }
    // evaluate the expression
    double resultValue=0.0;
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    qEqsExpRCP_[i]->evaluateFunction( resultValue, qEqsEqusVarValues_[i] );
#else
    qEqsExpRCP_[i]->evaluateFunction( resultValue);
#endif

    // add it to the Q vector
    // cew not sure about my change here; do we know that vars are in order?
    //  instead of adding i to index, should use indepVarOffset_ map
    //  but how do I get ver name at this level?
    int lidIndexVectorIndex = index + i + 1;	// +1 for V
    (*daeQVecPtr)[lidIndexVector[lidIndexVectorIndex]] += resultValue;
    Xyce::dout() << "loadDAEQVector:  entry for LID index " << lidIndexVectorIndex
    	<< ", varname " << offsetToIndepVar_[lidIndexVectorIndex-index] << " : " << resultValue << std::endl;
  }
}



//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
void MembraneUserDefined::loadDAEFVector( int segmentNumber,
                                                     std::vector< int > & lidIndexVector,
                                                     Linear::Vector * solnVecPtr,
                                                     Linear::Vector * daeFVecPtr,
                                                     double segArea)
{
  Xyce::dout() << "loadDAEFVector" << std::endl;

  int index = segmentNumber * numIndependentVars_;

  // leak contribution
  (*daeFVecPtr)[lidIndexVector[index]] += gMem_ * segArea * ((*solnVecPtr)[lidIndexVector[index]] - vRest_ );

  // to get at the independent vars of this segment we'll need to use the passed in lidIndexVector
  // as in (*solnVecPtr)[ lidIndexVector[index + some_offset ] ]

  // add contribution from current equation
  int numCurrentExp = currentEqusExpRCP_.size();
  for( int i=0; i<numCurrentExp; i++)
  {
    Xyce::dout() << "membrane V equation contribution from current equation # " << i << std::endl;
    // get the number of vars in this expression and update their values
    int numvars = currentEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      currentEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ currentEqusVarNames_[i][j] ] ] ];
      Xyce::dout() << "Segment " << segmentNumber << " current load  variable = " << currentEqusVarNames_[i][j] << " value = " << currentEqusVarValues_[i][j] << std::endl;
    }
    // evaluate the expression
    double resultValue=0.0;
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    currentEqusExpRCP_[i]->evaluateFunction( resultValue, currentEqusVarValues_[i] );
#else
    currentEqusExpRCP_[i]->evaluateFunction( resultValue );
#endif
    Xyce::dout() << "Segment " << segmentNumber << " current equ F = " << resultValue << std::endl;

    // add it to the F vector
    // cew - need to multiply by segArea, assuming current specified as current density
    (*daeFVecPtr)[lidIndexVector[index]] += segArea*resultValue;
  }

  // add contribution from extra vars
  int numExp = fEqsExpRCP_.size();
  for( int i=0; i<numExp; i++)
  {
    Xyce::dout() << "F terms from f equation # " << i << std::endl;
    // get the number of vars in this expression and update their values
    int numvars = fEqsEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      fEqsEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ fEqsEqusVarNames_[i][j] ] ] ];
      Xyce::dout() << "Segment " << segmentNumber << " extra var load variable = " << fEqsEqusVarNames_[i][j] << " value = " << fEqsEqusVarValues_[i][j] << std::endl;
    }
    // evaluate the expression
    double resultValue=0.0;
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    fEqsExpRCP_[i]->evaluateFunction( resultValue, fEqsEqusVarValues_[i] );
#else
    fEqsExpRCP_[i]->evaluateFunction( resultValue );
#endif

    // add it to the F vector
    //int lidIndexVectorIndex = index + numCurrentExp + i;
    int lidIndexVectorIndex = index + i + 1;	// +1 for V
    (*daeFVecPtr)[lidIndexVector[lidIndexVectorIndex]] += resultValue;
    Xyce::dout() << "Segment " << segmentNumber << " LID index " << lidIndexVectorIndex
    	<< ", extra vars equ F = " << resultValue << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
void MembraneUserDefined::loadDAEdQdx( int segmentNumber, int vOffset,
                                                std::vector< int > & lidIndexVector,
                                                std::vector< std::vector< int > > & jacobianOffsets,
                                                Linear::Vector * solnVecPtr,
                                                Linear::Matrix * dQdxMatPtr,
                                                double segArea)
{
  Xyce::dout() << "loadDAEdQdx" << std::endl;

  // while lidIndexVector lists LID's for just the segment variables (just V in the case
  // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
  // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

  // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
  // in the case of the passive cable numIndependentVars_=1.

  int index = segmentNumber * numIndependentVars_;
  int row = numExternalVars_ + index;               // numExternalVars_ a contant of 2 assumed in MembraneModel base class

  // Vseg equation
  (*dQdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset]] += cMem_ * segArea;
  int numCurrentExp = currentEqusExpRCP_.size();

  // cew - in theory, currents might contribute to dQdx, although they don't in our test case

  // add contribution for extra vars
  int numExp = qEqsExpRCP_.size();
  for( int i=0; i<numExp; i++)
  {
    std::string targetVarName = offsetToIndepVar_[i+1];	// +1 because offsetToIndepVar includes V, but numExp doesn't
    Xyce::dout() << "Q terms for q equation " << i << "; assumed to be for variable " << targetVarName << std::endl;
    // get the number of vars in this expression and update their values
    int numvars = qEqsEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      qEqsEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ qEqsEqusVarNames_[i][j] ] ]];
      Xyce::dout() << "contribution from " << qEqsEqusVarNames_[i][j] << ": " << qEqsEqusVarValues_[i][j] << std::endl;
    }
    // evaluate the expression
    double resultValue=0.0;
    std::vector<double> derivValues;
    derivValues.resize( numvars );
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    qEqsExpRCP_[i]->evaluate( resultValue, derivValues, qEqsEqusVarValues_[i] );
#else
    qEqsExpRCP_[i]->evaluate( resultValue, derivValues );
#endif
    Xyce::dout() << "expression result:  " << resultValue << std::endl;

    // add it to the dQdx matrix
    for( int j=0; j<numvars; j++)
    {
      std::string varName = qEqsEqusVarNames_[i][j];
      Xyce::dout() << "contribution of variable " << j << "(" << varName << ")" << std::endl;
      int targetVarOffset = indepVarOffset_[targetVarName];
      int targetVarIndex = index + targetVarOffset;
      int targetVarRow = row + targetVarOffset;
      int jacOffsetsCol = systemJacOffset_[targetVarOffset][varName];
      Xyce::dout() << "targerVarOffset: " << targetVarOffset << ", targetVarIndex:  " << targetVarIndex
          << ", jacobianRow: " << targetVarRow << ", jacobianCol: " << jacOffsetsCol << std::endl;
      (*dQdxMatPtr)[lidIndexVector[targetVarIndex]][jacobianOffsets[targetVarRow][jacOffsetsCol]] += derivValues[j];
      Xyce::dout() << "adding (*dQdxMatPtr)[ " << lidIndexVector[targetVarIndex] << " ][ " << jacobianOffsets[targetVarRow][jacOffsetsCol] << " ] = " << derivValues[j] << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::loadDAEdFdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, Cognitive Modeling & Richard Schiek Elec. Sys. Modeling
// Creation Date : 12/14/2010
//-----------------------------------------------------------------------------
void MembraneUserDefined::loadDAEdFdx( int segmentNumber, int vOffset,
                                             std::vector< int > & lidIndexVector,
                                             std::vector< std::vector< int > > & jacobianOffsets,
                                             Linear::Vector * solnVecPtr,
                                             Linear::Matrix * dFdxMatPtr,
                                             double segArea)
{
  Xyce::dout() << "loadDAEdFdx" << std::endl;

  // while lidIndexVector lists LID's for just the segment variables (just V in the case
  // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
  // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

  // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
  // in the case of the passive cable numIndependentVars_=1.

  int index = segmentNumber * numIndependentVars_;
  int row = numExternalVars_ + index;       // numExternalVars_ a contant of 2 assumed in MembraneModel base class

  // leak contribution
  (*dFdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset]] += gMem_ * segArea;
  Xyce::dout() << "leak current dFdx:  " << "vOffset = " << vOffset << " row = " << row << std::endl;
  Xyce::dout() << "adding (*dFdxMatPtr)[ " << lidIndexVector[index] << " ][ " << jacobianOffsets[row][vOffset ] << " ] = " << gMem_ * segArea << std::endl;

  // add contribution from current equation
  int numCurrentEq = currentEqusExpRCP_.size();
  for( int i=0; i<numCurrentEq; i++)
  {
    // get the number of vars in this expression and update their values
    int numvars = currentEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      currentEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ currentEqusVarNames_[i][j] ] ] ];
    }
    // evaluate the expression
    double resultValue=0.0;
    std::vector<double> derivValues;
    derivValues.resize( numvars );
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    currentEqusExpRCP_[i]->evaluate( resultValue, derivValues, currentEqusVarValues_[i] );
#else
    currentEqusExpRCP_[i]->evaluate( resultValue, derivValues );
#endif

    // now we have the derivValues[] but the order is determined by currentEqusVarNames_.  We
    // need to convert variable names to the appropriate indices into jacobianOffsets
    for( int j=0; j<numvars; j++)	// this loops through the variables AS THEY APPEAR in the current equation
    {
      std::string varName = currentEqusVarNames_[i][j];
      Xyce::dout() << "contribution of variable " << j << "(" << varName << ") to V" << std::endl;
      // cew jacOffsetsCol comes from the first row of systemJacOffset_ (for V) and the name of this variable
      // here, we have to allow for the fact that nonzero entries may be shifted by vOffset
      int jacOffsetsCol = systemJacOffset_[0][varName] + vOffset;
      Xyce::dout() << " jacOffsetsRow: " << row << " jacOffsetsCol: " << jacOffsetsCol << std::endl;
      // cew need to multiply these terms by segArea
      (*dFdxMatPtr)[lidIndexVector[index]][ jacobianOffsets[row][jacOffsetsCol] ] += segArea*derivValues[j];
      Xyce::dout() << "adding (*dFdxMatPtr)[ " << lidIndexVector[index] << " ][ " << jacobianOffsets[row][jacOffsetsCol] << " ] = " << segArea*derivValues[j] << std::endl;
    }
  }

  // add contribution from extra vars
  int numExp = fEqsExpRCP_.size();
  for( int i=0; i<numExp; i++)
  {
    std::string targetVarName = offsetToIndepVar_[i+1];
    Xyce::dout() << "F terms for f equation " << i << "; assumed to be for variable " << targetVarName << std::endl;
    // get the number of vars in this expression and update their values
    int numvars = fEqsEqusVarNames_[i].size();
    for( int j=0; j<numvars; j++)
    {
      fEqsEqusVarValues_[i][j]=(*solnVecPtr)[ lidIndexVector[index + indepVarOffset_[ fEqsEqusVarNames_[i][j] ] ] ];
    }
    // evaluate the expression
    double resultValue=0.0;
    std::vector<double> derivValues;
    derivValues.resize( numvars );
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    fEqsExpRCP_[i]->evaluate( resultValue, derivValues, fEqsEqusVarValues_[i] );
#else
    fEqsExpRCP_[i]->evaluate( resultValue, derivValues );
#endif

    // add it to the F vector
    for( int j=0; j<numvars; j++)	// this loops through variables AS THEY APPEAR in F eqn i
    {
      std::string varName = fEqsEqusVarNames_[i][j];
      Xyce::dout() << "contribution of variable " << j << "(" << varName << ")" << std::endl;
      int targetVarOffset = indepVarOffset_[targetVarName];
      int targetVarIndex = index + targetVarOffset;
      int targetVarRow = row + targetVarOffset;
      int jacOffsetsCol = systemJacOffset_[targetVarOffset][varName];
      Xyce::dout() << "targetVarOffset: " << targetVarOffset << " targetVarIndex: " << targetVarIndex
          << ", jacOffsetsRow: " << targetVarRow << " jacOffsetsCol: " << jacOffsetsCol << std::endl;
      (*dFdxMatPtr)[lidIndexVector[targetVarIndex]][ jacobianOffsets[targetVarRow][jacOffsetsCol] ] += derivValues[j];
      Xyce::dout() << "Adding (*dFdxMatPtr)[ " << lidIndexVector[targetVarIndex] << " ][ " <<  jacobianOffsets[targetVarRow][jacOffsetsCol] << " ] = " << derivValues[j] << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::convertStringsToExpression
// Purpose       : convert vectors of strings to expressions
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::convertStringsToExpression( std::vector< std::string > & stringInput, std::vector<RCP<Util::Expression> > & expRCPOut )
{
  int numStrings = stringInput.size();
  for( int i=0; i<numStrings; i++ )
  {
    //Xyce::dout() << "Making expression from :" << stringInput.at(i) << std::endl;
#if 0
    expRCPOut.push_back( rcp( new Util::Expression( solState.expressionGroup_, stringInput.at(i) ) ) );
#else
    expRCPOut.push_back( rcp( new Util::Expression( solState.getGroupWrapper()->expressionGroup_, stringInput.at(i) ) ) );
#endif
    /*
    int type=0;
    std::vector<string> names;
    expRCPOut.at(i)->get_names(type, names);

    Xyce::dout() << " type = " << type << std::endl;
    for( int j=0; j<names.size(); j++)
    {
      Xyce::dout() << "name[ " << j << " ] = " << names.at(j) << std::endl;
    }
    */

  }
  return;
}


//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::consolidateExpressions
// Purpose       : simplify parameters and functions in membrane equations
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::consolidateExpressions()
{
  // loop through parameters and try to set their value in the other expressions
  int numParams = extraParametersExpRCP_.size();
  for( int i=0; i<numParams; i+=2)
  {
    std::vector<std::string> names;
#if 0
    int type=0;
    extraParametersExpRCP_.at(i)->get_names(type, names);
#else
    // ERK. not clear what type of param is needed.  type=0 implies "give me everything".
    // whatever it is, it needs to be handled via group anyway.
    //extraParametersExpRCP_.at(i)->
#endif
    if( names.size() > 1 )
    {
      Xyce::dout() << "Warning MembraneUserDefined::consolidateExpressions() parameter name vec longer than expected." << names.size() << std::endl;
    }
    paramNames_.push_back( names[0] );
    double value=0;
    extraParametersExpRCP_.at(i+1)->evaluateFunction( value );
    paramValues_.push_back( value );
  }

  // loop over functions to get the function name to function value association
  int numFuncs = extraFunctionsExpRCP_.size();
  for( int i=0; i<numFuncs; i+=2 )
  {
    std::vector<std::string> names;
#if 0
    int type=0;
    extraFunctionsExpRCP_.at(i)->get_names(type, names);
#else
    // ERK. not clear what type of param is needed.  type=0 implies "give me everything".
    // whatever it is, it needs to be handled via group anyway.
    //extraFunctionsExpRCP_.at(i)->get
#endif
    funcNames_.push_back(names[0]);
    funcExpRCP_.push_back(extraFunctionsExpRCP_.at(i+1));
#if 0
    funcNumArgs_.push_back( extraFunctionsExpRCP_.at(i+1)->num_vars() );
#endif
  }

  // substitute parameters into current expression, extra equations and functions
  substituteParameters( currentEqusExpRCP_ );
  substituteParameters( fEqsExpRCP_ );
  substituteParameters( qEqsExpRCP_ );
  substituteParameters( funcExpRCP_ );

  // now substitute functions
  substituteFunctions( currentEqusExpRCP_ );
  substituteFunctions( fEqsExpRCP_ );
  substituteFunctions( qEqsExpRCP_ );

  // convert remaining symbols into variables
  makeSymbolSet();
  convertSymbolsToVars( currentEqusExpRCP_, currentEqusVarNames_, currentEqusVarValues_ );
  convertSymbolsToVars( fEqsExpRCP_, fEqsEqusVarNames_, fEqsEqusVarValues_ );
  convertSymbolsToVars( qEqsExpRCP_, qEqsEqusVarNames_, qEqsEqusVarValues_ );


  // debugging output
  for( int i=0; i<paramNames_.size(); i++ )
  {
    Xyce::dout() << "Param \"" << paramNames_.at(i) << "\" = " << paramValues_.at(i) << std::endl;
  }

  for( int i=0; i<funcNames_.size(); i++ )
  {
    Xyce::dout() << "Func \"" << funcNames_.at(i) << "\" = " << funcExpRCP_.at(i)->get_expression() << " num_args = " << funcNumArgs_.at(i)  << std::endl;
  }

  for( int i=0; i<currentEqusExpRCP_.size(); i++ )
  {
    Xyce::dout() << "Currentequ = " << currentEqusExpRCP_.at(i)->get_expression() << std::endl;
  }

  for( int i=0; i<fEqsExpRCP_.size(); i++ )
  {
    Xyce::dout() << "F: " << fEqsExpRCP_.at(i)->get_expression() << " Q: " << qEqsExpRCP_.at(i)->get_expression() << std::endl;
    double fval=0, qval=0;
    std::vector<double> fargs, qargs;
    fargs.push_back(0.2); fargs.push_back(-0.015);  qargs.push_back(0.2);
#if 0
    // ERK.  FIX THIS!   commenting out so this will compile
    fEqsExpRCP_.at(i)->evaluateFunction(fval, fargs);
    qEqsExpRCP_.at(i)->evaluateFunction(qval, qargs);
#else
    fEqsExpRCP_.at(i)->evaluateFunction(fval);
    qEqsExpRCP_.at(i)->evaluateFunction(qval);
#endif
    Xyce::dout() << "F(V=-0.015,n/m/h=0.2) = " << fval << " Q(V=-0.015,n/m/h=0.2) = " << qval << std::endl;
  }



  return;
}


//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::substituteParameters
// Purpose       : substitute parameter values in a vector of expressions
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::substituteParameters( std::vector<RCP<Util::Expression> > & expRCP_ )
{
  int numParams = paramNames_.size();
  for( int i=0; i<numParams; i++)
  {
    int numExp = expRCP_.size();
    for( int j=0; j<numExp; j++ )
    {
      expRCP_.at(j)->make_constant(paramNames_[i], paramValues_[i]);
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::substituteFunctions
// Purpose       : substitute user supplied functions in a vector of expressions
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::substituteFunctions( std::vector<RCP<Util::Expression> > & expRCP_ )
{
  int numFuncs = funcNames_.size();
  for( int i=0; i<numFuncs; i++)
  {
    int numExp = expRCP_.size();
    for( int j=0; j<numExp; j++ )
    {
      //expRCP_.at(j)->replace_func(funcNames_[i], *funcExpRCP_[i], funcNumArgs_[i]);
      //
      // ERK.  This is a neccessary but not sufficient change for the new expression library.
      //
      // replace_func is for the old expresison library, and it just does a string replacement.
      //
      expRCP_.at(j)->attachFunctionNode(funcNames_[i], *funcExpRCP_[i]);
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::makeSymbolSet
// Purpose       : Scan the user supplied equations to get the variables
//                 and make a map of variable name to offset that will
//                 be used in loads
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::makeSymbolSet()
{
  // get the user defined names
  // this assumes that the user hasn't repeated a name.  I should trap for this.
  // I used a set<string> for this initially, but the ordering isn't unique so the
  // vars would end up in a strange order.
  int type=0;
  int numIndependentVars = indepVarsExpRCP_.size();		// cew 1 less than numIndependentVars_, which includes V
  for( int i=0; i<numIndependentVars; i++ )
  {
    std::vector<std::string> expNames;
#if 0
    indepVarsExpRCP_.at(i)->get_names( type, expNames );
#else
    // ERK. not clear what type of param is needed.  type=0 implies "give me everything".
    // whatever it is, it needs to be handled via group anyway.
    //indepVarsExpRCP_.at(i)->get
#endif
    int numExpNames=expNames.size();
    for( int j=0; j<numExpNames; j++ )					// cew ? isn't there just one name per independent var?
    {
      Xyce::dout() << "makeSymbolSet, i=" <<i << ", j="<<j << ":  adding " << expNames.at(j) << " to userDefinedNames" << std::endl;
      userDefinedNames_.push_back( expNames.at(j) );
    }
  }

  // TODO:  check through the F and Q and current equations and make sure there aren't any variables that aren't
  // in userDefinedNames_
  // actually, methods called from consolidateExpressions kind of does this;
  // currently just prints a warning, sets value to 0, and continues

  // now userDefinedNames_ contains all the vars the user supplied
  // assign them unique offsets so that we can consistently index them in
  // the solution, F, Q and jacobian.

  std::vector<std::string>::iterator currName = userDefinedNames_.begin();
  std::vector<std::string>::iterator endName = userDefinedNames_.end();
  int offsetCount=0;
  while( currName != endName )
  {
    indepVarOffset_[ *currName ] = ++offsetCount;   // set up map
    offsetToIndepVar_[offsetCount] = *currName;     // and inverse map
    currName++;
  }

  // now add "V" to userDefinedNames_ and give it an offset of zero in
  // the indepVarOffset_ map: "V" must be offset zero -- it's first.
  userDefinedNames_.push_back("V");  // V is automatically part of the system of equations.
  indepVarOffset_["V"] = 0;
  offsetToIndepVar_[0]="V";

  // debugging output
  Xyce::dout() << "MembraneUserDefined::makeSymbolSet() Independent var offset map : " << std::endl;
  std::map<std::string,int>::iterator currVarOffset = indepVarOffset_.begin();
  std::map<std::string,int>::iterator endVarOffset = indepVarOffset_.end();
  while( currVarOffset != endVarOffset )
  {
    Xyce::dout() << "map[ " << currVarOffset->first << " ] = " << currVarOffset->second << std::endl;
    currVarOffset++;
  }
  std::map<int,std::string>::iterator currOffset = offsetToIndepVar_.begin();
  std::map<int,std::string>::iterator endOffset = offsetToIndepVar_.end();
  while( currOffset != endOffset )
  {
    Xyce::dout() << "inv-map[ " << currOffset->first << " ] = " << currOffset->second << std::endl;
    currOffset++;
  }

}


//-----------------------------------------------------------------------------
// Function      : MembraneUserDefined::convertSymbolsToVars
// Purpose       : Using user supplied list and what's left in an expression,
//                 convert remaining symbols to vars so that we can set them
//                 and evaluate the expression
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 09/08/2011
//-----------------------------------------------------------------------------
void MembraneUserDefined::convertSymbolsToVars( std::vector<RCP<Util::Expression> > & expRCP, std::vector< std::vector<std::string> > & expNamesVec, std::vector< std::vector<double> > & expValsVec )
{
  int type=0;
  std::vector<std::string>::iterator endName = userDefinedNames_.end();
  int numExp = expRCP.size();
  if( numExp > 0 )
  {
    expNamesVec.resize( numExp );
    expValsVec.resize( numExp );
  }
  for( int i=0; i<numExp; i++ )
  {
    // get the remaining symbol names in the expression
    std::vector<std::string> expNamesTmp;
#if 0
    expRCP.at(i)->get_names( type, expNamesTmp );
#else
    // ERK. not clear what type of param is needed.  type=0 implies "give me everything".
    // whatever it is, it needs to be handled via group anyway.
    //expRCP.at(i)->get
#endif

    int numExpNames=expNamesTmp.size();
    for( int j=0; j<numExpNames; j++ )
    {
      // if expNames[j] is in our set of known symbols that
      // are variables, then set its type to variable.  Otherwise
      // issue a warning.
      std::vector<std::string>::iterator findResult = find( userDefinedNames_.begin(), userDefinedNames_.end(), expNamesTmp.at(j) );
      if( findResult != endName )
      {
#if 0
        // found this one in the set so make it a variable
        expRCP.at(i)->make_var( expNamesTmp.at(j), 0.0 );
#else
        expRCP.at(i)->make_constant( expNamesTmp.at(j), 0.0 );  // ERK will this work?
#endif
        // save it in a std::vector<std::string> so we can quickly load values into a std::vector<double> to evaluate the expression.
        expNamesVec[i].push_back( expNamesTmp.at(j) );
      }
      else
      {
        Xyce::dout() << "Warning symbol \"" << expNamesTmp.at(j) << "\" not listed in independent variables.  Ignoring for now." << std::endl;
      }
    }
    // make the values vec the same size as the names vec.
    expValsVec[i].resize( expNamesVec[i].size() );

    Xyce::dout() << "MembraneUserDefined::convertSymbolsToVars: expression " << expRCP.at(i)->get_expression() << " Has vars: ";
    std::vector<std::string>::iterator cni = expNamesVec[i].begin();
    std::vector<std::string>::iterator eni = expNamesVec[i].end();
    while( cni != eni )
    {
      Xyce::dout() << *cni << ", ";
      cni++;
    }
    Xyce::dout() << std::endl;

  }
}

} // namespace Device
} // namespace Xyce
