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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/03/06
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_NLS_ConductanceExtractor.h>
#include <N_NLS_Manager.h>
#include <N_NLS_TwoLevelPrintJac.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_TOP_Topology.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

// ----------   Static Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::ConductanceExtractor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
ConductanceExtractor::ConductanceExtractor(
  NonLinearSolver & nls,
  Topo::Topology & topTmp)
  : gidsSetUpFlag_(false),
    nls_(nls),
    top_(topTmp),
    lasSysPtr_(0),
    rhsVectorPtr_(0),
    dfdvVectorPtr_(0),
    NewtonVectorPtr_(0),
    dxdvVectorPtr_(0),
    lasSolverPtr_(0),
    matrixDiagonalPtr_(0),
    jacobianMatrixPtr_(0),
    savedRHSVectorPtr_(0),
    savedNewtonVectorPtr_(0),
    gradVectorPtr_(0),
    columnVectorPtr_(0),
    columnMapPtr_(0)
{
  lasSysPtr_    = nls_.lasSysPtr_;
  rhsVectorPtr_ = nls_.rhsVectorPtr_;
  dfdvVectorPtr_ = nls_.rhsVectorPtr_; // using the same vector for dfdv.

  NewtonVectorPtr_      = nls_.NewtonVectorPtr_;
  dxdvVectorPtr_        = nls_.NewtonVectorPtr_; // using same vector for dxdv.
  lasSolverPtr_         = nls_.lasSolverRCPtr_.get();
  jacobianMatrixPtr_    = nls_.jacobianMatrixPtr_;

  // creations
  savedRHSVectorPtr_    = lasSysPtr_->builder().createVector();
  savedNewtonVectorPtr_ = lasSysPtr_->builder().createVector();
  matrixDiagonalPtr_    = lasSysPtr_->builder().createVector();

#ifdef Xyce_PARALLEL_MPI
  // construct column vector used for parallel construction of RHS's
  columnMapPtr_ = jacobianMatrixPtr_->getColMap( savedRHSVectorPtr_->pmap()->pdsComm() );
  columnVectorPtr_ = new Linear::Vector( *(savedRHSVectorPtr_->pmap()), *columnMapPtr_ );
#endif
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::~ConductanceExtractor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
ConductanceExtractor::~ConductanceExtractor()
{
  delete savedRHSVectorPtr_;
  delete savedNewtonVectorPtr_;
  delete matrixDiagonalPtr_;
  for (  std::vector<Linear::Vector *>::iterator it = dIdxPtrVector_.begin(), end = dIdxPtrVector_.end(); it != end; ++it)
    delete *it;
  delete columnVectorPtr_;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::setOptions
// Purpose       :
// Special Notes : ERK. 2/10/2019.  This is apparently never called. FIX?
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool ConductanceExtractor::setOptions(const Util::OptionBlock & option_block)
{
  bool bsuccess = true;

  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    if ((*it).uTag() == "DEBUGLEVEL")
    {
      setNonlinearConductanceDebugLevel((*it).getImmutableValue<int>());
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::setupIDs_
//
// Purpose       : This function sets up the various GIDs and matrix rows,
//                 that are needed in the extract function.
//
// Special Notes : There are 2 types of IDs needed, which are both associated
//                 with independent voltage sources.  One type is the GID of
//                 the current variable, which is needed for a "get" operation.
//                 The other is the LID of the matrix row for the voltage drop
//                 equation.  In serial, this is the same as the current GID,
//                 but in parallel they will generally be different.  For the
//                 matrix row, the LID is used for a "put" operation, so it
//                 makes more sense for it to be an LID.
//
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool ConductanceExtractor::setupIDs_(const std::map<std::string, double> &  inputMap)
{
  bool bsuccess = true;

  // Setup vectors for each terminal current
  int idSize = inputMap.size();
  currentGIDs_.resize(idSize);
  currentLIDs_.resize(idSize);
  vsrcPosGIDs_.resize(idSize);
  vsrcPosLIDs_.resize(idSize);
  for (int i=0; i<idSize ;++i)
  {
    dIdxPtrVector_.push_back(lasSysPtr_->builder().createVector());
    currentGIDs_[i] = -1;
    currentLIDs_[i] = -1;
    vsrcPosGIDs_[i] = -1;
    vsrcPosLIDs_[i] = -1;
  }

  // Loop over the map, which should contain, in the first
  // arguments, the names of the voltage sources we need.
  std::map<std::string,double>::const_iterator iterM = inputMap.begin();
  std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
  int i=0;
  for (; iterM != endM; ++i, ++iterM)
  {
    // Note: When passing the sourceName to topology, it must
    // be all CAPS, or topology won't find it.
    ExtendedString src = iterM->first;
    src.toUpper();
    std::string sourceName = src;
    char type;

    // This is to get the IDs for the currents through the
    // voltage sources specified in the map.
    std::vector<int> GIDList, extGIDList;
    top_.getNodeSVarGIDs(NodeID(sourceName, Xyce::_DNODE), GIDList, extGIDList, type);

    std::vector<int>::iterator iterI;
    if (!(GIDList.empty ()))
    {
      iterI = GIDList.begin();
      currentGIDs_[i] = *iterI;
      currentLIDs_[i] = dIdxPtrVector_[0]->pmap()->globalToLocalIndex(currentGIDs_[i]);

      iterI = extGIDList.begin();
      vsrcPosGIDs_[i] = *iterI;
      vsrcPosLIDs_[i] = dIdxPtrVector_[0]->pmap()->globalToLocalIndex(vsrcPosGIDs_[i]);

      if( vsrcPosLIDs_[i] == -1 )
      {
        Xyce::Report::DevelFatal().in("ConductanceExtractor::setupIDs_")
          << " The " << sourceName << " source has the positive node"
          << " owned by another processor.  The 2-level solve can't handle that.";
      }

      // check that vneg is connected to gnd
      ++iterI;
      int vnegGID = *iterI;
      if (vnegGID != -1)
      {
        Xyce::Report::DevelFatal().in("ConductanceExtractor::setupIDs_")
          << " The " + sourceName + " source has the negative node"
          << " connected to something other than ground!  The 2-level solve can't handle that.";
      }
    }

  }

  if (DEBUG_CONDUCTANCE && isActive(Diag::NONLINEAR_CONDUCTANCE))
  {
    Xyce::dout() << "current GIDs: " << std::endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      Xyce::dout() << "  currentGIDs_["<<i1<<"] = " << currentGIDs_[i1] << ", currentLIDs_["<<i1<<"] = " << currentLIDs_[i1] << std::endl;
    }

    Xyce::dout() << "Vsrc pos equation rows: " << std::endl;
    for( int i1=0; i1 < idSize; ++i1 )
    {
      Xyce::dout() << "  vsrcPosGIDs_["<<i1<<"] = " << vsrcPosGIDs_[i1] << ", vsrcPosLIDs_["<<i1<<"] = " << vsrcPosLIDs_[i1] << std::endl;
    }
  }

  gidsSetUpFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::setup_dIdX_Vectors_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
bool ConductanceExtractor::setup_dIdX_Vectors_ ()
{
  bool bsuccess = true;

  // Set up dIdx's.  These correspond to rows in the Jacobian.
  // There will be one for each Vsrc current.   In general, we are assuming
  // that each "connecting" voltage source is connected to ground at the
  // negative node, and to the rest of the circuit at the positive node,
  // so the node of interest (and thus KCL of interest) is the
  // positive node.  dIdx is a vector containing all the derivatives
  // from the corresponding KCL row, minus the derivative with respect to
  // I.  (I is a solution variable, but needs to be excluded from dIdx,
  // or the various dot products will cancel some terms that they
  // should not).
  int idSize = currentGIDs_.size();

  Linear::Vector * currentVec;

  for( int iC_row=0; iC_row < idSize; ++iC_row )
  {
#ifdef Xyce_PARALLEL_MPI
    currentVec = columnVectorPtr_;
#else
    currentVec = dIdxPtrVector_[iC_row];
#endif

    currentVec->putScalar(0.0);

    if( currentGIDs_[iC_row] != -1 )
    {
      int iRow = vsrcPosGIDs_[iC_row];
      int rowLength = jacobianMatrixPtr_->getRowLength(iRow);
      int numEntries = rowLength;
      std::vector<double> coeffs(rowLength, 0.0);
      std::vector<int> colIndices(rowLength, -1);

      jacobianMatrixPtr_->getRowCopy
        (iRow, rowLength, numEntries, &coeffs[0], &colIndices[0]);

      for (int ic=0;ic<rowLength;++ic)
      {
        // need to exclude entries that are with respect to
        // the the variable, 'I'.
        int gid = colIndices[ic];
        if (gid ==  currentGIDs_[iC_row]) coeffs[ic] = 0.0;
      }

      for (int icol=0;icol<rowLength;++icol)
      {
        double val = coeffs[icol];
        int gid = colIndices[icol];
        if (gid != -1)
          currentVec->setElementByGlobalIndex(gid, val, 0);
      }
    }
    currentVec->fillComplete();

#ifdef Xyce_PARALLEL_MPI
    *(dIdxPtrVector_[iC_row]) = *columnVectorPtr_;
#endif

    if (DEBUG_CONDUCTANCE && isActive(Diag::NONLINEAR_CONDUCTANCE))
    {
      Xyce::dout() << "\ndIdx[" << iC_row << "]:" << std::endl;
      dIdxPtrVector_[iC_row]->print(Xyce::dout());
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::extract
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/06
//-----------------------------------------------------------------------------
bool ConductanceExtractor::extract(
  const std::map<std::string,double> &  inputMap,
  std::vector<double> &                 outputVector,
  std::vector< std::vector<double> > &  jacobian)
{
  bool bsuccess = true;

  if (DEBUG_CONDUCTANCE && isActive(Diag::NONLINEAR_CONDUCTANCE))
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "ConductanceExtractor::extract" << std::endl;
    Xyce::dout() << subsection_divider << std::endl;
  }

  if (inputMap.empty() ||
      outputVector.empty() ||
      jacobian.empty() )
  {
    return false;
  }

  // Need the "next" vector.
  Linear::Vector * solnVecPtr = (nls_.dsPtr_)->nextSolutionPtr;

  if (!gidsSetUpFlag_)
  {
    setupIDs_(inputMap);
  }

  // Now that the solve is complete, obtain the currents and conductances
  // that are needed for the upper level nonlinear solve.

  // The currents are owned by the individual sources in the device package:
  // Note: this may need some parallel refactoring.

  // First obtain the currents.
  int idSize = currentGIDs_.size();
  for (unsigned int i=0;i<currentGIDs_.size();++i)
  {
    int index1 = currentGIDs_[i];
    if( index1 > -1 )
      outputVector[i] = solnVecPtr->getElementByGlobalIndex(index1);
    else
      outputVector[i] = 0.0;
  }

#ifdef Xyce_PARALLEL_MPI
  //sumAll to get all currents locally
  N_PDS_Comm &comm = dfdvVectorPtr_->pmap()->pdsComm();
  std::vector<double> tmpVector(idSize,0.0);
  for( int i = 0; i < idSize; ++i )
  {
    tmpVector[i] = outputVector[i];
    outputVector[i] = 0.0;
  }

  comm.sumAll( &(tmpVector[0]), &(outputVector[0]), idSize );
#endif

  if (DEBUG_CONDUCTANCE && isActive(Diag::NONLINEAR_CONDUCTANCE))
  {
    unsigned int itmp;
    for (itmp=0;itmp < outputVector.size();++itmp)
    {
      Xyce::dout() << "currentVector["<< itmp <<"] = " << outputVector[itmp]<<std::endl;
    }
  }

  // This function needs to solve one or more linear systems.  The linear
  // systems, (for a direct sensitivity) is:
  // J * dx/dv = df/dv
  //
  // where "v" is an applied voltage, not part of the system of unknowns,
  // which is represented by "x".
  //
  // J is the traditional Jacobian matrix, df/dx.
  //
  // For the direct approach, a different dx/dv is needed for each v, so
  // a different linear system, with a different df/dv will be solved for
  // each applied voltage, v.
  //
  // For an adjoint approach, it could be possible to do a single linear
  // solve, and get the sensitivities to all the v's all at once.
  //
  // Note:  In general, this extraction is done in the context of a
  // multi-level Newton solve, in which the extracted conductances are
  // used one "level" up in the simulation.  The inputs are all voltages,
  // which are applied via voltage sources.  This means that the equation
  // used to complete the system is of the form:
  //
  //   f = V(n) - Vexternal = 0
  //
  //   where V(n) is one of the solution variables in x, and Vexternal
  //   is the voltage applied from above, represented in the above
  //   equations as "v".
  //
  //   Given that this single, simple equation is the entire coupling,
  //   the df/dv derivative, for any Vexternal, is -1.0.

  // Save the old rhs and newton vectors.  (This might not be neccessary).
  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->putScalar(0.0);
  savedRHSVectorPtr_->addVec(1.0, *(rhsVectorPtr_));

  savedNewtonVectorPtr_->putScalar(0.0);
  savedNewtonVectorPtr_->addVec(1.0, *(NewtonVectorPtr_));

  // Before we try to do any linear solves, check that the Jacobian
  // actually has been loaded.  Sometimes, early in the run, the inner
  // solve will not load the Jacobian, b/c the norm of f is zero, or very
  // small. Typically this will happen if the inner problem is linear,
  // and is entirely driven by the upper level circuit via the connected
  // Vsrc devices.

  jacobianMatrixPtr_->getDiagonal ( *matrixDiagonalPtr_ );
  double diagInfNorm = 0.0;
  matrixDiagonalPtr_->infNorm(&diagInfNorm);

  // if the infinite norm of the diagonal is zero, then the matrix
  // probably hasn't been set up yet.  If that is the case, we should
  // call loadJacobian now.  It should be safe to call this, given
  // that the RHS vector was loaded, no matter what, as part of the
  // nonlinear solve.
  if (diagInfNorm < 1.0e-30)
  {
    if (DEBUG_CONDUCTANCE)
      Xyce::Report::UserWarning0()
        << "\n\tJacobian for inner problem not loaded.  Forcing a load.";

    nls_.jacobian_();
  }

  // Get dIdx vectors, one for each I-row.
  bool b1 = setup_dIdX_Vectors_(); bsuccess = bsuccess && b1;

  // This loop is over the different applied V's.
  // This loop is also over columns of the small (output) Jacobian.
  std::map<std::string,double>::const_iterator iterM = inputMap.begin();
  std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
  int iV_col=0;
  for (;iV_col<idSize;++iV_col,++iterM)
  {
    // set up dfdv:
    dfdvVectorPtr_->putScalar(0.0);
    int irow = currentLIDs_[iV_col];
    // note: dfdv = -1.0, but -dfdv needs to go in here.
    if( irow != -1 )
      (*dfdvVectorPtr_)[irow] = 1.0;   // = -dfdv.

    int solutionStatus;
    if (iV_col==0)
    {
      // ERK. this isn't set up to track statistics such as # of linear solves, etc.
      solutionStatus=lasSolverPtr_->solve(false); // don't reuse LU factors
    }
    else
    {
      solutionStatus=lasSolverPtr_->solve(true); // DO reuse LU factors
    }

    int iC_row=0;
    for (;iC_row<idSize;++iC_row)
    {
      // Get the dot product of dIdx and dxdv, to get dIdv.
      double dIdv = dIdxPtrVector_[iC_row]->dotProduct(*(dxdvVectorPtr_));

      if (DEBUG_CONDUCTANCE && isActive(Diag::NONLINEAR_CONDUCTANCE))
      {
        std::string vsrcName = iterM->first;
        print_(vsrcName);
        Xyce::dout() << "dIdv = " << dIdv << std::endl;
      }

      // put dIdV's into the small matrix:
      jacobian[iC_row][iV_col] = dIdv;
    }

  } // cols of output Jacobian (iV_col)

  // Restore the RHS and Newton vectors. (again, this may not be necessary).
  rhsVectorPtr_->putScalar(0.0);
  rhsVectorPtr_->addVec(1.0, *(savedRHSVectorPtr_));

  NewtonVectorPtr_->putScalar(0.0);
  NewtonVectorPtr_->addVec(1.0, *(savedNewtonVectorPtr_));

  if (VERBOSE_CONDUCTANCE)
    printJacobian_ (inputMap,jacobian);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::printJacobian_
// Purpose       : Prints the small STL jacobian to the screen.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
void ConductanceExtractor::printJacobian_
    (const std::map<std::string,double> & inputMap,
     std::vector< std::vector<double> > & jacobian)
{

  // set up a names vector, then call function
  std::vector<std::string> names;
  int numElectrodes = jacobian.size();
  std::map<std::string,double>::const_iterator iterM = inputMap.begin();
  std::map<std::string,double>::const_iterator  endM = inputMap.end  ();
  for (int iE1 = 0; iE1 < numElectrodes; ++iE1,++iterM) { names.push_back(iterM->first); }

  Xyce::Nonlinear::printJacobian(Xyce::dout(),std::string(""),names,jacobian);

  return;
}

//-----------------------------------------------------------------------------
// Function      : ConductanceExtractor::print_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/08/06
//-----------------------------------------------------------------------------
void ConductanceExtractor::print_(const std::string & varName)
{
  Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
  std::string srcName = varName;
  Xyce::dout() << "Info for input voltage: " << srcName << std::endl;
  Xyce::dout() << "Jacobian:" << std::endl;
  jacobianMatrixPtr_->print(Xyce::dout());

  // now print out the dxdv vector:
  Xyce::dout() << "dxdv:" << std::endl;
  dxdvVectorPtr_->print(Xyce::dout());

  Xyce::dout() << "dfdv:" << std::endl;
  dfdvVectorPtr_->print(Xyce::dout());
}

} // namespace Nonlinear
} // namespace Xyce
