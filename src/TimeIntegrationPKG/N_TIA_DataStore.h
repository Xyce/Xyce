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
// Purpose       : This file handles the class that defines the data arrays
//                 needed for the time integration algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_DataStore_h
#define Xyce_N_TIA_DataStore_h

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_TIA_fwd.h>

#include <N_TIA_TwoLevelError.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : DataStore
// Purpose       : This is the class for defining data arrays needed in the
//                 time integration algorithms.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class DataStore
{
  public:
    DataStore(int max_order, const Linear::Builder &linear_builder);
    ~DataStore();

    void allocateSensitivityArrays(int nP, 
        bool includeTransientDirect,
        bool includeTransientAdjoint);

    void deleteSensitivityArrays();

    void allocateHBVectors();
    void allocateWaMPDEVectors();

  private:
    DataStore(const DataStore& rhs);
    DataStore &operator=(const DataStore& rhs);

  // DataStore Functions
  public:

    void updateSolDataArrays();
    bool updateStateDataArrays();
    void updateSolDataArraysAdjoint (int timeIndex);

    void setConstantHistory();
    void setConstantSensitivityHistory();
    void setZeroHistory();
    void setConstantHistoryAdjoint ();

    void setErrorWtVector(const TIAParams &tia_params, const std::vector<char> &     variable_type);
    double WRMS_errorNorm();

    bool equateTmpVectors ();
    bool usePreviousSolAsPredictor ();

    void stepLinearCombo ();

    double partialErrorNormSum();
    double partialQErrorNormSum();

    double partialSum_m1(int currentOrder);
    double partialSum_p1(int currentOrder, int maxOrder);
    double partialSum_q1();

    double delta_x_errorNorm_q1();
  
    bool getSolnVarData( const int & gid, std::vector<double> & varData );
    bool getStateVarData( const int & gid, std::vector<double> & varData );
    bool setSolnVarData( const int & gid, const std::vector<double> & varData );
    bool setStateVarData( const int & gid, const std::vector<double> & varData );
    bool getStoreVarData( const int & gid, std::vector<double> & varData );
    bool setStoreVarData( const int & gid, const std::vector<double> & varData );

    int getNumSolnVarData() const { return 11; }
    int getNumStateVarData() const { return 7; }
    int getNumStoreVarData() const { return 3; }
  
    bool setNextSolVectorPtr (Linear::Vector * solVecPtr);
    bool setNextSolVectorPtr (Linear::Vector & solVecPtr);// TT: added
    bool unsetNextSolVectorPtr ();
  
    bool resetAll(double absolute_error_tolerance, double relative_error_tolerance);

    bool resetFastTimeData ();

  public:

    const Linear::Builder&        builder_;
 
    // limiter flag:
    bool limiterFlag;

    // TIA Arrays (pointers) for  Integration Solution Process:
    unsigned int maxOrder;
    unsigned int solutionSize;
    unsigned int stateSize;
    unsigned int storeSize;
    unsigned int leadCurrentSize;

    // temporary vectors:
    Linear::Vector * tmpSolVectorPtr;
    Linear::Vector * tmpStaVectorPtr;
    Linear::Vector * tmpStaDerivPtr;
    Linear::Vector * tmpStoVectorPtr;
    Linear::Vector * tmpLeadCurrentVectorPtr;
    Linear::Vector * tmpLeadDeltaVPtr;

    Linear::Vector * tmpLeadCurrentQDerivVectorPtr;

    // Predictors
    Linear::Vector * xn0Ptr;

    // Solutions:
    Linear::Vector * currSolutionPtr;
    Linear::Vector * lastSolutionPtr;
    Linear::Vector * nextSolutionPtr;

    // Used for pointer switching, related to nextSolPtrSwitched_
    Linear::Vector * savedNextSolutionPtr;

    // States:
    Linear::Vector * currStatePtr;
    Linear::Vector * lastStatePtr;
    Linear::Vector * nextStatePtr;

    // Storage:
    Linear::Vector * currStorePtr;
    Linear::Vector * nextStorePtr;
    
    // Lead current and power vectors 
    Linear::Vector * currLeadCurrentPtr;
    Linear::Vector * nextLeadCurrentPtr;
    
    Linear::Vector * currLeadDeltaVPtr;
    Linear::Vector * nextLeadDeltaVPtr;
    
    // for lead current calculations.  F component is
    // held in the store vector, Q component is here
    Linear::Vector * currLeadCurrentQPtr;
    Linear::Vector * nextLeadCurrentQPtr;    

    // Number of sensitivity parameters
    int numParams;

    // sensitivity vectors:
    Linear::MultiVector* sensRHSPtrVector;

    // adjoint sparse storage experiment
    Linear::FilteredMultiVector * sparseSensRHSMV;

    Linear::MultiVector* nextDfdpPtrVector;

    Linear::MultiVector* currDqdpPtrVector;
    Linear::MultiVector* nextDqdpPtrVector;

    Linear::MultiVector* nextDbdpPtrVector;

    Linear::MultiVector* currDXdpPtrVector;
    Linear::MultiVector* nextDXdpPtrVector;

    // adjoint sparse storage experiment
    std::vector<std::vector<int> > masterIndexVector;
    std::vector<int> masterIndexVectorSize;

    // matvecs:
    Linear::MultiVector* currDQdxDXdpPtrVector;
    Linear::MultiVector* lastDQdxDXdpPtrVector;
    Linear::MultiVector* nextDQdxDXdpPtrVector;

    Linear::MultiVector* currDFdxDXdpPtrVector;
    Linear::MultiVector* lastDFdxDXdpPtrVector;
    Linear::MultiVector* nextDFdxDXdpPtrVector;

    // Adjoint sensitivity solutions:
    Linear::Vector* nextLambdaPtr;
    Linear::Vector* currLambdaPtr;
    Linear::Vector* lastLambdaPtr;

    Linear::Vector* nextDQdxLambdaPtr;
    Linear::Vector* currDQdxLambdaPtr;
    Linear::Vector* lastDQdxLambdaPtr;

    Linear::Vector* nextDFdxLambdaPtr;
    Linear::Vector* currDFdxLambdaPtr;
    Linear::Vector* lastDFdxLambdaPtr;

    // Derivatives of States:
    Linear::Vector * currStateDerivPtr;
    Linear::Vector * nextStateDerivPtr;

    // Derivatives of Store for lead curent calculations
    Linear::Vector * currLeadCurrentQDerivPtr;
    Linear::Vector * nextLeadCurrentQDerivPtr;

    // Derivatives of dq/dp for sensitivity calculations
    Linear::MultiVector* nextDqdpDerivPtrVector;

    // Error Vectors
    Linear::Vector * errWtVecPtr;

    // Jacobian and RHS (pointers to objects in linear system)
    Linear::Matrix * JMatrixPtr;
    Linear::Vector * RHSVectorPtr;

    // NonLinear Solution Vectors
    Linear::Vector * newtonCorrectionPtr;
    Linear::Vector * qNewtonCorrectionPtr;

    // Mask for error norms (to allow some equations not to take part in
    // weighted norms)
    Linear::Vector * deviceErrorWeightMask_;

    // 2-level information:
    std::vector<TwoLevelError> innerErrorInfoVec;

    // new-DAE data (originally from the new-DAE derrived class)
    // Error Vectors
    Linear::Vector * qErrWtVecPtr;

    // DAE formulation vectors
    Linear::Vector * daeQVectorPtr;
    Linear::Vector * daeFVectorPtr;
    Linear::Vector * daeBVectorPtr;

    // DAE formulation matrices
    Linear::Matrix * dQdxMatrixPtr;
    Linear::Matrix * dFdxMatrixPtr;

    // HB temporary Matvec storage vectors
    Linear::Vector * dQdxVecVectorPtr;
    Linear::Vector * dFdxVecVectorPtr;

    // voltage limiting vectors
    Linear::Vector * dFdxdVpVectorPtr;
    Linear::Vector * dQdxdVpVectorPtr;

    // History arrays
    std::vector<Linear::Vector*> xHistory;
    std::vector<Linear::Vector*> qHistory;
    std::vector<Linear::Vector*> sHistory;    // state history
    std::vector<Linear::Vector*> stoHistory;  // store history
    std::vector<Linear::Vector*> leadCurrentHistory;  // history for lead current Q component.
    std::vector<Linear::Vector*> leadCurrentQHistory;  // history for lead current Q component.
    std::vector<Linear::Vector*> leadCurrentQDerivHistory;  // history for lead current dQ/dt component.
    std::vector<Linear::Vector*> leadDeltaVHistory;  // history for junction voltage
    
    // sensitivity histories
    std::vector< Linear::MultiVector* > dbdpHistory;
    std::vector< Linear::MultiVector* > dfdpHistory;
    std::vector< Linear::MultiVector* > dqdpHistory;

    // histories used for transient adjoints:
    int itAdjointIndex;
    bool adjointDcop;
    std::vector< int > orderHistory;
    std::vector< double > dtHistory;
    std::vector< double > timeHistory;
    std::vector< Linear::Vector *> solutionHistory;
    std::vector< Linear::Vector *> stateHistory;
    std::vector< Linear::Vector *> storeHistory;

    // outer loop over parameter list, inner is the history
    std::vector< Linear::MultiVector * > functionSensitivityHistory;

    // adjoint sparse storage experiment
    std::vector< Linear::FilteredMultiVector* > sparseFunctionSensitivityHistory;

    // adjoint tmp matrix storage
    Linear::Matrix * tmpMatrixPtr;

    // Predictors
    Linear::Vector * qn0Ptr;

    // Step-size selection temporary vectors for two-level Newton
    Linear::Vector * delta_x;

    // Temporary vectors for WaMPDE interpolation
    Linear::Vector * tmpXn0APtr;
    Linear::Vector * tmpXn0BPtr;

    // These are for MPDE fast time scale points
    std::vector<double> timeSteps;
    std::vector<bool> timeStepsBreakpointFlag;
    std::vector<Linear::Vector*> fastTimeSolutionVec;
    std::vector<Linear::Vector*> fastTimeStateVec;
    std::vector<Linear::Vector*> fastTimeQVec;
    std::vector<Linear::Vector*> fastTimeStoreVec;

    std::vector<double> objectiveVec_;
    std::vector<double> dOdpVec_;
    std::vector<double> dOdpAdjVec_;
    std::vector<double> scaled_dOdpVec_;
    std::vector<double> scaled_dOdpAdjVec_;
    std::vector<double> paramOrigVals_; 

  private:
    bool nextSolPtrSwitched_;

    std::vector<int> indexIVars;
    std::vector<int> indexVVars;
    std::vector<int> indexMaskedVars;

    double absErrTol_, relErrTol_;
    double solsMaxValue;

    // Vectors for new LTE strategies
    Linear::Vector * maxSolutionPtr;
    Linear::Vector * relSolutionPtr;
    int index; 

    bool allocateSensitivityArraysComplete_;
    bool includeTransientAdjoint_;
    bool includeTransientDirect_;
};

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_DataStore_h

