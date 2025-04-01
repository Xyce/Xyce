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

//-----------------------------------------------------------------------------
//
// Purpose        : This file contains the device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/30/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceInstance_h
#define Xyce_N_DEV_DeviceInstance_h

#include <list>
#include <map>
#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>

#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSupport.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_NodeSymbols.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceInstance  : public DeviceEntity
{
private:
  DeviceInstance();

public:
  DeviceInstance(
     const InstanceBlock &     instance_block,
     ParametricData<void> &    parametric_data,
     const FactoryBlock &      factory_block);

  virtual ~DeviceInstance();

private:
  DeviceInstance(const DeviceInstance &);
  DeviceInstance &operator=(const DeviceInstance &);

public:
  // Is the instance of this device linear.
  // NOTE: Only linear devices need to override this virtual function.
  virtual bool isLinearDevice() const { return false; }

  // NOTE: Only PDE devices need to override this virtual function.
  virtual bool isPDEDevice() const { return false; }

  virtual std::ostream &printName(std::ostream &os) const;

  // This function configures the device to request space in the store
  // vector for lead current calculations.  It must be called soon
  // after the constructor call before the store vector is allocated.
  virtual void enableLeadCurrentCalc();

  // This function returns an empty set if it's called by any device other
  // than a mutual inductor.  For the mutual inductor devices (linear and
  // non-linear), it is re-defined by those classes, and used to pass the 
  // set of inductor names to the addDeviceInstance() function, in the 
  // Device Manager.  This helps enable lead current and power calculations 
  // for mutual inductor devices.
  virtual std::vector< std::string > getInductorNames() const {
    std::vector< std::string > emptySet;
    return emptySet;
  }

  virtual std::vector< double > getInductorInductances() const {
    std::vector< double > emptySet;
    return emptySet;
  }

  virtual void setInductorInductances(std::vector< double > & set) {}

  virtual void registerLIDs( const LocalIdVector & intLIDVecRef, const LocalIdVector & extLIDVecRef )
  {}

  virtual void registerStateLIDs( const LocalIdVector & staLIDVecRef )
  {}

  virtual void registerStoreLIDs( const LocalIdVector & stoLIDVecRef )
  {}

  virtual void registerBranchDataLIDs( const LocalIdVector & branchDataLIDVecRef )
  {}

  virtual const std::vector<std::string> & getDepSolnVars();
  virtual const std::vector<int> & getDepSolnTypes();
  virtual void registerDepSolnGIDs( const std::vector< IdVector > & varList );

  virtual void registerDepSolnLIDs(const std::vector< IdVector > & depSolnLIDVecRef);
  
  virtual const JacobianStamp & jacobianStamp() const 
  {
    static JacobianStamp dummy;
    return dummy;
  }

  virtual void registerJacLIDs( const JacobianStamp & jacLIDVec );

  virtual void setupPointers()
  {}

  virtual const IdVector & getDepSolnGIDVec()
  {
    return expVarGIDs;
  }

  virtual void setupBreakPoints() {return;}

  virtual bool getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);

  virtual bool updateSource ();

  virtual bool applyScale () { return true; }
  virtual bool processParams () { return true; }
  virtual bool processInstanceParams () { return true; }

  virtual bool updateTemperature(const double & temp_tmp);

  virtual bool isConverged();

  virtual bool testDAEMatrices(const std::vector<const std::string *> &nameVec);

  virtual bool loadTrivialDAE_FMatrixStamp ();
  bool trivialStampLoader (Linear::Matrix * matPtr);

  virtual bool updateIntermediateVars() = 0;
  virtual bool updatePrimaryState() = 0;
  virtual bool updateSecondaryState ();
  virtual bool setIC ();

  // This indicates if the device has functions that can output plot files for
  // internal variables.
  virtual bool plotfileFlag () {return false;}

  // load zeros into mask for equations that should not be used
  // to compute error estimates.  Return true if any zeros set.
  // Default implementation just does nothing (leaves everything 1.0)
  virtual void loadErrorWeightMask() {}

  // tell device instance that current solution has been accepted at
  // current time.  Most devices don't care, but the transmission line
  // does.
  virtual void acceptStep() {}

  // new DAE functions:
  virtual bool loadDAEQVector ()=0;
  virtual bool loadDAEFVector ()=0;
  virtual bool loadDAEBVector (){return true;}

  virtual bool loadDAEdQdx ()=0;
  virtual bool loadDAEdFdx ()=0;

  virtual int getNumNoiseSources () const
  {
    return 0;
  }

  virtual void setupNoiseSources (Xyce::Analysis::NoiseData & noiseDataVec)
  {
    return;
  }

  virtual void getNoiseSources (Xyce::Analysis::NoiseData & noiseDataVec)
  {
    return;
  }

  const InstanceName &getName() const
  {
    return name_;
  }

  int getNumIntVars() const 
  {
    return numIntVars;
  }

  int getNumExtVars() const 
  {
    return numExtVars;
  }

  int getNumStateVars() const 
  {
    return numStateVars;
  }

  int getNumStoreVars() const 
  {
    return numStoreVars;
  }

  int getNumBranchDataVars() const
  {
    return numBranchDataVars;
  }
  
  void setNumStoreVars(int num_store_vars) 
  {
    numStoreVars = num_store_vars;
  }

  void setNumBranchDataVars(int num_branch_data_vars) 
  {
    numBranchDataVars = num_branch_data_vars;
  }


  virtual const std::vector<int> & getDevConMap();

  virtual DeviceState * getInternalState();
  virtual bool setInternalState( const DeviceState & state );

  virtual bool loadDFDV(int iElectrode, Linear::Vector * dfdvPtr);
  virtual bool calcConductance (int iElectrode, const Linear::Vector * dxdvPtr);

  /// Populates and returns the store name map.
  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const = 0;

  virtual bool outputPlotFiles(bool force_final_output) {return true;}

  // two level newton and PDE-continuation
  virtual bool enablePDEContinuation(int &max_PDE_continuation_steps);
  virtual bool disablePDEContinuation();
  virtual void setPDEContinuationAlpha (double alpha);
  virtual void setPDEContinuationBeta  (double beta );

  virtual bool setInitialGuess ();
  virtual double getMaxTimeStepSize  ();
  virtual bool maxTimeStepSupported () {return false;};
  virtual bool getFastSourceFlag() const {return false;};

  virtual void varTypes( std::vector<char> & varTypeVec ) {}

  bool getNumericalSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );


  bool getNumericalMatrixSensitivity ( const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs );

  bool getNumericalBSensVectorsforAC ( const std::string & paramName,
    std::vector< std::complex<double> > & dbdp,
    std::vector<int> &        BindicesVec);

protected:
  void jacStampMap(
    const JacobianStamp &       stamp_parent,
    IdVector &                  map_parent,
    JacobianStamp &             map2_parent,
    JacobianStamp &             stamp,
    IdVector &                  map,
    JacobianStamp &             map2,
    int                         from,
    int                         to,
    int                         original_size);

  void jacStampMap_fixOrder(
    const JacobianStamp &       stamp_parent,
    JacobianStamp &             map2_parent,
    JacobianStamp &             stamp,
    JacobianStamp &             map2);

  void computeJacStampAndMaps(
     const PairVector & jacobianElements,
     const PairVector & collapsedNodes,
     JacobianStamp & jacStamp,
     IdVector & nodeMap,
     PairMap & pairToJacStampMap,
     const int numberOfNodes);
  
  void outputJacStamp(const JacobianStamp & jac);
  void outputJacMaps(const std::vector<int>  & jacMap, const JacobianStamp & jacMap2);

public:
  bool getOrigFlag() const 
  {
    return origFlag;
  }

  void setOrigFlag(bool origFlag_local) 
  {
    origFlag = origFlag_local;
  }

  void consolidateDevLIDs() 
  {
    if (devLIDs.empty())
    {
      devLIDs = extLIDVec;
      devLIDs.insert(devLIDs.end(), intLIDVec.begin(), intLIDVec.end());
      devLIDs.insert(devLIDs.end(), expVarLIDs.begin(), expVarLIDs.end());
    }
  }

  const IdVector &getDevLIDs() const 
  {
    return devLIDs;
  }

  const std::vector<IdVector > &getDevJacLIDs() const 
  {
    return devJacLIDs;
  }

  const IdVector &getStaLIDVec() const 
  {
    return staLIDVec;
  }

  bool getMergeRowColChecked() const 
  {
    return mergeRowColChecked;
  }
  void setMergeRowColChecked(bool mergeRowColChecked_local) 
  {
    mergeRowColChecked = mergeRowColChecked_local;
  }

  const MatrixLoadData &getMatrixLoadData() const 
  {
    return mlData;
  }

  MatrixLoadData &getMatrixLoadData() 
  {
    return mlData;
  }

  const ExternData &getExternData() const {
    return extData;
  }

private:
  InstanceName          name_;
  MatrixLoadData &      mlData;

protected:
  const ExternData &    extData;

  IdVector              intLIDVec;
  IdVector              extLIDVec;
  IdVector              staLIDVec;
  IdVector              stoLIDVec;
  IdVector              devLIDs;                ///< devLIDs is a combined LID vector, containing int, ext, and expVar ID's.
  JacobianStamp         devJacLIDs;

  // device support class: (limiter functions, etc.)
  DeviceSupport devSupport;

private:
  bool configuredForLeadCurrent;  

public:
  std::vector<int> & cols;
  std::vector<double> & vals;

  NumericalJacobian * numJacPtr;

  bool origFlag;

  int numIntVars;
  int numExtVars;
  int numStateVars;
  int numStoreVars;

  bool loadLeadCurrent;           // flag indicating that we want to load lead current data during F & Q load

  int numBranchDataVars;            // number of spaces to allocate in lead current and junction voltage arrays for lead current and power calculations
                                    // this is initially set to zero but the next variable is added to this if the parser find a lead current or power
                                    // calculation called for by the simulation.
  int numBranchDataVarsIfAllocated; // number of spaces to allocate in lead current and junction voltage arrays for lead current and power calculations

  std::vector<int> devConMap;

  bool mergeRowColChecked;
};

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
inline bool DeviceInstance::updateSecondaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/22/03
//-----------------------------------------------------------------------------
inline bool DeviceInstance::setIC ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getInternalState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
inline DeviceState * DeviceInstance::getInternalState()
{
  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDevConMap
// Purpose       : Get connectivity map for leads.  Zero means a lead is
//                 connected to ground.  Other values indicate subsets of
//                 leads that have connection to each other.  Example would
//                 be a mosfet which would have 1 for gate and 2 for drain
//                 and source and zero for bulk, assuming bulk is grounded
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/20/05
//-----------------------------------------------------------------------------
inline const std::vector<int> & DeviceInstance::getDevConMap()
{
  return devConMap;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::loadDFDV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool DeviceInstance::loadDFDV(int iElectrode, Linear::Vector * dfdvPtr)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::calcConductance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool DeviceInstance::calcConductance (int iElectrode, const Linear::Vector * dxdvPtr)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance:isConverged ()
// Purpose       : Return whether a device has done something that should
//                  be interpreted as invalidating other convergence tests
//                  (i.e. that means this step should not be considered
//                   converged even if norms are good)
//                 Since origFlag is set to true by the DeviceInstance
//                 constructor, this is a suitable base class method for
//                 almost all devices.  Devices with more complex convergence
//                 issues can override.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
inline bool DeviceInstance::isConverged()
{
  return origFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getInstanceBreakPoints
// Purpose       : virtual function for obtaining breakpoints from a device.
//
// Special Notes : No-op for the base class version.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
inline bool DeviceInstance::getInstanceBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes)
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateSource
// Purpose       : virtual function for obtaining breakpoints from a device.
//
// Special Notes : No-op for the base class version.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
inline bool DeviceInstance::updateSource ()
{
  return true;
}

inline void addInternalNode(Util::SymbolTable &symbol_table, int index, const InstanceName &instance_name, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::SOLUTION_SYMBOL, index, spiceInternalName(instance_name, lead_name));
}

inline void addInternalNode(Util::SymbolTable &symbol_table, int index, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::SOLUTION_SYMBOL, index, lead_name);
}

inline void addStoreNode(Util::SymbolTable &symbol_table, int index, const InstanceName &instance_name, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::STORE_SYMBOL, index, spiceStoreName(instance_name, lead_name));
}

inline void addStoreNode(Util::SymbolTable &symbol_table, int index, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::STORE_SYMBOL, index, lead_name);
}

inline void addStateNode(Util::SymbolTable &symbol_table, int index, const InstanceName &instance_name, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::STATE_SYMBOL, index, spiceInternalName(instance_name, lead_name));
}

inline void addBranchDataNode(Util::SymbolTable &symbol_table, int index, const InstanceName &instance_name, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::BRANCH_SYMBOL, index, spiceStoreName(instance_name, lead_name));
}

inline void addBranchDataNode(Util::SymbolTable &symbol_table, int index, const std::string &lead_name) {
  Util::addSymbol(symbol_table, Util::BRANCH_SYMBOL, index, lead_name);
}

} // namespace Device
} // namespace Xyce

#endif

