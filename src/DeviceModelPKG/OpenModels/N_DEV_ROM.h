//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 12/11/09
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ROM_h
#define Xyce_N_DEV_ROM_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace ROM {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "ROM";}
  static const char *deviceTypeName() {return "ROM level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 1000;}
  static const char *instanceDefaultParameter() {return "BASE_FILENAME";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Citer,
     const FactoryBlock &      factory_block);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:

  bool isLinearDevice()
  { 
    if (usePortDesc)
      return false;
    else
      return true;
  }

  // Additional Public Declarations
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp) { return true; }

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  bool setIC ();

  void loadErrorWeightMask ();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // Data Members for Class Attributes
  bool isCSparse;
  bool isGSparse;

  // User-specified parameters:
  bool maskROMVars;
  int usePortDesc;
  int numROMVars;
  std::string baseFileName;   // base file name for reduced-order model files
  std::vector<double> Chat;  // Reduced-order model V'*C*V
  std::vector<int> Chat_colIdx, Chat_rowPtr;  // Chat structures if stored in CSR format
  std::vector<double> Ghat;  // Reduced-order model V'*G*V
  std::vector<int> Ghat_colIdx, Ghat_rowPtr;  // Ghat structures if stored in CSR format
  std::vector<int> CG_colIdx, CG_rowPtr;  // Union of Chat and Ghat maps stored in CSR format
  std::vector<double> Bhat;  // Reduced-order model V'*B
  std::vector<double> Lhat;  // Reduced-order model L'*V
  std::vector<double> Qhat;  // Workspace Qhat = Chat * xhat
  std::vector<double> Fhat;  // Workspace Fhat = [Iq - Lhat'* xhat; Ghat*xhat - Bhat*up]
  std::vector<double> i_ip;  // Storage for Iq

  // Two-level stamps (BNB)
  std::vector<double> Jstamp; // stamp for Jacobian
  std::vector<double> Fstamp; // stamp for F
  std::vector<double> G2;     // intermediate variable
  std::vector<double> C2;     // intermediate varaible
  std::vector<double> A2;     // intermediate varaible
  std::vector<double> A2last;
  std::vector<double> G2p;    // intermediate varaible
  std::vector<double> Gp2;    // intermediate varaible
  std::vector<double> A2sol;  // intermediate variable
  double dt, dt_last, alph, alph_last, coef, coefLast;
  double currentOrder, usedOrder;
  int lastTimeStepNumber;
  std::vector<int> ipiv_A2;  // for LAPACK math

  //local id's (offsets)
  std::vector<int> li_ROM;  // Interior variables
  std::vector<int> li_state; // Internal state

  // Offsets for Jacobian
  std::vector<int> AEqu_up_NodeOffset;
  std::vector<int> AEqu_ip_NodeOffset;
  std::vector< std::vector<int> > AEqu_NodeOffset;
  std::vector<int> ROMEqu_Lt_NodeOffset;
  std::vector<int> ROMEqu_B_NodeOffset;
  std::vector<int> ROMEqu_GpC_NodeOffset;
  // Offsets for sparse C and C in Jacobian
  std::vector<int> ROMEqu_C_NodeOffset;
  std::vector<int> ROMEqu_G_NodeOffset;

  // Pointers for Jacobian
  std::vector<double *> fEqu_up_NodePtr;
  std::vector<double *> fEqu_ip_NodePtr;
  std::vector<double *> fEqu_un_NodePtr; // BNB

  std::vector<double *> qROMEqu_Chat_VarsPtrs;
  std::vector<double *> fROMEqu_Ghat_VarsPtrs;
  std::vector<double *> fROMEqu_Lhat_VarsPtrs;
  std::vector<double *> fROMEqu_Bhat_VarsPtrs;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Model  : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;

public:

  Model(
     const Configuration &       configuration,
     const ModelBlock &      MB,
     const FactoryBlock &    factory_block);
  ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;
  bool processParams ();
  bool processInstanceParams ();

public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  void setupBaseInstanceContainer()
  {
    std::vector<Instance*>::iterator iter = instanceContainer.begin();
    std::vector<Instance*>::iterator end   = instanceContainer.end();
    for ( ; iter!=end; ++iter)
    {
      Xyce::Device::DeviceModel::baseInstanceContainer.push_back( static_cast<Xyce::Device::DeviceInstance *>(*iter) );
    }
  }

private:
  std::vector<Instance*> instanceContainer;

private:

};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1),
    separateInstances_(false)
  {}

  void printMatrix (std::string vname, double * Matrix, int Nrows, int Ncols); // BNB

  virtual bool updateState (double * solVec, double * staVec, double * stoVec)
  { return updateState( solVec, staVec, stoVec, ALL ); }
  virtual bool updateState (double * solVec, double * staVec, double * stoVec, int loadType);
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV)
  { return loadDAEVectors( solVec, fVec, qVec, bVec, leadF, leadQ, junctionV, ALL ); }
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV, int loadType);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  { return loadDAEMatrices( dFdx, dQdx, ALL ); }
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType);

private:
  bool separateInstances_;
  InstanceVector      linearInstances_;            ///< List of owned linear rom instances
  InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear rom instances
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace ROM
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_ROM_h

