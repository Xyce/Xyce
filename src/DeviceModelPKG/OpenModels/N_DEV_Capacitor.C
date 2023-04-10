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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Capacitor.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Math.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {
namespace Capacitor {

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Capacitor::Traits::loadInstanceParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the instance constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : David Baur
// Creation Date : 1/27/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the instance parameter map.
///
/// @param p     instance parameter map
///
/// @see Xyce::Device::Resistor::Traits::loadInstanceParameters
///
void Traits::loadInstanceParameters(ParametricData<Capacitor::Instance> &p)
{
  p.addPar("C", 1.e-6, &Capacitor::Instance::C)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_FARAD)
    .setDescription("Capacitance")
    .setAnalyticSensitivityAvailable(true)
    .setSensitivityFunctor(&capSens)
    .setAnalyticMatrixSensitivityAvailable(true)
    .setMatrixSensitivityFunctor(&capMatrixSens);

  p.addPar("Q", 0.0, &Capacitor::Instance::Q)
    .setExpressionAccess(ParameterType::SOLN_DEP)
    .setUnit(U_COULOMB)
    .setDescription("Charge");
  p.addPar("M", 1.0, &Capacitor::Instance::multiplicityFactor)
    .setUnit(U_NONE)
    .setDescription("Multiplicity Factor");
  p.addPar("IC", 0.0, &Capacitor::Instance::IC)
    .setGivenMember(&Capacitor::Instance::ICGiven)
    .setUnit(STANDARD);
  p.addPar("L", 1.0, &Capacitor::Instance::length)
    .setUnit(U_METER)
    .setDescription("Semiconductor capacitor width");
  p.addPar("W", 1.e-6, &Capacitor::Instance::width)
    .setUnit(U_METER)
    .setDescription("Semiconductor capacitor length");
  p.addPar("AGE", 0.0, &Capacitor::Instance::age)
    .setUnit(U_HOUR)
    .setDescription("Age of capacitor");
  p.addPar("D", 0.0233, &Capacitor::Instance::ageCoef)
    .setDescription("Age degradation coefficient");
  p.addPar("TEMP", 0.0, &Capacitor::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(STANDARD)
    .setDescription("Device temperature");

  p.addPar("TC1", 0.0, &Capacitor::Instance::tempCoeff1)
    .setGivenMember(&Capacitor::Instance::tempCoeff1Given)
    .setUnit(U_DEGCM1)
    .setDescription("Linear Temperature Coefficient");
  p.addPar("TC2", 0.0, &Capacitor::Instance::tempCoeff2)
    .setGivenMember(&Capacitor::Instance::tempCoeff2Given)
    .setUnit(U_DEGCM2)
    .setDescription("Quadratic Temperature Coefficient");
  p.makeVector("TC", 2);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Capacitor::Traits::loadModelParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the model constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : David Baur
// Creation Date : 1/27/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the model parameter map.
///
/// @param p     model parameter map
///
/// @see Xyce::Device::Resistor::Traits::loadInstanceParameters
///
void Traits::loadModelParameters(ParametricData<Capacitor::Model> &p)
{
  p.addPar("C", 1.0, &Capacitor::Model::capacitanceMultiplier)
    .setUnit(U_NONE)
    .setDescription("Capacitance multiplier");
  p.addPar("CJ", 0.0, &Capacitor::Model::cj)
    .setUnit(U_FARADMM2)
    .setDescription("Junction bottom capacitance");
  p.addPar("CJSW", 0.0, &Capacitor::Model::cjsw)
    .setUnit(U_FARADMM1)
    .setDescription("Junction sidewall capacitance");
  p.addPar("DEFW", 1.e-6, &Capacitor::Model::defWidth)
    .setUnit(U_METER)
    .setDescription("Default device width");
  p.addPar("NARROW", 0.0, &Capacitor::Model::narrow)
    .setUnit(U_METER)
    .setDescription("Narrowing due to side etching");
  p.addPar("TC1", 0.0, &Capacitor::Model::tempCoeff1)
    .setUnit(STANDARD);
  p.addPar("TC2", 0.0, &Capacitor::Model::tempCoeff2)
    .setUnit(STANDARD);
  p.addPar("TNOM", 0.0, &Capacitor::Model::tnom)
    .setGivenMember(&Capacitor::Model::tnomGiven)
    .setUnit(STANDARD);
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
/// Process parameters.
///
/// @return true on success
///
/// In general, the processParams method is intended as a place for complex
/// manipulation of parameters that must happen if temperature or other
/// external changes happen.  
///
/// The Capacitor device supports an "AGE" parameter and a degradation
/// rate parameter that together determine how to modify the
/// capacitance given on the instance line.  Further, Xyce supports a
/// "semiconductor capacitor" model which allows the user to specify
/// the capacitance through a combination of model parameters (junction
/// capacitance and junction sidewall capacitance) and instance parameters
/// (length and width).
///
/// Both of these methods of capacitance value determination need to
/// be done after the normal processing of netlist parameters.  That
/// processing is done here.
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   6/03/02
bool Instance::processParams ()
{
  //  Set any non-constant parameter defaults:

  if (!given("W"))
    width = model_.defWidth;
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();

  if (!tempCoeff1Given)
    tempCoeff1=model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2=model_.tempCoeff2;

  
  baseCap = C;

  // test for various error conditions on the instance line
  if (!given("C") && !given("L") && !given("Q"))
  {
    UserError(*this) << "Could find neither C, Q, nor L parameters in instance.";
  }

  // Now we know we have either capacitance or length (or both) specified.
  // Must specify C if AGE is specified.
  if (!given("C") && given("AGE"))
  {
    UserError(*this) << "Age (A) defined, but no C instance parameter given. Can't use age with semiconductor capacitor options.";
  }

  // the age aware capacitor simply modifies the base capacitance.
  if (given("AGE") && age >= 1)
  {
    ageFactor = (1-ageCoef*log10(age));
    baseCap *= ageFactor;
  }
 
  // If both C and L are specified then C will be used.  If only L is specified
  // then it will be used.
  if (!given("C") && !given("Q"))
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      dout() << "Semiconductor capacitor " <<  getName() 
        //<< Util::push << std::endl
        << std::endl
             << "cj = " << model_.cj << std::endl
             << "cjsw = " << model_.cjsw << std::endl
             << "width = " << width << std::endl
             << "length = " << length << std::endl
             << "narrow = " << model_.narrow << 
             //Util::pop << std::endl;
             std::endl;
    }

    baseCap = C =
              model_.cj*(length-model_.narrow)*(width-model_.narrow) +
              2*model_.cjsw*(length+width-2*model_.narrow);
  }
  
  // M must be non-negative
  if (multiplicityFactor <= 0)
  {
    UserError(*this) << "Multiplicity Factor (M) must be non-negative" << std::endl;
  }

  // If there are any time dependent parameters, set their values for
  // the current time.

  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
///
/// Update the parameters that depend on the temperature of the device
///
/// @param temp_tmp temperature
///
/// Xyce has a number of mechanisms that allow temperature to be changed
/// after a device has been instantiated.  These include .STEP loops over
/// temperature.  When temperature is changed, any device that has parameters
/// that depend on temperature must be updated.  That updating happens here.
///
/// The capacitor device  supports temperature-dependent resistance through its
/// TC1 (linear dependence) and TC2 (quadratic dependence) parameters.
/// If these parameters are specified, the capacitance must be updated.
///
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  //double difference, factor;
  double difference;

  difference = temp - model_.tnom;
  temperatureFactor = model_.capacitanceMultiplier*(1.0 + tempCoeff1*difference +  
                                        tempCoeff2*difference*difference);
  C = baseCap*temperatureFactor;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    dout() << "Capacitor " << getName() << " updateTemperature()" 
      //<< Util::push << std::endl
      << std::endl
           << "C = " << C << std::endl
           << "temp = " << temp << std::endl
           << "temp_tmp = " << temp_tmp << std::endl
           << "tnom = " << model_.tnom << std::endl
           << "difference = " << difference << std::endl
           << "tempCoeff1 = " << tempCoeff1 << std::endl
           << "tempCoeff2 = " << tempCoeff2 << std::endl
           << "baseCap = " << baseCap << std::endl
           << "temperatureFactor = " << temperatureFactor  << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Construct a Capacitor model from a "model block" that was created
/// by the netlist parser.
///
/// @param configuration
/// @param model_block
/// @param factory_block
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   3/16/00
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    expPtr(0),
    expNumVars(0),
    multiplicityFactor(1.0),
    IC(0),
    temperatureFactor(1.0),
    ageFactor(1.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    baseCap(0.0),
    tempGiven(0),
    ICGiven(false),
    solVarDepC(false),
    solVarDepQ(false),
    UIC_(false),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_branch_data(0),
    li_QState(-1),
    li_vcapState(-1),
    li_capState(-1),
    APosEquPosNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquBraNodeOffset(-1),
    APosEquBraNodeOffset(-1),
    ANegEquBraNodeOffset(-1)

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  , qPosEquPosNodePtr(0),
    qNegEquPosNodePtr(0),
    qPosEquNegNodePtr(0),
    qNegEquNegNodePtr(0),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fBraEquBraNodePtr(0),
    fPosEquBraNodePtr(0),
    fNegEquBraNodePtr(0)
#endif

{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 1;
  setNumStoreVars(0);
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  devConMap[0] = 1;
  devConMap[1] = 2;

  /// Unlike the resistor, the capacitor's jacobian stamp is set up directly
  /// in the constructor, and is not static.  This is because the capacitor
  /// supports some options that the resistor does not:
  ///
  ///  - The capacitor instance line may be given an IC=value that
  ///    will be used as the initial voltage drop across the
  ///    capacitance at the operating point.
  ///
  ///  - The capacitance on the instance line may be an expression
  ///    that is permitted to be a function of other solution variables (such 
  ///    as voltages elsewhere in the circuit).
  ///
  ///  - Instead of specifying the capacitance, the user may specify the
  ///    expression for charge directly.
  ///
  ///  Both of these require that the Jacobian stamp for the device be modified.
  ///  Use of a static Jacobian stamp would prevent this flexibility, because
  ///  the static stamp would be used by all capacitor devices, even those
  ///  that do not make use of the options.
  ///
  /// Since the setting of the Jacobian stamp in this otherwise simple device
  /// is so complex, we will document how it is set here.
  ///
  /// In its simplest form, the charge \f$q_0\f$ on the capacitor is
  /// \f$C(V_+-V_-)\f$.  Thus, the current out of the positive node is
  /// \f$d(q_0)/dt\f$, and the current "out" of the negative node is
  /// \f$-d(q0)/dt\f$.  In the Xyce formulation, we load \f$q_0\f$ into
  /// the Q vector for the positive node, and \f$-q_0\f$ into the Q
  /// vector for the negative node; the time integator will later
  /// differentiate this to obtain the currents.  Thus, the contribution
  /// to the Jacobian from the capacitor (with constant capacitance and
  /// no initial condition) will require loading the dQdx matrix with the
  /// derivatives of \f$q_0\f$ with respect to the voltage nodes:
  /// \f[
  /// \left[\begin{array}{rr}
  /// C & -C \\
  /// -C & C  
  /// \end{array} \right] \f]
  ///
  /// The jacobian stamp in this case is the same as the one defined by
  /// the resistor: it has two rows, one for the positive node equation
  /// and one for negative node equation, and two columns, one for the
  /// positive node and one for the negative.  Column 0's value in each
  /// row is 0 to reflect that the first nonzero value of the jacobian
  /// row is the one corresponding to the positive node, and column 1's value
  /// is 1 to reflect that the second nonzero corresponds to the dependence  on
  /// the negative node.
  /// 
  /// If an initial condition is present, however, the circuit at DC is
  /// the same as if there were only a voltage source across the
  /// positive and negative nodes, and in transient it is the same as
  /// the capacitor without the voltage source present.  Thus, at DC the
  /// current out of the positive node would be equal to the voltage
  /// source branch current, and the current out of the negative node
  /// would be the negative of that.  Since these quantities are not
  /// differentiated they would be placed in the F vector, no the Q
  /// vector.  An extra equation, called the "branch equation"
  /// stipulates that the voltage drop between the positive and negative
  /// node be equal to the initial condition, so this element of 
  /// the F vector would be loaded with \f$(V_+-V_-)-V_{ic}\f$.
  /// Therefore in DC the dFdx matrix would be loaded with:
  /// \f[
  /// \left[ \begin{array}{rrr}
  /// 0 & 0 & 1 \\
  /// 0 & 0 & -1 \\
  /// 1 & -1 & 0
  /// \end{array}\right]  \f]
  /// The the first two rows are for the positive and negative node
  /// equations, and the third row is for the branch equation.  The
  /// third column is for the branch current variable. At DC nothing
  /// would be loaded into either the Q vector or its derivative.
  ///
  /// In transient, the loads into the Q vector are the same as
  /// without the initial condition.  Xyce requires that a single jacobian
  /// stamp be used for both the dFdx and dQdx matrices, and does not allow
  /// this matrix to vary between DC and transient.  Thus the dQdx matrix would
  /// be loaded with:
  /// \f[
  /// \left[ \begin{array}{rrr}
  /// C  & -C & 0 \\
  /// -C &  C & 0 \\
  /// 0  &  0 & 0
  /// \end{array}\right]  \f]
  /// The dFdx matrix must be loaded with a single value to turn off the
  /// branch equation and prevent a singular Jacobian:
  /// \f[
  /// \left[ \begin{array}{rrr}
  /// 0 & 0 & 0 \\
  /// 0 & 0 & 0 \\
  /// 0 & 0 & 1
  /// \end{array}\right]  \f]
  /// The net result of this modification is that now, irrespective of whether
  /// we are at DC or in transient, every element of the 3x3 matrix is potentially
  /// nonzero in one or the other of dFdx or dQdx, and therefore our jacobian
  /// stamp is also a dense 3x3 matrix, with each column's value equal to that
  /// column's number.
  ///
  /// Finally, if the capacitance is solution-variable dependent, each of the
  /// rows for positive and negative nodes must be augmented with an additional
  /// column for each variable that the capacitance depends on.  These rows
  /// are similarly dense, and each value of the jacobian stamp for each column
  /// is equal to its column number.
  if( jacStamp.empty() )
  {
    jacStamp_IC.resize(3);
    jacStamp_IC[0].resize(3);
    jacStamp_IC[1].resize(3);
    jacStamp_IC[2].resize(3);
    jacStamp_IC[0][0] = 0;
    jacStamp_IC[0][1] = 1;
    jacStamp_IC[0][2] = 2;
    jacStamp_IC[1][0] = 0;
    jacStamp_IC[1][1] = 1;
    jacStamp_IC[1][2] = 2;
    jacStamp_IC[2][0] = 0;
    jacStamp_IC[2][1] = 1;
    jacStamp_IC[2][2] = 2;

    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);


  // Handle case where capacitance is solution-variable dependent:
  if (getDependentParams().size()>0)
  {
    std::vector<Depend>::const_iterator d;
    std::vector<Depend>::const_iterator begin=getDependentParams().begin();
    std::vector<Depend>::const_iterator end=getDependentParams().end();

    for (d=begin; d!=end; ++d)
    {
      if (d->expr->getNumDdt() != 0)
      {
        UserError(*this) << "Dependent expression " << d->expr->get_expression() << " for parameter " << d->name << " contains time derivatives";
      }

      if (d->numVars > 0)
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          dout() << "Capacitor " << getName()
                 << ": Found solution-dependent parameter " << d->name << " depending on " << expNumVars << " variables" << std::endl;
          if (d->expr->isTimeDependent())
          {
            dout() << "     " << "Expression is time-dependent."  << std::endl;
          }
          else
          {
            dout() << "     " << "Expression is not time-dependent."  << std::endl;
          }
        }

        if (d->name == "C")
        {
          expNumVars = d->numVars;
          solVarDepC = true;
          // To do the proper integration of the charge, we need to save the
          // voltage drop, the old capacitance and
          // the derivatives of Q and C from the last step.
          numStateVars += 2+2*expNumVars;
          expPtr = d->expr;
          dependentParamExcludeMap_[d->name] = 1;
        }

        if (d->name == "Q")
        {
          expNumVars = d->numVars;
          solVarDepQ = true;
          expPtr = d->expr;
          dependentParamExcludeMap_[d->name] = 1;
        }

        if (solVarDepC || solVarDepQ)
        {
          // We now need to extend the pos and neg rows of the jacstamps
          // to account for the additional dependencies:
          jacStamp[0].resize(2+expNumVars);
          jacStamp[1].resize(2+expNumVars);
          jacStamp_IC[0].resize(3+expNumVars);
          jacStamp_IC[1].resize(3+expNumVars);
          for (int i=0; i<expNumVars; ++i)
          {
            jacStamp[0][2+i]=2+i;
            jacStamp[1][2+i]=2+i;

            jacStamp_IC[0][3+i]=3+i;
            jacStamp_IC[1][3+i]=3+i;
          }

          // finally, allocate space to hold the derivatives of C or Q w.r.t.
          // the variables it depends on:
          expVarDerivs.resize(expNumVars);
          // and LIDs for state vector if doing C, but not Q
          if (solVarDepC)
          {
            li_dQdXState.resize(expNumVars);
            li_dCdXState.resize(expNumVars);
          }
        }
        else
        {
          UserError(*this) << d->name << " cannot depend on solution variables. This is only allowed for the C and Q parameters" ;
        }
      }
    }
  }

  if (solVarDepQ && solVarDepC)
    UserError(*this) << "Both C and Q have been specified as expression parameters.  Only one may be specified at a time. ";

  if (solVarDepQ && ICGiven)
    UserError(*this) << "Q has been specified as an expression parameter and an IC given.  IC with Q specified is not implemented.";

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params:

  processParams ();

  // we're gonna have to fake a voltage source at the operating point
  if (ICGiven ) numIntVars = 1;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations


//-----------------------------------------------------------------------------
// Function      : Instance::isLinearDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Instance::isLinearDevice() const
{
  if( loadLeadCurrent )
  {
    return false;
  }

  const std::vector<Depend> & depVec = const_cast<Xyce::Device::Capacitor::Instance*>(this)->getDependentParams();
  if ( depVec.size() )
  {
    std::vector<Depend>::const_iterator d;
    std::vector<Depend>::const_iterator begin=depVec.begin();
    std::vector<Depend>::const_iterator end=depVec.end();

    for (d=begin; d!=end; ++d)
    {
      int expNumVars = d->numVars;
      int expNumGlobal = d->numGlobals;
      Util::Expression* expPtr = d->expr;

      if (expNumVars > 0 || expPtr->isTimeDependent() || expNumGlobal > 0 )
      {   
        return false;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
///
/// Register local IDs
///
/// Register the local internal and external node IDs.
///
/// @param intLIDVecRef internal local IDs from topology package
/// @param extLIDVecRef external local IDs from topology package
/// 
/// Instantiation (calling the device constructor) of the device
/// sets up variables numIntVars and numExtVars, the numbers of internal and
/// external variables associated with the device.  This information is then
/// used by the Topology package to assign locations in the solution vector
/// (and all other vectors of the same shape) for those variables.
/// The "Local IDs" (LIDs) of these locations are provided by Topology
/// so the device can know where to load its data.
///
/// This method saves the LIDs from Topology and associates each one with
/// a particular local name for the internal or external variable.  They 
/// are then used when we load the F and Q vectors.
///
/// The Capacitor device has no internal variables, so this method makes no use 
/// of the intLIDVecRef array.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   6/20/02
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)

{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.  For the matrix indices, first do the
  // rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  // For fake voltage source at operating point
  if (ICGiven)
  {
    li_Bra = intLIDVec[0];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    dout() << "Capacitor " << getName() << " Instance::registerLIDs" 
     // << Util::push << std::endl
      << std::endl
           << "li_Pos_ = " << li_Pos << std::endl
           << "li_Neg_ = " << li_Neg << std::endl;

    if (ICGiven)
      dout() << "li_Bra = "<< li_Bra<< std::endl;

    //dout() << Util::pop << std::endl;
    dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
///
/// Register the local state IDs
///
/// @param staLIDVecRef State variable local IDs
///
/// In general, devices may declare at construction that they require storage
/// locations in the "state vector."  Topology assigns locations in the 
/// state vector and returns "Local IDs" (LIDs) for devices to use for their
/// state vector entries.  If a device uses state vector entries, it
/// uses the registerStateLIDs method to save the local IDs for later use.
/// 
/// The capacitor has at least one state variable (the charge) and as many
/// as three (the charge plus state variables used to support the "voltage
/// dependent capacitance" feature) plus the number of variables on whihc
/// the capacitance depends.
///
/// @note The storage of the charge as a state variable when the capacitance
/// is a constant is a holdover from older implementations, and in the future
/// may be saved only as part of the voltage-dependent capacitance feature.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   06/20/02
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  int i=0;

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  li_QState = staLIDVec[i++];

  // If the capacitance is voltage dependent, we have additional state  vars
  // (does not apply to voltage-dependent Q expressions
  if (solVarDepC)
  {
    li_vcapState = staLIDVec[i++];
    li_capState = staLIDVec[i++];

    for (int j = 0; j<expNumVars; ++j)
    {
      li_dQdXState[j] = staLIDVec[i++];
    }

    for (int j = 0; j<expNumVars; ++j)
    {
      li_dCdXState[j] = staLIDVec[i++];
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 01/23/2013
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "store" vector.  As with other such vectors, the device
/// declares at construction time how many store vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   1/23/2013
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Capacitor::Instance::registerBranchDataLIDs
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
/// The Capacitor device uses exactly one "branch data vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(C1)" for the current through C1. and a junction voltage.
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
  if (ICGiven)
    addInternalNode(symbol_table, li_Bra, getName(), "branch");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
///
/// Return Jacobian stamp that informs topology of the layout of the
/// capacitor jacobian.
///
/// @return const reference to a std::vector of std::vector of
/// integers describing Jacobian stamp shape
//
/// The Jacobian stamp describes the shape of the Jacobian to the
/// Topology subsystem.  The Topology subsystem, in turn, returns
/// the offsets into the matrix and solution vectors where this
/// instance data is located.
///
/// The Jacobian stamp of the capacitor depends on whether an initial
/// condition is given or not.  
///
/// @author Robert Hoekstra
/// @date 8/20/2001
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (ICGiven)
    return jacStamp_IC;
  else
    return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
///
/// Register the Jacobian local IDs
///
/// @param jacLIDVec Jacobian local Ids
///
/// @see Xyce::Device::Capacitor::Instance::Capacitor
///
/// Having established local IDs for the solution variables, Topology must
/// also assign local IDs for the elements of the Jacobian matrix.
///
/// For each non-zero element that was identified in the jacobianStamp,
/// Topology will assign a Jacobian LID.  The jacLIDVec will have the 
/// same structure as the JacStamp, but the values in the jacLIDVec will
/// be offsets into the row of the sparse Jacobian matrix corresponding
/// to the non-zero identified in the stamp.
/// 
/// These offsets are stored in this method for use later when we load
/// the Jacobian.
///
/// @note Because the capacitor's Jacobian stamp depends on whether an
/// initial condition was given or not, this method is slightly more
/// complex than the corresponding method in the Resistor.  The method
/// is further complicated by the possibilty that the capacitance may
/// be given as an expression depending on solution variables; in this case,
/// the capacitor jacstamp has extra columns for the device's dependence on
/// those other variables. 
///
/// @see Xyce::Device::Resistor::Instance::registerJacLIDs
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   08/27/02
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];

  if (ICGiven)
  {
    APosEquBraNodeOffset = jacLIDVec[0][2];
    ANegEquBraNodeOffset = jacLIDVec[1][2];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
    ABraEquBraNodeOffset = jacLIDVec[2][2];
  }
  // set additional offsets if we have a solution-variable dependent C or Q
  if (solVarDepC || solVarDepQ)
  {
    int depVarsBaseIndex = 2;
    if (ICGiven)
    {
      depVarsBaseIndex=3;
    }

    APosEquDepVarOffsets.resize(expNumVars);
    ANegEquDepVarOffsets.resize(expNumVars);

    for ( int i=0; i<expNumVars; ++i)
    {
      APosEquDepVarOffsets[i] = jacLIDVec[0][depVarsBaseIndex+i];
      ANegEquDepVarOffsets[i] = jacLIDVec[1][depVarsBaseIndex+i];
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    dout() << "Capacitor " << getName() << " Instance::registerJacLIDs" 
      //<< Util::push << std::endl
      << std::endl
           << "APosEquPosNodeOffset: " << APosEquPosNodeOffset << std::endl
           << "APosEquNegNodeOffset: " << APosEquNegNodeOffset << std::endl
           << "ANegEquPosNodeOffset: " << ANegEquPosNodeOffset << std::endl
           << "ANegEquNegNodeOffset: " << ANegEquNegNodeOffset << std::endl
           << "APosEquBraNodeOffset: " << APosEquBraNodeOffset << std::endl
           << "ANegEquBraNodeOffset: " << ANegEquBraNodeOffset << std::endl
           << "ABraEquPosNodeOffset: " << ABraEquPosNodeOffset << std::endl
           << "ABraEquNegNodeOffset: " << ABraEquNegNodeOffset << std::endl
           << "ABraEquBraNodeOffset: " << ABraEquBraNodeOffset << std::endl;
    if (solVarDepC || solVarDepQ)
    {
      for ( int i=0; i<expNumVars; ++i)
      {
        dout() << "APosEquDepVarOffsets["<<i<<"]: " << APosEquDepVarOffsets[i] << std::endl
               << "ANegEquDepVarOffsets["<<i<<"]: " << ANegEquDepVarOffsets[i] << std::endl;
      }
    }
    //dout() << Util::pop << std::endl;
    dout() << std::endl;
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
///
/// Setup direct access pointer to solution matrix and vectors.
///
/// @see Xyce::Device::Capacitor::Instance::registerJacLIDs
///
/// As an alternative to the row offsets defined in registerJacLIDs, it 
/// is also possible to obtain direct pointers of the Jacobian elements.
///
/// This method uses the offsets obtained in registerJacLIDs to retrieve
/// the pointers.
///
/// In the capacitor device the pointers to the matrix are only saved
/// (and are only used for matrix access) if
/// Xyce_NONPOINTER_MATRIX_LOAD is NOT defined at compile time.  For
/// some devices the use of pointers instead of array indexing can be
/// a performance enhancement.
///
/// Use of pointers in this device is disabled by defining
/// Xyce_NONPOINTER_MATRIX_LOAD at compile time.  When disabled, array
/// indexing with the offsets from registerJacLIDs is used in
/// the matrix load methods.
///
/// @author Eric Keiter, SNL
/// @date   11/30/08
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  qPosEquPosNodePtr = &(dQdx[li_Pos][APosEquPosNodeOffset]);
  qPosEquNegNodePtr = &(dQdx[li_Pos][APosEquNegNodeOffset]);
  qNegEquPosNodePtr = &(dQdx[li_Neg][ANegEquPosNodeOffset]);
  qNegEquNegNodePtr = &(dQdx[li_Neg][ANegEquNegNodeOffset]);

  if (solVarDepC || solVarDepQ)
  {
    qPosEquDepVarsPtrs.resize(expNumVars);
    qNegEquDepVarsPtrs.resize(expNumVars);

    for (int i=0; i<expNumVars; ++i)
    {
      qPosEquDepVarsPtrs[i]=&(dQdx[li_Pos][APosEquDepVarOffsets[i]]);
      qNegEquDepVarsPtrs[i]=&(dQdx[li_Neg][ANegEquDepVarOffsets[i]]);
    }
  }

  // there are no contributions to the dFdx matrix from dependent C's, so
  // we don't bother with those pointers.

  if (ICGiven)
  {
    fPosEquBraNodePtr = &(dFdx[li_Pos][APosEquBraNodeOffset]);
    fNegEquBraNodePtr = &(dFdx[li_Neg][ANegEquBraNodeOffset]);
    fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
    fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
    fBraEquBraNodePtr = &(dFdx[li_Bra][ABraEquBraNodeOffset]);
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
///
/// Update the state variables.
///
/// @return true on success
///
/// The capacitor's state variables are used to store the charge on
/// the capacitor.  In the case of a constant capacitance the charge
/// is \f$q_0=CV\f$, but the computation is much more complex if the
/// capacitance is variable.
///
/// @note This method is called by the default implementation of the
/// loadState master function. While the Capacitor class reimplements
/// the "Master" "loadState" function that loads the contributions
/// from all capacitor devices in a single loop, because this method
/// is so complicated for the case of solution-dependent capacitances,
/// the overloaded method falls back on this function in that case.
///
/// @note Even though this method IS called in some cases, note that it does
/// NOT call updateIntermediateVars as other device updatePrimaryState methods
/// do.  The capacitor device doesn't actually *HAVE* an updateIntermediateVars
/// method, and all the hard work is done either in the Master class
/// updateState function directly for the simple, constant-capacitance
/// case, or here, for the solution-variable-dependent-capacitance case.
/// 
/// @see Xyce::Device::Capacitor::Master::updateState
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/29/01
///
bool Instance::updatePrimaryState ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * staVec = extData.nextStaVectorRawPtr;
  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    dout() << " ----------------------------------" << std::endl;
    dout() << "Instance::updatePrimaryState:" << std::endl;
  }

  vcap = v_pos-v_neg;

  q0_Jdxp = 0.0;

  if( getSolverState().dcopFlag && ICGiven ) vcap = IC;
  if (ICGiven && UIC_)
  {
    q0_Jdxp = C*(IC-vcap);
    vcap = IC;
    UIC_ = false;
  }

  if (!solVarDepC && !solVarDepQ)
  {
    // Obtain the "current"  value for the charge stored in the capacitor.
    q0 = C*vcap;
    //staVec[li_QState] = q0;
  }
  else if (solVarDepQ)
  {
    expPtr->evaluate(Q,expVarDerivs);
    q0 = Q;
  }
  else
  {
    // The capacitance depends on solution variables, the work load just
    // went up.
    expPtr->evaluate(C,expVarDerivs);

    // Redo the age and temperature modifications, if necessary,
    // and apply to both C and derivatives.
    // This uses factors previously computed in updateTemperature and processParams.
    baseCap = C*ageFactor; 
    C *= ageFactor*temperatureFactor;
    for (int ii=0;ii<expNumVars; ++ii) { expVarDerivs[ii] *= ageFactor*temperatureFactor; }

    // At DC, the charge really still is V*C.
    if (getSolverState().dcopFlag)
    {
      q0 = vcap*C;
      // dQ/dX is just dC/dX*vcap
      for (int i=0;i<expNumVars; ++i)
      {
        expVarDerivs[i] *= vcap;
      }
    }
    else
    {
      ///
      /// For solution-variable dependent cap in transient, we can't
      /// use the expression \f$q_0=CV\f$ because the capacitance is
      /// actually dQ/dV, not Q/V.  We must integrate CdV to get the
      /// charge.  We approximate this by incrementally adding
      /// C'*deltaV as V changes.  C' is the average capacitance
      /// between this and the previous step.  Using the average
      /// assures charge conservation, at least when C is a function
      /// of vcap alone.
      ///
      double * oldstaVector = extData.currStaVectorRawPtr;
      double oldC;
      double oldVcap;
      q0=oldstaVector[li_QState];
      oldC=oldstaVector[li_capState];
      oldVcap=oldstaVector[li_vcapState];

      q0 += 0.5*(oldC+C)*(vcap-oldVcap);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        dout() << "  Derivatives of C w.r.t variables: " << std::endl;
        for (int i=0;i<expNumVars;++i)
        {
          dout() << " expVarDerivs[ "<< i << " ] = " << expVarDerivs[i] << std::endl;
        }
      }
      ///
      /// When C is not a function of vcap alone, we have additional derivative
      /// terms that must also be integrated for proper computation of
      /// all dQdx entries.
      ///
      // If expressions contain the "ddt" (time differentiation)
      // function, everything becomes harder.  So for now, we're going
      // to disallow use of ddt in C, so we don't have to deal with
      // that added complexity.  This restriction is enforced in the
      // constructor, where we throw a fatal error if the user gives
      // us a ddt-dependent C.  This code is NOT sufficient if ddt is allowed
      // in capacitance expressions.
      //
      // Now we have some trickiness because we are integrating CdV to get Q.
      // In order to get dQ/dX we have two cases:
      //   X is one of our capacitor's voltage nodes:  dQ/dX = C or -C depending
      //       on whether X is the pos or negative node
      //  X is NOT one of our nodes: dQ/dX = integral( dC/dX *dV)

      // For now, because we don't have an easy way to tell which of our
      // expression nodes is which, we'll just calculate all the dC/dX and dQ/dX
      // the same way.  When it comes time to assemble the jacobian, we'll use
      // the node/equation offsets to know when to skip adding in this component.

      // The logic here is similar to the computation of the charge itself.
      // We'll use "expVarDerivs" to hold the final dQ/dX values.


      // Need to save the dC/dX value for next step.
      for (int i=0;i<expNumVars; ++i)
      {
        staVec[li_dCdXState[i]] = expVarDerivs[i];
      }

      // we have to integrate

      // dQ/dx = olddQdX + .5*(olddCdX+newdCdX)*(vcap-oldvcap)

      for (int i=0; i< expNumVars; ++i)
      {
        expVarDerivs[i] = oldstaVector[li_dQdXState[i]]
          + 0.5*(oldstaVector[li_dCdXState[i]]+expVarDerivs[i])*
          (vcap-oldstaVector[li_vcapState]);

      }
    }
    // Regardless of whether it's dcop or not, expVarDerivs now contains
    // dQ/dX for all the X's.  It's WRONG if X is one of our voltage nodes,
    // so we have to be careful not to use it in that case.  This logic
    // is handled down in loadDAEdQdx
    // Save to state:
    for (int i=0;i<expNumVars; ++i)
    {
      staVec[li_dQdXState[i]] = expVarDerivs[i];
    }



    staVec[li_QState] = q0;
    staVec[li_vcapState]=vcap;
    staVec[li_capState]=C;
  }
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
///
/// Load the DAE Q vector.
///
/// @return true on success
///
/// The Xyce DAE formulation solves the differential-algebraic
/// equations \f$F(x)+dQ(x)/dt=0\f$ These are vector equations
/// resulting from the modified nodal analysis of the circuit.
/// 
/// This method loads the Q-vector contributions for a single capacitor
/// instance.
///
/// In this method, the offsets defined in registerLIDs are used to
/// store the device's Q contributions into the Q vector.
///
/// The Q vector is used for devices that store charge or magnetic
/// energy.  The capacitor is such a device
///
/// @note This method is called by the default implementation of the
/// loadDAEVectors master function. Since the capacitor class
/// reimplements the "Master" "loadDAEVectors" function that loads the
/// contributions from all capacitor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEQVector method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the Capacitor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Capacitor::Master::loadDAEVectors
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/24/03
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  qVec[li_Pos] += (q0 * multiplicityFactor);
  qVec[li_Neg] += (-q0 * multiplicityFactor);
  if( loadLeadCurrent )
  {
    double * leadQ = extData.nextLeadCurrQCompRawPtr;
    leadQ[li_branch_data] = q0 * multiplicityFactor;
  }
  if (q0_Jdxp != 0)
  {
    double *dQdxdVp = extData.dQdxdVpVectorRawPtr;
    dQdxdVp[li_Pos] += q0_Jdxp*multiplicityFactor;
    dQdxdVp[li_Neg] -= q0_Jdxp*multiplicityFactor;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 For the capacitor this doesn't do anything, except in
//                 the case of IC= being specified.  In that case, then
//                 some extra stuff is contributed that doesn't have time
//                 derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
///
/// Load the DAE F vector.
///
/// @return true on success
///
/// The Xyce DAE formulation solves the differential-algebraic
/// equations \f$F(x)+dQ(x)/dt=0\f$ These are vector equations
/// resulting from the modified nodal analysis of the circuit.
/// 
/// This method loads the F-vector contributions for a single capacitor
/// instance.
///
/// In this method, the offsets defined in registerLIDs are used to
/// store the device's F contributions into the F vector.
///
/// The only time a capacitor adds anything to the F vector is in the 
/// DC phase of a computation if and only if an initial condition is given
/// on the capacitor instance line. 
///
/// @note This method is called by the default implementation of the
/// loadDAEVectors master function. Since the Capacitor class
/// reimplements the "Master" "loadDAEVectors" function that loads the
/// contributions from all capacitor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEFVector method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the Capacitor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Capacitor::Master::loadDAEVectors
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/24/03
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  double Vpos = 0.0;
  double Vneg = 0.0;
  double v_tmp = 0.0;
  double * fVec = extData.daeFVectorRawPtr;
  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    double * solVec = extData.nextSolVectorRawPtr;

    if ( ICGiven && getSolverState().dcopFlag )
    {
      leadF[li_branch_data] = solVec[li_Bra];
    }
    else
    {
      leadF[li_branch_data] = 0;
    }

    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }


  if (ICGiven && getSolverState().dcopFlag)
  {
    // If we're doing the operating point and we have an initial condition,
    // get the current from the branch equation
    Vpos = (*extData.nextSolVectorPtr)[li_Pos];
    Vneg = (*extData.nextSolVectorPtr)[li_Neg];

    // load current into the F vector
    fVec[li_Pos] += (*extData.nextSolVectorPtr)[li_Bra];
    fVec[li_Neg] += -(*extData.nextSolVectorPtr)[li_Bra];
  }

  // Initial condition stuff.
  v_tmp=0;
  if (ICGiven && getSolverState().dcopFlag)
  {
    v_tmp= (Vpos-Vneg-IC);
  }

  // Do this whenever there's a Branch equation, but only if there is one.
  // We'll be using 0 if we're not the OP.
  if (ICGiven)
  {
    fVec[li_Bra] += v_tmp;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 capacitor instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
///
/// Load the DAE the derivative of the Q vector with respect to the
/// solution vector x, dFdx
///
/// Loads the contributions for a single capacitor instance to the 
/// dFdx matrix (the Q contribution to the Jacobian).
///
/// This method uses the Jacobian LIDs (row offsets) that were stored by
/// registerJacLIDs.
///
/// @see Xyce::Device::Capacitor::Instance::registerJacLIDs
///
/// @note This method is called by the default implementation of the
/// loadDAEMatrices master function. Since the Capacitor class
/// reimplements the "Master" "loadDAEMatrices" function that loads the
/// contributions from all capacitor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the capacitor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Capacitor::Master::loadDAEMatrices
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   03/05/04
bool Instance::loadDAEdQdx ()
{
  if (!(ICGiven&& getSolverState().dcopFlag))
  {
    Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
    if (!solVarDepQ)
    {
      dQdx[li_Pos][APosEquPosNodeOffset] += (C * multiplicityFactor);
      dQdx[li_Pos][APosEquNegNodeOffset] -= (C * multiplicityFactor);
      dQdx[li_Neg][ANegEquPosNodeOffset] -= (C * multiplicityFactor);
      dQdx[li_Neg][ANegEquNegNodeOffset] += (C * multiplicityFactor);


      // Remember the comments in updatePrimaryState: expVarDerivs
      // contains dQ/dX, but they are only correct when X is not one of
      // our nodal voltages.  If X *IS* one of our nodal voltages, dQ/dX
      // is either C or -C and is already handled above.  We need only
      // do the stuff below for the dependencies on voltages that are
      // NOT our nodal voltages.
      if (solVarDepC)
      {
        for (int i=0; i< expNumVars; ++i)
        {
          if ( (APosEquDepVarOffsets[i] != APosEquPosNodeOffset)
               && (APosEquDepVarOffsets[i] != APosEquNegNodeOffset))
          {
            dQdx[li_Pos][APosEquDepVarOffsets[i]] += (expVarDerivs[i] * multiplicityFactor);
          }
          if ( (ANegEquDepVarOffsets[i] != ANegEquPosNodeOffset)
               && (ANegEquDepVarOffsets[i] != ANegEquNegNodeOffset))
          {
            dQdx[li_Neg][ANegEquDepVarOffsets[i]] -= (expVarDerivs[i] * multiplicityFactor);
          }
        }
      }
    }
    else  // solVarDepQ
    {
      // For the voltage-dependent Q, all of the correct derivatives are
      // in expVarDerivs already, just copy them out.
      for (int i=0; i< expNumVars; ++i)
      {
        dQdx[li_Pos][APosEquDepVarOffsets[i]] += (expVarDerivs[i] * multiplicityFactor);
        dQdx[li_Neg][ANegEquDepVarOffsets[i]] -= (expVarDerivs[i] * multiplicityFactor);
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 For the capacitor this doesn't do anything, unless IC=
//                 has been specified for an initial condition.  Then,
//                 there are extra equations that do not contain time
//                 derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
///
/// Load the DAE the derivative of the F vector with respect to the
/// solution vector x, dFdx
///
/// Loads the contributions for a single capacitor instance to the 
/// dFdx matrix (the F contribution to the Jacobian).
///
/// This method uses the Jacobian LIDs (row offsets) that were stored by
/// registerJacLIDs.
///
/// @see Xyce::Device::Capacitor::Instance::registerJacLIDs
///
/// The capacitor only loads the dFdx matrix when an initial condition is
/// given on the instance line for the device.
///
/// @note This method is called by the default implementation of the
/// loadDAEMatrices master function. Since the Capacitor class
/// reimplements the "Master" "loadDAEMatrices" function that loads the
/// contributions from all capacitor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the Capacitor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Capacitor::Master::loadDAEMatrices
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   03/05/04
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (ICGiven && getSolverState().dcopFlag)
  {
    // Special Jacobian if we're doing operating point when IC given
    dFdx[li_Pos][APosEquBraNodeOffset] +=  1.0;
    dFdx[li_Neg][ANegEquBraNodeOffset] += -1.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] +=  1.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += -1.0;
  }
  else
  {
    if (ICGiven)
      dFdx[li_Bra][ABraEquBraNodeOffset] += 1.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       : Make sure initial conditions are applied with NOOP or UIC
//                 transient runs.
// Special Notes : This function is only called when Xyce is given
//                 a UIC or NOOP keyword on the tran line, once, at the very
//                 start of transient.
//                 All it does now is set a boolean to tell the computational
//                 bits of the capacitor code that they have something special
//                 to do.  It replaces an earlier version that used to try
//                 to write values into the solution vector.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 09/06/18
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  if (ICGiven)
    UIC_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  if (ICGiven)
  {
    varTypeVec.resize(1);
    varTypeVec[0] = 'I';
  }
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
///
/// Process model parameters
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   6/03/02
bool Model::processParams ()
{

  if (!tnomGiven)
    tnom = getDeviceOptions().tnom;

  // If there are any time dependent parameters, set their values for
  // the current time.

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
///
/// Process the instance parameters of instance owned by this model
///
/// This method simply loops over all instances associated with this
/// model and calls their processParams method.
///
/// @return true
///
/// @author Dave Shirely, PSSI
/// @date   03/23/06
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/17/00
//-----------------------------------------------------------------------------
///
/// Construct a capacitor model from a "model block" that was created
/// by the netlist parser.
///
/// @param configuration
/// @param model_block
/// @param factory_block
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   5/17/00
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
    capacitanceMultiplier(1.0),
    cj(0.0),
    cjsw(0.0),
    defWidth(10e-6),
    narrow(0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    tnomGiven(0)
{

  // Set params to constant default values :
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams(model_block.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

Model::~Model()
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
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------

std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of capacitor instances: " << isize << std::endl;
  os << "    name\t\tmodelName\tParameters" << std::endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << "\t\tC = " << (*iter)->C;
    os << "\tIC = " << (*iter)->IC;
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


// Capacitor Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Update state for all capacitor instances, regardless of model.
///
/// @param solVec solution vector
/// @param staVec state vector
/// @param stoVec store vector
///
/// @return true on success
///
/// @note Because the capacitor device re-implements the base-class
/// Master::updateState, the Instance::updatePrimaryState method is never
/// called, nor is the Instance::updateIntermediateVars method.  This method
/// replaces those, and does the same work but inside a loop over all
/// capacitor instances.
///
/// Because the computation of state variables is so complex in the
/// event that the capacitance is given by an expression that depends
/// on solution variables, this method falls back on calling the
/// instance's updatePrimaryState method instead of reimplementing the
/// computation here.
///
/// @see Xyce::Device::Capacitor::Instance::updatePrimaryState
/// @author Eric Keiter, SNL
/// @date   11/26/08
bool Master::updateState (double * solVec, double * staVec, double * stoVec, int loadType)
{
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);

    double v_pos = solVec[ci.li_Pos];
    double v_neg = solVec[ci.li_Neg];
    ci.vcap = v_pos-v_neg;

    ci.q0_Jdxp = 0.0;

    if( getSolverState().dcopFlag && ci.ICGiven )
    {
      ci.vcap = ci.IC;
    }
    if (ci.ICGiven && ci.UIC_)
    {
      ci.q0_Jdxp = ci.C*(ci.IC-ci.vcap);
      ci.vcap = ci.IC;
      ci.UIC_ = false;
    }

    if (!ci.solVarDepC && !ci.solVarDepQ)
    {
      // Obtain the "current"  value for the charge stored in the capacitor.
      ci.q0 = ci.C * ci.vcap;
      //staVec[ci.li_QState] = ci.q0;
    }
    else
    {
      // fall back on old pre-turbo scheme
      bool tmpBool=true;
      tmpBool = ci.updatePrimaryState ();
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
///
/// Load DAE vectors of all capacitor instances, regardless of model
///
/// @param solVec solution vector
/// @param fVec f vector
/// @param qVec q vector
/// @param leadF store lead current f vector
/// @param leadQ store lead current q vector
///
/// @return true on success
///
/// @note Because the capacitor device re-implements the base-class
/// Master::loadDAEVectors, the Instance::loadDAEFVector method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all capacitor instances.
///
/// @see Xyce::Device::Capacitor::Instance::loadDAEFVector
///
/// @author Eric Keiter, SNL
/// @date   11/26/08
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, 
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    dout() << " ----------------------------------" << std::endl;
    dout() << " Master::loadDAEVectors: " << std::endl;
  }

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  InstanceVector::const_iterator it, end;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);
    if (ci.ICGiven)
    {
      double Vpos (0.0), Vneg (0.0), v_tmp (0.0);

      // Initial condition
      if (getSolverState().dcopFlag)
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          dout() << " loading dcop F vector for cap " << ci.getName() << ":" << std::endl;
        }
        // If doing the DCOP and have IC=, get current from branch equation
        // ci.i0   = solVec[ci.li_Bra]; moved to CapacitorMaster::updateState() where stovec is passed in.
        Vpos    = solVec[ci.li_Pos];
        Vneg    = solVec[ci.li_Neg];
        fVec [ci.li_Pos] += solVec[ci.li_Bra];
        fVec [ci.li_Neg] += -solVec[ci.li_Bra];

        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          dout() << " f[ " << ci.li_Pos << " ] += " << solVec[ci.li_Bra]<< std::endl;
          dout() << " f[ " << ci.li_Neg << " ] += " << -solVec[ci.li_Bra] << std::endl;
        }

        v_tmp= (Vpos-Vneg-ci.IC);

        // Do this only if there's a Branch equation
        fVec[ci.li_Bra] += v_tmp;
      }
      else
      {
        solVec[ci.li_Bra] = 0.0;
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        dout() << " f[ " << ci.li_Bra << " ] += " << v_tmp << std::endl;
      }
    }

    qVec[ci.li_Pos] += (ci.q0 * ci.multiplicityFactor);
    qVec[ci.li_Neg] += (-ci.q0 * ci.multiplicityFactor);

    if (ci.q0_Jdxp != 0.0)
    {
      double *dQdxdVp = ci.extData.dQdxdVpVectorRawPtr;
      dQdxdVp[ci.li_Pos] += ci.q0_Jdxp*ci.multiplicityFactor;
      dQdxdVp[ci.li_Neg] -= ci.q0_Jdxp*ci.multiplicityFactor;
    }

    if( ci.loadLeadCurrent )
    {
      // If an IC is given then include the current through the voltage
      // source used to enforce that IC at the DCOP.  Otherwise, leadF
      // is zero.
      if ( (ci.ICGiven) && (getSolverState().dcopFlag) )
      {
        leadF[ci.li_branch_data] = solVec[ci.li_Bra];
      }
      else
      {
        leadF[ci.li_branch_data] = 0;
      }
      leadQ[ci.li_branch_data] = (ci.q0 * ci.multiplicityFactor);
      junctionV[ci.li_branch_data] = solVec[ci.li_Pos] - solVec[ci.li_Neg];
 
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      dout() << " loading Q vector for cap " << ci.getName() << ":" << std::endl;
      dout() << " q[ " << ci.li_Pos << " ] += " << ci.q0 << std::endl;
      dout() << " q[ " << ci.li_Neg <<  " ] += " << -ci.q0 << std::endl;

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
///
/// Load DAE matrices for all capacitor instances, regardless of model
///
/// @param dFdx matrix of derivatives of F vector with respect to solution
/// @param dQdx matrix of derivatives of Q vector with respect to solution
///
/// @return true on success
///
/// @note Because the capacitor device re-implements the base-class
/// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all capacitor instances.
///
/// @see Xyce::Device::Capacitor::Instance::loadDAEdFdx
///
/// @author Eric Keiter, SNL
/// @date   11/26/08
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    dout() << " ----------------------------------" << std::endl;
    dout() << " Master::loadDAEMatrices: " << std::endl;
  }

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  InstanceVector::const_iterator it, end;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      dout() << " loads for capacitor " << ci.getName() << std::endl;
    }

    if (ci.ICGiven && getSolverState().dcopFlag)
    {
      // Special Jacobian if we're doing operating point when IC given
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *(ci.fPosEquBraNodePtr) +=  1.0;
      *(ci.fNegEquBraNodePtr) += -1.0;
      *(ci.fBraEquPosNodePtr) +=  1.0;
      *(ci.fBraEquNegNodePtr) += -1.0;
#else
      dFdx[ci.li_Pos][ci.APosEquBraNodeOffset] +=  1.0;
      dFdx[ci.li_Neg][ci.ANegEquBraNodeOffset] += -1.0;
      dFdx[ci.li_Bra][ci.ABraEquPosNodeOffset] +=  1.0;
      dFdx[ci.li_Bra][ci.ABraEquNegNodeOffset] += -1.0;
#endif
    }
    else
    {
      if (ci.ICGiven)
      {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
        *(ci.fBraEquBraNodePtr) += 1.0;
#else
        dFdx[ci.li_Bra][ci.ABraEquBraNodeOffset] += 1.0;
#endif
      }
    }

    if (!(ci.ICGiven&& getSolverState().dcopFlag))
    {
      if (!ci.solVarDepQ)
      {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
        *(ci.qPosEquPosNodePtr) += (ci.C * ci.multiplicityFactor);
        *(ci.qPosEquNegNodePtr) -= (ci.C * ci.multiplicityFactor);
        *(ci.qNegEquPosNodePtr) -= (ci.C * ci.multiplicityFactor);
        *(ci.qNegEquNegNodePtr) += (ci.C * ci.multiplicityFactor);

        if (ci.solVarDepC)
        {
          for (int i=0; i< ci.expNumVars; ++i)
          {
            // Similar to logic in loadDAEdQdX:
            if ((ci.qPosEquDepVarsPtrs[i] != ci.qPosEquPosNodePtr)
                && (ci.qPosEquDepVarsPtrs[i] != ci.qPosEquNegNodePtr))
            {
              *(ci.qPosEquDepVarsPtrs[i]) += (ci.expVarDerivs[i] * ci.multiplicityFactor);
            }

            if ((ci.qNegEquDepVarsPtrs[i] != ci.qNegEquPosNodePtr)
                && (ci.qNegEquDepVarsPtrs[i] != ci.qNegEquNegNodePtr))
            {
              *(ci.qNegEquDepVarsPtrs[i]) -= (ci.expVarDerivs[i] * ci.multiplicityFactor);
            }
            if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
            {
              if ((ci.qPosEquDepVarsPtrs[i] != ci.qPosEquPosNodePtr)
                  && (ci.qPosEquDepVarsPtrs[i] != ci.qPosEquNegNodePtr))
              {
                dout() << " q[pos][ " << ci.APosEquDepVarOffsets[i] << " ] += "
                       << (ci.expVarDerivs[i]  * ci.multiplicityFactor) << std::endl;
              }
              if ((ci.qNegEquDepVarsPtrs[i] != ci.qNegEquPosNodePtr)
                  && (ci.qNegEquDepVarsPtrs[i] != ci.qNegEquNegNodePtr))
              {
                dout() << " q[neg][ " << ci.ANegEquDepVarOffsets[i] << " ] += "
                       << (ci.expVarDerivs[i]  * ci.multiplicityFactor) << std::endl;
              }
            }
          }
        }
#else
	dQdx[ci.li_Pos][ci.APosEquPosNodeOffset] += (ci.C * ci.multiplicityFactor);
	dQdx[ci.li_Pos][ci.APosEquNegNodeOffset] -= (ci.C * ci.multiplicityFactor);
	dQdx[ci.li_Neg][ci.ANegEquPosNodeOffset] -= (ci.C * ci.multiplicityFactor);
	dQdx[ci.li_Neg][ci.ANegEquNegNodeOffset] += (ci.C * ci.multiplicityFactor);

        if (ci.solVarDepC)
        {
          for (int i=0; i< ci.expNumVars; ++i)
          {
            if ( (ci.APosEquDepVarOffsets[i] != ci.APosEquPosNodeOffset)
                 && (ci.APosEquDepVarOffsets[i] != ci.APosEquNegNodeOffset))
            {
              dQdx[ci.li_Pos][ci.APosEquDepVarOffsets[i]] +=
                (ci.expVarDerivs[i] * ci.multiplicityFactor);
            }
            if ( (ci.ANegEquDepVarOffsets[i] != ci.ANegEquPosNodeOffset)
                 && (ci.ANegEquDepVarOffsets[i] != ci.ANegEquNegNodeOffset))
            {
              dQdx[ci.li_Neg][ci.ANegEquDepVarOffsets[i]] -=
                (ci.expVarDerivs[i] * ci.multiplicityFactor);
            }

            if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
            {
              if ( (ci.APosEquDepVarOffsets[i] != ci.APosEquPosNodeOffset)
                   && (ci.APosEquDepVarOffsets[i] != ci.APosEquNegNodeOffset))
              {
                dout() << " q[pos][ " << ci.APosEquDepVarOffsets[i] << " ] += "
                       << (ci.expVarDerivs[i] * ci.multiplicityFactor)<< std::endl;
              }
              if ( (ci.ANegEquDepVarOffsets[i] != ci.ANegEquPosNodeOffset)
                   && (ci.ANegEquDepVarOffsets[i] != ci.ANegEquNegNodeOffset))
              {
                dout() << " q[neg][ " << ci.ANegEquDepVarOffsets[i] << " ] += "
                       << (ci.expVarDerivs[i] * ci.multiplicityFactor) << std::endl;
              }
            }
          }
        }
#endif
      }
      else // solVarDepQ
      {
        // Just as in loadDAEdQdX, we just copy expVarDerivs into dQdX if
        // we are doing an expression-based Q.
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
          for (int i=0; i< ci.expNumVars; ++i)
          {
            *(ci.qPosEquDepVarsPtrs[i]) +=
              (ci.expVarDerivs[i] * ci.multiplicityFactor);
            *(ci.qNegEquDepVarsPtrs[i]) -=
              (ci.expVarDerivs[i] * ci.multiplicityFactor);
          }
#else
          for (int i=0; i< ci.expNumVars; ++i)
          {
              dQdx[ci.li_Pos][ci.APosEquDepVarOffsets[i]] +=
                (ci.expVarDerivs[i] * ci.multiplicityFactor);
              dQdx[ci.li_Neg][ci.ANegEquDepVarOffsets[i]] -=
                (ci.expVarDerivs[i] * ci.multiplicityFactor);
          }
#endif
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Capacitor::Traits::factory
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// Create a new instance of the Capacitor device.
///
/// @param configuration
/// @param factory_block
///
Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Capacitor::registerDevice
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// Define how to use the device in a netlist.
///
/// This method is called from the Xyce::Device::registerOpenDevices
/// function, which in turn is called by the device manager.
///
/// The device is declared here to be an "C" device, which may optionally
/// take a model card of type "C".  This device will correspond to model
/// level 1 of capacitor models.
void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() || (deviceMap.find("C")!=deviceMap.end())))
  {
    initialized = true;
    Config<Traits>::addConfiguration()
      .registerDevice("c", 1)
      .registerModelType("c", 1)
      .registerModelType("cap", 1);
  }
}

//-----------------------------------------------------------------------------
// Function      : capSensitivity::operator()
// Purpose       : produces df/dp and dq/dp, where p=C.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/31/2014
//-----------------------------------------------------------------------------
void capSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp,
    std::vector<double> & dqdp,
    std::vector<double> & dbdp,
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  double * solVec = in->extData.nextSolVectorRawPtr;
  double v_pos = solVec[in->li_Pos];
  double v_neg = solVec[in->li_Neg];
  double vcap = v_pos-v_neg;

  double dqdpLoc = vcap;

  dqdp.resize(2);
  dqdp[0] = +dqdpLoc;
  dqdp[1] = -dqdpLoc;

  Qindices.resize(2);
  Qindices[0] = in->li_Pos;
  Qindices[1] = in->li_Neg;
}


//-----------------------------------------------------------------------------
// Function      : capMatrixSensitivity::operator()
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void capMatrixSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  d_dqdx_dp.clear();
  d_dqdx_dp.resize(2);
  d_dqdx_dp[0].resize(2);
  d_dqdx_dp[1].resize(2);
  d_dqdx_dp[0][0] = +(in->multiplicityFactor);
  d_dqdx_dp[0][1] = -(in->multiplicityFactor);
  d_dqdx_dp[1][0] = -(in->multiplicityFactor);
  d_dqdx_dp[1][1] = +(in->multiplicityFactor);

  Q_lids.resize(2);
  Q_lids[0] = in->li_Pos;
  Q_lids[1] = in->li_Neg;

  Q_jacLIDs.clear();
  Q_jacLIDs.resize(2);
  Q_jacLIDs[0].resize(2);
  Q_jacLIDs[1].resize(2);
  Q_jacLIDs[0][0] = in->APosEquPosNodeOffset;
  Q_jacLIDs[0][1] = in->APosEquNegNodeOffset;
  Q_jacLIDs[1][0] = in->ANegEquPosNodeOffset;
  Q_jacLIDs[1][1] = in->ANegEquNegNodeOffset;
}

} // namespace Capacitor
} // namespace Device
} // namespace Xyce
