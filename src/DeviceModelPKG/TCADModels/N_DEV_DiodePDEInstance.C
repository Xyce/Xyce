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
// Purpose        : One dimensional PDE device, instance class
//                  implementation.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/06/03
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>
#include <cstdio>

#include <N_DEV_fwd.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DiodePDE.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_RegionData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {
namespace DiodePDE {

namespace {
typedef unsigned int    UINT;
}

// default number of mesh points:
static const int NUM_MESH_POINTS = 11;

// default maximum number of nonzero entries in a matrix row
static const int MAX_COLS_PER_ROW = 40;

void Traits::loadInstanceParameters(ParametricData<DiodePDE::Instance> &p)
{
  p.addPar("BASE.LOC", 0.5e-3, &DiodePDE::Instance::baseLocation)
    .setGivenMember(&DiodePDE::Instance::baseLocationGiven)
    .setUnit(U_CM)
    .setDescription("Location of base contact (necessary if running with three terminals).")
    .setCategory(CAT_GEOMETRY);

  p.addPar("AREA", 1.0, &DiodePDE::Instance::area)
    .setUnit(U_CMM2)
    .setDescription("Cross sectional area of the device.")
    .setCategory(CAT_GEOMETRY);

  // user-specified scaling vars:
  p.addPar("X0", 1.0e-7, &DiodePDE::Instance::x0_user)
    .setUnit(U_CM)
    .setDescription("Length scalar; adjust to mitigate convergence problems."
                    "The model will do all of its scaling automatically, so it is generally not "
                    "necessary to specify it manually.")
    .setCategory(CAT_SCALING);

  p.addPar("C0", 1.0e+15, &DiodePDE::Instance::C0_user)
    .setUnit(U_CMM3)
    .setDescription("Density scalar; adjust to mitigate convergence problems."
                    "The model will do all of its scaling automatically, so it is generally not "
                    "necessary to specify it manually.")
    .setCategory(CAT_SCALING);

  p.addPar("t0", 1.0e-6, &DiodePDE::Instance::t0_user)
    .setUnit(U_SECOND)
    .setDescription("Time scalar; adjust to mitigate convergence problems."
                    "The model will do all of its scaling automatically, so it is generally not "
                    "necessary to specify it manually.")
    .setCategory(CAT_SCALING);

  p.addPar("SCALEDENSITYTOMAXDOPING", true, &DiodePDE::Instance::scaleDensityToMaxDoping_)
    .setUnit(U_LOGIC)
    .setDescription("If set the density will be scaled by a fraction of the maximum doping."
                    "The model will do all of its scaling automatically, so it is generally not "
                    "necessary to specify it manually.")
    .setCategory(CAT_SCALING);

  p.addPar("DENSITYSCALARFRACTION", 1.0e-1, &DiodePDE::Instance::densityScalarFraction_)
    .setUnit(U_LOGIC)
    .setDescription("Fraction of the maximum doping by which density will be scaled."
                    "The model will do all of its scaling automatically, so it is generally not "
                    "necessary to specify it manually.")
    .setCategory(CAT_SCALING);

  p.addPar("NA", 1.0e+15, &DiodePDE::Instance::Na)
    .setUnit(U_CMM3)
    .setDescription("Acceptor doping level")
    .setCategory(CAT_DOPING);

  p.addPar("ND", 1.0e+15, &DiodePDE::Instance::Nd)
    .setUnit(U_CMM3)
    .setDescription("Donor doping level")
    .setCategory(CAT_DOPING);

  p.addPar("WJ", 1.0e-4, &DiodePDE::Instance::WJ)
    .setUnit(U_CM)
    .setDescription("Junction width, if graded junction enabled.")
    .setCategory(CAT_DOPING);

  p.addPar("TEMP", 300.15, &DiodePDE::Instance::Temp)
    .setUnit(STANDARD)
    .setDescription("Temperature");

  p.addPar("ANODE.AREA", 0.0, &DiodePDE::Instance::anodeArea)
    .setUnit(U_CMM2)
    .setDescription("Anode area (used for two-terminal devices)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("CATHODE.AREA", 0.0, &DiodePDE::Instance::cathodeArea)
    .setUnit(U_CMM2)
    .setDescription("Cathode area (used for two-terminal devices)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("EMITTER.AREA", 0.0, &DiodePDE::Instance::emitterArea)
    .setUnit(U_CMM2)
    .setDescription("Emitter area (used for three-terminal (BJT) devices)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("BASE.AREA", 0.0, &DiodePDE::Instance::baseArea)
    .setUnit(U_CMM2)
    .setDescription("Base area (used for three-terminal (BJT) devices)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("COLLECTOR.AREA", 0.0, &DiodePDE::Instance::collectorArea)
    .setUnit(U_CMM2)
    .setDescription("Collector area (used for three-terminal (BJT) devices)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("L", 1.0e-3, &DiodePDE::Instance::length)
    .setGivenMember(&DiodePDE::Instance::lengthGiven)
    .setUnit(U_CM)
    .setDescription("Device width. (Synonym with W parameter)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("W", 1.0e-3, &DiodePDE::Instance::width)
    .setGivenMember(&DiodePDE::Instance::widthGiven)
    .setUnit(U_CM)
    .setDescription("Device width. (Synonym with L parameter)")
    .setCategory(CAT_GEOMETRY);

  p.addPar("OUTPUTINTERVAL", 0.0, &DiodePDE::Instance::outputInterval)
    .setGivenMember(&DiodePDE::Instance::outputIntervalGiven)
    .setUnit(U_SECOND)
    .setDescription("Time interval for tecplot output (if tecplot is enabled).")
    .setCategory(CAT_OUTPUT);

  if (DEBUG_DEVICE)
  {
    p.addPar("ANODEINDEX", 1, &DiodePDE::Instance::anodeIndex_user)
      .setGivenMember(&DiodePDE::Instance::anodeIndex_userGiven);

    p.addPar("CATHODEINDEX", 0, &DiodePDE::Instance::cathodeIndex_user)
      .setGivenMember(&DiodePDE::Instance::cathodeIndex_userGiven);
  }

  // Set up map for non-double precision variables:
  p.addPar("GRADED", false, &DiodePDE::Instance::gradedJunctionFlag)
    .setUnit(U_LOGIC)
    .setDescription("Flag for graded junction vs. abrupt junction. – (1/true=graded, 0/false=abrupt)")
    .setCategory(CAT_DOPING);

  p.addPar("MOBMODEL", std::string("ARORA"), &DiodePDE::Instance::mobModelName)
    .setDescription("Mobility model.");

  p.addPar("FIELDDEP", false, &DiodePDE::Instance::fieldDependentMobility)
    .setGivenMember(&DiodePDE::Instance::fieldDependentMobilityGiven)
    .setUnit(U_LOGIC)
    .setDescription("If true, use field dependent mobility.");

  p.addPar("BULKMATERIAL", std::string("SI"), &DiodePDE::Instance::bulkMaterial)
    .setDescription("Bulk semiconductor material");

  p.addPar("OFFSETOUTPUTVOLTAGE", true, &DiodePDE::Instance::useVoltageOutputOffset_)
    .setUnit(U_LOGIC)
    .setDescription("This is an output parameter that determines the ``zero'' of the potential at output.  If OFFSETOUTPUTVOLTAGE=true (default) it will adjust the voltages at output so that the minimum voltage is zero. If true and also FIRSTELECTRODEOFFSET=true, then the voltage of the first electrode is the zero point.  If OFFSETOUTPUTVOLTAGE=false, the output voltage sets the intrisic Fermi level to zero.  Depending on circumstances each of these may be more or less convenient for plotting.")
    .setCategory(CAT_OUTPUT);

  p.addPar("FIRSTELECTRODEOFFSET", false, &DiodePDE::Instance::offsetWithFirstElectrode_)
    .setUnit(U_LOGIC)
    .setDescription("This is an output parameter.  It is only used if OFFSETOUTPUTVOLTAGE=true. (see description of that paramaeter")
    .setCategory(CAT_OUTPUT);

  //p.addPar("DISPLCUR", false, &DiodePDE::Instance::displCurrentFlag)
    //.setUnit(U_LOGIC)
    //.setDescription("If true, displacement current is computed and output");

  p.addPar("OUTPUTNLPOISSON", false, &DiodePDE::Instance::outputNLPoisson)
    .setUnit(U_LOGIC)
    .setDescription("Flag to determine if the results of the nonlinear Poisson "
                    "calculation is included in the output files.  Normally, this calculation"
                    " is used to initialize a drift-diffusion calculation and isn't of interest.")
    .setCategory(CAT_OUTPUT);

  p.addPar("AUGER", true, &DiodePDE::Instance::includeAugerRecomb)
    .setUnit(U_LOGIC)
    .setDescription("Flag to turn on/off Auger recombination");

  p.addPar("SRH", true, &DiodePDE::Instance::includeSRHRecomb)
    .setUnit(U_LOGIC)
    .setDescription("Flag to turn on/off Shockley-Read-Hall (SRH) recombination.");

  p.addPar("GNUPLOTLEVEL", 1, &DiodePDE::Instance::gnuplotLevel)
    .setDescription("Flag for gnuplot output.\n"
                    "0 - no gnuplot files.\n"
                    "1 - gnuplot files.\n"
                    "gnuplot is an open source plotting program that is usually installed on Linux "
                    "systems. gnuplot files will have the *Gnu.dat suffix, and the prefix will be the"
                    "name of the device instance.")
    .setCategory(CAT_OUTPUT);

  p.addPar("TECPLOTLEVEL", 1, &DiodePDE::Instance::tecplotLevel)
    .setDescription("Setting for Tecplot output:\n"
                    "0 - no Tecplot files\n"
                    "1 - Tecplot files, each output in a separate file. 2 - Tecplot file, each output"
                    "appended to a single file.\n"
                    "Tecplot files will have the .dat suffix, and the prefix will be the name of the device instance")
    .setCategory(CAT_OUTPUT);

  p.addPar("SGPLOTLEVEL", 0, &DiodePDE::Instance::sgplotLevel)
    .setDescription("Flag for sgplot output.\n"
                    "0 - no sgplot files.\n"
                    "1 - sgplot files.\n"
                    "sgplot is a plotting program that comes as part of the SG Framework. sgplot "
                    "files will have the *.res suffix, and the prefix will be the name of the "
                    "device instance")
    .setCategory(CAT_OUTPUT);

  // Doping file params:
  p.addPar("DOPING_FILE", std::string("NOFILE"), &DiodePDE::Instance::dopingFileName)
    .setDescription("File containing doping profile.")
    .setCategory(CAT_DOPING);

  p.addPar("PDOPE_FILE", std::string("NOFILE"), &DiodePDE::Instance::pdopeFileName)
    .setCategory(CAT_DOPING)
    .setDescription("File containing doping profile for P-type dopants.");

  p.addPar("NDOPE_FILE", std::string("NOFILE"), &DiodePDE::Instance::ndopeFileName)
    .setCategory(CAT_DOPING)
    .setDescription("File containing doping profile for N-type dopants.");

  p.addPar("NX", 11, &DiodePDE::Instance::NX)
    .setGivenMember(&DiodePDE::Instance::NXGiven)
    .setDescription("Number of mesh points");

  // Beginning of undocumented parameters section.
  // parameters that should not be included in the guides for various reasons:
  //
  p.addPar("MESHFILE", std::string("internal.msh"), &DiodePDE::Instance::meshFileName);

  // PN diode voltage BC's, if running "uncoupled"
  p.addPar("ANODE.BC", 0.5, &DiodePDE::Instance::anodebc)
    .setUnit(U_VOLT)
    .setDescription("Anode voltage boundary condition.  Only used if device is uncoupled from circuit, and running in diode mode.\n")
    .setCategory(CAT_BOUNDARYCONDITIONS);

  p.addPar("CATHODE.BC", 0.0, &DiodePDE::Instance::cathodebc)
    .setUnit(U_VOLT)
    .setDescription("Cathode voltage boundary condition.  Only used if device is uncoupled from circuit, and running in diode mode.\n")
    .setCategory(CAT_BOUNDARYCONDITIONS);

  // BJT voltage BC's, if running "uncoupled"
  p.addPar("EMITTER.BC", 0.5, &DiodePDE::Instance::emitterbc)
    .setUnit(U_VOLT)
    .setDescription("Emitter voltage boundary condition.  Only used if device is uncoupled from circuit, and running in BJT mode.\n")
    .setCategory(CAT_BOUNDARYCONDITIONS);

  p.addPar("COLLECTOR.BC", 0.0, &DiodePDE::Instance::collectorbc)
    .setUnit(U_VOLT)
    .setDescription("Collector voltage boundary condition.  Only used if device is uncoupled from circuit, and running in BJT mode.\n")
    .setCategory(CAT_BOUNDARYCONDITIONS);

  p.addPar("BASE.BC", 0.0, &DiodePDE::Instance::basebc)
    .setUnit(U_VOLT)
    .setDescription("Base voltage boundary condition.  Only used if device is uncoupled from circuit, and running in BJT mode.\n")
    .setCategory(CAT_BOUNDARYCONDITIONS);

  p.addPar("MASKVARSTIA", false, &DiodePDE::Instance::maskVarsTIAFlag_)
    .setUnit(U_LOGIC)
    .setDescription("If set to true, then some variables are excluded from the time integration error control calculation.");

  p.addPar("VOLTLIM", false, &DiodePDE::Instance::voltLimFlag)
    .setUnit(U_LOGIC)
   .setDescription("Flag to apply voltage limiting.  This is only relevant for an experimental two-level Newton solver.");

  p.addPar("MAXVOLTDELTA", 0.025, &DiodePDE::Instance::maxVoltDelta)
    .setUnit(U_VOLT)
    .setDescription("Maximum voltage change used by two-level Newton algorithm.");

  p.addPar("USEOLDNI", false, &DiodePDE::Instance::useOldNi)
    .setUnit(U_LOGIC)
    .setGivenMember(&DiodePDE::Instance::useOldNiGiven)
    .setDescription("Flag for using old(inaccurate) intrinsic carrier calculation.");

  p.addPar ("FERMIDIRAC", false, &DiodePDE::Instance::fermiDiracFlag)
    .setUnit(U_LOGIC)
    .setDescription("Use Fermi-Dirac statistics.");

  p.addPar("THERMIONICEMISSION", false, &DiodePDE::Instance::thermionicEmissionFlag)
    .setUnit(U_LOGIC);

  p.addPar("TUNNELING", std::string("none"), &DiodePDE::Instance::tunnelingModelName);
  // End of undocumented parameters section.

  p.addComposite("NODE", PDE_1DElectrode::getParametricData(), &DiodePDE::Instance::electrodeMap);
  p.addComposite("DOPINGPROFILES", DopeInfo::getParametricData(), &DiodePDE::Instance::dopeInfoMap);
  p.addComposite("REGION", DopeInfo::getParametricData(), &DiodePDE::Instance::dopeInfoMap);
  p.addComposite("LAYER", MaterialLayer::getParametricData(), &DiodePDE::Instance::materialVec);
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
//
// Purpose       : This function contains much of the initialization for
//                 the Instance class.  Most of this was
//                 originally in the constructor.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  updateTemperature(Temp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       model,
  const FactoryBlock &factory_block)
  : DevicePDEInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(model),
    indicesSetup_(false),
    includeBaseNode_(false),
    useElectrodeSpec_(false),
    maskVarsTIAFlag_(false),
    scaleDensityToMaxDoping_(true),
    densityScalarFraction_(1.0e-1),
    useVoltageOutputOffset_(true),
    offsetWithFirstElectrode_(false),
    VoltageOffset_(0.0),
    useLayerCompositeDoping_(false),
    Emax(0.0),
    VminExp(0.0),
    VmaxExp(0.0),
    diodeCap(0.0),
    Na(1.0e15),
    Nd(1.0e15),
    WJ(1.0e-4),
    XC(0.0),
    XL(0.0),
    XR(0.0),
    NnMax(1.0e15),
    NpMax(1.0e15),
    NnMin(1.0e5),  // approx...
    NpMin(1.0e5),
    NX(NUM_MESH_POINTS),
    LX(NX-1),
    NXGiven(false),
    maxVoltDelta(0.025), // thermal voltage.
    enableContinuationCalled(false),
    useOldNi(false),
    useOldNiGiven(false),
    meshFileName(""),
    dopingFileName("NOFILE"),
    ndopeFileName("NOFILE"),
    pdopeFileName("NOFILE"),
    width(1.0e-3),
    length(1.0e-3),
    widthGiven(false),
    lengthGiven(false),
    area(1.0),
    anodebc(0.0),
    cathodebc(0.0),

    emitterbc(0.0),
    collectorbc(0.0),
    basebc(0.0),

    anodeArea(0.0),
    cathodeArea(0.0),

    emitterArea(0.0),
    collectorArea(0.0),
    baseArea(0.0),

    baseLocation(0.5e-3),
    baseLocationGiven(false),

    gradedJunctionFlag(false),
    calledBeforeUIVB(false),
    callsOTEC(0),
    callsOSG(0),
    displCurrentFlag(false),
    equationSet(0),
    outputInterval(0.0),
    outputIntervalGiven(false),
    outputIndex(0),
    outputNLPoisson(false),
    lastOutputTime(-10.0),
    tecplotLevel(0),
    gnuplotLevel(0),
    sgplotLevel(0),
    voltLimFlag(false),
    includeAugerRecomb(true),
    includeSRHRecomb(true),
    fermiDiracFlag(false),
    thermionicEmissionFlag(false),
    tunnelingModelName("none"),

    anodeIndex_userGiven(false),
    cathodeIndex_user(0),
    cathodeIndex_userGiven(false),
    NUMRC(NX*3),

    maxColsPerRow(MAX_COLS_PER_ROW),
    numElectrodes(2),
    columnReorderingFlag(false),
    layerCompositeSpecified(false)
{
  bcVec.clear();

  // these 4 mesh things change later.
  numIntVars   = 3*NUM_MESH_POINTS;
  numExtVars   = 2;
  if (IB.numExtVars != 0)
  {
    numExtVars   = IB.numExtVars;
  }
  numStateVars = 2;

  if (numExtVars < 3)
  {
    includeBaseNode_ = false;
  }
  else if (numExtVars == 3)
  {
    includeBaseNode_ = true;
  }
  else if (numExtVars > 3)
  {
    UserFatal(*this) << "Too many external nodes are set!  Set no more than 3.";
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // check doping files...
  if ( given("PDOPE_FILE") && !given("NDOPE_FILE") )
  {
    UserFatal(*this) << "Ndope file specified with no Pdope file.  Exiting.";
  }

  if ( !given("PDOPE_FILE") && given("NDOPE_FILE") )
  {
    UserFatal(*this) << "Pdope file specified with no Ndope file.  Exiting.";
  }

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    Temp = getDeviceOptions().temp.getImmutableValue<double>();

  if (given("MESHFILE"))
  {
    UserFatal(*this) << "Mesh file was specified.  The 1D device doesn't need a mesh file."
                     << " Either add a model statement of level=2, or get rid of the mesh"
                     << " file specification.";
  }

  if (lengthGiven && !widthGiven)
  {
    width = length;
  }

  if (given("GNUPLOTLEVEL") && !given("TECPLOTLEVEL"))
  {
    tecplotLevel = gnuplotLevel;
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  processParams ();

  // calculate dependent (ie computed) params and check for errors:
  ExtendedString tmpName = mobModelName;
  tmpName.toLower();
  mobModelName = tmpName;

  bool bsuccess = true;
  bool bs1 = true;

  bs1 = setupDefaultLayer ();  bsuccess = bsuccess && bs1;
  bs1 = setupNumVars ();       bsuccess = bsuccess && bs1;
  bs1 = doAllocations ();      bsuccess = bsuccess && bs1;
  bs1 = setupMesh ();          bsuccess = bsuccess && bs1;
  bs1 = setupMaterialArrays ();bsuccess = bsuccess && bs1;
  bs1 = setupNodes ();         bsuccess = bsuccess && bs1;
  bs1 = setupDopingProfile (); bsuccess = bsuccess && bs1;
  bs1 = setupMiscConstants (); bsuccess = bsuccess && bs1;
  bs1 = setupScalingVars ();   bsuccess = bsuccess && bs1;

  bs1 = setupJacStamp ();      bsuccess = bsuccess && bs1;
  bs1 = cleanupJacStamp ();    bsuccess = bsuccess && bs1;

  if (!given("AREA")) area = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  for (std::map<std::string, DopeInfo *>::iterator it = dopeInfoMap.begin(); 
      it != dopeInfoMap.end(); ++it)
  {
    delete (*it).second;
  }

  for (std::map<std::string, PDE_1DElectrode *>::iterator it = electrodeMap.begin(); 
      it != electrodeMap.end(); ++it)
  {
    delete (*it).second;
  }

  int size = materialVec.size();
  for (int i=0;i<size;++i)
  {
    MaterialLayer *matPtr = materialVec[i];
    if (matPtr != 0)
    {
      delete matPtr;
      materialVec[i] = 0;
    }
  }
  materialVec.clear();
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/14/05
//-----------------------------------------------------------------------------
CompositeParam *
Instance::constructComposite(const std::string & compositeName, const std::string & paramName)
{
  if (compositeName == "DOPINGPROFILES" || compositeName == "REGION")
  {
    DopeInfo *n = new DopeInfo();
    dopeInfoMap[paramName] = n;
    return static_cast<CompositeParam *> (n);
  }
  else if (compositeName == "NODE" || compositeName == "ELECTRODE")
  {
    bcData bc;
    ExtendedString electrodeName = paramName;
    electrodeName.toUpper ();

    bc.eName = electrodeName;
    bc.nName = paramName;
    bc.given = true;
    bc.index = 0;

    if (electrodeName =="ANODE")
    {
      bc.meshIndex = 0;
      bc.neighborNode = 1;
    }
    else
    {
      bc.meshIndex = NUM_MESH_POINTS-1;
      bc.neighborNode = NUM_MESH_POINTS-2;
    }

    if (bc.given) ++numElectrodes;
    if (bc.given) bcVec.push_back(bc);

    PDE_1DElectrode *n = new PDE_1DElectrode();
    electrodeMap[paramName] = n;
    return static_cast<CompositeParam *> (n);
  }
  else if (compositeName == "LAYER")
  {
    layerCompositeSpecified = true;
    MaterialLayer *matPtr = new MaterialLayer();
    materialVec.push_back(matPtr);
    return (static_cast<CompositeParam *> (matPtr));
  }
  else
  {
    DevelFatal(*this).in("Instance::constructComposite")
      << "Unrecognized composite name: "
      <<  compositeName;
  }

  return NULL;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::doAllocations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/29/10
//-----------------------------------------------------------------------------
bool Instance::doAllocations ()
{
  // Set up a bunch of mesh-based arrays:
  dxVec.resize        (NX,0.0);
  xVec.resize         (NX,0.0);
  CVec.resize         (NX,0.0);
  CdonorVec.resize    (NX,0.0);
  CacceptorVec.resize (NX,0.0);
  VVec.resize         (NX,0.0);
  ExVec.resize        (NX,0.0);
  JnxVec.resize       (NX,0.0);
  JpxVec.resize       (NX,0.0);
  RVec.resize         (NX,0.0);
  SVec.resize         (NX,0.0);
  nnVec.resize        (NX,0.0);
  npVec.resize        (NX,0.0);

  NcVec.resize        (NX,0.0);
  NvVec.resize        (NX,0.0);
  EcVec.resize       (NX,0.0);
  EvVec.resize       (NX,0.0);
  EcEffVec.resize    (NX,0.0);
  EvEffVec.resize    (NX,0.0);
  bgnCVec.resize      (NX,0.0);
  bgnVVec.resize      (NX,0.0);
  NiVec.resize        (NX,0.0);
  NiEffVec.resize     (NX,0.0);
  EiVec.resize        (NX,0.0);
  EiEffVec.resize     (NX,0.0);
  EfVec.resize        (NX,0.0);
  EfEffVec.resize     (NX,0.0);
  relPermVec.resize   (NX,0.0);
  bulkMaterialVec.resize(NX);

  dRdpVec.resize      (NX,0.0);
  dRdnVec.resize      (NX,0.0);

  dJndn1Vec.resize   (NX,0.0);
  dJndn2Vec.resize   (NX,0.0);
  dJndV1Vec.resize   (NX,0.0);
  dJndV2Vec.resize   (NX,0.0);
  dJndp1Vec.resize   (NX,0.0);
  dJndp2Vec.resize   (NX,0.0);

  dJpdn1Vec.resize   (NX,0.0);
  dJpdn2Vec.resize   (NX,0.0);
  dJpdV1Vec.resize   (NX,0.0);
  dJpdV2Vec.resize   (NX,0.0);
  dJpdp1Vec.resize   (NX,0.0);
  dJpdp2Vec.resize   (NX,0.0);

  tnVec.resize(NX,0.0);
  tpVec.resize(NX,0.0);
  unE_Vec.resize (NX-1,0.0);
  upE_Vec.resize (NX-1,0.0);

  // indexing arrays, local.  jacStamp is resized elsewhere
  li_Vrowarray.resize(NX,0);
  li_Nrowarray.resize(NX,0);
  li_Prowarray.resize(NX,0);

  // displacement current stuff
  stateDispl.resize(NX,0);
  stateDispl_owned.resize(NX,0);
  displCurrent.resize(NX,0.0);
  li_stateDispl.resize(NX,0);

  // set up the boundary stencil:
  boundarySten.resize(NX,0);
  edgeBoundarySten.resize(NX,0);
  internalBoundarySten.resize(NX,0);
  heterojunctionSten.resize(NX,0);
  matIndex.resize(NX,0);

  // these will always be set.
  edgeBoundarySten[0]=1;
  edgeBoundarySten[LX]=1;
  boundarySten[0]=1;
  boundarySten[LX]=1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNodes
//
// Purpose       : This sets up the bcVec container.  bcVec is a vector of
//                 bcData classes, which contain bondary condition
//                 related data.
//
//                 The key issues for a boundary are:
//                    - determine circuit node
//                    - determine mesh boundary location
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/29/10
//-----------------------------------------------------------------------------
bool Instance::setupNodes ()
{
  // If the user did not use the ELECTRODE/NODE vector-composite
  // specification, then set up a default (implicit) set of electrodes.
  //
  // If 2 terminals, assume a diode, with a specification consistent with the
  // SPICE diode noder order:  D cathode anode
  //
  // If 3 terminals, assume a BJT with the node order being the same as
  // the SPICE Gummel-Poon specification:  Q col bas emit
  //
  if ( bcVec.empty() )
  {
    bcVec.clear();
    bcVec.resize(numExtVars);

    useElectrodeSpec_ = false;

    if (includeBaseNode_)
    {
      // collector:
      int collectorIndex=0;
      bcIndexMap["collector"] = collectorIndex;
      bcVec[collectorIndex].eName = "collector";
      bcVec[collectorIndex].Vequ = collectorbc;
      bcVec[collectorIndex].VequGiven = given("COLLECTOR.BC");
      bcVec[collectorIndex].area = collectorArea;
      bcVec[collectorIndex].areaGiven = given("COLLECTOR.AREA");
      bcVec[collectorIndex].meshIndex = LX;
      bcVec[collectorIndex].neighborNode = LX-1;
      if (!given("COLLECTOR.AREA")) bcVec[collectorIndex].area = area;

      // base:
      int baseIndex=1;
      bcIndexMap["base"] = baseIndex;
      bcVec[baseIndex].eName = "base";
      bool found=false;
      int bIndex=0;
      double minDelta = length;
      for (int i=0;i<NX;++i)
      {
        double deltaX=fabs(baseLocation-xVec[i]);
        if (deltaX < minDelta)
        {
          bIndex=i;
          minDelta=deltaX;
        }
      }

      bcVec[baseIndex].Vequ = basebc;
      bcVec[baseIndex].VequGiven = given("BASE.BC");
      bcVec[baseIndex].area = baseArea;
      bcVec[baseIndex].areaGiven = given("BASE.AREA");
      //bcVec[baseIndex].meshIndex = static_cast<int> (LX/2);
      bcVec[baseIndex].meshIndex = bIndex;
      bcVec[baseIndex].neighborNode = bcVec[baseIndex].meshIndex-1;
      if (!given("BASE.AREA")) bcVec[baseIndex].area = area;

      // emitter:
      int emitterIndex=2;
      bcIndexMap["emitter"] = emitterIndex;
      bcVec[emitterIndex].eName = "emitter";
      bcVec[emitterIndex].Vequ = emitterbc;
      bcVec[emitterIndex].VequGiven = given("EMITTER.BC");
      bcVec[emitterIndex].area = emitterArea;
      bcVec[emitterIndex].areaGiven = given("EMITTER.AREA");
      bcVec[emitterIndex].meshIndex = 0;
      bcVec[emitterIndex].neighborNode = 1;
      if (!given("EMITTER.AREA")) bcVec[emitterIndex].area = area;
    }
    else
    {
      // anode:
      int anodeIndex=1;
      bcIndexMap["anode"] = anodeIndex;
      bcVec[anodeIndex].eName = "anode";
      bcVec[anodeIndex].Vequ = anodebc;
      bcVec[anodeIndex].VequGiven = given("ANODE.BC");
      bcVec[anodeIndex].area = anodeArea;
      bcVec[anodeIndex].areaGiven = given("ANODE.AREA");
      bcVec[anodeIndex].meshIndex = 0;
      bcVec[anodeIndex].neighborNode = 1;
      if (!given("ANODE.AREA")) bcVec[anodeIndex].area = area;

      // cathode:
      int cathodeIndex=0;
      bcIndexMap["cathode"] = cathodeIndex;
      bcVec[cathodeIndex].eName = "cathode";
      bcVec[cathodeIndex].Vequ = cathodebc;
      bcVec[cathodeIndex].VequGiven = given("CATHODE.BC");
      bcVec[cathodeIndex].area = cathodeArea;
      bcVec[cathodeIndex].areaGiven = given("CATHODE.AREA");
      bcVec[cathodeIndex].meshIndex = LX;
      bcVec[cathodeIndex].neighborNode = LX-1;
      if (!given("CATHODE.AREA")) bcVec[cathodeIndex].area = area;
    }
  }
  else  // user did use the ELECTRODE/NODE specification.
  {
    useElectrodeSpec_ = true;

    std::vector<int> tmpMeshSten(NX,0);

    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      PDE_1DElectrode & electrode = *(electrodeMap[bcVec[iBC].nName]);

      bcIndexMap[ bcVec[iBC].eName ] = iBC;

      if (electrode.sideGiven)
      {
        ExtendedString side = electrode.side;
        side.toLower();
        if (side == "left")
        {
          bcVec[iBC].meshIndex = 0;
          bcVec[iBC].neighborNode = 1;
          tmpMeshSten[0] = 1;
        }
        else if (side == "right")
        {
          bcVec[iBC].meshIndex = LX;
          bcVec[iBC].neighborNode = LX-1;
          tmpMeshSten[LX] = 1;
        }
        else if (side == "middle" || side == "mid")
        {

          double location = electrode.location;
          bool found=false;
          int bIndex=0;
          double minDelta = length;
          for (int imesh=0;imesh<NX;++imesh)
          {
            double deltaX=fabs(location-xVec[imesh]);
            if (deltaX < minDelta)
            {
              bIndex=imesh;
              minDelta=deltaX;
            }
          }

          bcVec[iBC].meshIndex = bIndex;
          bcVec[iBC].neighborNode = bIndex-1;
          // assuming current coming from the emitter direciton

          // check to make sure that bIndex isn't already used.
          if (tmpMeshSten[bIndex] == 1)
          {
            DevelFatal(*this).in("Instance::setupNodes")
              << "Failed to find mesh index for " << bcVec[iBC].eName;
          }
        }
        else
        {
          DevelFatal(*this).in("Instance::setupNodes")
            << "  Unrecognized side specified.";
        }
      }
      else
      {
        //DevelFatal(*this).in("Instance::setupNodes") << "Side NOT specified.";
      }

      bcVec[iBC].areaGiven = electrode.areaGiven;
      if (electrode.areaGiven)
      {
        bcVec[iBC].area = electrode.area;
      }
    }
  }

  indicesSetup_=true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << " area               = " << area << std::endl;
    Xyce::dout() << " areaGiven          = " << given("AREA") << std::endl;
    int isize=bcVec.size();
    for (int i=0;i<isize;++i)
    {
      Xyce::dout() << " bcVec["<<i<<"].area      = " << bcVec[i].area << std::endl;
      Xyce::dout() << " bcVec["<<i<<"].areaGiven = " << bcVec[i].areaGiven << std::endl;
      Xyce::dout() << " bcVec["<<i<<"].meshIndex = " << bcVec[i].meshIndex << std::endl;
    }
  }

  for (int i=0;i<bcVec.size();++i) { bcVec[i].dFdVckt.resize(1,0.0); } // one entry each for V 

  // allocate conductance array:
  numElectrodes = bcVec.size(); // should be n x n,
  // where n=number of terminals.
  condVec.resize(numElectrodes);
  for (int iE=0;iE<numElectrodes;++iE)
  {
    condVec[iE].resize(numElectrodes,0.0);
  }

  // initialize the boundary stencils.
  // Note: two of the points will be at meshIndex=0 and meshIndex=LX.  If there
  // is a 3rd terminal (for the base of a BJT) it will be somewhere in the middle.

  for (int i=0;i<bcVec.size();++i)
  {
    int meshIndex=bcVec[i].meshIndex;

    if (meshIndex==0 || meshIndex==LX)
    {
      edgeBoundarySten[meshIndex]=1;
    }
    else
    {
      internalBoundarySten[meshIndex]=1;
    }
    boundarySten[meshIndex]=1;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNumVars
//
// Purpose       : mostly sets up numIntVars.   numExtVars was set earlier
//                 and is easy.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/11/03
//-----------------------------------------------------------------------------
bool Instance::setupNumVars ()
{
  // Determine the proper size of NX, LX, and get the givens.
  bool compositeGiven=true;
  if(layerCompositeSpecified)
  {
    int matVecSize=materialVec.size();
    NX=0;
    width=0.0;
    for (int i=0;i<matVecSize;++i)
    {
      MaterialLayer & matLay = *(materialVec[i]);
      NX += matLay.NX;
      width += matLay.width;

      matLay.LX = matLay.NX-1;
    }

    LX = NX-1;
  }

  if (NXGiven)
  {
    LX = NX-1;
    numIntVars    = 3*NX;
    numStateVars  = numExtVars + NX - 1; // the NX-1 is for the displacement current.
    maxColsPerRow = MAX_COLS_PER_ROW;
  }
  else
  {
    UserFatal(*this) << "NX parameter was not specified.";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/11/03
//-----------------------------------------------------------------------------
bool Instance::setupJacStamp ()
{
  int iMeshNode;
  int numVars   = 3;
  int Voffset   = 0;
  int Noffset   = 1;
  int Poffset   = 2;
  int baseIndex = 0;
  int baseIndex_m1 = 0;
  int baseIndex_p1 = 0;

  int Vindex    = 0;
  int Nindex    = 0;
  int Pindex    = 0;

  // "minus one" indices
  int Vindex_m1    = 0;
  int Nindex_m1    = 0;
  int Pindex_m1    = 0;

  // "plus one" indices
  int Vindex_p1    = 0;
  int Nindex_p1    = 0;
  int Pindex_p1    = 0;

  int extVarOffset = numExtVars;
  int iBC;

  int jacSize = numIntVars + extVarOffset;

  meshToLID.clear(); meshToLID.resize(NX,-1);
  jacStamp.clear();   jacStamp.resize(jacSize);

  // set up the meshToLID converter first.
  int bcSize=bcVec.size();
  int lid=0;
  for (int iBC=0;iBC<bcSize;++iBC)
  {
    int meshIndex=bcVec[iBC].meshIndex;
    meshToLID[meshIndex] = lid;
    lid++;
  }

  for (int i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    meshToLID[i] = lid;
    lid++;
  }

  // external vars (from circuit) first:
  
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    iMeshNode = bcVec[iBC].meshIndex;

    if (edgeBoundarySten[iMeshNode]!=1 &&  internalBoundarySten[iMeshNode]!=1)
    {
      DevelFatal(*this).in("Instance::setupJacStamp")
        << "Boundary point not in the stencil.";
    }

    int iNN   = bcVec[iBC].neighborNode;
    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    int numConnectedMeshPoints=6;
    jacStamp[iBC].resize(numConnectedMeshPoints+bcVec.size(),-1);

    int col=0;
    if (iMeshNode > iNN) // i=0, or right-looking
    {
      baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;

      Vindex_m1 = baseIndex_m1 + Voffset;
      Nindex_m1 = baseIndex_m1 + Noffset;
      Pindex_m1 = baseIndex_m1 + Poffset;

      jacStamp[iBC][col++] = iBC;
      jacStamp[iBC][col++] = Vindex;
      jacStamp[iBC][col++] = Vindex_m1;
      jacStamp[iBC][col++] = Nindex;
      jacStamp[iBC][col++] = Nindex_m1;
      jacStamp[iBC][col++] = Pindex;
      jacStamp[iBC][col++] = Pindex_m1;
    }
    else // i=LX, or left-looking
    {
      baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

      Vindex_p1 = baseIndex_p1 + Voffset;
      Nindex_p1 = baseIndex_p1 + Noffset;
      Pindex_p1 = baseIndex_p1 + Poffset;

      jacStamp[iBC][col++] = iBC;
      jacStamp[iBC][col++] = Vindex;
      jacStamp[iBC][col++] = Vindex_p1;
      jacStamp[iBC][col++] = Nindex;
      jacStamp[iBC][col++] = Nindex_p1;
      jacStamp[iBC][col++] = Pindex;
      jacStamp[iBC][col++] = Pindex_p1;
    }

    // add columns for the rest of the connected circuit nodes.
    // This isn't needed for conventional Newton solves, but is 
    // needed for 2-level solves.
    for (int iBC2=0;iBC2<bcVec.size();++iBC2)
    {
      if (iBC2==iBC) continue;
      jacStamp[iBC][col++] = iBC2;
    }

    // these will be replaced later, but they are needed for that replacement
    bcVec[iBC].dIdXcols.clear();
    for (int ii=1;ii<7;ii++)
    {
      bcVec[iBC].dIdXcols.push_back(jacStamp[iBC][ii]);
      bcVec[iBC].dIdX.push_back(0.0);
    }
  }

  // Variables associated with the mesh.
  // First do the mesh points that are boundary condition points.
  // Do the anode (BC) mesh point:

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    iMeshNode = bcVec[iBC].meshIndex;
    int iNN   = bcVec[iBC].neighborNode;
    int NodeIndex = bcIndexMap[bcVec[iBC].eName]; // iBC?

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    jacStamp[Vindex].resize(3,-1);
    jacStamp[Nindex].resize(4,-1);
    jacStamp[Pindex].resize(4,-1);

    if (edgeBoundarySten[iMeshNode]==1)
    {
      if (iMeshNode < iNN) // i=0
      {
        baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;
        Vindex_p1 = baseIndex_p1 + Voffset;
        Nindex_p1 = baseIndex_p1 + Noffset;
        Pindex_p1 = baseIndex_p1 + Poffset;

        int col=0;
        jacStamp[Vindex][col++] = NodeIndex;
        jacStamp[Vindex][col++] = Vindex;
        jacStamp[Vindex][col++] = Vindex_p1;

        col=0;
        jacStamp[Nindex][col++] = Nindex;
        jacStamp[Nindex][col++] = Nindex_p1;
        jacStamp[Nindex][col++] = Pindex;
        jacStamp[Nindex][col++] = Pindex_p1;

        col=0;
        jacStamp[Pindex][col++] = Pindex;
        jacStamp[Pindex][col++] = Pindex_p1;
        jacStamp[Pindex][col++] = Nindex;
        jacStamp[Pindex][col++] = Nindex_p1;
      }
      else // i=LX
      {
        baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
        Vindex_m1 = baseIndex_m1 + Voffset;
        Nindex_m1 = baseIndex_m1 + Noffset;
        Pindex_m1 = baseIndex_m1 + Poffset;

        int col=0;
        jacStamp[Vindex][col++] = Vindex_m1;
        jacStamp[Vindex][col++] = Vindex;
        jacStamp[Vindex][col++] = NodeIndex;

        col=0;
        jacStamp[Nindex][col++] = Nindex_m1;
        jacStamp[Nindex][col++] = Nindex;
        jacStamp[Nindex][col++] = Pindex_m1;
        jacStamp[Nindex][col++] = Pindex;

        col=0;
        jacStamp[Pindex][col++] = Pindex_m1;
        jacStamp[Pindex][col++] = Pindex;
        jacStamp[Pindex][col++] = Nindex_m1;
        jacStamp[Pindex][col++] = Nindex;
      }
    }
    else if (internalBoundarySten[iMeshNode]==1) // probably base node
    {
      // base node applies a BC to the potential and majority carrier.
      // The minority carrier does not get a BC, and is treated like an
      // internal point.  So the stamp here is the same as an interior
      // point plus a little extra.  Not all of these will be used.

      baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
      baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
      baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

      Vindex_m1 = baseIndex_m1 + Voffset;
      Nindex_m1 = baseIndex_m1 + Noffset;
      Pindex_m1 = baseIndex_m1 + Poffset;
      Vindex    = baseIndex    + Voffset;
      Nindex    = baseIndex    + Noffset;
      Pindex    = baseIndex    + Poffset;
      Vindex_p1 = baseIndex_p1 + Voffset;
      Nindex_p1 = baseIndex_p1 + Noffset;
      Pindex_p1 = baseIndex_p1 + Poffset;

      // voltage col arrays:
      int col=0;
      jacStamp[Vindex].resize(6,-1);
      jacStamp[Vindex][col++] = NodeIndex;
      jacStamp[Vindex][col++] = Vindex_m1;
      jacStamp[Vindex][col++] = Vindex;
      jacStamp[Vindex][col++] = Vindex_p1;
      jacStamp[Vindex][col++] = Nindex;
      jacStamp[Vindex][col++] = Pindex;

      // electron col arrays:
      col=0;
      jacStamp[Nindex].resize(9,-1);
      jacStamp[Nindex][col++] = Nindex_m1;
      jacStamp[Nindex][col++] = Nindex;
      jacStamp[Nindex][col++] = Nindex_p1;
      jacStamp[Nindex][col++] = Vindex_m1;
      jacStamp[Nindex][col++] = Vindex;
      jacStamp[Nindex][col++] = Vindex_p1;
      jacStamp[Nindex][col++] = Pindex_m1;
      jacStamp[Nindex][col++] = Pindex;
      jacStamp[Nindex][col++] = Pindex_p1;

      // hole col arrays:
      col=0;
      jacStamp[Pindex].resize(9,-1);
      jacStamp[Pindex][col++] = Pindex_m1;
      jacStamp[Pindex][col++] = Pindex;
      jacStamp[Pindex][col++] = Pindex_p1;
      jacStamp[Pindex][col++] = Vindex_m1;
      jacStamp[Pindex][col++] = Vindex;
      jacStamp[Pindex][col++] = Vindex_p1;
      jacStamp[Pindex][col++] = Nindex_m1;
      jacStamp[Pindex][col++] = Nindex;
      jacStamp[Pindex][col++] = Nindex_p1;

    }
    else// not a boundary.  oops!
    {
      DevelFatal(*this).in("Instance::setupJacStamp")
        << "Boundary point not in the stencil.";
    }
  }

  // Now do the non-BC mesh points.
  for (iMeshNode=0;iMeshNode<NX;++iMeshNode)
  {
    if (boundarySten[iMeshNode]==1) continue;

    baseIndex_m1 = numVars*meshToLID[iMeshNode-1] + extVarOffset;
    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    baseIndex_p1 = numVars*meshToLID[iMeshNode+1] + extVarOffset;

    Vindex_m1 = baseIndex_m1 + Voffset;
    Nindex_m1 = baseIndex_m1 + Noffset;
    Pindex_m1 = baseIndex_m1 + Poffset;
    Vindex    = baseIndex    + Voffset;
    Nindex    = baseIndex    + Noffset;
    Pindex    = baseIndex    + Poffset;
    Vindex_p1 = baseIndex_p1 + Voffset;
    Nindex_p1 = baseIndex_p1 + Noffset;
    Pindex_p1 = baseIndex_p1 + Poffset;

    // voltage col arrays:
    int col=0;
    jacStamp[Vindex].resize(5,-1);
    jacStamp[Vindex][col++] = Vindex_m1;
    jacStamp[Vindex][col++] = Vindex;
    jacStamp[Vindex][col++] = Vindex_p1;
    jacStamp[Vindex][col++] = Nindex;
    jacStamp[Vindex][col++] = Pindex;

    // electron col arrays:
    col=0;
    jacStamp[Nindex].resize(9,-1);
    jacStamp[Nindex][col++] = Nindex_m1;
    jacStamp[Nindex][col++] = Nindex;
    jacStamp[Nindex][col++] = Nindex_p1;
    jacStamp[Nindex][col++] = Vindex_m1;
    jacStamp[Nindex][col++] = Vindex;
    jacStamp[Nindex][col++] = Vindex_p1;
    jacStamp[Nindex][col++] = Pindex_m1;
    jacStamp[Nindex][col++] = Pindex;
    jacStamp[Nindex][col++] = Pindex_p1;

    // hole col arrays:
    col=0;
    jacStamp[Pindex].resize(9,-1);
    jacStamp[Pindex][col++] = Pindex_m1;
    jacStamp[Pindex][col++] = Pindex;
    jacStamp[Pindex][col++] = Pindex_p1;
    jacStamp[Pindex][col++] = Vindex_m1;
    jacStamp[Pindex][col++] = Vindex;
    jacStamp[Pindex][col++] = Vindex_p1;
    jacStamp[Pindex][col++] = Nindex_m1;
    jacStamp[Pindex][col++] = Nindex;
    jacStamp[Pindex][col++] = Nindex_p1;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    // dump the jacStamp to stdout:
    int jacSize = jacStamp.size ();
    Xyce::dout() << "jacStamp size = " << jacSize << std::endl;

    for(int i=0;i<jacSize;++i)
    {
      int colSize = jacStamp[i].size();
      for (int j=0;j<colSize;++j)
      {
        Xyce::dout() << "  jacStamp["<<i<<"]["<<j<<"] = " << jacStamp[i][j] << std::endl;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::cleanupJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/12/10
//-----------------------------------------------------------------------------
bool Instance::cleanupJacStamp ()
{
#if 1
  // set up normal jacMap for when all resistances nonzero
  // If nothing is remapped, this amounts to a null operation when the
  // map is used later.  The maps become important when we start
  // remapping nodes because of zero lead resistances
  jacMap.clear();
  jacMap2.clear();
  jacMap.resize(jacStamp.size());
  jacMap2.resize(jacStamp.size());

  int mapSize = jacMap.size();
  for (int i=0;i<mapSize;++i)
  {
    jacMap[i]=i;
    jacMap2[i].resize(jacStamp[i].size());
    for (int j=0;j<jacStamp[i].size();++j)
    {
      jacMap2[i][j] = j;
    }
  }

  // Now fix the ordering of the columns in the jacStamp.  If the columns in each row
  // are not in ascending order, then the jacStampMap calls below (for removing
  // variables) will not work correctly.
  //
  // NOTE:  This is probably not safe to do, for the PDE devices, so column reordering is off
  // by default.
  if (columnReorderingFlag)
  {
    std::vector< std::vector<int> > tempStamp_eric;
    std::vector< std::vector<int> > tempMap2_eric;
    jacStampMap_fixOrder(jacStamp, jacMap2, tempStamp_eric, tempMap2_eric);
    jacStamp = tempStamp_eric;
    jacMap2 = tempMap2_eric;
  }

#endif // if 1

  return true;
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
    // set up the name map for this device.
    for (int i = 0; i < NX; ++i)
    {
      if (li_Vrowarray[i] != -1)
      {
        std::ostringstream oss;
        oss << "_V_" << i;
        addInternalNode(symbol_table, li_Vrowarray[i], getName(), oss.str());
      }

      if (li_Nrowarray[i] != -1)
      {
        std::ostringstream oss;
        oss << "_N_" << i;
        addInternalNode(symbol_table, li_Nrowarray[i], getName(), oss.str());
      }

      if (li_Prowarray[i] != -1)
      {
        std::ostringstream oss;
        oss << "_P_" << i;
        addInternalNode(symbol_table, li_Prowarray[i], getName(), oss.str());
      }
    }
}

//-----------------------------------------------------------------------------
// Function      : DiodeDPEInstance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::registerLIDs:\n";
    Xyce::dout() << "          name = " << getName() << std::endl;
    Xyce::dout() << "        numInt = " << numIntVars << std::endl;
    Xyce::dout() << "        numEXt = " << numExtVars << std::endl;
    Xyce::dout() << "        NX     = " << NX << std::endl;
  }

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    bcVec[iBC].lid = extLIDVec[iBC];
  }

  // First do the boundary condition mesh points.  There will be imposed
  // boundary conditions at each boundary
  // for potential, hole density and electron density.

  // The electrostatic potential, from the perspective of the device
  // simulation, is not neccessarily the same as the voltage used in the
  // circuit part of the code.  Also, obviously, the densities are not used by
  // the circuit sim., so these boundary conditions are considered internal
  // variables.

  int meshIndex = 0;
  int intIndex = 0;

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    meshIndex = bcVec[iBC].meshIndex;
    li_Vrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Nrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Prowarray[meshIndex] = intLIDVec[intIndex++];
  }

  // now do the interior points.  These will be blocked (V,N,P) together.
  //meshIndex=1;
  //while (meshIndex < LX)
  for (meshIndex=0;meshIndex<NX;++meshIndex)
  {
    if (boundarySten[meshIndex]==1) continue;

    li_Vrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Nrowarray[meshIndex] = intLIDVec[intIndex++];
    li_Prowarray[meshIndex] = intLIDVec[intIndex++];
  }

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int size=bcVec[iBC].dIdXcols.size();
    int extSize = extLIDVec.size();
    for (int ii=0;ii<size;ii++)
    {
      int stampRowIndex = bcVec[iBC].dIdXcols[ii];
      bcVec[iBC].dIdXcols[ii] = intLIDVec[stampRowIndex-extSize];
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "\n  solution indices:\n";

    for (int i=0;i<NX;++i)
    {
      Xyce::dout() << "     li_Vrowarray["<<i<<"] = " << li_Vrowarray[i];
      Xyce::dout() << "\tli_Nrowarray["<<i<<"] = " << li_Nrowarray[i];
      Xyce::dout() << "\tli_Prowarray["<<i<<"] = " << li_Prowarray[i] << std::endl;
    }
    Xyce::dout() << section_divider << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "  In Instance::registerStateLIDs\n\n";
    Xyce::dout() << "  name             = " << getName() << std::endl;
    Xyce::dout() << "  Number of State LIDs: " << numStateVars << std::endl;
  }

  // Copy over the local ID lists:
  staLIDVec = staLIDVecRef;

  int i;
  int j;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].li_stateC = staLIDVec[i];
  }

  for (i=0,j=2;i<NX-1;++i,++j)
  {
    li_stateDispl[i] = staLIDVec[j];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  State indices:" << std::endl;
    Xyce::dout() << std::endl;
    for (i=0;i<bcVec.size();++i)
    {
      Xyce::dout() << "bcVec["<<i<<"].li_stateC = "<<bcVec[i].li_stateC<< std::endl;
    }
    Xyce::dout() << std::endl;

    Xyce::dout() << "  Displacement current state variable local indices:" << std::endl;
    for (i=0;i<NX-1;++i)
    {
      Xyce::dout() << "  li_stateDispl["<<i<<"] = " << li_stateDispl[i] << std::endl;
    }
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/11/03
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
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/12/03
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs
( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs ( jacLIDVec );

  int numVars   = 3;
  int baseIndex;
  int Vindex;
  int Nindex;
  int Pindex;
  int Voffset = 0;
  int Noffset = 1;
  int Poffset = 2;
  int i,j;

  int extVarOffset = numExtVars;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::registerJacLIDs" << std::endl;

    int jacLIDSize = jacLIDVec.size();
    Xyce::dout() << "jacLIDSize = " << jacLIDSize << std::endl;
    for (i=0;i<jacLIDSize;++i)
    {
      int jacLIDcolSize = jacLIDVec[i].size();
      Xyce::dout() << std::endl;
      Xyce::dout() << "jacLIDVec["<<i<<"].size = " << jacLIDcolSize << std::endl;
      for (j=0;j<jacLIDcolSize;++j)
      {
        Xyce::dout() << "jacLIDVec["<<i<<"]["<<j<<"] = ";
        Xyce::dout() << jacLIDVec[i][j] << std::endl;
      }
    }
  }

  li_Vcolarray.resize(NX);
  li_Ncolarray.resize(NX);
  li_Pcolarray.resize(NX);


  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i1=0;
    int numCols = jacLIDVec[iBC].size();
    int crossOffsetSize = bcVec.size()-1;
    int colArraySize = numCols-crossOffsetSize;

    bcVec[iBC].li_colArray.resize(colArraySize,-1);
    bcVec[iBC].crossOffsets.clear();

    for (i1=0;i1<colArraySize;++i1)
    {
      bcVec[iBC].li_colArray[i1] = jacLIDVec[iBC][i1];
    }
    for (;i1<numCols;++i1)
    {
      bcVec[iBC].crossOffsets.push_back(jacLIDVec[iBC][i1]);
    }
    if (bcVec[iBC].crossOffsets.size() != crossOffsetSize)
    {
        DevelFatal(*this).in("Instance::registerJacLIDs")
          << "Wrong number of columns in jac LID array";
    }

    bcVec[iBC].lidOffset = jacLIDVec[iBC][0];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl;
      for(i=0;i<bcVec[iBC].li_colArray.size();++i)
      {
        Xyce::dout() << bcVec[iBC].eName << ": li_colArray["<<i<<"] = "<<bcVec[iBC].li_colArray[i]<< std::endl;
      }
      Xyce::dout() << std::endl;
      
      for(i=0;i<bcVec[iBC].crossOffsets.size();++i)
      {
        Xyce::dout() << bcVec[iBC].eName << ": crossOffsets["<<i<<"] = "<<bcVec[iBC].crossOffsets[i]<< std::endl;
      }
      Xyce::dout() << std::endl;
    }
  }

  // Do the non-BC mesh points.
  for (int iMeshNode=0;iMeshNode<NX;++iMeshNode)
  {
    if (boundarySten[iMeshNode]==1) continue;

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;

    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    // voltage col arrays:
    int i1=0;
    int j1=0;
    li_Vcolarray[iMeshNode].resize(5,-1);
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
    li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

    // electron col arrays:
    i1=0; j1=0;
    li_Ncolarray[iMeshNode].resize(9,-1);
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
    li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];


    // hole col arrays:
    i1=0; j1=0;
    li_Pcolarray[iMeshNode].resize(9,-1);
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "registerJacLIDs: iMeshNode = " << iMeshNode << std::endl;
      Xyce::dout() << "jacLIDVec[Vindex].size = " << jacLIDVec[Vindex].size()<<std::endl;
      Xyce::dout() << "jacLIDVec[Nindex].size = " << jacLIDVec[Nindex].size()<<std::endl;
      Xyce::dout() << "jacLIDVec[Pindex].size = " << jacLIDVec[Pindex].size()<<std::endl;

      for (i=0;i<5;++i)
      {
        Xyce::dout() << " li_Vcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Vcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<7;++i)
      {
        Xyce::dout() << " li_Ncolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Ncolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<7;++i)
      {
        Xyce::dout() << " li_Pcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Pcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
    }
  }

  // Do the BC mesh points.  These are actually "first" in the LID array.
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int iMeshNode = bcVec[iBC].meshIndex;
    int iNN   = bcVec[iBC].neighborNode;

    baseIndex    = numVars*meshToLID[iMeshNode  ] + extVarOffset;
    Vindex = baseIndex + Voffset;
    Nindex = baseIndex + Noffset;
    Pindex = baseIndex + Poffset;

    if (edgeBoundarySten[iMeshNode]==1)
    {
      int i1=0;
      int j1=0;
      li_Vcolarray[iMeshNode].resize(3,-1);
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

      // The stamp for N and P is dependent on the direction.
      if (iMeshNode < iNN)// i=0
      {
        i1=0; j1=0;
        li_Ncolarray[iMeshNode].resize(3,-1);
        li_Ncolarray[iMeshNode][i1++] = -1;
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

        i1=0; j1=0;
        li_Pcolarray[iMeshNode].resize(3,-1);
        li_Pcolarray[iMeshNode][i1++] = -1;
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      }
      else
      {
        i1=0; j1=0;
        li_Ncolarray[iMeshNode].resize(3,-1);
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
        li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

        i1=0; j1=0;
        li_Pcolarray[iMeshNode].resize(3,-1);
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
        li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      }
    }
    else if (internalBoundarySten[iMeshNode]==1) // probably base node
    {
      // voltage col arrays:
      int i1=0;
      int j1=0;
      li_Vcolarray[iMeshNode].resize(6,-1);
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];
      li_Vcolarray[iMeshNode][i1++] = jacLIDVec[Vindex][j1++];

      // electron col arrays:
      i1=0; j1=0;
      li_Ncolarray[iMeshNode].resize(9,-1);
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];
      li_Ncolarray[iMeshNode][i1++] = jacLIDVec[Nindex][j1++];

      // hole col arrays:
      i1=0; j1=0;
      li_Pcolarray[iMeshNode].resize(9,-1);
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
      li_Pcolarray[iMeshNode][i1++] = jacLIDVec[Pindex][j1++];
    }
    else// not a boundary.  oops!
    {
      DevelFatal(*this).in("Instance::registerJacLIDs")
        << "Boundary point not in the stencil.";
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "registerJacLIDs: ("<<bcVec[iBC].eName<<") iMeshNode = ";
      Xyce::dout() << iMeshNode << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Vcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Vcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Ncolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Ncolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
      for (i=0;i<3;++i)
      {
        Xyce::dout() << " li_Pcolarray["<<iMeshNode<<"]["<<i<<"] = ";
        Xyce::dout() << li_Pcolarray[iMeshNode][i] << std::endl;
      }
      Xyce::dout() << std::endl;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers()
//
// Purpose       : Sets up raw pointers for optimized matrix loads.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/10
//-----------------------------------------------------------------------------
void Instance::setupPointers()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  fVmatPtr.resize(NX);
  fNmatPtr.resize(NX);
  fPmatPtr.resize(NX);
  qVmatPtr.resize(NX);
  qNmatPtr.resize(NX);
  qPmatPtr.resize(NX);

  for (int i=0;i<NX;++i)
  {
    int Vrow = li_Vrowarray[i];
    int Nrow = li_Nrowarray[i];
    int Prow = li_Prowarray[i];

    int vSize = li_Vcolarray[i].size();
    fVmatPtr[i].resize(vSize);
    qVmatPtr[i].resize(vSize);
    for (int j=0;j<vSize;++j)
    {
      fVmatPtr[i][j] = &(dFdx[Vrow][li_Vcolarray[i][j]]);
      qVmatPtr[i][j] = &(dQdx[Vrow][li_Vcolarray[i][j]]);
    }

    int nSize = li_Ncolarray[i].size();
    fNmatPtr[i].resize(nSize);
    qNmatPtr[i].resize(nSize);
    for (int j=0;j<nSize;++j)
    {
      fNmatPtr[i][j] = &(dFdx[Nrow][li_Ncolarray[i][j]]);
      qNmatPtr[i][j] = &(dQdx[Nrow][li_Ncolarray[i][j]]);
    }

    int pSize = li_Pcolarray[i].size();
    fPmatPtr[i].resize(pSize);
    qPmatPtr[i].resize(pSize);
    for (int j=0;j<pSize;++j)
    {
      fPmatPtr[i][j] = &(dFdx[Prow][li_Pcolarray[i][j]]);
      qPmatPtr[i][j] = &(dQdx[Prow][li_Pcolarray[i][j]]);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermdiateVars
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
  bool bs1 = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "updateIntermediateVars.  name = " << getName() << std::endl;
  }

  bs1 = obtainSolution ();        bsuccess = bsuccess && bs1;
  bs1 = calcEfield ();            bsuccess = bsuccess && bs1;
  bs1 = calcMobilities   ();      bsuccess = bsuccess && bs1;
  bs1 = calcRecombination ();     bsuccess = bsuccess && bs1;
  bs1 = calcElectronCurrent ();   bsuccess = bsuccess && bs1;
  bs1 = calcHoleCurrent ();       bsuccess = bsuccess && bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcTerminalCurrents
// Purpose       : Calculates total diode current(s) to be used in the
//                 circuit KCL equations.
//
// Special Notes : Two options:
//
//  1) use the fluxes from the PDE calculation
//  2) use the integrated emission/capture rates from the rxn network.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/04/01
//-----------------------------------------------------------------------------
bool Instance::calcTerminalCurrents ()
{
  bool bsuccess = true;
  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  // Calculate the diode current using DD fluxes.
  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int index = bcVec[iBC].meshIndex;
    int iNN=bcVec[iBC].neighborNode;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;

    double sign =  ((iNN > index)?1.0:-1.0);
    int edgeIndex= ((iNN > index)?index:iNN);

    bcVec[iBC].elecCurrent = sign*JnxVec[edgeIndex]*A0;
    bcVec[iBC].holeCurrent = sign*JpxVec[edgeIndex]*A0;

    if (edgeBoundarySten[index]==1)
    {
      bcVec[iBC].currentSum = bcVec[iBC].elecCurrent + bcVec[iBC].holeCurrent;
    }
    else if (internalBoundarySten[index]==1)
    {
      std::string & type = bcVec[iBC].type;

      // only majority carrier goes to the boundary
      if (type=="ntype")
      {
        bcVec[iBC].currentSum = bcVec[iBC].elecCurrent;
      }
      else if (type=="ptype")
      {
        bcVec[iBC].currentSum = bcVec[iBC].holeCurrent;
      }
      else // oops.
      {
        DevelFatal(*this).in("Instance::calcTerminalCurrents")
          << "Unrecognized type on boundary.";
      }
    }
    else
    {
      DevelFatal(*this).in("Instance::calcTerminalCurrents")
        << "Unrecognized boundary.";
    }

    if (displCurrentFlag)
    {
      bcVec[iBC].currentSum += bcVec[iBC].displCurrent;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::subsection_divider << std::endl
                 << "Calculated currents, etc., coming from the DD calculation:" << std::endl
                 << "  scalingVars.J0 = " << J0<<std::endl
                 << "  scalingVars.a0 = " << a0<<std::endl
                 << Xyce::subsection_divider << std::endl;
    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      Xyce::dout() << bcVec[iBC];
    }
    Xyce::dout() << Xyce::subsection_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdTerminalCurrents
//
// Purpose       : This function calculates partial derivatives associated
//                 with the Jacobian loads for the KCL equation rows.
//
// Special Notes : Originally, this work was performed in the loadJacDD
//                 function.  However, some of this information is also needed
//                 for the decoupled 2-level Newton, to calculate the
//                 terminal conductances.
//
//                 To calculate the terminal conductances, the following is
//                 needed for each electrode.
//
//                  dIdVckt - derivative of terminal current w.r.t. Vckt.
//                            This is the also the Jacobian contribution
//                            for the (KCL row, KCL col) entry of the matrix.
//
//                  dFdVckt - derivative of the RHS vector w.r.t. Vckt.
//                            This is a vector quantity, and corresponds to
//                            the (*, Vckt) column of the PDE matrix
//                            sub-block.
//
//                  dIdX  - derivative of the terminal current w.r.t. the
//                          vector of PDE solution variables. (ie not
//                          including Vckt, as that is not part of the PDE
//                          domain).  This is a vector quantity.
//                          This corresponds to the (KCL row, *) entry of
//                          the matrix, modulo dIdVckt.
//
//                  With the  "new" boundary conditions, the dFdVckt vector
//                  should only have one nonzero element, which corresponds
//                  to dIdVckt.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/03
//-----------------------------------------------------------------------------
bool Instance::pdTerminalCurrents ()
{
  std::string msg;

  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  // calculate dIdVckt.----------------------------------------
  // inspired from the loadMatKCLDDForm function.
  int iBC;
  int bcSize = bcVec.size();
  for (iBC=0; iBC < bcSize; ++iBC)
  {
    // This is always zero.
    bcVec[iBC].dIdVckt = 0.0;
  } 

  // calculate dFdVckt.----------------------------------------
  // dFdVckt is a column of the matrix, which corresponds to the variable
  // Vckt, but only includes rows for internal PDE device variables.
  // There is a single entry in this column.
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int dFdVcktIndex = 0;
    bcVec[iBC].dFdVckt[dFdVcktIndex] = -scalingVars.rV0; // V
  }

  // calculate dIdX.-------------------------------------------------
  // This is adapted from the function loadJacKCLDDFormulation.
  // dIdX is a KCL row from the matrix, but with only the
  // columns of internal PDE device variables included.
  for (iBC=0;iBC<bcSize;++iBC)
  {
    int i=bcVec[iBC].meshIndex; //i is the mesh point

    if (edgeBoundarySten[i]!=1 &&  internalBoundarySten[i]!=1)
    {
      DevelFatal(*this).in("Instance::pdTerminalCurrents")
        << "Unrecognized boundary.";
    }

    int iNN=bcVec[iBC].neighborNode;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;

    double sign = ((iNN > i)?1.0:-1.0);
    double dJndV = 0.0;
    double dJpdV = 0.0;
    double dJndn = 0.0;
    double dJpdp = 0.0;

    double dJndV_nn = 0.0;
    double dJpdV_nn = 0.0;
    double dJndn_nn = 0.0;
    double dJpdp_nn = 0.0;

    // if looking right, then edge=i, if looking left, then edge=i-1.
    if (i<iNN)
    {
      int edgeIndex=i;
      dJndV = dJndV1Vec[edgeIndex];
      dJpdV = dJpdV1Vec[edgeIndex];
      dJndn = dJndn1Vec[edgeIndex];
      dJpdp = dJpdp1Vec[edgeIndex];

      dJndV_nn = dJndV2Vec[edgeIndex];
      dJpdV_nn = dJpdV2Vec[edgeIndex];
      dJndn_nn = dJndn2Vec[edgeIndex];
      dJpdp_nn = dJpdp2Vec[edgeIndex];
    }
    else
    {
      int edgeIndex=iNN;
      dJndV = dJndV2Vec[edgeIndex];
      dJpdV = dJpdV2Vec[edgeIndex];
      dJndn = dJndn2Vec[edgeIndex];
      dJpdp = dJpdp2Vec[edgeIndex];

      dJndV_nn = dJndV1Vec[edgeIndex];
      dJpdV_nn = dJpdV1Vec[edgeIndex];
      dJndn_nn = dJndn1Vec[edgeIndex];
      dJpdp_nn = dJpdp1Vec[edgeIndex];
    }

    // center (boundary) point:
    double Vcoef = (sign* dJndV + sign* dJpdV)*A0;
    double Ncoef = (sign* dJndn)*A0;
    double Pcoef = (sign* dJpdp)*A0;

    // neighbor point:
    double Vcoef_nn = (sign* dJndV_nn + sign* dJpdV_nn)*A0;
    double Ncoef_nn = (sign* dJndn_nn)*A0;
    double Pcoef_nn = (sign* dJpdp_nn)*A0;

    if (internalBoundarySten[i]==1)
    {
      // only majority carrier goes to the boundary
      if      (bcVec[iBC].type=="ntype") { Pcoef=0.0; Pcoef_nn=0.0; }
      else if (bcVec[iBC].type=="ptype") { Ncoef=0.0; Ncoef_nn=0.0; }
      else // oops.
      {
        DevelFatal(*this).in("Instance::pdTerminalCurrents")
          << "Unrecognized type " << bcVec[iBC].type << " on boundary.";
      }
    }

    // nn=i+1 for right-looking or nn=i-1 for left-looking
    // The "dIdXcols" arrays were set up elsewhere, in a combination of the 
    // setupJacStamp and registerLID functions.  They should not contain offsets 
    // and don't need to be set up here.
    int count = 0;
    bcVec[iBC].dIdX[count++] = Vcoef; // derivative w.r.t. V[i]:
    bcVec[iBC].dIdX[count++] = Vcoef_nn; // derivative w.r.t. V[nn].
    bcVec[iBC].dIdX[count++] = Ncoef; // derivative w.r.t. n[i]:
    bcVec[iBC].dIdX[count++] = Ncoef_nn; // derivative w.r.t. n[nn].
    bcVec[iBC].dIdX[count++] = Pcoef; // derivative w.r.t. p[i]:
    bcVec[iBC].dIdX[count++] = Pcoef_nn; // derivative w.r.t. p[nn].
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDXDV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::calcDXDV ()
{

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDFDV
// Purpose       : Load -dfdv into the RHS vector for the specified
//                 electrode.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::loadDFDV (int ielectrode, Linear::Vector * dfdvPtr)
{
  Linear::Vector & dfdv = *(dfdvPtr);
  bcData & bc = bcVec[ielectrode];

  // load only the V term:
  int Vrow = li_Vrowarray[bc.meshIndex]; 
  int dFdVindex = 0;
  dfdv[Vrow] = - (bc.dFdVckt[dFdVindex]);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcConductance
//
// Purpose       : Calculates device conductances for a single electrode of
//                 the PDE device.  This function is for
//                 calculating conductances between extern circuit nodes.
//
//                 The point of calculating these quantities is to provide
//                 a lumped parameter substitute for the full device, when
//                 running 2-level Newton.
//
// Special Notes : This function is (ultimately) invoked from the nonlinear
//                 solver, as that part of the code, when running in
//                 2-level mode, knows when this information is needed.
//
//                 It is assumed that when this function is called, the
//                 "deltaX" vector contains the information:  dXdVckt, where
//                 X is solution vector variables associated with the PDE
//                 device, while Vckt is the voltage on the attached
//                 circuit node.  The reason it is in the "deltaX" vector
//                 is that it was obtained via a linear solver of the
//                 problem:
//
//                 dXdVckt = J^-1 . dFdVckt
//
//                 dFdVckt was calculated previously in pdTerminalCurrents,
//                 and J is the Jacobian.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/07/03
//-----------------------------------------------------------------------------
bool Instance::calcConductance (int iElectrode, const Linear::Vector * dxdvPtr)
{
  bool bsuccess = true;
  const Linear::Vector & dxdv = *dxdvPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "calcConductances  name = " << getName() << std::endl;
    Xyce::dout() << "electrode = " << bcVec[iElectrode].eName;
    Xyce::dout() << "  dIdVckt["<<iElectrode<<"] = " << bcVec[iElectrode].dIdVckt;
    Xyce::dout() << "  currentSum["<<iElectrode<<"] = " << bcVec[iElectrode].currentSum;
    Xyce::dout() << std::endl;
    Xyce::dout() << std::endl;
  }

  if (!(bcVec[iElectrode].dxdvAllocated))
  {
    bcVec[iElectrode].dxdvPtr = extData.lasSysPtr->builder().createVector();
    bcVec[iElectrode].dxdvAllocated = true;
  }

  // A linear solve should have just been performed up in the Newton
  // solver.  The result of that solve, dxdv was placed in the RHS vector.
  // dxdv is needed later, so save a copy.
  *(bcVec[iElectrode].dxdvPtr) = *(dxdvPtr);

  // doing the iElectrode Column of the condVec array.
  // This should correspond to the bcVec[iElectrode].lid column of the
  // Jacobian.

  double Gij = 0.0;
  for (int iEqu=0;iEqu< numElectrodes; ++iEqu)
  {
    // conductance Gij .
    //
    // subscript i = variable, which is one of the electrode voltages.
    // subscript j = electrode
    //
    //   Gij = dot( dIi/dX , dX/dVj ) + dIi/dVj
    //
    //  if i != j, then the last term is zero.  
    //  For the 1D device, unlike the 2D device, the dIdVckt term is ALWAYS zero.
    //  So for the 1D device, Gij = dot (dIi/dX, dX/dVj)

    // load dIdX:
    Linear::Vector & dIdX = *(extData.tmpdIdXPtr);
    dIdX.putScalar (0.0);
    int DIDXSize = bcVec[iEqu].dIdX.size();

    for (int iDIDX=0;iDIDX<DIDXSize;++iDIDX)
    {
      dIdX[bcVec[iEqu].dIdXcols[iDIDX]] = bcVec[iEqu].dIdX[iDIDX];
    }

    if (DEBUG_DEVICE)
    {
      std::ostringstream oss;
      oss << "dIdX" << std::setw(2) << std::setfill('0') << iEqu << ".txt";
      dIdX.writeToFile(oss.str().c_str());
    }

    // get dot product:
    Gij = dxdv.dotProduct( dIdX );
    condVec[iEqu][iElectrode] = Gij;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      if (false)
      //if (iEqu==2 && iElectrode==2)
      {
        // ERK. 2/214/2019, added during work on bug 1155.
        //
        // This debug output is here b/c when I try to use 2-level on the 3-terminal 
        // BJT test, the 2,2 entry in the extracted Jacobian winds up being zero.  Also 
        // that extracted Jacobian isn't symmetric.  So there is a problem.  This problem 
        // doesn't happen with the 2-terminal tests I've tried.
        //
        // To the precision of what is ouptut here, the dot product should not be evaluating 
        // to zero.  But it is, which suggests that at full precision all the terms cancel.  
        //
        // As of this writing, this problems isn't resolved and I need to move on.
        Xyce::dout() << "(2,2) dIdX vector:"<<std::endl;
        dIdX.printPetraObject(Xyce::dout());

        Xyce::dout() << "(2,2) dxdv vector:"<<std::endl;
        dxdv.printPetraObject(Xyce::dout());

        double testDotProd = dxdv.dotProduct( dIdX );
        Xyce::dout() << "test dot product = " << testDotProd <<std::endl;
      }

      char outstring[128];
      double Itmp = bcVec[iEqu].currentSum;
      double Vtmp = bcVec[iEqu].Vckt - bcVec[iElectrode].Vckt;
      Vtmp *= scalingVars.V0;
      double GV = Gij*Vtmp;
      for(int i=0;i<128;++i) outstring[i] = static_cast<char>(0);
      sprintf(outstring, "(%2d,%2d): G=%12.4e", iEqu,iElectrode,Gij);
      Xyce::dout() << std::string(outstring) << std::endl;
      sprintf(outstring, "(%2d,%2d): G=%12.4e G*V=%12.4e I=%12.4e V=%12.4e",
              iEqu,iElectrode,Gij,GV,Itmp,Vtmp);
      Xyce::dout() << std::string(outstring) << std::endl;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  updateIntermediateVars ();
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);

  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int li_state=bcVec[iBC].li_stateC;
    staVector[li_state] = bcVec[iBC].currentSum;
  }

  // Now store the dielectric displacement in the state vector for
  // displacement current.
  int i;
  for (i = 0; i< NX-1; ++i)
  {
    double D = eSi * e0 * scalingVars.E0 * ExVec[i];
    staVector[li_stateDispl[i]] = D;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & staDeriv =  *(extData.nextStaDerivVectorPtr);

  // Now get displacement current.
  for( int i = 0; i< NX-1; ++i)
  {
    displCurrent[i] = staDeriv[li_stateDispl[i]];
  }

  displCurrent[LX] = displCurrent[LX-1];

  calcTerminalCurrents ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadErrorWeightMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
//                 For this device, the electrostatic potential is
//                 optionally excluded from time integration norm
//                 calculations.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/07/10
//-----------------------------------------------------------------------------
void Instance::loadErrorWeightMask()
{
  if (maskVarsTIAFlag_)
  {
    Linear::Vector * maskVectorPtr = extData.deviceErrorWeightMask_;

    for (int i=0;i<NX;++i)
    {
      int Vrow = li_Vrowarray[i];
      (*maskVectorPtr)[Vrow] = 0.0;
      (*maskVectorPtr)[Vrow] = 0.0;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool Instance::setInitialGuess ()
{
  bool bsuccess = true;
  bool bs1 = true;

  if (variablesScaled)
  {
    bs1 = unScaleVariables ();    bsuccess = bsuccess && bs1;
  }
  bs1 = calcDensityBCs      (); bsuccess = bsuccess && bs1;
  bs1 = calcVequBCs         (); bsuccess = bsuccess && bs1;
  bs1 = calcInitialGuess    (); bsuccess = bsuccess && bs1;
  bs1 = calcMobilities      (); bsuccess = bsuccess && bs1;
  bs1 = calcLifetimes       (); bsuccess = bsuccess && bs1;
  bs1 = scaleVariables      (); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecNLPoisson
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
bool Instance::loadVecNLPoisson (double * rhs)
{
  bool bsuccess = true;
  int i(-1);
  int Vrow(-1), Nrow(-1), Prow(-1);
  double coef(0); 
  double coef2(0);

  Ut = Vt/scalingVars.V0;

  // KCL equations for the two connecting terminals:
  // For the NL poisson, there is no coupling to the circuit,
  // so nothing to do here.


  // boundary conditions:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    // no coupling to ckt, so effectively Vckt is hardwired to zero.
    //rhs[Vrow] += VVec[i] - (bcVec[iBC].Vckt + bcVec[iBC].Vequ);
    rhs[Vrow] += VVec[i] - (bcVec[iBC].Vequ);
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;
  }

  // mesh points associated with heterojunction boundaries:
  int hetSize=heterojunctionBCs.size();
  for (int ihet=0;ihet<hetSize;++ihet)
  {
    int i1=heterojunctionBCs[ihet].first;
    int i2=heterojunctionBCs[ihet].second;

    // recall ExVec[i] = -(VVec[i+1] - VVec[i])/dxVec[i];
    coef = (relPermVec[i1-1]*ExVec[i1-1]-relPermVec[i2]*ExVec[i2]);
    Vrow = li_Vrowarray[i1];
    Nrow = li_Nrowarray[i1];
    Prow = li_Prowarray[i1];
    rhs[Vrow] += coef;
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;

    // only need 1 equation to enforce this.
    Vrow = li_Vrowarray[i2];
    Nrow = li_Nrowarray[i2];
    Prow = li_Prowarray[i2];
    rhs[Vrow] = 0.0;
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;
  }

  // mesh points for the PDE problem:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    if ( heterojunctionSten[i]!=0 ) continue;

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);
    coef  = (relPermVec[i]*ExVec[i]-relPermVec[i-1]*ExVec[i-1])/aveDx;
    coef  *= scalingVars.L0;

    double holeDens = getVoltDepHoleDens ( VminExp, VVec[i], Na);
    double elecDens = getVoltDepElecDens ( VmaxExp, VVec[i], Nd);
    coef2 = -(holeDens-elecDens+CVec[i]);

    coef += coef2;
    rhs[Vrow] += coef;

    // Now do electron, hole continuity
    rhs[Nrow] = 0.0;
    rhs[Prow] = 0.0;
  } // row loop...

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadVecDDForm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
bool Instance::loadVecDDForm (double * rhs)
{
  bool bsuccess = true;
  int i;
  int Vrow, Nrow, Prow;
  double coef, coef2;

  // KCL equations for the two connecting terminals:
  // if this is the inner loop of a multilevel Newton solve, don't do the
  // KCL-related loads.
  if ( !(getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM))
  {
    for (int iBC=0;iBC<bcVec.size();++iBC)
    {
      rhs[bcVec[iBC].lid] += bcVec[iBC].currentSum;
    }
  } // end of twoLevelNewtonCouplingMode if statement.

  // mesh points for the PDE problem:

  // boundary conditions:
  // all of these take a Dirchlet BC on voltage.
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    rhs[Vrow] += VVec[i]-bcVec[iBC].Vbc;

    if (edgeBoundarySten[i])
    {
      rhs[Nrow] = nnVec[i]-bcVec[iBC].nnbc;
      rhs[Prow] = npVec[i]-bcVec[iBC].npbc;
    }
    else if (internalBoundarySten[i])
    {
      std::string & type = bcVec[iBC].type;
      double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);

      if (type=="ntype")  // boundary condition on e-, let h+ flow
      {
        rhs[Nrow] = nnVec[i]-bcVec[iBC].nnbc;
        rhs[Prow] = -(JpxVec[i]-JpxVec[i-1])/aveDx - RVec[i];
      }
      else if (type=="ptype")  // boundary condition on h+, let e- flow
      {
        rhs[Nrow] = (JnxVec[i]-JnxVec[i-1])/aveDx - RVec[i];
        rhs[Prow] = npVec[i]-bcVec[iBC].npbc;
      }
      else
      {
        DevelFatal(*this).in("Instance::loadVecDDForm")
          << "Unrecognized type on boundary.";
      }
    }
    else
    {
        DevelFatal(*this).in("Instance::loadVecDDForm")
          << "Unrecognized stencil on boundary.";
    }
  }

  // interior mesh points:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    if ( heterojunctionSten[i]!=0 ) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);
    coef  = (relPermVec[i]*ExVec[i]-relPermVec[i-1]*ExVec[i-1])/aveDx;
    coef  *= scalingVars.L0;
    coef2 = -(npVec[i]-nnVec[i]+CVec[i]);

    coef += coef2;
    rhs[Vrow] += coef;

    // Now do electron continuity
    // get electron time derivative and scale.
    rhs[Nrow] = (JnxVec[i]-JnxVec[i-1])/aveDx - RVec[i];

    // Now do hole continuity
    // get hole time derivative and scale.
    rhs[Prow] = -(JpxVec[i]-JpxVec[i-1])/aveDx - RVec[i];

  } // row loop...

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
  std::vector<Util::BreakPoint> &breakPointTimes)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatNLPoisson
// Purpose       : This function performs an analytic Jacobian matrix load for
//                 the diode-pde class, for the case of solving a nonlinear
//                 poisson equation.
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/18/01
//-----------------------------------------------------------------------------
bool Instance::loadMatNLPoisson (Linear::Matrix & mat)
{
  bool bsuccess = true;
  int Vrow, Nrow, Prow;
  int i,j;

  Ut = Vt/scalingVars.V0;
  double rUt = 1.0/Ut;

  double elecDens, holeDens;
  double dx1, dx2;
  int iBC;

  // Load the jacobian, row by row.

  // For the NL poisson option, device is not coupled to the ckt,
  // so put 1's on the diagonal.
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    mat[bcVec[iBC].lid][bcVec[iBC].lidOffset] = 1.0;
  }

  // boundary conditions on the mesh:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    int offset1 = li_Vcolarray[i][0];
    int offset2 = li_Vcolarray[i][1];
    int offset3 = li_Vcolarray[i][2];

    if (i==0)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 1.0;
      mat[Vrow][offset3] = 0.0;
    }
    else if (i==LX)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 1.0;
      mat[Vrow][offset3] = 0.0;
    }
    else if (internalBoundarySten[i]==1)
    {
      mat[Vrow][offset1] = 0.0;
      mat[Vrow][offset2] = 0.0;
      mat[Vrow][offset3] = 1.0;
    }

    mat[Nrow][li_Ncolarray[i][1]] = 1.0;
    mat[Prow][li_Pcolarray[i][1]] = 1.0;
  }


  // mesh points associated with heterojunction boundaries:
  int hetSize=heterojunctionBCs.size();
  for (int ihet=0;ihet<hetSize;++ihet)
  {
    int i1=heterojunctionBCs[ihet].first;
    int i2=heterojunctionBCs[ihet].second;

    // recall ExVec[i] = -(VVec[i+1] - VVec[i])/dxVec[i];
    //coef = (relPermVec[i1-1]*ExVec[i1-1]-relPermVec[i2]*ExVec[i2]);
    Vrow = li_Vrowarray[i1];
    Nrow = li_Nrowarray[i1];
    Prow = li_Prowarray[i1];

    int offset1 = li_Vcolarray[i1][0];
    int offset2 = li_Vcolarray[i1][1];
    int offset3 = li_Vcolarray[i1][2];
    //int offset4 = li_Vcolarray[i][3];  // not available yet

    dx1 = dxVec[i1-1];
    dx2 = dxVec[i2];

    double perm1=relPermVec[i1-1];
    double perm2=relPermVec[i2];

    mat[Vrow][offset1] = perm1/dx1;
    mat[Vrow][offset2] = -perm1/dx1;
    mat[Vrow][offset3] = -perm2/dx2;
    //mat[Vrow][offset4] = perm2/dx2; // not available yet

    mat[Nrow][li_Ncolarray[i1][1]] = 1.0;
    mat[Prow][li_Pcolarray[i1][1]] = 1.0;

    // only need 1 equation to enforce this.
    Vrow = li_Vrowarray[i2];
    Nrow = li_Nrowarray[i2];
    Prow = li_Prowarray[i2];
    mat[Vrow][li_Vcolarray[i2][1]] = 1.0;
    mat[Nrow][li_Ncolarray[i2][1]] = 1.0;
    mat[Prow][li_Pcolarray[i2][1]] = 1.0;
  }

  // rows associated with the PDE mesh:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    if ( heterojunctionSten[i]!=0 ) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    holeDens = getVoltDepHoleDens ( VminExp   , VVec[i], Na);
    elecDens = getVoltDepElecDens ( VmaxExp   , VVec[i], Nd);

    dx1 = dxVec[i-1];
    dx2 = dxVec[i];
    double L0 = scalingVars.L0 * MaterialSupport::getRelPerm(semi);
   // double L0 = scalingVars.L0 * relPermVec[i];

    int Vrow = li_Vrowarray[i];
    int Nrow = li_Nrowarray[i];
    int Prow = li_Prowarray[i];

    int offset1 = li_Vcolarray[i][0];
    int offset2 = li_Vcolarray[i][1];
    int offset3 = li_Vcolarray[i][2];

    mat[Vrow][offset1] = -L0/(dx1*dx2);
    mat[Vrow][offset2] = 2.0*L0/(dx1*dx2) + rUt*holeDens + rUt*elecDens;
    mat[Vrow][offset3] = -L0/(dx1*dx2);

    mat[Nrow][li_Ncolarray[i][1]] = 1.0;
    mat[Prow][li_Pcolarray[i][1]] = 1.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatKCLDDForm
// Purpose       : Loads drift-diffusion-KCL equations into a matrix.
// Special Notes : This function is used for both old and new DAE.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool Instance::loadMatKCLDDForm (Linear::Matrix & mat)
{
  int Vrow;
  int iBC;

  double & J0 = scalingVars.J0;
  double & a0 = scalingVars.a0;

  int i,j;
  int liOffIndex;
  int count  = 0;

  // rows associated with the connecting terminal KCL's:
  int li_row; // circuit node lid

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    liOffIndex = 1;
    i=bcVec[iBC].meshIndex; //i is the mesh point

    if (edgeBoundarySten[i]!=1 &&  internalBoundarySten[i]!=1)
    {
      DevelFatal(*this).in("Instance::loadMatKCLForm")
        << "Unrecognized boundary.";
    }

    int iNN=bcVec[iBC].neighborNode;
    li_row = bcVec[iBC].lid;
    double & area = bcVec[iBC].area;
    double A0=J0*a0*area;
    std::vector<int> & colA = bcVec[iBC].li_colArray;

    double sign = ((iNN > i)?1.0:-1.0);
    double dJndV = 0.0;
    double dJpdV = 0.0;
    double dJndn = 0.0;
    double dJndp = 0.0;
    double dJpdp = 0.0;
    double dJpdn = 0.0;

    double dJndV_nn = 0.0;
    double dJpdV_nn = 0.0;
    double dJndn_nn = 0.0;
    double dJndp_nn = 0.0;
    double dJpdp_nn = 0.0;
    double dJpdn_nn = 0.0;

    // if looking right, then edge=i, if looking left, then edge=i-1.
    int edgeIndex=i;

    if (i<iNN)
    {
      edgeIndex=i;
      dJndV = dJndV1Vec[edgeIndex];
      dJpdV = dJpdV1Vec[edgeIndex];
      dJndn = dJndn1Vec[edgeIndex];
      dJndp = dJndp1Vec[edgeIndex];
      dJpdp = dJpdp1Vec[edgeIndex];
      dJpdn = dJpdn1Vec[edgeIndex];

      dJndV_nn = dJndV2Vec[edgeIndex];
      dJpdV_nn = dJpdV2Vec[edgeIndex];
      dJndn_nn = dJndn2Vec[edgeIndex];
      dJndp_nn = dJndp2Vec[edgeIndex];
      dJpdp_nn = dJpdp2Vec[edgeIndex];
      dJpdn_nn = dJpdn2Vec[edgeIndex];
    }
    else
    {
      edgeIndex=iNN;
      dJndV = dJndV2Vec[edgeIndex];
      dJpdV = dJpdV2Vec[edgeIndex];
      dJndn = dJndn2Vec[edgeIndex];
      dJndp = dJndp2Vec[edgeIndex];
      dJpdp = dJpdp2Vec[edgeIndex];
      dJpdn = dJpdn2Vec[edgeIndex];

      dJndV_nn = dJndV1Vec[edgeIndex];
      dJpdV_nn = dJpdV1Vec[edgeIndex];
      dJndn_nn = dJndn1Vec[edgeIndex];
      dJndp_nn = dJndp1Vec[edgeIndex];
      dJpdp_nn = dJpdp1Vec[edgeIndex];
      dJpdn_nn = dJpdn1Vec[edgeIndex];
    }

    // center (boundary) point:
    double Vcoef = sign*(dJndV + dJpdV)*A0;
    double Ncoef = sign*(dJndn + dJpdn)*A0;
    double Pcoef = sign*(dJndp + dJpdp)*A0;

    // neighbor point:
    double Vcoef_nn = sign*(dJndV_nn + dJpdV_nn)*A0;
    double Ncoef_nn = sign*(dJndn_nn + dJpdn_nn)*A0;
    double Pcoef_nn = sign*(dJndp_nn + dJpdp_nn)*A0;

    if (internalBoundarySten[i]==1)
    {
      std::string & type = bcVec[iBC].type;

      // only majority carrier goes to the boundary
      if (type=="ntype")
      {
       // ERK.  2/20/2019 not correcting Vcoef seems wrong
        Vcoef = sign*(dJndV)*A0;
        Vcoef_nn = sign*(dJndV_nn)*A0;
        Pcoef=0.0;
        Pcoef_nn=0.0;
      }
      else if (type=="ptype")
      {
        // ERK.  2/20/2019 not correcting Vcoef seems wrong
        Vcoef = sign*(dJpdV)*A0;
        Vcoef_nn = sign*(dJpdV_nn)*A0;
        Ncoef=0.0;
        Ncoef_nn=0.0;
      }
      else // oops.
      {
        DevelFatal(*this).in("Instance::loadMatKCLForm")
          << "Unrecognized type on boundary.";
      }
    }

    // nn=i+1 for right-looking or nn=i-1 for left-looking
    // derivative w.r.t. V[i]:
    mat[li_row][colA[liOffIndex++]] += Vcoef;

    // derivative w.r.t. V[nn].
    mat[li_row][colA[liOffIndex++]] += Vcoef_nn;

    // derivative w.r.t. n[i]:
    mat[li_row][colA[liOffIndex++]] += Ncoef;

    // derivative w.r.t. n[nn].
    mat[li_row][colA[liOffIndex++]] += Ncoef_nn;

    // derivative w.r.t. p[i]:
    mat[li_row][colA[liOffIndex++]] += Pcoef;

    // derivative w.r.t. p[nn].
    mat[li_row][colA[liOffIndex++]] += Pcoef_nn;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatCktTrivial
// Purpose       : This function handles rows associated with the connecting
//                 terminal KCL's:
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool Instance::loadMatCktTrivial (Linear::Matrix & mat)
{
  bool bsuccess = true;
  bool bs1 = true;
  int j;
  int iBC;

  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    mat[bcVec[iBC].lid][bcVec[iBC].lidOffset] = 1.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadMatDDForm
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
bool Instance::loadMatDDForm (Linear::Matrix & mat)
{
  bool bsuccess = true;
  bool bs1 = true;

  int Vrow, Nrow, Prow;

  int i,j;
  double coef;

  // set up some of the partial derivative arrays:
  bs1 = pdRecombination ();    bsuccess = bsuccess && bs1;
  bs1 = pdElectronCurrent ();  bsuccess = bsuccess && bs1;
  bs1 = pdHoleCurrent ();      bsuccess = bsuccess && bs1;
  bs1 = pdTerminalCurrents (); bsuccess = bsuccess && bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    for (int i=1;i<LX;++i)
    {
      double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);
      Xyce::dout()<<"\t" <<dJndp1Vec[i]/aveDx<<"\t"<<dJndp2Vec[i]/aveDx
               <<"\t" <<dJpdn1Vec[i]/aveDx<<"\t"<<dJpdn2Vec[i]/aveDx<<std::endl;
    }
  }

  if ( !(getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM))
  {
    bs1 = loadMatKCLDDForm ( mat );
    bsuccess = bsuccess && bs1;
  }
  else
  {
    bs1 = loadMatCktTrivial ( mat );
    bsuccess = bsuccess && bs1;
  }

  // boundary conditions on the mesh:
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    if (i==0)
    {
      mat[Vrow][li_Vcolarray[i][0]] = -scalingVars.rV0;
      mat[Vrow][li_Vcolarray[i][1]] =  1.0;
      mat[Vrow][li_Vcolarray[i][2]] =  0.0;

      mat[Nrow][li_Ncolarray[i][1]] = +1.0;
      mat[Prow][li_Pcolarray[i][1]] = 1.0;
    }
    else if (i==LX)
    {
      mat[Vrow][li_Vcolarray[i][0]] =  0.0;
      mat[Vrow][li_Vcolarray[i][1]] =  1.0;
      mat[Vrow][li_Vcolarray[i][2]] = -scalingVars.rV0;

      mat[Nrow][li_Ncolarray[i][1]] = +1.0;
      mat[Prow][li_Pcolarray[i][1]] = 1.0;
    }
    else if (internalBoundarySten[i]==1)
    {
      mat[Vrow][li_Vcolarray[i][0]] = -scalingVars.rV0;
      mat[Vrow][li_Vcolarray[i][2]] = 1.0;

      std::string & type = bcVec[iBC].type;

      if (type=="ntype")  // boundary condition on e-, let h+ flow
      {
        mat[Nrow][li_Ncolarray[i][1]] = +1.0;

        double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

        // derivative w.r.t. npVec[i-1]:
        mat[Prow][li_Pcolarray[i][0]] = dJpdp1Vec[i-1]/aveDx;

        // derivative w.r.t. npVec[i  ]:
        mat[Prow][li_Pcolarray[i][1]] =
          -(dJpdp1Vec[i] - dJpdp2Vec[i-1])/aveDx - dRdpVec[i];

        // derivative w.r.t. npVec[i+1]:
        mat[Prow][li_Pcolarray[i][2]] = -dJpdp2Vec[i]/aveDx;

        // derivative w.r.t.  VVec[i-1]:
        mat[Prow][li_Pcolarray[i][3]] = (dJpdV1Vec[i-1]/aveDx);

        // derivative w.r.t.  VVec[i  ]:
        mat[Prow][li_Pcolarray[i][4]] =
          -(dJpdV1Vec[i] - dJpdV2Vec[i-1])/aveDx;

        // derivative w.r.t.  VVec[i+1]:
        mat[Prow][li_Pcolarray[i][5]] = (-dJpdV2Vec[i]/aveDx);

        // derivative w.r.t. nnVec[i  ]:
        //mat[Prow][li_Pcolarray[i][6]] = -dRdnVec[i];

        // derivative w.r.t. nnVec[i-1]:
        mat[Prow][li_Pcolarray[i][6]] = dJpdn1Vec[i-1]/aveDx;

        // derivative w.r.t. nnVec[i  ]:
        mat[Prow][li_Pcolarray[i][7]] =
          -(dJpdn1Vec[i] - dJpdn2Vec[i-1])/aveDx -dRdnVec[i];

        // derivative w.r.t. nnVec[i+1]:
        mat[Prow][li_Pcolarray[i][8]] = -dJpdn2Vec[i]/aveDx;
      }
      else if (type=="ptype")  // boundary condition on h+, let e- flow
      {
        mat[Prow][li_Pcolarray[i][1]] = 1.0;

        double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

        // derivative w.r.t. nnVec[i-1]:
        mat[Nrow][li_Ncolarray[i][0]] = -dJndn1Vec[i-1]/aveDx;

        // derivative w.r.t. nnVec[i  ]:
        mat[Nrow][li_Ncolarray[i][1]] =
          (dJndn1Vec[i] - dJndn2Vec[i-1])/aveDx - dRdnVec[i];

        // derivative w.r.t. nnVec[i+1]:
        mat[Nrow][li_Ncolarray[i][2]] = dJndn2Vec[i]/aveDx;

        // derivative w.r.t.  VVec[i-1]:
        mat[Nrow][li_Ncolarray[i][3]] = (-dJndV1Vec[i-1]/aveDx);

        // derivative w.r.t.  VVec[i  ]:
        mat[Nrow][li_Ncolarray[i][4]] =
          (dJndV1Vec[i] - dJndV2Vec[i-1])/aveDx;

        // derivative w.r.t.  VVec[i+1]:
        mat[Nrow][li_Ncolarray[i][5]] = (dJndV2Vec[i]/aveDx);

        // derivative w.r.t. npVec[i-1]:
        mat[Nrow][li_Ncolarray[i][6]] = dJndp1Vec[i-1]/aveDx;

        // derivative w.r.t. npVec[i  ]:
        mat[Nrow][li_Ncolarray[i][7]] =
          -(dJndp1Vec[i] - dJndp2Vec[i-1])/aveDx - dRdpVec[i];

        // derivative w.r.t. npVec[i+1]:
        mat[Nrow][li_Ncolarray[i][8]] = -dJndp2Vec[i]/aveDx;
      }
      else
      {
        DevelFatal(*this).in("Instance::loadMatDDForm")
          << "Unrecognized type on boundary.";
      }

    }
  }

  // load the rows associated with the PDE mesh:
  for (i=0;i<NX;++i)
  {
    if (boundarySten[i]==1) continue;
    if ( heterojunctionSten[i]!=0 ) continue;

    ExtendedString semi = bulkMaterial; semi.toLower();

    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    // poisson equation row -------------------------------------
    double dx1 = dxVec[i-1];
    double dx2 = dxVec[i];

    double L0 = scalingVars.L0 * MaterialSupport::getRelPerm(semi);
    //double L0 = scalingVars.L0 * relPermVec[i];

    // del^2 V elements:
    *(fVmatPtr[i][0])=(-L0/(dx1*dx2));
    *(fVmatPtr[i][1])=(2.0*L0/(dx1*dx2));
    *(fVmatPtr[i][2])=(-L0/(dx1*dx2));

    //double dfdn = q/eps;
    // electron density dependency:
    *(fVmatPtr[i][3]) = +1.0;

    // hole density dependency:
    *(fVmatPtr[i][4]) = -1.0;

    // electron continuity row -------------------------------------
    double aveDx = 0.5*(dxVec[i-1]+dxVec[i]);

    // derivative w.r.t. nnVec[i-1]:
    *(fNmatPtr[i][0]) = -dJndn1Vec[i-1]/aveDx;

    // derivative w.r.t. nnVec[i  ]:
    *(fNmatPtr[i][1]) =
      (dJndn1Vec[i] - dJndn2Vec[i-1])/aveDx - dRdnVec[i];


    // derivative w.r.t. nnVec[i+1]:
    *(fNmatPtr[i][2]) = dJndn2Vec[i]/aveDx;

    // derivative w.r.t.  VVec[i-1]:
    *(fNmatPtr[i][3]) = (-dJndV1Vec[i-1]/aveDx);

    // derivative w.r.t.  VVec[i  ]:
    *(fNmatPtr[i][4]) =
      (dJndV1Vec[i] - dJndV2Vec[i-1])/aveDx;

    // derivative w.r.t.  VVec[i+1]:
    *(fNmatPtr[i][5]) = (dJndV2Vec[i]/aveDx);

    // derivative w.r.t. npVec[i-1]:
    *(fNmatPtr[i][6]) = -dJndp1Vec[i-1]/aveDx;

    // derivative w.r.t. npVec[i  ]:
    *(fNmatPtr[i][7]) =
      (dJndp1Vec[i] - dJndp2Vec[i-1])/aveDx -dRdpVec[i];

    // derivative w.r.t. npVec[i+1]:
    *(fNmatPtr[i][8]) = dJndp2Vec[i]/aveDx;

    // hole continuity row -------------------------------------

    // derivative w.r.t. npVec[i-1]:
    *(fPmatPtr[i][0]) = dJpdp1Vec[i-1]/aveDx;

    // derivative w.r.t. npVec[i  ]:
    *(fPmatPtr[i][1]) =
      -(dJpdp1Vec[i] - dJpdp2Vec[i-1])/aveDx - dRdpVec[i];

    // derivative w.r.t. npVec[i+1]:
    *(fPmatPtr[i][2]) = -dJpdp2Vec[i]/aveDx;

    // derivative w.r.t.  VVec[i-1]:
    *(fPmatPtr[i][3]) = (dJpdV1Vec[i-1]/aveDx);

    // derivative w.r.t.  VVec[i  ]:
    *(fPmatPtr[i][4]) =
      -(dJpdV1Vec[i] - dJpdV2Vec[i-1])/aveDx;

    // derivative w.r.t.  VVec[i+1]:
    *(fPmatPtr[i][5]) = (-dJpdV2Vec[i]/aveDx);

    // derivative w.r.t. nnVec[i-1]:
    *(fPmatPtr[i][6]) = -dJpdn1Vec[i-1]/aveDx;

    // derivative w.r.t. nnVec[i  ]:
    *(fPmatPtr[i][7]) =
      (dJpdn1Vec[i] - dJpdn2Vec[i-1])/aveDx -dRdnVec[i];

    // derivative w.r.t. nnVec[i+1]:
    *(fPmatPtr[i][8]) = dJpdn2Vec[i]/aveDx;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcLifetimes
// Purpose       : This function calculates the electron and hole lifetimes
//                 and places them into the tn and tp arrays.
// Special Notes : This function assumes scaling off.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcLifetimes ()
{
  for (int i=0;i<NX;++i)
  {
    tnVec[i] = MaterialSupport::calcLt (false, fabs(CVec[i]));
    tpVec[i] = MaterialSupport::calcLt (true , fabs(CVec[i]));
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcMobilities
// Purpose       : This function calculates the electron and hole mobilities
//                 and places them into the un and up arrays.
//
// Special Notes : The mobility functions assume that scaling is off, so
//                 this function unscales and rescales things as needed.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcMobilities ()
{
  bool bsuccess = true;
  int i;

  MobInfo<pdeFadType> mi;
  mi.mobModelName = mobModelName;
  mi.materialName = bulkMaterial;
  mi.fieldDependent = fieldDependentMobility;

  // Now do edge mobilities:
  for (i=0;i<LX;++i)
  {
    // possibly these doping related quantities should be determined
    // using splines instead.
    mi.N = (fabs(CVec[i+1])+fabs(CVec[i]))*0.5;
    mi.N *= ((variablesScaled)?scalingVars.C0:1.0);

    mi.Na = (fabs(CacceptorVec[i])+fabs(CacceptorVec[i+1]))*0.5
            *((variablesScaled)?scalingVars.C0:1.0);

    mi.Nd = (fabs(CdonorVec[i])+fabs(CdonorVec[i+1]))*0.5
            *((variablesScaled)?scalingVars.C0:1.0);

    if (mi.N == 0.0) mi.N = 1.0;

    pdeFadType v1=VVec[i];
    pdeFadType v2=VVec[i+1];
    pdeFadType n1=nnVec[i];
    pdeFadType n2=nnVec[i+1];
    pdeFadType p1=npVec[i];
    pdeFadType p2=npVec[i+1];

    v1.diff(0,6);
    v2.diff(1,6);
    n1.diff(2,6);
    n2.diff(3,6);
    p1.diff(4,6);
    p2.diff(5,6);

    pdeFadType Efield = (-(v2-v1)/dxVec[i]);
    // this is the most consistent, as it relies on the SG approximation for the midpoint density.
    mi.n = fabs(nMidpoint(n1,n2,Efield,dxVec[i],-1))*((variablesScaled)?scalingVars.C0:1.0);
    mi.p = fabs(nMidpoint(p1,p2,Efield,dxVec[i],+1))*((variablesScaled)?scalingVars.C0:1.0);
    mi.epar = fabs(Efield)*((variablesScaled)?scalingVars.E0:1.0);

    //electron mobility
    mi.holeFlag = false;
    unE_Vec[i] = MaterialSupport::calcMob(mi);
    unE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

    // hole mobility
    mi.holeFlag = true;
    upE_Vec[i] = MaterialSupport::calcMob(mi);
    upE_Vec[i] /= ((variablesScaled)?scalingVars.u0:1.0);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "\tunE["<<i<<"]="<<unE_Vec[i];
      Xyce::dout() << "\tupE["<<i<<"]="<<upE_Vec[i];
      Xyce::dout() << std::endl;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
//
// Special Notes : This won't quite work yet, because Ni and the bandgap
//                 functions (up in the material support class) are
//                 not yet really temperature dependent.
//
//     Things that change with temperature:
//
//        thermal voltage (Vt)
//        scaling variable, scalingVars.V0 = Vt
//        intrinsic concentration, Ni.
//        bandgap, Eg.
//        mobilities (updated in real time, so no need to change here).
//        density boundary conditions (depend on Ni)
//        Vequ (may depend on Ni and Eg).
//        defect reactions
//        other??
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  bool bs1 = true;

  if (indicesSetup_) // instance block constructor sets this flag,
    // but default constructor does not.  If it is set,
    // then the device is ready to process this function,
    // but not otherwise.
  {
    Temp = temp_tmp;

    // first un-scale everything, if neccessary:
    if (variablesScaled)
    {
      bs1 = unScaleVariables ();    bsuccess = bsuccess && bs1;
    }

    bs1 = setupMiscConstants ();    bsuccess = bsuccess && bs1;
    bs1 = setupScalingVars ();      bsuccess = bsuccess && bs1;

    bs1 = calcDensityBCs       ();  bsuccess = bsuccess && bs1;
    bs1 = calcVequBCs          ();  bsuccess = bsuccess && bs1;
    bs1 = calcMobilities       ();  bsuccess = bsuccess && bs1;

    if (!variablesScaled)
    {
      bs1 = scaleVariables ();      bsuccess = bsuccess && bs1;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVoltDepDensities
// Purpose       : This function calculates electron and hole densities,
//                 based on the electrostatic potential.  It is only to be
//                 called during the initialization phase.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcVoltDepDensities ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<NX;++i)
  {
    npVec[i] = getVoltDepHoleDens(VminExp , VVec[i], Na);
    nnVec[i] = getVoltDepElecDens(VmaxExp , VVec[i], Nd);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDopingProfile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/28/07
//-----------------------------------------------------------------------------
bool Instance::setupDopingProfile ()
{
  bool bsuccess (false);
  bool fromFile (false);
  int i;

  if (useLayerCompositeDoping_)
  {
    return true;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Instance::setupDopingProfile\n";
  }

  if ( dopingFileName != "NOFILE" )
  {
    DopeInfo::readDopingFile (dopingFileName, xloc_ndope_vec, ndope_vec, pdope_vec);
    xloc_pdope_vec.clear();
    xloc_pdope_vec.resize( xloc_ndope_vec.size(), 0.0);
    xloc_pdope_vec = xloc_ndope_vec ;
    ndopeInterpolator.clear(); ndopeInterpolator.init(xloc_ndope_vec, ndope_vec);
    pdopeInterpolator.clear(); pdopeInterpolator.init(xloc_pdope_vec, pdope_vec);
    bsuccess=true;
    fromFile=true;
  }
  else if ( ( ( ndopeFileName != "NOFILE" ) && ( pdopeFileName != "NOFILE") ) )
  {
    DopeInfo::readDopingFile (ndopeFileName, xloc_ndope_vec, ndope_vec);
    DopeInfo::readDopingFile (pdopeFileName, xloc_pdope_vec, pdope_vec);
    ndopeInterpolator.clear(); ndopeInterpolator.init(xloc_ndope_vec, ndope_vec);
    pdopeInterpolator.clear(); pdopeInterpolator.init(xloc_pdope_vec, pdope_vec);
    bsuccess=true;
    fromFile=true;
  }
  else
  {
    bsuccess = calcDopingProfile ();
  }


  // use the N and P dopants to create the C vector.
  if (fromFile)
  {
    Na = 0.0;
    Nd = 0.0;
    for (i=0;i<NX;++i)
    {
      double xtmp = xVec[i];
      double ndopeDopeValue(0.0), pdopeDopeValue(0.0);

      ndopeInterpolator.eval(xloc_ndope_vec, ndope_vec, xtmp, ndopeDopeValue);
      pdopeInterpolator.eval(xloc_pdope_vec, pdope_vec, xtmp, pdopeDopeValue);

      CVec[i] = ndopeDopeValue-pdopeDopeValue;

      if (Na > CVec[i]) Na = CVec[i];
      if (Nd < CVec[i]) Nd = CVec[i];
    }
    Na = fabs(Na);
    Nd = fabs(Nd);
  }

  // now that we have the C vector, loop over the boundary
  // conditions and determine if n-type or p-type.

  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int i = bcVec[iBC].meshIndex;
    if (CVec[i] > 0.0)
    {
      bcVec[iBC].type = "ntype";
    }
    else
    {
      bcVec[iBC].type = "ptype";
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Na = " << Na << std::endl;
    Xyce::dout() << "Nd = " << Nd << std::endl;
    for (i=0;i<NX;++i)
    {
      Xyce::dout() << "x[" << i << "] = " << xVec[i] << "\t";
      Xyce::dout() << "C[" << i << "] = " << CVec[i] << std::endl;
    }

    Xyce::dout() << section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDopingProfile
//
// Purpose       : This function sets up the initial doping profile.
//
// Special Notes : 03/31/03. This function is being modified to handle a
//                 more general doping specification.  The old way of
//                 specifying doping, which assumes a PN junction, with Na
//                 and Nd may eventually be phased out.  For now, both
//                 methods of specification are supported, with the new
//                 style over-riding the old, in the event that both are
//                 specified.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcDopingProfile ()
{
  bool bsuccess = true;
  double midpoint;
  int i;

  // Check which style of doping specification to use.  If the dopeInfoMap
  // is empty, then assume the old method.  If not, use the dopeInfoMap.

  if (dopeInfoMap.empty ())
  {
    // first setup, check the graded junction parameters:
    if (gradedJunctionFlag)
    {
      // if junction width was not specified, set to 1/10 the diode width.
      if (!given("WJ"))
        WJ = 0.1 * width;

      midpoint = width/2.0;
      XL = midpoint - WJ/2.0;
      XR = midpoint + WJ/2.0;
    }

    for (i=0;i<NX;++i)
    {
      if (gradedJunctionFlag)
      {
        if (xVec[i] <= XL) CVec[i] = +Nd;
        else if (xVec[i]>XL && xVec[i]<XR)
          CVec[i] = Nd-(Na+Nd)/(XR-XL)*(xVec[i]-XL);
        else CVec[i] = -Na;
      }
      else
      {
        if (xVec[i] < xVec[LX]/2.0) 
        {
          CVec[i] = Nd;
        }
        else 
        {
          CVec[i] = -Na;
        }
      }
    }
  }
  else
  {
    // loop over the dope info map, and sum contributions from each
    // doping entity into the total doping array, CVec.
    std::map<std::string, DopeInfo *>::iterator iter;
    std::map<std::string, DopeInfo *>::iterator start = dopeInfoMap.begin();
    std::map<std::string, DopeInfo *>::iterator end   = dopeInfoMap.end();

    for ( iter = start; iter != end; ++iter )
    {
      DopeInfo & di = *(iter->second);
      di.setupInfo(CVec,CdonorVec,CacceptorVec,xVec);
    }

    Na = 0.0;
    Nd = 0.0;
    for (i=0;i<NX;++i)
    {
      if (Na > CVec[i]) Na = CVec[i];
      if (Nd < CVec[i]) Nd = CVec[i];
    }
    Na = fabs(Na);
    Nd = fabs(Nd);

  } // if statement

  if (Na == 0.0 || Nd == 0.0)
  {
    UserError(*this) << "Mistake in doping. Na=" << Na << " Nd=" << Nd;
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDefaultLayer
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/8/2014
//-----------------------------------------------------------------------------
bool Instance::setupDefaultLayer ()
{
  // if the layer specification was not used, then create one.
  if (!layerCompositeSpecified)
  {
    MaterialLayer *matPtr = new MaterialLayer(bulkMaterial);

    matPtr->materialGiven=true;
    matPtr->name="FULLDOMAIN";
    matPtr->nameGiven=true;
    matPtr->begin=0;
    matPtr->end=NX;
    matPtr->NX=NX;
    matPtr->NXGiven=NXGiven;
    matPtr->width=width;
    matPtr->widthGiven=widthGiven;

    matPtr->processParams();

    materialVec.resize(1, matPtr);
    layerCompositeSpecified=true;
    useLayerCompositeDoping_=false;
  }
  else
  {
    useLayerCompositeDoping_=true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMesh
// Purpose       : This function sets up the mesh.  Should only be called once.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::setupMesh ()
{
  if(layerCompositeSpecified)
  {
    int matVecSize=materialVec.size();
    int totalMeshIndex = 1;
    for (int imat=0;imat<matVecSize;++imat)
    {
      MaterialLayer & matLay = *(materialVec[imat]);

      matLay.LX = matLay.NX-1;
      matLay.begin = totalMeshIndex-1;
      matLay.end = matLay.begin + matLay.LX;

      heterojunctionSten[matLay.begin] = 1;
      heterojunctionSten[matLay.end] = 1;

      double dx = matLay.width /(static_cast<double>(matLay.LX));

      if (imat>0)
      {
        MaterialLayer & matLayPrev = *(materialVec[imat-1]);
        xVec[matLay.begin] = xVec[matLayPrev.end];
        Xyce::dout() << "Setting xVec["<<matLay.begin<<"] to "<<xVec[matLay.begin] <<"."<<std::endl; 

        std::pair<int,int> hetPair = std::make_pair (matLayPrev.end, matLay.begin);
        heterojunctionBCs.push_back(hetPair);
      }

      double base=xVec[matLay.begin];
      int iloc=0;
      for (int ix=matLay.begin;ix<=matLay.end;++ix,++totalMeshIndex,++iloc)
      {
        double extra=static_cast<double>(iloc)*dx;
        xVec[ix] = base + extra;
        matIndex[ix] = imat;
      }
      for (int ix=matLay.begin;ix<matLay.end;++ix)
      {
        dxVec[ix] = xVec[ix+1]-xVec[ix];
      }
    }
    dxVec[LX] = dxVec[LX-1];

    // un-set the heterojunction stencil at the boundaries, as these are more
    // appropriately covered by the boundarySten stencil.
    heterojunctionSten[0] = 0;
    heterojunctionSten[LX] = 0;
  }
  else
  {
    // set up an evenly spaced mesh:
    double dx_tmp = width/(static_cast<double>(LX));

    for (int i=0;i<NX;++i)
    {
      xVec[i] = static_cast<double>(i)*dx_tmp;
    }

    for (int i=0;i<LX;++i)
    {
      dxVec[i] = xVec[i+1]-xVec[i];
    }
    dxVec[LX] = dxVec[LX-1];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for (int i=0;i<NX;++i)
    {
      Xyce::dout() << "x["<<i<<"] = " << xVec[i];
      Xyce::dout() << "\tdx["<<i<<"] = " << dxVec[i];
      Xyce::dout() << std::endl;
    }

    Xyce::dout() << "heterojunction boundary points:" <<std::endl;
    for (int i=0;i<heterojunctionBCs.size();++i)
    {
      Xyce::dout() << "(" 
        << heterojunctionBCs[i].first << "," 
        << heterojunctionBCs[i].second << ")"  <<std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMaterialArrays
// Purpose       : sets up material arrays that are owned by the instance.
//
// Special Notes : If it is temperature-dependent, it is in the instance.
//                 (temperature is an instance variable).  If not, it
//                 is in the model.
//
// This function calculates the density of states arrays.
// For reference, see p32-33 of "Fundamentals of III-V Devices"
// by William Liu.
//
// Equation 1-23 from Liu (conduction band DOS):
//
//    N_c = 2 \left(\frac{2 \pi k T}{h^2}\right)^{3/2} m_{de}^{*{3/2}}
//
//    where m_{de} is the DOS effective mass for electrons.
//
// Equation 1-27 from Liu (valance band DOS):
//
//    N_v = 2 \left(\frac{2 \pi k T}{h^2}\right)^{3/2} (m_{lh}^{*{3/2}} + m_{hh}^{*{3/2}})
//
//    where m_{lh} is the DOS effective mass for light holes
//    where m_{hh} is the DOS effective mass for light holes
//
//    According to Liu, the valance band DOS is dependent on two different
//    degenerate energy bands, each of which has its own hole effective mass.
//
// Scope         : public
// Creation Date :
// Creator       : Eric R. Keiter, SNL
//-----------------------------------------------------------------------------
bool Instance::setupMaterialArrays ()
{
  double dnbnd0 = 2.0*M_PI*e_mass*kb*Temp/pow(h_planck,2.0);
  dnbnd0 = 2.0*pow(dnbnd0,1.5)/1.0e6;
  double kbq = 8.6173324e-5; // boltzmann's constant in eV K^-1


  int size = materialVec.size();
  for (int im=0;im<size; ++im)
  {
    MaterialLayer & matLayer = *(materialVec[im]);

    int begin = matLayer.begin;
    int end = matLayer.end;

    double Ec = matLayer.Ec;
    double Ev = matLayer.Ev;
    double bg = fabs(Ec - Ev);

    double Nc = dnbnd0*(matLayer.dnco);
    double Nv = dnbnd0*(matLayer.dnva);

    double EcEff = matLayer.Ec-matLayer.narco;
    double EvEff = matLayer.Ev+matLayer.narva;
    double bgEff = fabs(EcEff-EvEff);

    // Should this Ni include band-gap narrowing correction?  Generally YES.
    double Ni = sqrt( Nc * Nv ) * exp (-bg/(2.0*kbq*Temp));
    double NiEff = sqrt( Nc * Nv ) * exp (-bgEff/(2.0*kbq*Temp));

    matLayer.Ni = Ni;
    matLayer.NiEff = NiEff;
    matLayer.bg = bg;
    matLayer.bgEff = bgEff;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "layer="<<matLayer.name<< "\t"<<matLayer.material;
      Xyce::dout() << "\tNc="<<Nc<<"\tNv="<<Nv<<"\tbg="<<bg<<"\tNi="<<Ni<<std::endl;
    }

    for (int i=begin;i<=end;++i)
    {
      relPermVec[i] = matLayer.diel;
      EcVec[i]     = matLayer.Ec;
      EvVec[i]     = matLayer.Ev;
      EcEffVec[i]  = matLayer.EcEff;
      EvEffVec[i]  = matLayer.EvEff;

      bgnCVec[i]    = matLayer.narco;
      bgnVVec[i]    = matLayer.narva;

      if (useLayerCompositeDoping_)
      {
        CdonorVec[i]  = matLayer.Cdonor;
        CacceptorVec[i] = matLayer.Cacceptor;
        CVec[i] = matLayer.Cdonor - matLayer.Cacceptor;
      }

      NiVec[i]      = matLayer.Ni;
      NiEffVec[i]   = matLayer.NiEff;

      EiVec[i]      = 0.5*(Ec+Ev)+0.5*kb*Temp*log(Nv/Nc);
      EiEffVec[i]   = 0.5*(EcEff+EvEff)+0.5*kb*Temp*log(Nv/Nc);
      EfVec[i]      = 1.0;
      EfEffVec[i]   = 1.0;

      bulkMaterialVec[i] = matLayer.name;
    }
  }

  if (useLayerCompositeDoping_)
  {
    Na = 0.0;
    Nd = 0.0;
    for (int i=0;i<NX;++i)
    {
      if (Na > CVec[i]) Na = CVec[i];
      if (Nd < CVec[i]) Nd = CVec[i];
    }
    Na = fabs(Na);
    Nd = fabs(Nd);
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for (int im=0;im<size; ++im)
    {
      MaterialLayer & matLayer = *(materialVec[im]);

      int begin = matLayer.begin;
      int end = matLayer.end;

      Xyce::dout() << matLayer.name << "\tbegin="<<begin<<"\tend="<<end<<std::endl;
    }

    for (int i=0;i<NiVec.size();++i)
    {
      Xyce::dout() << i <<"\t" 
                   << EcVec[i] << "\t"
                   << EvVec[i] << "\t"
                   << bgnCVec[i] << "\t"
                   << bgnVVec[i] << "\t"
                   << NcVec[i] << "\t"
                   << NvVec[i] << "\t"
                   << NiVec[i] << "\n";
    }
    Xyce::dout() << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMiscConstants
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/20/03
//-----------------------------------------------------------------------------
bool Instance::setupMiscConstants ()
{
  if (useOldNi)
  {
    Ni = MaterialSupport::getNi_old (bulkMaterial, Temp); // this is not accurate
  }
  else
  {
    Ni = MaterialSupport::getNi (bulkMaterial, Temp);
  }
  Vt = kb*Temp/charge;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupScalingVars
// Purpose       : This function sets up scaling variables.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/01
//-----------------------------------------------------------------------------
bool Instance::setupScalingVars ()
{
  bool bsuccess = true;

  // just to make sure:
  Vt = kb*Temp/charge;

  if (given("X0"))
  {
    scalingVars.x0 = x0_user;
  }
  else
  {
    scalingVars.x0  = width;// distance scaling (cm)
  }

  // For the 1D device, cross-sectional area is kind of a weird concept.
  // For the equations within the device, area really doesn't factor into
  // the discretization.  It is only important at the electrodes, for
  // calculating the current coming out of the device.
  scalingVars.a0 = scalingVars.x0*scalingVars.x0;

  scalingVars.T0  = Temp;// temperature scaling (K)  (not really used)

  // electrostatic potential scaling (V)
  scalingVars.V0 = Vt;
  scalingVars.rV0 = 1.0/scalingVars.V0;

  // concentration scaling (cm^-3);
  if (given("C0"))
  {
    scalingVars.C0 = C0_user;
  }
  else if (scaleDensityToMaxDoping_)  // this is the default
  {
    if (Na >= Nd) scalingVars.C0  = Na;
    else          scalingVars.C0  = Nd;

    scalingVars.C0 *= densityScalarFraction_; // 1.0e-2 by default
  }
  else
  {
    scalingVars.C0 = 1.0e+17;
  }

  if (given("t0"))
  {
    scalingVars.t0 = t0_user;
    scalingVars.D0  = (scalingVars.x0*scalingVars.x0)/scalingVars.t0;
  }
  else
  {
    // diffusion coefficient scaling (cm^2/s)
    scalingVars.D0  = 35.0;

    // time scaling (s)
    scalingVars.t0  = (scalingVars.x0*scalingVars.x0)/scalingVars.D0;
  }

  // mobility coefficient scaling (cm^2/V/s)
  scalingVars.u0  = scalingVars.D0/scalingVars.V0;

  // recombination rate scaling (cm^-3/s)
  scalingVars.R0  = scalingVars.D0*scalingVars.C0/(scalingVars.x0*scalingVars.x0);
  scalingVars.rR0 = 1.0/scalingVars.R0;

  // electric field scaling (V/cm)
  scalingVars.E0  = scalingVars.V0/scalingVars.x0;

  // particle flux scaling (cm^-2/s)
  scalingVars.F0  = scalingVars.D0*scalingVars.C0/scalingVars.x0;

  // current density scaling (A/cm^2)
  scalingVars.J0  = charge*scalingVars.D0*scalingVars.C0/scalingVars.x0;

  // Laplacian scaling constant
  //scalingVars.L0  = scalingVars.V0*eps/(charge*scalingVars.x0*scalingVars.x0*scalingVars.C0);
  scalingVars.L0  = scalingVars.V0*e0/(charge*scalingVars.x0*scalingVars.x0*scalingVars.C0);  // Laplacian scaling constant

  // rate constant scaling.   k0 = 1/(C0*t0) = cm^3/sec
  scalingVars.rk0 = scalingVars.C0*scalingVars.t0;
  scalingVars.rt0 = 1.0/scalingVars.t0;
  scalingVars.k0 = 1.0/scalingVars.rk0;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "scalingVars.x0 = " << scalingVars.x0 << std::endl;
    Xyce::dout() << "scalingVars.a0 = " << scalingVars.a0 << std::endl;
    Xyce::dout() << "scalingVars.T0 = " << scalingVars.T0 << std::endl;
    Xyce::dout() << "scalingVars.V0 = " << scalingVars.V0 << std::endl;
    Xyce::dout() << "scalingVars.C0 = " << scalingVars.C0 << std::endl;
    Xyce::dout() << "scalingVars.D0 = " << scalingVars.D0 << std::endl;
    Xyce::dout() << "scalingVars.u0 = " << scalingVars.u0 << std::endl;
    Xyce::dout() << "scalingVars.R0 = " << scalingVars.R0 << std::endl;
    Xyce::dout() << "scalingVars.t0 = " << scalingVars.t0 << std::endl;
    Xyce::dout() << "scalingVars.E0 = " << scalingVars.E0 << std::endl;
    Xyce::dout() << "scalingVars.F0 = " << scalingVars.F0 << std::endl;
    Xyce::dout() << "scalingVars.J0 = " << scalingVars.J0 << std::endl;
    Xyce::dout() << "scalingVars.L0 = " << scalingVars.L0 << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::scaleVariables
//
// Purpose       : This function performs scaling on all the relevant variables.
//
// Special Notes : It should only be called at the end of the initial setup.
//                 Calculations done during the course of the calculation are
//                 performed with the assumption of scaling.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/28/01
//-----------------------------------------------------------------------------
bool Instance::scaleVariables ()
{
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na /= scalingVars.C0;
  Nd /= scalingVars.C0;
  Ni /= scalingVars.C0;

  int i;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].Vbc  /= scalingVars.V0;
    bcVec[i].Vequ /= scalingVars.V0;
    bcVec[i].nnbc /= scalingVars.C0;
    bcVec[i].npbc /= scalingVars.C0;

    bcVec[i].area /= scalingVars.a0;
  }
  area /= scalingVars.a0;

  VminExp    /= scalingVars.V0;
  VmaxExp    /= scalingVars.V0;

  maxVoltDelta /= scalingVars.V0;

  for (i=0;i<NX;++i)
  {
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
    CVec[i]  /= scalingVars.C0;
    CdonorVec[i] /= scalingVars.C0;
    CacceptorVec[i] /= scalingVars.C0;
    VVec[i]  /= scalingVars.V0;
    tnVec[i] /= scalingVars.t0;
    tpVec[i] /= scalingVars.t0;
    xVec[i]  /= scalingVars.x0;
    dxVec[i] /= scalingVars.x0;

    (*solVectorPtr)[li_Vrowarray[i]] /= scalingVars.V0;
    (*solVectorPtr)[li_Nrowarray[i]] /= scalingVars.C0;
    (*solVectorPtr)[li_Prowarray[i]] /= scalingVars.C0;
  }

  variablesScaled = true;

  return true;

}

//-----------------------------------------------------------------------------
// Function      : Instance::unScaleVariables
// Purpose       : This function is the inverse of scaleVariables.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/20/03
//-----------------------------------------------------------------------------
bool Instance::unScaleVariables ()
{
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  Na *= scalingVars.C0;
  Nd *= scalingVars.C0;
  Ni *= scalingVars.C0;

  int i;
  for (i=0;i<bcVec.size();++i)
  {
    bcVec[i].Vbc  *= scalingVars.V0;
    bcVec[i].Vequ *= scalingVars.V0;
    bcVec[i].nnbc *= scalingVars.C0;
    bcVec[i].npbc *= scalingVars.C0;

    bcVec[i].area *= scalingVars.a0;
  }
  area *= scalingVars.a0;

  VminExp    *= scalingVars.V0;
  VmaxExp    *= scalingVars.V0;

  maxVoltDelta *= scalingVars.V0;

  for (i=0;i<NX;++i)
  {
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
    CVec[i]  *= scalingVars.C0;
    CdonorVec[i] *= scalingVars.C0;
    CacceptorVec[i] *= scalingVars.C0;
    VVec[i]  *= scalingVars.V0;
    tnVec[i] *= scalingVars.t0;
    tpVec[i] *= scalingVars.t0;
    xVec[i]  *= scalingVars.x0;
    dxVec[i] *= scalingVars.x0;

    (*solVectorPtr)[li_Vrowarray[i]] *= scalingVars.V0;
    (*solVectorPtr)[li_Nrowarray[i]] *= scalingVars.C0;
    (*solVectorPtr)[li_Prowarray[i]] *= scalingVars.C0;
  }

  variablesScaled = false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcInitialGuess
// Purpose       : This function calculates the initial e-, h+ densties
//                 and the intial voltage.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::calcInitialGuess ()
{
  bool bsuccess = true;
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  // set up an initial guess for nn and np, based on the doping profile,
  // and the equilibrium density expressions.  Place these in the
  // solution vector.
  if (!fermiDiracFlag)
  {
    double tmp;
    double Ci;
    double Cisq, Nisq;
    for (int i=0;i<NX;++i)
    {
      Ci = CVec[i];
      Cisq = Ci*Ci;
      Nisq = Ni*Ni;  // Ni is the intrinsic concentration

      // equilibrium electron concentration:
      tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
      nnVec[i] = ((Ci>=0)?(tmp):(0.0)) + ((Ci<0)?(Nisq/tmp):(0.0));

      // equilibrium hole concentration:
      tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
      npVec[i] = ((Ci<=0)?(tmp):(0.0)) + ((Ci>0)?(Nisq/tmp):(0.0));
    }
  }
  else
  {
    for (int i=0;i<NX;++i)
    {
      // using Fermi-Dirac statistics
      double kbq_ = 8.6173324e-5; // boltzmann's constant in eV K^-1
      double dope_us_ = CVec[i];
      double temp_us_ = Temp;

      double Nc = NcVec[i];
      double Nv = NvVec[i];

      double cond_band = EcVec[i];
      double vale_band = EvVec[i];
      double bandgap = cond_band-vale_band;

      // N-type
      if (dope_us_ >= 0.0)
      {
        // Assume n approximately equal to Nd
        nnVec[i] = CdonorVec[i];

        // Get Ef - Ec from the inverse fermi function
        double ef_m_ec_ = kbq_*temp_us_* fdinvObj (dope_us_/Nc);

        // BGN value.  FIX THIS!
        double bgn_ = 3.23e-8 * std::pow(dope_us_, 1.0/3.0);

        // Calculate Ef - Ev via (Ef - Ec + Eg) = (Ef - Ec + Ec - Ev) = (Ef - Ev)
        double ef_m_ev_ = ef_m_ec_ + bandgap-bgn_;

        // Calculate the minority carrier concentration using Boltzmann-style approx.
        npVec[i] = Nv*std::exp(-ef_m_ev_/(kbq_*temp_us_));
      }
      // P-type
      else
      {
        dope_us_ = std::fabs(dope_us_);
        npVec[i] = CacceptorVec[i];

        // Get Ev - Ef from the inverse fermi function
        double ev_m_ef_ = kbq_*temp_us_* fdinvObj (dope_us_/Nv);

        // BGN value
        double bgn_ = 2.55e-8 * std::pow(dope_us_, 1.0/3.0);

        // Calculate Ec - Ef via (Ev - Ef + Eg) = (Ev - Ef + Ec - Ev) = (Ec - Ef)
        double ec_m_ef_ = ev_m_ef_ + (bandgap-bgn_);

        // Calculate the minority carrier concentration using Boltzmann-style approx.
        nnVec[i] = Nc*std::exp(-ec_m_ef_/(kbq_*temp_us_));
      }
    }
  }

  // set up initial guess for V, place in solution vector
  double Vmax = -1.0e99;
  double Vmin = +1.0e99;
  for (int i=0;i<NX;++i)
  {
    // the doping is n-type.
    if (nnVec[i]>=npVec[i])
    {
      VVec[i] = + Vt * log(nnVec[i]/Ni);
    }
    // the doping is p-type.
    else
    {
      VVec[i] = - Vt * log(npVec[i]/Ni);
    }

    if (Vmax < VVec[i]) Vmax = VVec[i];
    if (Vmin > VVec[i]) Vmin = VVec[i];
  }

  // get the maximum and minimum potentials.
  VmaxExp = -1.0e99;
  VminExp = +1.0e99;

  for (int i=0;i<NX;++i)
  {
    if (VmaxExp < VVec[i]) VmaxExp = VVec[i];
    if (VminExp > VVec[i]) VminExp = VVec[i];
  }

  for (int i=0;i<NX;++i)
  {
    (*solVectorPtr)[li_Vrowarray[i]] = VVec[i];
    (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
    (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcVequBCs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/17/03
//-----------------------------------------------------------------------------
bool Instance::calcVequBCs ()
{
  bool bsuccess = true;
  Vt = kb*Temp/charge;

  double VminBC =+1.0e+99;
  double VmaxBC =-1.0e+99;

  int bcSize=bcVec.size();
  for (int i=0;i<bcSize;++i)
  {
    int mIndex = bcVec[i].meshIndex;
    double Ci = CVec[mIndex];
    double Cisq = Ci*Ci;
    double Nisq = Ni*Ni;  // Ni is the intrinsic concentration
    double tmp, nnTmp, npTmp;

    // equilibrium electron concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    nnTmp = ((Ci>=0)?(tmp):(0.0)) + ((Ci<0)?(Nisq/tmp):(0.0));

    // equilibrium hole concentration:
    tmp = (fabs(Ci)+sqrt(Cisq+4*Nisq))/2.0;
    npTmp = ((Ci<=0)?(tmp):(0.0)) + ((Ci>0)?(Nisq/tmp):(0.0));

    //ExtendedString bulkMat = bulkMaterialVec[mIndex];
    ExtendedString bulkMat = bulkMaterial;
    bulkMat.toLower();
    ExtendedString mater = bcVec[i].material;
    mater.toLower();

    if (bcVec[i].VequGiven != 1)
    {
      if (mater=="neutral")
      {
        // the doping is n-type.
        if (Ci>0)
        {
          bcVec[i].Vequ = + Vt * log(nnTmp/Ni);
        }
        else        // the doping is p-type.
        {
          bcVec[i].Vequ = - Vt * log(npTmp/Ni);
        }
      }
      else // this electrode is a schottky barrier.
      {
        // the doping is n-type.
        if (Ci>0)
        {
          bcVec[i].Vequ = + MaterialSupport::workfunc(mater)
                          - MaterialSupport::affin(bulkMat)
                          - 0.5 * MaterialSupport::bandgap(bulkMat, Temp)
                          + 2.0 * Vt * log(nnTmp/Ni);
        }
        else        // the doping is p-type.
        {
          bcVec[i].Vequ = + MaterialSupport::workfunc(mater)
                          - MaterialSupport::affin(bulkMat)
                          - 0.5 * MaterialSupport::bandgap(bulkMat, Temp)
                          - 2.0 * Vt * log(npTmp/Ni);
        }
      }
    }

    if (VminBC > bcVec[i].Vequ) VminBC = bcVec[i].Vequ;
    if (VmaxBC < bcVec[i].Vequ) VmaxBC = bcVec[i].Vequ;
  }

  VoltageOffset_ = -VminBC;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcDensityBCs
// Purpose       : This function sets up the boundary condition variables
//                 for each electrode.
//
// Special Notes : This function is similar to calcBoundaryConditions, but
//                 this one only calculates BC's on N and P.  Since these
//                 never change, they only need to be calculated once.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/03/03
//-----------------------------------------------------------------------------
bool Instance::calcDensityBCs ()
{
  bool bsuccess = true;

  NnMax = -1.0e+99;
  NpMax = -1.0e+99;

  NnMin = +1.0e+99;
  NpMin = +1.0e+99;

  // This density boundary condition is from Selberherr,
  // enforcing thermal equilibrium and
  // vanishing space charge at ohmic contacts
  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    int i1 = bcVec[iBC].meshIndex;
    bcVec[iBC].nnbc = 0.5*(sqrt(CVec[i1]*CVec[i1]+4*Ni*Ni)+CVec[i1]);
    bcVec[iBC].npbc = 0.5*(sqrt(CVec[i1]*CVec[i1]+4*Ni*Ni)-CVec[i1]);

    if (NnMax < bcVec[iBC].nnbc) NnMax = bcVec[iBC].nnbc;
    if (NpMax < bcVec[iBC].npbc) NpMax = bcVec[iBC].npbc;
  }

  if (NnMin <= 0) NnMin = 1.56269e-9;  // just a guess.
  if (NpMin <= 0) NpMin = 1.56269e-9;  // just a guess.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcBoundaryConditions
// Purpose       : This function sets up the boundary condition variables
//                 for each electrode.
//
// Special Notes : If a continuation method is being used, a good parameter
//                 to vary is the voltage boundary condition  on the
//                 electrodes.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Instance::calcBoundaryConditions ()
{
  bool bsuccess = true;
  int iBC;

  int bcSize=bcVec.size();
  if (getSolverState().PDEcontinuationFlag_)
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vbc = bcVec[iBC].Vckt_ramp + bcVec[iBC].Vequ;
    }
  }
  else
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vbc = (bcVec[iBC].Vckt + bcVec[iBC].Vequ);
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::obtainNodeVoltages.
//
// Purpose       : This function obtains the nodal voltages from the
//                 solution vector, which are applied as boundary
//                 conditions on the electrodes.
//
// Special Notes : This was originally part of function obtainSolution, but
//                 is needed also by function enableContinuation.  So I've
//                 put it in one place.
//
//                 If voltage limiting is turned on, this is the function
//                 in which to apply it.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/13/02
//-----------------------------------------------------------------------------
bool Instance::obtainNodeVoltages ()
{
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  int iBC;
  for (iBC=0;iBC<bcVec.size();++iBC)
  {
    bcVec[iBC].Vckt = (*solVectorPtr)[bcVec[iBC].lid];
    bcVec[iBC].Vckt /= scalingVars.V0;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::applyVoltageLimiting
//
// Purpose       : if voltage limiting is turned on, this function
//                 applies it to the Vckt values.
//
// Special Notes : This is only really set up to work when the 2-level
//                 Newton is being used.
//
//                 For now, this is just a test capability.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/02
//-----------------------------------------------------------------------------
bool Instance::applyVoltageLimiting ()
{
  for (int iBC=0;iBC<bcVec.size();++iBC)
  {
    double v1     = bcVec[iBC].Vckt * scalingVars.V0;
    double v1_old = bcVec[iBC].Vckt_old * scalingVars.V0;
    double delV1 = v1 - v1_old;

    if ( delV1 > 1.25 )  v1 = v1_old + 1.25;

    if ( delV1 < -0.75) v1 = v1_old - 0.75;

    bcVec[iBC].Vckt       = v1/scalingVars.V0;
    bcVec[iBC].Vckt_final = v1/scalingVars.V0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::obtainSolution
// Purpose       : This function extracts V, nn, and np from the solution
//                 vector and copies them into local arrays.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::obtainSolution ()
{
  bool bsuccess = true;
  bool bs1 = true;
  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;

  // First get the two circuit node voltages:
  bsuccess = obtainNodeVoltages ();

  // set up the V solution array:
  int i;
  for (i=0;i<NX;++i)
  {
    VVec[i] = (*solVectorPtr)[li_Vrowarray[i]];
  }

  // If the previous solution is from the nonlinear Poisson solution,
  // then calculate what the electron and hole densities must be, and
  // place them into the solution vector.

  // If we are past the nonlinear Poisson phase, then simply obtain
  // nn and np from the solution vector and move on.

  if (getSolverState().dcopFlag && getSolverState().doubleDCOPStep==0)
  {
    calcVoltDepDensities ();

    for (i=0;i<NX;++i)
    {
      (*solVectorPtr)[li_Nrowarray[i]] = nnVec[i];
      (*solVectorPtr)[li_Prowarray[i]] = npVec[i];
    }
  }
  else
  {
    for (i=0;i<NX;++i)
    {
      nnVec[i] = (*solVectorPtr)[li_Nrowarray[i]];

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      if (nnVec[i] < 0.0) nnVec[i] = 0.0;
#endif

      npVec[i] = (*solVectorPtr)[li_Prowarray[i]];

#ifdef Xyce_PDE_DENSITY_CONSTRAINT
      if (npVec[i] < 0.0) npVec[i] = 0.0;
#endif
    }

    // now set boundary conditions:
    // if the circuit is coupled to the PDE device, then bc's
    // must be updated everytime.
    //
    // If the circuit and PDE device are not coupled, then the
    // circuit node voltages can be considered constant, and the
    // BC's only need updating at the first Newton step.
    if ( !(getSolverState().twoLevelNewtonCouplingMode==Nonlinear::INNER_PROBLEM))
    {
      bs1 = calcBoundaryConditions (); bsuccess = bsuccess && bs1;
    }
    else
    {
      if (getSolverState().newtonIter == 0)
      {
        bs1 = calcBoundaryConditions (); bsuccess = bsuccess && bs1;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/22/03
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;
  bool bs1 = true;
  bool skipOutput = false;

  // usually, don't bother outputting nonlinear Poisson result.
  if (equationSet == 0 && !(outputNLPoisson))  return bsuccess;

  // If using output interval, check if enough time has passed to do
  // another output.  (only applies for transient - not DCOP).
  if ( !(getSolverState().dcopFlag) &&
       !(force_final_output) &&
       outputIntervalGiven)
  {
    double outMult = static_cast<double> (outputIndex);
    double nextOutputTime = outMult * outputInterval;

    if (nextOutputTime > getSolverState().currTime_)
    {
      skipOutput = true;
    }
  }

  // If this is a "forced" final output, make sure that it didn't already output.
  // This can happen if the output interval is an exact multiple of the
  // total simulation time.
  if (force_final_output &&
      getSolverState().currTime_==lastOutputTime) skipOutput=true;

  if (skipOutput) return bsuccess;
  ++outputIndex;
  lastOutputTime = getSolverState().currTime_;

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << std::endl << "Doing an output at time = " << getSolverState().currTime_ << std::endl;
  }

  if (tecplotLevel > 0) {bs1 = outputTecplot (); bsuccess = bsuccess && bs1;}
  if (sgplotLevel  > 0) {bs1 = outputSgplot  (); bsuccess = bsuccess && bs1;}

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTecplot
// Purpose       : This function outputs a file which is easily plottable
//                 by tecplot.  Simply run tecplot "filename.dat" <return>
//
//
// Special Notes : This file can also be plotted using gnuplot.  If the
//                 name of the file is "Z1DIODE_000.dat", plot inside of
//                 gnuplot using:
//
//                 plot "Z1DIODE_000.dat" using 1:3 w l
//
//                 or, if you want a log plot, for the doping:
//
//                 plot "Z1DIODE_000.dat" using 1:(log($6)) w l
//
//                 The "$" and both pairs of parens are needed for some
//                 reason.
//
// Special Notes : If tecplot level is set to 1, then output each dataset
//                 in a separate file.  If not, then append to a single file.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::outputTecplot ()
{
  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  if (tecplotLevel == 1)
  {
    sprintf(filename,"%s_%03d.dat",outputName.c_str(),callsOTEC);
  }
  else
  {
    sprintf(filename,"%s.dat",outputName.c_str());
  }

  double time = getSolverState().currTime_;
  FILE *fp1;

  if (tecplotLevel == 1)
  {
    fp1 = fopen(filename,"w");
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fp1 = fopen(filename,"w");
    }
    else
    {
      fp1 = fopen(filename,"a");
    }
  }

  if (tecplotLevel == 1)
  {
    if (equationSet == 0)
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds. equation set = nonlinear Poisson\",\n",
              outputName.c_str(),time);
    }
    else
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds. equation set = drift diffusion\",\n",
              outputName.c_str(),time);
    }
  }
  else
  {
    if (callsOTEC <= 0)
    {
      fprintf(fp1,
              " TITLE = \"Spatially Dependent data for PDE diode: %s  time = %20.12e seconds.\",\n",
              outputName.c_str(),time);
    }
  }

  int rSize=0;
  int cSize=0;

  if (callsOTEC <= 0 || tecplotLevel == 1)
  {
    fprintf(fp1,"%s","\tVARIABLES = \"X \",\n");

    fprintf(fp1,"%s","\t    \"V \",\n");
    fprintf(fp1,"%s","\t    \"nn (electron dens.) \",\n");
    fprintf(fp1,"%s","\t    \"np (hole dens.) \",\n");
    fprintf(fp1,"%s","\t    \"nn*np (carrier product) \",\n");
    fprintf(fp1,"%s","\t    \"Dopant density \",\n");
    fprintf(fp1,"%s","\t    \"fabs(Dopant density)\",\n");
    fprintf(fp1,"%s","\t    \"electron lifetime \",\n");
    fprintf(fp1,"%s","\t    \"hole lifetime \",\n");
    //fprintf(fp1,"%s","\t    \"electron mobility \",\n");
    //fprintf(fp1,"%s","\t    \"hole mobility \",\n");
    fprintf(fp1,"%s","\t    \"Jn \",\n");
    fprintf(fp1,"%s","\t    \"Jp \",\n");
    fprintf(fp1,"%s","\t    \"R  \",\n");
    fprintf(fp1,"%s","\t    \"Ex \",\n");
    fprintf(fp1,"%s","\t    \"Idispl \", \n");

#if 0
    fprintf(fp1,"%s","\t    \"Conduction Band, uncorrected \", \n");
    fprintf(fp1,"%s","\t    \"Valance Band, uncorrected \", \n");

    fprintf(fp1,"%s","\t    \"Band-gap narrowing, Conduction Band \", \n");
    fprintf(fp1,"%s","\t    \"Band-gap narrowing, Valance Band \", \n");

    fprintf(fp1,"%s","\t    \"Conduction Band, corrected for BGN \", \n");
    fprintf(fp1,"%s","\t    \"Valance Band, corrected for BGN \", \n");
    fprintf(fp1,"%s","\t    \"Fermi Level\", \n");

    fprintf(fp1,"%s","\t    \"conduction band DOS\", \n");
    fprintf(fp1,"%s","\t    \"valance band DOS\", \n");

    fprintf(fp1,"\t    \"n0, Fermi-Dirac \",\n");
    fprintf(fp1,"\t    \"p0, Fermi-Dirac \",\n");
    fprintf(fp1,"\t    \"n0, Boltzmann\",\n");
    fprintf(fp1,"\t    \"p0, Boltzmann\",\n");
    fprintf(fp1,"\t    \"np0 Fermi-Dirac\",\n");
    fprintf(fp1,"\t    \"Ni^2 (Boltzmann np0)\",\n");
    fprintf(fp1,"%s","\t    \"Ni (intrinsic concentration) \", \n");
#endif
  }

  fprintf(fp1,"\tZONE F=POINT,I=%d", NX);

  if (getSolverState().dcopFlag)
  {
    fprintf(fp1,"  T = \"DCOP step = %d\" \n", callsOTEC);
  }
  else
  {
    fprintf(fp1,"  T = \"time step = %d time = %20.12e\" AUXDATA time = \"%20.12e seconds\" \n", callsOTEC , time, time);
  }

  double vcorrection = 0.0;
  if (useVoltageOutputOffset_)
  {
    if (offsetWithFirstElectrode_) // not the default.  This is here to match Wampler's 1D code.
    {
      vcorrection = -VVec[0]*scalingVars.V0;
    }
    else
    {
      vcorrection = VoltageOffset_;
    }
  }

  for (i=0;i<NX;++i)
  {
    fprintf(fp1,"  %20.12e",xVec[i]*scalingVars.x0);
    fprintf(fp1,"  %20.12e", (VVec[i]*scalingVars.V0 + vcorrection) );
    fprintf(fp1,"  %20.12e",nnVec[i]*scalingVars.C0);
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",npVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",nnVec[i]*scalingVars.C0*npVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",CVec[i]*scalingVars.C0);
    fprintf(fp1,"  %20.12e",fabs(CVec[i]*scalingVars.C0));
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",tnVec[i]*scalingVars.t0);
    fprintf(fp1,"  %20.12e",tpVec[i]*scalingVars.t0);
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",JnxVec[i]*scalingVars.J0);
    fprintf(fp1,"  %20.12e",JpxVec[i]*scalingVars.J0);
    fprintf(fp1,"%s","\n");
    fprintf(fp1,"  %20.12e",RVec[i]*scalingVars.R0);
    fprintf(fp1,"  %20.12e",ExVec[i]*scalingVars.E0);
    fprintf(fp1,"  %20.12e",displCurrent[i]*scalingVars.J0);
    fprintf(fp1,"%s","\n");

#if 0
    fprintf(fp1,"  %20.12e", EcVec[i]);
    fprintf(fp1,"  %20.12e", EvVec[i]);
    fprintf(fp1,"  %20.12e", bgnCVec[i]);
    fprintf(fp1,"  %20.12e", bgnVVec[i]);

    double con = EcVec[i]-bgnCVec[i];
    double val = EvVec[i]+bgnVVec[i];
    fprintf(fp1,"  %20.12e", con);
    fprintf(fp1,"  %20.12e", val);

    fprintf(fp1,"%s","\n");

    fprintf(fp1,"  %20.12e", EfVec[i]);
    fprintf(fp1,"  %20.12e", NcVec[i]);
    fprintf(fp1,"  %20.12e", NvVec[i]);

    double Ni=NiVec[i];
    double n0,p0;
    n0_and_p0(
        (nnVec[i]*scalingVars.C0), (npVec[i]*scalingVars.C0), 
        Ni, con, val, NcVec[i], NvVec[i], Temp, n0, p0);

    fprintf(fp1,"  %20.12e",n0);
    fprintf(fp1,"  %20.12e",p0);

    if (CdonorVec[i] > CacceptorVec[i])
    {
      n0 = (CdonorVec[i]-CacceptorVec[i]); p0=1.0;
      if (n0 != 0.0) { p0 = Ni*Ni/n0; }
    }
    else
    {
      p0 = (CacceptorVec[i]-CdonorVec[i]); n0=1.0;
      if (p0 != 0.0) { n0 = Ni*Ni/p0; }
    }
    fprintf(fp1,"  %20.12e",n0);
    fprintf(fp1,"  %20.12e",p0);

    double np0 = np0_calculation(
        (nnVec[i]*scalingVars.C0), (npVec[i]*scalingVars.C0), 
        Ni, con, val, NcVec[i], NvVec[i], Temp);

    fprintf(fp1,"  %20.12e",np0);
    fprintf(fp1,"  %20.12e",Ni*Ni);
    fprintf(fp1,"  %20.12e",Ni);
    fprintf(fp1,"\n");

    fprintf(fp1,"%s","\n");
#endif
  }

  ++callsOTEC;
  fclose(fp1);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputSgplot
// Purpose       : Outputs data in a format readable by the simgen plotting
//                 tools.
//
// Special Notes : Despite the name, the output file should be plotted using
//                 the program oneplot(a 1D plotter), not sgplot
//                 (a 2d plotter).
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool Instance::outputSgplot ()
{
  int i;
  char fileName[32];

  static const int LEN_IDENT2 = 31;

  for (i = 0 ; i < 32; ++i)
    fileName[i] = static_cast<char>(0);

  sprintf(fileName,"%s_%03d.res",outputName.c_str(),callsOSG);
  ++callsOSG;

  FILE * handle1 = fopen(fileName, "w");

  UINT numArrays  = 3;
  double timeVar = 0.0;

  UINT inx = NX;

  char title[64];
  sprintf(title, "%s", "Xyce diodePDE 1D output");

  fwrite(&inx      , sizeof(  UINT), 1, handle1);  // array size.
  fwrite(&numArrays, sizeof(  UINT), 1, handle1);  // number of arrays, besides x.
  fwrite( title    , sizeof(  char),64, handle1);  // title
  fwrite(&timeVar  , sizeof(double), 1, handle1);  // time.

  char names[3][LEN_IDENT2];
  sprintf(names[0], "%s", "V");
  sprintf(names[1], "%s", "Ne");
  sprintf(names[2], "%s", "Np");

  // output the variable names, other than x:
  for(i=0;i<numArrays;++i)
  {
    fwrite(names[i], sizeof(char),(LEN_IDENT2), handle1);
  }

  double vcorrection = 0.0;
  if (useVoltageOutputOffset_)
  {
    if (offsetWithFirstElectrode_) // not the default.  This is here to match Wampler's 1D code.
    {
      vcorrection = -VVec[0]*scalingVars.V0;
    }
    else
    {
      vcorrection = VoltageOffset_;
    }
  }

  for (i=0;i<inx;++i)
  {
    xVec[i] *= scalingVars.x0;
    VVec[i] *= scalingVars.V0 + vcorrection;
    nnVec[i] *= scalingVars.C0;
    npVec[i] *= scalingVars.C0;
  }

  // output x-axis:
  fwrite( &xVec[0], sizeof(double),inx , handle1 );

  // output V
  fwrite( &VVec[0], sizeof(double),inx , handle1 );

  // output nn
  fwrite( &nnVec[0], sizeof(double),inx , handle1 );

  // output np
  fwrite( &npVec[0], sizeof(double),inx , handle1 );

  for (i=0;i<inx;++i)
  {
    xVec[i] /= scalingVars.x0;
    VVec[i] /= scalingVars.V0;
    nnVec[i] /= scalingVars.C0;
    npVec[i] /= scalingVars.C0;
  }

  fclose(handle1);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcRecombination
// Purpose       :
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcRecombination ()
{
  if (!includeAugerRecomb && !includeSRHRecomb) return true;

  for (int i=0;i<NX;++i)
  {
    double Rsrh=0.0;
    double Raug=0.0;

    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    if (includeSRHRecomb)
    {
      Rsrh = MaterialSupport::calcRsrh (bulkMaterial, Ni,n,p,tn,tp);
    }

    if (includeAugerRecomb)
    {
      Raug = MaterialSupport::calcRaug (bulkMaterial, Ni,n,p);
    }

    RVec[i] = (Rsrh + Raug);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdRecombination
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the recombination term.
// Special Notes : This function assumes scaling is turned on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdRecombination ()
{
  if (!includeAugerRecomb && !includeSRHRecomb) return true;

  int i;
  double A, B, C;
  double dAdn, dAdp;
  double dBdn, dBdp;

  for (i=0;i<NX;++i)
  {
    double dRsrhdn=0.0;
    double dRsrhdp=0.0;
    double dRaugdn=0.0;
    double dRaugdp=0.0;

    // Rsrh derivatives: checklater.
    // (Rsrch = A*B)

    double n  = nnVec[i];
    double p  = npVec[i];
    double tn = tnVec[i];
    double tp = tpVec[i];

    if (includeSRHRecomb)
    {
      dRsrhdn = MaterialSupport::pdRsrhN(bulkMaterial,Ni,n,p,tn,tp);
      dRsrhdp = MaterialSupport::pdRsrhP(bulkMaterial,Ni,n,p,tn,tp);
    }

    if (includeAugerRecomb)
    {
      dRaugdn = MaterialSupport::pdRaugN(bulkMaterial,Ni,n,p);
      dRaugdp = MaterialSupport::pdRaugP(bulkMaterial,Ni,n,p);
    }

    dRdnVec[i] = dRsrhdn + dRaugdn;
    dRdpVec[i] = dRsrhdp + dRaugdp;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcElectronCurrent
// Purpose       :
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcElectronCurrent ()
{
  Ut = Vt/scalingVars.V0;
  for (int i=0;i<LX;++i)
  {
    if (i>0 && i< LX && heterojunctionSten[i]!=0 && heterojunctionSten[i+1]!=0 )
    {
      JnxVec[i] = JnxVec[i-1]; // kludge for now
    }
    else
    {
      JnxVec[i] =
        -J_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);
    }
  }
  JnxVec[LX] = JnxVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdElectronCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with electron current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdElectronCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    dJndn1Vec[i] =
      -dJdn1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndn2Vec[i] =
      -dJdn2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndV1Vec[i] =
      -dJdV1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndV2Vec[i] =
      -dJdV2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndp1Vec[i] =
      -dJdp1_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);

    dJndp2Vec[i] =
      -dJdp2_qdep(nnVec[i], nnVec[i+1], ExVec[i], unE_Vec[i], dxVec[i],-1);
  }

  dJndn1Vec[LX] = dJndn1Vec[LX-1];
  dJndn2Vec[LX] = dJndn2Vec[LX-1];
  dJndV1Vec[LX] = dJndV1Vec[LX-1];
  dJndV2Vec[LX] = dJndV2Vec[LX-1];
  dJndp1Vec[LX] = dJndp1Vec[LX-1];
  dJndp2Vec[LX] = dJndp2Vec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcHoleCurrent
// Purpose       : This function assumes scaling is on.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcHoleCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    if (i>0 && i< LX && heterojunctionSten[i]!=0 && heterojunctionSten[i+1]!=0 )
    {
      JpxVec[i] = JpxVec[i-1]; //kludge for now
    }
    else
    {
      JpxVec[i] =
        J_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);
    }
  }

  JpxVec[LX] = JpxVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::pdHoleCurrent
// Purpose       : This function sets up the arrays of partial derivatives
//                 associated with the hole current.
// Special Notes : This function assumes scaling is on.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/01
//-----------------------------------------------------------------------------
bool Instance::pdHoleCurrent ()
{
  Ut = Vt/scalingVars.V0;

  for (int i=0;i<LX;++i)
  {
    dJpdp1Vec[i] =
      dJdn1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdp2Vec[i] =
      dJdn2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdV1Vec[i] =
      dJdV1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdV2Vec[i] =
      dJdV2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdn1Vec[i] =
      -dJdp1_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);

    dJpdn2Vec[i] =
      -dJdp2_qdep(npVec[i], npVec[i+1], ExVec[i], upE_Vec[i], dxVec[i],+1);
  }

  dJpdn1Vec[LX] = dJpdn1Vec[LX-1];
  dJpdn2Vec[LX] = dJpdn2Vec[LX-1];
  dJpdV1Vec[LX] = dJpdV1Vec[LX-1];
  dJpdV2Vec[LX] = dJpdV2Vec[LX-1];
  dJpdn1Vec[LX] = dJpdn1Vec[LX-1];
  dJpdn2Vec[LX] = dJpdn2Vec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::calcEfield
// Purpose       : This function works with or without scaling.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/15/01
//-----------------------------------------------------------------------------
bool Instance::calcEfield ()
{
  double absEx;
  Emax = 0.0;

  for (int i=0;i<LX;++i)
  {
    if (i>0 && i< LX && heterojunctionSten[i]!=0 && heterojunctionSten[i+1]!=0 )
    {
      ExVec[i] = ExVec[i-1]; // kludge for now (to avoid nan's or inf's in the output)
    }
    else
    {
      ExVec[i] = -(VVec[i+1] - VVec[i])/dxVec[i];
    }

    absEx = fabs(ExVec[i]);
    if (absEx > Emax) Emax = absEx;
  }
  Emax *= scalingVars.E0;

  ExVec[LX] = ExVec[LX-1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::enablePDEContinuation
// Purpose       : Sets up the various parameters neccessary for a continuation
//                 calculation.  Mainly, it sets up the voltage step size
//                 for all the voltage BC's.
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::enablePDEContinuation(int &max_PDE_continuation_steps)
{
  bool bnoChange = true;
  int iBC;
  int bcSize=bcVec.size();

  continuationAlpha = 1.0;

  if (!enableContinuationCalled)
  {
    for (iBC=0;iBC<bcSize;++iBC)
    {
      bcVec[iBC].Vckt_old = bcVec[iBC].Vckt;
    }
  }

  obtainNodeVoltages ();

  for (iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_final = bcVec[iBC].Vckt;
    bcVec[iBC].Vckt_orig  = bcVec[iBC].Vckt;
  }

  // This (voltlim) is a very new thing.  Use carefully...
  if (getDeviceOptions().voltageLimiterFlag && voltLimFlag)
  {
    applyVoltageLimiting ();
  }

  for (iBC=0;iBC<bcSize;++iBC)
  {
    double dV,tmp1V, tmp2V;
    tmp1V = bcVec[iBC].Vckt_final;
    tmp2V = bcVec[iBC].Vckt_old;
    dV    = tmp1V - tmp2V;

    bcVec[iBC].Vckt_delta = dV;

    bcVec[iBC].Vckt_deltaC = dV/static_cast<double>(max_PDE_continuation_steps);

    // if this deltaC is too big, then we need to change the
    // number of continuation steps.
    double maxDelta = maxVoltDelta;

    if (fabs(bcVec[iBC].Vckt_deltaC) > maxDelta)
    {
      int tmp_steps = static_cast<int>(fabs(dV)/maxDelta) + 1;
      max_PDE_continuation_steps = tmp_steps;

      bcVec[iBC].Vckt_deltaC = dV/static_cast<double>(max_PDE_continuation_steps);
    }

    if (fabs(dV) > 1.0e-3) bnoChange = false;

    bcVec[iBC].Vckt_ramp     = bcVec[iBC].Vckt_old;
    bcVec[iBC].Vckt_ramp_old = bcVec[iBC].Vckt_old;
  }

  if (!enableContinuationCalled) enableContinuationCalled = true;

  // if none of the boundary conditions  have changed, then
  // return a false.
  return (!bnoChange);
}

//-----------------------------------------------------------------------------
// Function      : Instance::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool Instance::disablePDEContinuation ()
{
  int iBC;

  int bcSize=bcVec.size();
  for (iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_old   = bcVec[iBC].Vckt_final;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setPDEContinuationAlpha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/18/03
//-----------------------------------------------------------------------------
void Instance::setPDEContinuationAlpha (double alpha)
{
  if (DEBUG_DEVICE)
  {
    Xyce::dout() << section_divider << std::endl;
  Xyce::dout() << "Instance::setPDEContinuationAlpha" << std::endl;
  }


  // First do the voltage boundary conditions:
  int bcSize=bcVec.size();
  for (int iBC=0;iBC<bcSize;++iBC)
  {
    bcVec[iBC].Vckt_ramp = bcVec[iBC].Vckt_old + (bcVec[iBC].Vckt_delta)*alpha;

    // make sure we haven't gone too far:
    if ((bcVec[iBC].Vckt_delta >  0 && bcVec[iBC].Vckt_ramp >  bcVec[iBC].Vckt_final) ||
        (bcVec[iBC].Vckt_delta <= 0 && bcVec[iBC].Vckt_ramp <= bcVec[iBC].Vckt_final) )
    {
      bcVec[iBC].Vckt_ramp = bcVec[iBC].Vckt_final;
    }

    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "  " << bcVec[iBC].eName << "  Vckt_ramp = " << bcVec[iBC].Vckt_ramp << std::endl;
    }
  }
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  Device *device = new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);

  return device;
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() ||
      ((deviceMap.find("PDE")!=deviceMap.end()) && (levelSet.find(1)!=levelSet.end()))))
  {
    initialized = true;

    Config<Traits>::addConfiguration()
      .registerDevice("pde", 1)
      .registerModelType("pde", 1)
      .registerModelType("zod", 1);
  }
}

} // namespace DiodePDE
} // namespace Device
} // namespace Xyce
