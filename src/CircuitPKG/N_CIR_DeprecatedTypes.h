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
// Purpose        : Deprecated type declaraions
//
// Special Notes  : These should NOT be used, but if you just have to use the old names, here they are.
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2015/05/12
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_DeprecatedTypes_h
#define Xyce_N_CIR_DeprecatedTypes_h

#include <N_ANP_fwd.h>
#include <N_CIR_fwd.h>
#include <N_DEV_fwd.h>
#include <N_ERH_fwd.h>
#include <N_IO_Measure_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

typedef Xyce::Analysis::OutputMgrAdapter N_ANP_OutputMgrAdapter;
typedef Xyce::Analysis::SweepParam N_ANP_SweepParam;
typedef Xyce::Analysis::AC N_ANP_AC;
typedef Xyce::Analysis::NOISE N_ANP_NOISE;
typedef Xyce::Analysis::NoiseData N_ANP_NoiseData;
typedef Xyce::Analysis::AnalysisBase N_ANP_AnalysisBase;
typedef Xyce::Analysis::AnalysisManager N_ANP_AnalysisManager;
typedef Xyce::Analysis::Transient N_ANP_Transient;
typedef Xyce::Analysis::Step N_ANP_Step;
typedef Xyce::Analysis::DCSweep N_ANP_DCSweep;
typedef Xyce::Analysis::MPDE N_ANP_MPDE;
typedef Xyce::Analysis::HB N_ANP_HB;
typedef Xyce::Analysis::Dakota N_ANP_Dakota;
typedef Xyce::Analysis::MOR N_ANP_MOR;
typedef Xyce::Circuit::Simulator N_CIR_Xyce;
typedef Xyce::Circuit::Simulator N_CIR_Xyce;
typedef Xyce::Device::Depend Depend;
typedef Xyce::Device::Depend Depend;
typedef Xyce::Device::DeviceBuilder N_DEV_DeviceBuilder;
typedef Xyce::Device::DeviceEntity N_DEV_DeviceEntity;
typedef Xyce::Device::DeviceInstance N_DEV_DeviceInstance;
typedef Xyce::Device::DeviceInterface N_DEV_DeviceInterface;
typedef Xyce::Device::DeviceMgr N_DEV_DeviceMgr;
typedef Xyce::Device::DeviceModel N_DEV_DeviceModel;
typedef Xyce::Device::DeviceOptions N_DEV_DeviceOptions;
typedef Xyce::Device::DeviceSensitivities N_DEV_DeviceSensitivities;
typedef Xyce::Device::DeviceState N_DEV_DeviceState;
typedef Xyce::Device::DeviceSupport N_DEV_DeviceSupport;
typedef Xyce::Device::ExternalSimulationData N_DEV_ExternalSimulationData;
typedef Xyce::Device::ExternCodeInterface N_DEV_ExternCodeInterface;
typedef Xyce::Device::ExternData N_DEV_ExternData;
typedef Xyce::Device::InstanceBlock N_DEV_InstanceBlock;
typedef Xyce::Device::MatrixLoadData N_DEV_MatrixLoadData;
typedef Xyce::Device::ModelBlock N_DEV_ModelBlock;
typedef Xyce::Device::NumericalJacobian N_DEV_NumericalJacobian;
typedef Xyce::Device::Param N_DEV_Param;
typedef Xyce::Device::Region N_DEV_Region;
typedef Xyce::Device::RegionData N_DEV_RegionData;
typedef Xyce::Device::RxnRegion N_DEV_RxnRegion;
typedef Xyce::Device::RxnRegion2 N_DEV_RxnRegion2;
typedef Xyce::Device::RxnRegionData N_DEV_RxnRegionData;
typedef Xyce::Device::SolverState N_DEV_SolverState;
typedef Xyce::Device::SourceInstance N_DEV_SourceInstance;
typedef Xyce::Device::XyceInterface N_DEV_XyceInterface;
typedef Xyce::Device::SourceData N_DEV_SourceData;
typedef Xyce::Device::SinData N_DEV_SinData;
typedef Xyce::Device::ExpData N_DEV_ExpData;
typedef Xyce::Device::ACData N_DEV_ACData;
typedef Xyce::Device::PulseData N_DEV_PulseData;
typedef Xyce::Device::PWLinData N_DEV_PWLinData;
typedef Xyce::Device::SFFMData N_DEV_SFFMData;
typedef Xyce::Device::ConstData N_DEV_ConstData;
typedef Xyce::Device::DevicePDEInstance N_DEV_DevicePDEInstance;
typedef Xyce::Device::DevicePDEModel N_DEV_DevicePDEModel;
typedef Xyce::Device::PDE_Electrode N_DEV_PDE_Electrode;
typedef Xyce::Device::PDE_1DElectrode N_DEV_PDE_1DElectrode;
typedef Xyce::Device::PDE_2DElectrode N_DEV_PDE_2DElectrode;
typedef Xyce::Device::ScalingVars N_DEV_ScalingVars;
typedef Xyce::Device::SpecieSource N_DEV_SpecieSource;
typedef Xyce::Device::ExternDevice::Instance N_DEV_ExternDeviceInstance;
typedef Xyce::Device::ExternDevice::Model N_DEV_ExternDeviceModel;
typedef Xyce::IO::Measure::Base N_IO_MeasureBase;
typedef Xyce::IO::Measure::Average N_IO_MeasureAverage;
typedef Xyce::IO::Measure::DerivativeEvaluation N_IO_MeasureDerivativeEvaluation;
typedef Xyce::IO::Measure::Duty N_IO_MeasureDuty;
typedef Xyce::IO::Measure::EquationEvaluation N_IO_MeasureEquationEvaluation;
typedef Xyce::IO::Measure::FindWhen N_IO_MeasureFindWhen;
typedef Xyce::IO::Measure::Fourier N_IO_MeasureFourier;
typedef Xyce::IO::Measure::Frequency N_IO_MeasureFrequency;
typedef Xyce::IO::Measure::IntegralEvaluation N_IO_MeasureIntegralEvaluation;
typedef Xyce::IO::Measure::Max N_IO_MeasureMax;
typedef Xyce::IO::Measure::Min N_IO_MeasureMin;
typedef Xyce::IO::Measure::OffTime N_IO_MeasureOffTime;
typedef Xyce::IO::Measure::OnTime N_IO_MeasureOnTime;
typedef Xyce::IO::Measure::PeakToPeak N_IO_MeasurePeakToPeak;
typedef Xyce::IO::Measure::RMS N_IO_MeasureRMS;
typedef Xyce::IO::Measure::RelativeError N_IO_MeasureRelativeError;
typedef Xyce::IO::Measure::RiseFallDelay N_IO_MeasureRiseFallDelay;
typedef Xyce::IO::CircuitBlock N_IO_CircuitBlock;
typedef Xyce::IO::CircuitContext N_IO_CircuitContext;
typedef Xyce::IO::CircuitMetadata N_IO_CircuitMetadata;
typedef Xyce::IO::CmdParse N_IO_CmdParse;
typedef Xyce::IO::DeviceBlock N_IO_DeviceBlock;
typedef Xyce::IO::DeviceMetadata N_IO_DeviceMetadata;
typedef Xyce::IO::FourierMgr N_IO_FourierMgr;
typedef Xyce::IO::FunctionBlock N_IO_FunctionBlock;
typedef Xyce::IO::OutputMgr N_IO_OutputMgr;
typedef Xyce::IO::Measure::Manager N_IO_MeasureMgr;
typedef Xyce::IO::ParameterBlock N_IO_ParameterBlock;
typedef Xyce::IO::PkgOptionsMgr N_IO_PkgOptionsMgr;
typedef Xyce::IO::PkgOptionsReg N_IO_PkgOptionsReg;
typedef Xyce::IO::RestartMgr N_IO_RestartMgr;
typedef Xyce::IO::RestartNode N_IO_RestartNode;
typedef Xyce::IO::SpiceSeparatedFieldTool N_IO_SpiceSeparatedFieldTool;
typedef Xyce::IO::OutputFileBase N_IO_OutputFileBase;
typedef Xyce::IO::FunctionBlock FunctionBlock;
typedef Xyce::IO::CircuitContext CircuitContext;
typedef Xyce::Linear::AmesosSolver N_LAS_AmesosSolver;
typedef Xyce::Linear::AztecOOSolver N_LAS_AztecOOSolver;
typedef Xyce::Linear::BlockMatrix N_LAS_BlockMatrix;
typedef Xyce::Linear::BlockVector N_LAS_BlockVector;
typedef Xyce::Linear::Builder N_LAS_Builder;
typedef Xyce::Linear::HBBuilder N_LAS_HBBuilder;
typedef Xyce::Linear::Matrix N_LAS_Matrix;
typedef Xyce::Linear::MultiVector N_LAS_MultiVector;
typedef Xyce::Linear::PrecondFactory N_LAS_PrecondFactory;
typedef Xyce::Linear::Preconditioner N_LAS_Preconditioner;
typedef Xyce::Linear::Problem N_LAS_Problem;
typedef Xyce::Linear::Solver N_LAS_Solver;
typedef Xyce::Linear::System N_LAS_System;
typedef Xyce::Linear::Transform N_LAS_Transform;
typedef Xyce::Linear::Vector N_LAS_Vector;
typedef Xyce::Linear::MatrixFreeEpetraOperator N_LAS_MatrixFreeEpetraOperator;
typedef Xyce::Linear::QueryUtil N_LAS_QueryUtil;
typedef Xyce::Loader::CktLoader N_LOA_CktLoader;
typedef Xyce::Loader::NonlinearEquationLoader N_LOA_NonlinearEquationLoader;
typedef Xyce::Loader::HBLoader N_LOA_HBLoader;
typedef Xyce::Loader::Loader N_LOA_Loader;
typedef Xyce::Nonlinear::DampedNewton N_NLS_DampedNewton;
typedef Xyce::Nonlinear::ConstraintBT N_NLS_ConstraintBT;
typedef Xyce::Nonlinear::ParamMgr N_NLS_ParamMgr;
typedef Xyce::Nonlinear::TwoLevelNewton N_NLS_TwoLevelNewton;
typedef Xyce::Nonlinear::ConductanceExtractor N_NLS_ConductanceExtractor;
typedef Xyce::Nonlinear::Sensitivity N_NLS_Sensitivity;
typedef Xyce::Nonlinear::NonLinearSolver N_NLS_NonLinearSolver;
typedef Xyce::Nonlinear::NLParams N_NLS_NLParams;
typedef Xyce::Nonlinear::NonLinInfo N_NLS_NonLinInfo;
typedef Xyce::Nonlinear::ReturnCodes N_NLS_ReturnCodes;
typedef Xyce::Nonlinear::Manager N_NLS_Manager;
typedef Xyce::Parallel::Communicator N_PDS_Comm;
typedef Xyce::Parallel::Manager N_PDS_Manager;
typedef Xyce::Parallel::Communicator N_PDS_Comm;
typedef Xyce::Parallel::Manager N_PDS_Manager;
typedef Xyce::TimeIntg::DataStore N_TIA_DataStore;
typedef Xyce::TimeIntg::TimeIntegrationMethod N_TIA_TimeIntegrationMethod;
typedef Xyce::TimeIntg::StepErrorControl N_TIA_StepErrorControl;
typedef Xyce::TimeIntg::OneStep N_TIA_OneStep;
typedef Xyce::TimeIntg::MPDEInterface N_TIA_MPDEInterface;
typedef Xyce::TimeIntg::TIAParams N_TIA_TIAParams;
typedef Xyce::TimeIntg::TwoLevelError N_TIA_TwoLevelError;
typedef Xyce::TimeIntg::WorkingIntegrationMethod N_TIA_WorkingIntegrationMethod;
typedef Xyce::Topo::CktGraph N_TOP_CktGraph;
typedef Xyce::Topo::CktGraphBasic N_TOP_CktGraphBasic;
typedef Xyce::Topo::CktGraphCreator N_TOP_CktGraphCreator;
typedef Xyce::Topo::CktGraphCreatorBasic N_TOP_CktGraphCreatorBasic;
typedef Xyce::Topo::CktGraphSupport N_TOP_CktGraphSupport;
typedef Xyce::Topo::CktNode N_TOP_CktNode;
typedef Xyce::Topo::CktNodeCreator N_TOP_CktNodeCreator;
typedef Xyce::Topo::CktNode_Dev N_TOP_CktNode_Dev;
typedef Xyce::Topo::CktNode_V N_TOP_CktNode_V;
typedef Xyce::Topo::Directory N_TOP_Directory;
typedef Xyce::Topo::Indexor N_TOP_Indexor;
typedef Xyce::Topo::Manager N_TOP_Manager;
typedef Xyce::Topo::Node N_TOP_Node;
typedef Xyce::Topo::NodeDevBlock N_TOP_NodeDevBlock;
typedef Xyce::Topo::TopoLSUtil N_TOP_TopoLSUtil;
typedef Xyce::Topo::Topology N_TOP_Topology;
typedef Xyce::NodeID NodeID;
typedef Xyce::Util::MachineDependentParams N_UTL_MachineDependentParams;
typedef Xyce::Util::OptionBlock N_UTL_OptionBlock;
typedef Xyce::Util::Param N_UTL_Param;
typedef Xyce::Util::Timer N_UTL_Timer;
typedef Xyce::Util::Expression N_UTL_Expression;
typedef Xyce::Util::ExpressionInternals N_UTL_ExpressionInternals;
typedef Xyce::Util::ExpressionData N_UTL_ExpressionData;
typedef Xyce::Util::BreakPoint N_UTL_BreakPoint;

#endif
