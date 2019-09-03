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
// Purpose        : Classes for remeasure objects
//
// Special Notes  :
//
// Creator        : Pete Sholander, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 2/15/18
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Remeasure_h
#define Xyce_N_IO_Remeasure_h

// ----------   Standard Includes   ----------
#include <string>
#include <list>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_OutputMgr.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {
namespace Measure{

class RemeasureBase;
class RemeasureAC;
class RemeasureDC;
class RemeasureTRAN;

//-------------------------------------------------------------------------
// Class         : RemeasureBase
// Purpose       : Base class for remeasure objects.
// Special Notes : 
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 2/15/18
//-------------------------------------------------------------------------
class RemeasureBase 
{
 
public:
  RemeasureBase(
    N_PDS_Comm &pds_comm, 
    Manager &measure_manager, 
    OutputMgr &output_manager,
    Analysis::AnalysisManager &analysis_manager,
    Analysis::AnalysisCreatorRegistry &analysis_registry);
  virtual ~RemeasureBase();

  // getters and setters for base object variables
  int getIndepVarCol() {return index;}
  Analysis::Mode getAnalysisMode() { return analysis_mode; }
  void setAnalysisMode(Analysis::Mode value) {analysis_mode = value;}
  
  // All of these functions must be defined in the derived classes
  virtual void setIndepVarCol(int rank, int i, std::string colName) = 0;
  virtual void checkIndepVarCol(int rank, int i) = 0;
  virtual void setIndepVar(double value) = 0;
  virtual double getIndepVar() = 0;
  virtual void updateMeasures(Linear::Vector& varValuesVec) = 0;

  // These functions are not pure virtual because they are only redefined for 
  // the RemeasureDC object.
  virtual void resetSweepVars() { return; }
  virtual void updateSweepVars() { return; }

private:
  RemeasureBase(const RemeasureBase &right);
  RemeasureBase &operator=(const RemeasureBase &right);

protected:
  N_PDS_Comm &pds_comm;  
  Manager &measure_manager;
  OutputMgr &output_manager;
  Analysis::AnalysisManager &analysis_manager;
  Analysis::AnalysisCreatorRegistry &analysis_registry;

  // this enum is used to set the file suffix (.ma,.ms or .mt) and also during checking of the
  // measure mode vs. the analysis mode. 
  Analysis::Mode analysis_mode; 

  // index of the column with the independent variable in the remeasured data file
  int index;  
};

//-------------------------------------------------------------------------
// Class         : RemeasureAC
// Purpose       : Class for AC remeasure objects.
// Special Notes : 
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/15/18
//-------------------------------------------------------------------------
 class RemeasureAC : public RemeasureBase
{
 
public:
  RemeasureAC(
    N_PDS_Comm &pds_comm, 
    Manager &measure_manager, 
    OutputMgr &output_manager,
    Analysis::AnalysisManager &analysis_manager,
    Analysis::AnalysisCreatorRegistry &analysis_registry);
  ~RemeasureAC();
  
  void setIndepVarCol(int rank, int i, std::string colName);
  void checkIndepVarCol(int rank, int i);
  void setIndepVar(double value) {frequency=value;}  
  double getIndepVar() { return frequency;}

  // Used to update the AC measure values.  This function call will need to 
  // be updated when/if lead currents are supported for AC analyses.
  void updateMeasures(Linear::Vector& varValuesVec)
  {
    measure_manager.updateACMeasures(pds_comm.comm(), frequency, &varValuesVec, 0, 0);
  }

private:
  double frequency;  // frequency is the independent variable for AC remeasure
};

//-------------------------------------------------------------------------
// Class         : RemeasureDC
// Purpose       : Class for DC remeasure objects.
// Special Notes : 
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
 class RemeasureDC : public RemeasureBase
{
 
public:
  RemeasureDC(
    N_PDS_Comm &pds_comm, 
    Manager &measure_manager, 
    OutputMgr &output_manager,
    Analysis::AnalysisManager &analysis_manager,
    Analysis::AnalysisCreatorRegistry &analysis_registry);
  ~RemeasureDC();
  
  void setIndepVarCol(int rank, int i, std::string colName);
  void checkIndepVarCol(int rank, int i);
  void setIndepVar(double value) {DCIndex = value;}
  double getIndepVar() { return DCIndex;}

  // used to update the DC measure values
  void updateMeasures(Linear::Vector& varValuesVec)
  {
    measure_manager.updateDCMeasures(pds_comm.comm(), dcParamsVec, &varValuesVec, 0, 0, &varValuesVec, 0, 0 );
  }

  // used to manage the DC Sweep Vector when re-measuring DC data
  void resetSweepVars();
  void updateSweepVars();

private:
  Analysis::DCSweep *SweepVecPtr;  // used during the partial instantiation of the Analysis Manager
  std::vector<Analysis::SweepParam> dcParamsVec;   // used for local recreation of DC Params Vector
  int DCIndex;      // the value of the INDEX column is the "independent variable" for DC remeasure
  int dcStepCount;  // used to sense steps caused by .DC line, and to track position within the current step
};

//-------------------------------------------------------------------------
// Class         : RemeasureTRAN
// Purpose       : Class for TRAN remeasure objects.
// Special Notes : 
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/15/18
//-------------------------------------------------------------------------
class RemeasureTRAN : public RemeasureBase
{
 
public:
  RemeasureTRAN(
    N_PDS_Comm &pds_comm, 
    Manager &measure_manager, 
    OutputMgr &output_manager,
    Analysis::AnalysisManager &analysis_manager,
    Analysis::AnalysisCreatorRegistry &analysis_registry);
  ~RemeasureTRAN();
  
  void setIndepVarCol(int rank, int i, std::string colName);
  void checkIndepVarCol(int rank, int i);

  // for TRAN remeasure, the TIME is the independent variable, and it is kept by the 
  // output_manager, rather than as private variable in the RemeasureTRAN object.  This 
  // is because many of the TRAN measures make calls to output_manager.getCircuitTime()
  void setIndepVar(double value) { output_manager.setCircuitTime(value); }
  double getIndepVar() { return output_manager.getCircuitTime(); }

  // used to update the TRAN measure values
  void updateMeasures(Linear::Vector& varValuesVec)
  {
    measure_manager.updateTranMeasures(pds_comm.comm(), output_manager.getCircuitTime(), &varValuesVec, 0, 0, &varValuesVec, 0, 0 );
  }

private:
 
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Remeasure_h
