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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Outputter_h
#define Xyce_N_IO_Outputter_h

#include <list>
#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_IO_OutputTypes.h>
#include <N_LAS_fwd.h>
#include <N_ANP_fwd.h>

#include <N_ANP_UQSupport.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Op.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {

// print formats
namespace Format {
enum Format {STD, TECPLOT, PROBE, CSV, RAW, RAW_ASCII, DAKOTA, TS1, TS2};
}

// data formats for Touchstone output
namespace DataFormat {
enum DataFormat {RI, MA, DB};
}

//-----------------------------------------------------------------------------
// Class         : Table
//
// Purpose       : This struct manages the header data
//
// Special Notes :
//
// Creator       : Todd Coffey, 1414
// Creation Date : 9/19/08
//-----------------------------------------------------------------------------
struct Table
{
  // Enum for column justification (left/center/right) in std header output
  enum Justification
    {
      JUSTIFICATION_LEFT,
      JUSTIFICATION_CENTER,
      JUSTIFICATION_RIGHT,
      JUSTIFICATION_NONE
    };

  struct Column
  {
    Column()
      : name_(),
        format_(std::ios_base::scientific),
        width_(17),
        precision_(9),
        justification_(JUSTIFICATION_LEFT)
    {}

    Column(const Column &column)
      : name_(column.name_),
        format_(column.format_),
        width_(column.width_),
        precision_(column.precision_),
        justification_(column.justification_)
    {}

    Column(const std::string &name, std::ios_base::fmtflags format, int width, int precision, Justification justification)
      : name_(name),
        format_(format),
        width_(width),
        precision_(precision),
        justification_(justification)
    {}

    std::string             name_;
    std::ios_base::fmtflags format_;
    int                     width_;
    int                     precision_;
    Justification           justification_;
  };

  typedef std::vector<Column> ColumnList;

  Table()
  {}

  Table(const Table &table)
    : columnList_(table.columnList_.begin(), table.columnList_.end())
  {}

  Table &operator=(const Table &table)
  {
    columnList_.assign(table.columnList_.begin(), table.columnList_.end());

    return *this;
  }

  virtual ~Table()
  {}

  void addColumn(std::string name, std::ios_base::fmtflags format, int width, int precision, Justification justification)
  {
    columnList_.push_back(Column(name, format, width, precision, justification));
  }

  void addColumn(std::string name, int width, int precision, Justification justification)
  {
    columnList_.push_back(Column(name, std::ios_base::scientific, width, precision, justification));
  }

  ColumnList          columnList_;
};


//-----------------------------------------------------------------------------
// Class         : Table
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 17 15:16:48 2014
//-----------------------------------------------------------------------------
///
/// A fallback print parameters are generated from a none specific print
/// type.  I.E. .PRINT TRAN creates a fallback .PRINT HOMOTOPY print
/// parameters, while .PRINT HOMOTOPY does not.  Fallback print
/// parameters are removed if a specific print parameters is provided.
///
struct PrintParameters
{
  PrintParameters()
    : fallback_(false),
      filename_(),
      suffix_(),
      defaultExtension_(),
      overrideRawFilename_(),
      dashoFilename_(),
      analysisMode_(Analysis::ANP_MODE_INVALID),
      overrideRaw_(false),
      asciiRaw_(false),
      formatSupportsOverrideRaw_(true),
      dashoRequested_(false),
      format_(Format::STD),
      dataFormat_(DataFormat::RI),
      RFparamType_("S"),
      printIndexColumn_(true),
      printStepNumColumn_(false),
      variableList_(),
      table_(),
      streamWidth_(17),
      streamPrecision_(9),
      timeWidth_(8),
      delimiter_(),
      outputTimeScaleFactor_(1.0),
      filter_(0.0),
      expandComplexTypes_(false),
      addGnuplotSpacing_(false),
      addSplotSpacing_(false),
      outputPCEsampleStats_(true),
      outputAllPCEsamples_(false),
      outputPCECoeffs_(false)
  {}

  PrintParameters(const PrintParameters &print_parameters)
    : fallback_(print_parameters.fallback_),
      filename_(print_parameters.filename_),
      suffix_(print_parameters.suffix_),
      defaultExtension_(print_parameters.defaultExtension_),
      overrideRawFilename_(print_parameters.overrideRawFilename_),
      dashoFilename_(print_parameters.dashoFilename_),
      netlistLocation_(print_parameters.netlistLocation_),
      analysisMode_(print_parameters.analysisMode_),
      overrideRaw_(print_parameters.overrideRaw_),
      asciiRaw_(print_parameters.asciiRaw_),
      formatSupportsOverrideRaw_(print_parameters.formatSupportsOverrideRaw_),
      dashoRequested_(print_parameters.dashoRequested_),
      format_(print_parameters.format_),
      dataFormat_(print_parameters.dataFormat_),
      RFparamType_(print_parameters.RFparamType_),
      printIndexColumn_(print_parameters.printIndexColumn_),
      printStepNumColumn_(print_parameters.printStepNumColumn_),
      variableList_(print_parameters.variableList_.begin(), print_parameters.variableList_.end()),
      table_(print_parameters.table_),
      streamWidth_(print_parameters.streamWidth_),
      streamPrecision_(print_parameters.streamPrecision_),
      timeWidth_(print_parameters.timeWidth_),
      delimiter_(print_parameters.delimiter_),
      outputTimeScaleFactor_(print_parameters.outputTimeScaleFactor_),
      filter_(print_parameters.filter_),
      expandComplexTypes_(print_parameters.expandComplexTypes_),
      addGnuplotSpacing_(print_parameters.addGnuplotSpacing_),
      addSplotSpacing_(print_parameters.addSplotSpacing_),
      outputPCEsampleStats_(print_parameters.outputPCEsampleStats_),
      outputAllPCEsamples_(print_parameters.outputAllPCEsamples_),
      outputPCECoeffs_(print_parameters.outputPCECoeffs_)
  {}

  PrintParameters &operator=(const PrintParameters &print_parameters)
  {
    fallback_ = print_parameters.fallback_;
    filename_ = print_parameters.filename_;
    suffix_ = print_parameters.suffix_;
    defaultExtension_ = print_parameters.defaultExtension_;
    overrideRawFilename_ = print_parameters.overrideRawFilename_,
    dashoFilename_ = print_parameters.dashoFilename_,
    netlistLocation_ = print_parameters.netlistLocation_;
    analysisMode_ = print_parameters.analysisMode_;
    overrideRaw_ = print_parameters.overrideRaw_;
    asciiRaw_ = print_parameters.asciiRaw_;
    formatSupportsOverrideRaw_ = print_parameters.formatSupportsOverrideRaw_;
    dashoRequested_ = print_parameters.dashoRequested_;
    format_ = print_parameters.format_;
    dataFormat_ = print_parameters.dataFormat_;
    RFparamType_ = print_parameters.RFparamType_;
    printIndexColumn_ = print_parameters.printIndexColumn_;
    printStepNumColumn_ = print_parameters.printStepNumColumn_;
    variableList_.assign(print_parameters.variableList_.begin(), print_parameters.variableList_.end());
    table_ = print_parameters.table_;
    streamWidth_ = print_parameters.streamWidth_;
    streamPrecision_ = print_parameters.streamPrecision_;
    timeWidth_ = print_parameters.timeWidth_;
    delimiter_ = print_parameters.delimiter_;
    outputTimeScaleFactor_ = print_parameters.outputTimeScaleFactor_;
    filter_ = print_parameters.filter_;
    expandComplexTypes_ = print_parameters.expandComplexTypes_;
    addGnuplotSpacing_ = print_parameters.addGnuplotSpacing_;
    addSplotSpacing_ = print_parameters.addSplotSpacing_;
    outputPCEsampleStats_ = print_parameters.outputPCEsampleStats_;
    outputAllPCEsamples_ = print_parameters.outputAllPCEsamples_;
    outputPCECoeffs_ = print_parameters.outputPCECoeffs_;

    return *this;
  }

public:
  virtual ~PrintParameters()
  {}

  bool                          fallback_;                      ///< True if the creating PrintType is generic
  std::string                   filename_;                      ///< Output filename
  std::string                   suffix_;                        ///< Suffix (placed between filename and extension)
  std::string                   defaultExtension_;              ///< Extension if none provided
  std::string                   overrideRawFilename_;           ///< Output filename from -r command line option
  std::string                   dashoFilename_;                 ///< Output filename from -o command line option

  NetlistLocation               netlistLocation_;               ///< Location in netlist creating this PrintParameters
  Analysis::Mode                analysisMode_;                  ///< Analysis mode that triggered the current print
  bool                          overrideRaw_;                   ///< True if -r specified on command line
  bool                          asciiRaw_;                      ///< True if -a specified on command line
  bool                          formatSupportsOverrideRaw_;     ///< should -r output be generated for this format_
  bool                          dashoRequested_;                ///< true if -o specified on command line, but not -r
  Format::Format                format_;                        ///< Print file format specified
  DataFormat::DataFormat        dataFormat_;                    ///< Data format specified for Touchstone output
  std::string                   RFparamType_;                   ///< Parameter type (e.g., S, Y or Z) output in Touchstone file
  bool                          printIndexColumn_;              ///< True if INDEX column is to be printed
  bool                          printStepNumColumn_;            ///< True if STEPNUM column is to be printed
  Util::ParamList               variableList_;                  ///< Description of variables to be printed
  Table                         table_;                         ///< Formatting of table for print
  int                           streamWidth_;                   ///< Column width
  int                           streamPrecision_;               ///< Column precision
  int                           timeWidth_;                     ///< Really?
  std::string                   delimiter_;                     ///< Delimiter
  double                        outputTimeScaleFactor_;         ///< Output time in something other than seconds (such as milli-seconds)
  double                        filter_;
  bool                          expandComplexTypes_;            ///< For example, output V(1) as two columns, Re(V(1)) and Im(V(1))
  bool                          addGnuplotSpacing_;             ///< True means add two blank lines after steps 1,2,...,N-1 if
                                                                ///< there are more than one step.  For Gnuplot compatibility.
  bool                          addSplotSpacing_;               ///< True means add one blank line after steps 1,2,...,N-1 if
                                                                ///< there are more than one step.  For Gnuplot compatibility.
  bool                          outputPCEsampleStats_;          ///< Used by EmbeddedSampling outputters
  bool                          outputAllPCEsamples_;           ///< Used by EmbeddedSampling outputters
  bool                          outputPCECoeffs_;               ///< Used by EmbeddedSampling outputters
};

namespace Outputter {

class Interface
{
public:
  Interface()
  {}

  virtual ~Interface()
  {}

  void setAnalysisMode(Analysis::Mode analysis_mode);

  void startStep(int step, int max_step);

  void resetIndex();

  void output(
    Parallel::Machine           comm,
    const Linear::Vector &        solution_vector,
    const Linear::Vector &        state_vector,
    const Linear::Vector &        store_vector,
    const Linear::Vector &        lead_current_vector,
    const Linear::Vector &        junction_voltage_vector);

  void outputAC(
    Parallel::Machine           comm,
    double                      frequency,
    double                      fStart,
    double                      fStop,
    const Linear::Vector &        real_solution_vector,
    const Linear::Vector &        imaginary_solution_vector,
    const Util::Op::RFparamsData & RFparams);

  void outputSensitivityAC(
    Parallel::Machine                 comm,
    double                            frequency,
    const Linear::Vector &            freqDomainSolnVecReal,
    const Linear::Vector &            freqDomainSolnVecImaginary,
    const std::vector<double> &       paramVals,
    const std::vector<std::string> &  paramNameVec,
    const std::vector<std::string> &  objFuncVars,
    const std::vector<double> &       objectiveVec,
    const std::vector<double> &       dOdpVec,
    const std::vector<double> &       dOdpAdjVec,
    const std::vector<double> &       scaled_dOdpVec,
    const std::vector<double> &       scaled_dOdpAdjVec);

  void outputSParams(
    Parallel::Machine           comm,
    double                      frequency,
    double                      numFreq,
    std::vector<double> &       Z0sVec,
    const Util::Op::RFparamsData & RFparams);

  virtual void outputNoise(
    Parallel::Machine   comm,
    double              frequency,
    const Linear::Vector &        real_solution_vector,
    const Linear::Vector &        imaginary_solution_vector,
    double              totalOutputNoiseDens_, 
    double              totalInputNoiseDens_, 
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_);

  virtual void outputEmbeddedSampling(
    Parallel::Machine comm,
    bool regressionPCEenable,
    bool projectionPCEenable,
    int  numSamples,
    const std::vector<std::string> & regressionPCEcoeffs,
    const std::vector<std::string> & projectionPCEcoeffs,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec_);

  // Used for HB time-domain output such as .PRINT HB_TD lines.  This is
  // not used for .PRINT HB_STARTUP or .PRINT HB_IC lines though.
  void outputHB_TD(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const Linear::BlockVector &   timeDomainSolutionVec,
    const Linear::BlockVector &   timeDomainLeadCurrentVec,
    const Linear::BlockVector &   timeDomainJunctionVoltageVec);

  // Used for HB frequency-domain output such as .PRINT HB_FD lines.
  void outputHB_FD(
    Parallel::Machine           comm,
    const std::vector<double> & freqPoints,
    const Linear::BlockVector &   freqDomainSolutionVecReal,
    const Linear::BlockVector &   freqDomainSolutionVecImaginary,
    const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
    const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary);

  void outputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const Linear::BlockVector &   solution_vector);

  void outputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         parameter_values,
    const Linear::Vector &                solution_vector);

  void outputSensitivity(
    Parallel::Machine                   comm,
    const std::vector<double> &         objective_values,
    const std::vector<double> &         direct_values,
    const std::vector<double> &         adjoint_values,
    const std::vector<double> &         scaled_direct_values,
    const std::vector<double> &         scaled_adjoint_values,
    const Linear::Vector &                solution_vector,
    const Linear::Vector &                state_vector,
    const Linear::Vector &                store_vector);

  void finishOutput();

  void steppingComplete();

private:
  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode) = 0;

  virtual void doOutputTime(
    Parallel::Machine           comm,
    const Linear::Vector &        solution_vector,
    const Linear::Vector &        state_vector,
    const Linear::Vector &        store_vector,
    const Linear::Vector &        lead_current_vector,
    const Linear::Vector &        junction_voltage_vector) {}

  virtual void doOutputFrequency(
    Parallel::Machine           comm,
    double                      frequency,
    double                      fStart,
    double                      fStop,
    const Linear::Vector &        real_solution_vector,
    const Linear::Vector &        imaginary_solution_vector,
    const std::map<std::string, Teuchos::SerialDenseMatrix<int, std::complex<double> > * > & RFparams) {}

  virtual void doOutputSensitivityAC(
    Parallel::Machine                   comm,
    double                              frequency,
    const Linear::Vector &              real_solution_vector,
    const Linear::Vector &              imaginary_solution_vector,
    const std::vector<double> &         paramVals,
    const std::vector<std::string> &    paramNameVec,
    const std::vector<std::string> &    objFuncVars,
    const std::vector<double> &         objectiveVec,
    const std::vector<double> &         dOdpVec,
    const std::vector<double> &         dOdpAdjVec,
    const std::vector<double> &         scaled_dOdpVec,
    const std::vector<double> &         scaled_dOdpAdjVec) {}

   virtual void doOutputSParams(
    Parallel::Machine           comm,
    double                      frequency,
    double                      numFreq,
    std::vector<double> & Z0sVec,
    const Util::Op::RFparamsData & RFparams) {}

  virtual void doOutputNoise(
    Parallel::Machine   comm,
    double              frequency,
    const Linear::Vector &        real_solution_vector,
    const Linear::Vector &        imaginary_solution_vector,
    double              totalOutputNoiseDens_, 
    double              totalInputNoiseDens_, 
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_) {}

  virtual void doOutputEmbeddedSampling(
    Parallel::Machine   comm,
    bool                regressionPCEenable,
    bool                projectionPCEenable,
    int                 numSamples,
    const std::vector<std::string> & regressionPCEcoeffs,
    const std::vector<std::string> & projectionPCEcoeffs,
    const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec) {}

  virtual void doOutputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const Linear::BlockVector &   timeDomainSolutionVec,
    const Linear::BlockVector &   freqDomainSolutionVecReal,
    const Linear::BlockVector &   freqDomainSolutionVecImaginary,
    const Linear::BlockVector &   timeDomainLeadCurrentVec,
    const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
    const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
    const Linear::BlockVector &   timeDomainJunctionVoltageVec,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary) {}

  virtual void doOutputHB_TD(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const Linear::BlockVector &   timeDomainSolutionVec,
    const Linear::BlockVector &   timeDomainLeadCurrentVec,
    const Linear::BlockVector &   timeDomainJunctionVoltageVec) {}

  virtual void doOutputHB_FD(
    Parallel::Machine           comm,
    const std::vector<double> & freqPoints,
    const Linear::BlockVector &   freqDomainSolutionVecReal,
    const Linear::BlockVector &   freqDomainSolutionVecImaginary,
    const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
    const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
    const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary) {}

  virtual void doOutputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const Linear::BlockVector &   solution_vector) {}

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
     const std::vector<std::string> &   parameter_names,
     const std::vector<double> &        param_values,
     const Linear::Vector &               solution_vector) {}

  virtual void doOutputSensitivity(
    Parallel::Machine                    comm,
     const std::vector<double> &         objective_values,
     const std::vector<double> &         direct_values,
     const std::vector<double> &         adjoint_values,
     const std::vector<double> &         scaled_direct_values,
     const std::vector<double> &         scaled_adjoint_values,
     const Linear::Vector &                solution_vector,
     const Linear::Vector &                state_vector,
     const Linear::Vector &                store_vector) {}

  virtual void doFinishOutput() {}

  virtual void doStartStep(int step, int max_step) {}

  virtual void doResetIndex() {}

  virtual void doSteppingComplete() {}
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Outputter_h
