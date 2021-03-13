//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Outputter for .LIN analysis.  This handles the output of
//                  S-, Y- and Z-parameters in Touchstone1 format.
// Special Notes  :
//
// Creator        : Pete Sholander
//
// Creation Date  : 3/25/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterSParamTS1_h
#define Xyce_N_IO_OutputterSParamTS1_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// S-Parameter outputters

class SParamTS1 : public Interface
{
public:
  SParamTS1(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~SParamTS1();

private:
  SParamTS1(const SParamTS1 &);
  SParamTS1 &operator=(const SParamTS1 &);

public:

  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doFinishOutput();

  virtual void doStartStep(int step, int max_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

  virtual void doOutputSParams(
    Parallel::Machine   comm,
    double              frequency,
    double              numFreq,
    std::vector<double> & Z0sVec,
    const Util::Op::RFparamsData & RFparams);

private:
  void sparamHeader(
    Parallel::Machine   comm,
    double              numFreq,
    std::vector<double> & Z0sVec,
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams);

private:
  OutputMgr &                   outputManager_;
  PrintParameters               printParameters_;
  std::string                   outFilename_;
  std::ostream *                os_;
  int                           printCount_;
  int                           index_;
  int                           currentStep_;
  int                           numberOfSteps_;

  double numPorts_;
  std::vector<std::string>      sensParFullNames_;
  Table::ColumnList             columnList_;
  int numColumns_;  // number of columns in the output, including the Freq column.
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterSParamsTS1_h
