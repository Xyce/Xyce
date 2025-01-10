//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_Measure_fwd_h
#define Xyce_N_IO_Measure_fwd_h

namespace Xyce {
namespace IO {
namespace Measure {

class Base;
class Extrema;
class FFT;
class Stats;
class TranStats;
class Average;
class DerivativeEvaluationBase;
class DerivativeEvaluation;
class DerivativeEvaluationCont;
class Duty;
class EquationEvaluation;
class ErrorFunctions;
class Err1;
class Err2;
class WhenAT;
class FindWhenBase;
class FindWhen;
class FindWhenCont;
class Fourier;
class Frequency;
class IntegralEvaluation;
class Max;
class Min;
class OffTime;
class OnTime;
class PeakToPeak;
class RMS;
class Error;
class RiseFallDelay;
class TrigTargBase;
class TrigTarg;
class TrigTargCont;
class FFTFind;
class ENOB;
class SFDR;
class SNDR;
class SNR;
class THD;

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Measure_fwd_h
