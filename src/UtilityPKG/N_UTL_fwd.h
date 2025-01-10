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

#ifndef Xyce_N_UTL_fwd_h
#define Xyce_N_UTL_fwd_h

#include <complex>
#include <iosfwd>

#include <unordered_map>
using std::unordered_map;

#include <list>
#include <string>
#include <utility>
#include <vector>

#include <N_UTL_NoCase.h>

namespace Xyce {

std::ostream &lout();
std::ostream &dout();
std::ostream &pout();

extern const char *section_divider;             // Defined in LogStream
extern const char *subsection_divider;          // Defined in LogStream

typedef unordered_map<std::string, int, HashNoCase, EqualNoCase> NodeNameMap;

template <class T>
struct DataTypeTrait;

typedef std::complex<double> complex;

class ExtendedString;

namespace Util {

class BreakPoint;
class Expression;
class ExpressionData;
class ExpressionInternals;
class baseExpressionGroup;
class MachineDependentParams;
class Marshal;
class OptionBlock;
class Param;
class Timer;
struct FreqVecEntry; 
struct FreqMatEntry; 

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &push(std::basic_ostream<Ch, Tr> &os);

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &pop(std::basic_ostream<Ch, Tr> &os);

namespace Op {

typedef size_t Identifier;      ///< Identifies the operator with a single value

class Operator;
class BuilderManager;

template<class T, class R, class E = T>
class Op;
struct OpData;

typedef std::vector<Operator *> OpList;

} // namespace Op
} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_fwd_h
