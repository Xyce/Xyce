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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 12/21/2014
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_NoiseData_h
#define Xyce_N_ANP_NoiseData_h

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations -------
#include <iosfwd>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Class         : NoiseData
// Purpose       : 
// Special Notes : 
// Creator       : Eric Keiter
// Creation Date : 12/21/2014
//-----------------------------------------------------------------------------
class NoiseData
{
private:
protected:
public:
  NoiseData ():
    deviceName(""),
    noiseNames(0),
    omega(0.0),
    freq(0.0),
    inputNoiseTotal(0),
    outputNoiseTotal(0),
    inputNoiseDens(0),
    outputNoiseDens(0),
    noiseDens(0),
    lnNoiseDens(0),
    lastLnNoiseDens(0),
    li_Pos(0),
    li_Neg(0),
    li_PosCorl(1,-1),
    li_NegCorl(1,-1),
    totalNoise(0.0),
    totalOutputNoise(0.0),
    totalInputNoise(0.0),
    numSources(0), 
    T0(0.0),
    T2(0.0),
    T3(0.0)
  {};

  virtual ~NoiseData ();

  void resize(int size);

private:
protected:
public:
  // device name
  std::string deviceName;

  // noise source names
  std::vector<std::string> noiseNames;

  double omega;
  double freq;

  std::vector<double> inputNoiseTotal;
  std::vector<double> outputNoiseTotal;

  std::vector<double> inputNoiseDens;
  std::vector<double> outputNoiseDens;
  std::vector<double> noiseDens;
  std::vector<double> lnNoiseDens;
  std::vector<double> lastLnNoiseDens;
  std::vector<int> li_Pos;
  std::vector<int> li_Neg;
  std::vector<int> li_PosCorl;
  std::vector<int> li_NegCorl;

  double totalNoise;
  double totalOutputNoise;
  double totalInputNoise;

  int numSources;

  double T0;
  double T2;
  double T3;
};

// inline functions
//-----------------------------------------------------------------------------
// Function      : NoiseData::resize
// Purpose       : resizes everything
// Creator       : Eric R. Keiter
// Creation Date : 1/6/2014
//-----------------------------------------------------------------------------
inline void NoiseData::resize(int size)
{
  noiseNames.clear();
  inputNoiseDens.clear();
  outputNoiseDens.clear();
  noiseDens.clear();
  lnNoiseDens.clear();
  lastLnNoiseDens.clear();
  li_Pos.clear();
  li_Neg.clear();
  li_PosCorl.clear();
  li_NegCorl.clear();

  noiseNames.resize(size);

  inputNoiseTotal.resize(size,0.0);
  outputNoiseTotal.resize(size,0.0);

  inputNoiseDens.resize(size,0.0);
  outputNoiseDens.resize(size,0.0);
  noiseDens.resize(size,0.0);
  lnNoiseDens.resize(size,0.0);
  lastLnNoiseDens.resize(size,0.0);

  li_Pos.resize(size);
  li_Neg.resize(size);
  li_PosCorl.resize(size,-1);
  li_NegCorl.resize(size,-1);
}

//-----------------------------------------------------------------------------
// Function      : NoiseData::~NoiseData
// Purpose       : destructor
// Creator       : Eric R. Keiter
// Creation Date : 1/6/2014
//-----------------------------------------------------------------------------
inline NoiseData::~NoiseData()
{
  noiseNames.clear();
  inputNoiseTotal.clear();
  outputNoiseTotal.clear();
  inputNoiseDens.clear();
  outputNoiseDens.clear();
  noiseDens.clear();
  lnNoiseDens.clear();
  lastLnNoiseDens.clear();
  li_Pos.clear();
  li_Neg.clear();
  li_PosCorl.clear();
  li_NegCorl.clear();
}

//-----------------------------------------------------------------------------
// Function      : NoiseData::operator<<
// Purpose       : "<<" operator
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/21/2014
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const NoiseData & nd)
{
  os << Xyce::section_divider << std::endl;
  os << "noise device name = " << nd.deviceName << "\n";

  os << Xyce::section_divider << std::endl;
  os << std::endl;
  return os;
}

} // namespace Device
} // namespace Xyce

#endif
