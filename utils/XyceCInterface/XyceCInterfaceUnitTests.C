//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Unit tests for top-level Xyce::Simulator class
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 2024/2022
//
//
//
//
//-------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "Xyce_config.h"
#include <N_CIR_XyceCInterface.h>

//
// Xyce::Circuit::Simulator functions that need to be tested
//
// x constructor
// x initialize()
// x finalize()
// x getTime()
// x getDACDeviceNames()
// x getADCMap()
// getTimeVoltagePairs()
// updateTimeVoltagePairs()
// x simulationComplete()
// x getFinalTime()
// x setCircuitParameter()
// x getCircuitValue()
// x simulateUntil() in stead of provisionalStep() & acceptProvisionalStep()
// 


TEST ( XyceCInterface, Open)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( *((long *)xycePtr) != 0 );
  
}

TEST ( XyceCInterface, OpenAndClose)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  
}


//-------------------------------------------------------------------------------
int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

