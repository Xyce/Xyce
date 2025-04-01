//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include <Xyce_config.h>
#include <N_IO_ParsingHelpers.h>
#include <N_IO_CmdParse.h>

namespace {
TEST(IO_CmdParse, removeExtensionsFromDashoFilename)
{
  // the default extensions should be removed
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.HB.FD.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.FD.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.ES.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.NOISE.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.PCE.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.prn"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.s1p"));
  ASSERT_EQ("dasho", Xyce::IO::removeExtensionsFromDashoFilename("dasho.s12p"));

  // the input names should not be modified
  ASSERT_EQ("dasho.prnn", Xyce::IO::removeExtensionsFromDashoFilename("dasho.prnn"));
  ASSERT_EQ("dasho.ss1p", Xyce::IO::removeExtensionsFromDashoFilename("dasho.ss1p"));
  ASSERT_EQ("dasho.sp", Xyce::IO::removeExtensionsFromDashoFilename("dasho.sp"));
  ASSERT_EQ("dasho.sap", Xyce::IO::removeExtensionsFromDashoFilename("dasho.sap"));
  ASSERT_EQ(".prn", Xyce::IO::removeExtensionsFromDashoFilename(".prn"));
  ASSERT_EQ(".s1p", Xyce::IO::removeExtensionsFromDashoFilename(".s1p"));
}

} // namespace

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
