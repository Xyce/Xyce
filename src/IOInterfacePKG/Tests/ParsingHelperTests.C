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

#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include <Xyce_config.h>
#include <N_IO_ParsingHelpers.h>

namespace {
TEST(IO_Parsing_Helpers, isAbsolutePath)
{
  ASSERT_EQ(true, Xyce::IO::isAbsolutePath("/path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("./path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("../path"));

  #ifdef HAVE_WINDOWS_H
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("\\\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("C:\\path"));
    ASSERT_EQ(true, Xyce::IO::isAbsolutePath("c:\\path"));
  #else
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("\\\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("C:\\path"));
    ASSERT_EQ(false, Xyce::IO::isAbsolutePath("c:\\path"));
  #endif

  // These are relative paths even on Windows
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("C:path"));
  ASSERT_EQ(false, Xyce::IO::isAbsolutePath("c:path"));
}

TEST(IO_Parsing_Helpers, hasWinDriveLetter)
{
  #ifdef HAVE_WINDOWS_H
  ASSERT_EQ(true,Xyce::IO::hasWinDriveLetter("c:"));
  ASSERT_EQ(true,Xyce::IO::hasWinDriveLetter("C:"));
  #else
  ASSERT_EQ(false,Xyce::IO::hasWinDriveLetter("c:"));
  ASSERT_EQ(false,Xyce::IO::hasWinDriveLetter("C:"));
  #endif
}

} // namespace

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
