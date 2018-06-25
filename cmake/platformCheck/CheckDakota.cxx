//
// Dakota Test Program
//
// Author: Aaron Gibson
//
// This file tests whether we can compile with Dakota or not.
//

#include <DakotaInterface.hpp>
#include <DakotaModel.hpp>
#include <DakotaResponse.hpp>
#include <DakotaVariables.hpp>

#include <LibraryEnvironment.hpp>
#include <DirectApplicInterface.hpp>
#include <ProblemDescDB.hpp>

int main(int argc, char** argv)
{
  ::Dakota::ProgramOptions opts;

  // Defaults constructs the MPIManager, which assumes COMM_WORLD
  ::Dakota::LibraryEnvironment* dakotaEnv = new ::Dakota::LibraryEnvironment(opts);
  delete dakotaEnv;

  return 0;
}
