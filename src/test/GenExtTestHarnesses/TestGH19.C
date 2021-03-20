#include <iostream>
#include <string>
#include <cstring>
#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_ERH_Message.h>

void report_handler(const char *message, unsigned type) {
  std::cout << message << std::endl;
}

// This very cheesy test harness is designed to test the fix to
// Github issue 19/Gitlab-ex issue 150.
// In that issue, it was discovered that static data in the ErrorHandlingPKG
// was not being cleared after a failed use of an API simulator object,
// so that a subsequent allocation and use of a new object broke.
//
// In this program we allocate a GenCouplingSimulator object and invoke it
// to initialize itself using a known-bad netlist that will cause an error.
// This simulator object will go out of scope and be deleted and a new one
// allocated and initialized with a known-good netlist.
//
// Prior to the fix, the second would exit with a fatal error because the
// static data from the previous run still indicated pending fatal errors.
// As of the fix to the issue, this data is cleared out upon allocation of
// the simulator object.

int main(int argc, char *argv[])
{
  char *fakeargv[2];
  int fakeargc=2;
  fakeargv[0]=argv[0];
  {
    std::string netlist= "netlist_bad.txt";
    char *netlistptr = (char *)malloc(netlist.length()+1);
    strcpy(netlistptr,netlist.c_str());
    fakeargv[1]= netlistptr;
    Xyce::Circuit::GenCouplingSimulator xyce1;
    Xyce::set_report_handler(report_handler);
    try {
      xyce1.initialize(fakeargc, fakeargv);
    } catch (std::runtime_error &e) {
      std::cout << std::endl << "first instance crashed" << std::endl;
    }
  }
  {
    std::string netlist= "netlist_good.txt";
    char *netlistptr = (char *)malloc(netlist.length()+1);
    strcpy(netlistptr,netlist.c_str());
    fakeargv[1]= netlistptr;
    Xyce::Circuit::GenCouplingSimulator xyce2;
    Xyce::set_report_handler(report_handler);
    try {
      xyce2.initialize(fakeargc, fakeargv);
    } catch (std::runtime_error &e) {
      std::cout << std::endl << "second instance crashed" << std::endl;
    }
  }
}
