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

//-------------------------------------------------------------------------
//
// Purpose        : External code output abstract interface
//
// Special Notes  : 
//
// Creator        : Tom Russo, SNL, Electrical Models and Simulation
//
// Creation Date  : 2/2/2018
//
//-------------------------------------------------------------------------
#ifndef N_IO_ExtOutInterface_H
#define N_IO_ExtOutInterface_H

#include<string>
#include<vector>
#include <N_UTL_fwd.h>
#include <N_IO_OutputTypes.h>

namespace Xyce {
namespace IO {

/// External code coupling output interface class
///
/// This "abstract" interface is intended to be a base class to allow
/// external, coupled codes to bypass netlist ".print" output formats,
/// and use instead a set of API calls to establish the desired
/// outputs, and provide a means to deliver them from Xyce to the
/// external code.
///
/// Because the methods here are not pure virtual functions, but
/// rather empty stubs, it isn't really an abstract interface in the
/// usual sense.  The intent here is that the user of this type of
/// interface should only have to implement the handful of methods
/// actually needed to get the job done.
///
/// Objects of classes inherited from this interface will be used by Xyce to
/// call back to the external code.  These callbacks will either be to
/// query the external code for information about requested output,
/// or to provide the output and output metadata (such as error conditions,
/// header information, and information about step analyses.
///
class ExternalOutputInterface
{

public:
  ExternalOutputInterface(){};
  virtual ~ExternalOutputInterface(){};

  /// Inform Xyce what type of "print" output we want to handle with this
  /// interface.
  ///
  /// This corresponds to the second word in a .print statement, as in
  /// .print TRAN ...
  ///
  /// Because we expect that most use of this interface will be for transient
  /// output, we default to that.  Users need not reimplement this method
  /// unless they're using something else.
  virtual OutputType::OutputType getOutputType()
  {
    return OutputType::TRAN;
  };

  /// Give Xyce a name for this output interface
  ///
  /// This function will be called when a "location" is needed for
  /// outputting error messages.  It will be used in place of a netlist
  /// name.
  ///
  /// It is not technically required that one implement this and make
  /// every output interface have a unique name, but it could make debugging
  /// difficult should one have more than one outputter active *and* one
  /// of them is getting syntax or other errors.  So it is recommended that
  /// every object of classes derived from this interface implement this
  /// function, and make each object have a unique name.
  virtual std::string getName()
  {
    return "ExternalOutput";
  }
  
  /// Inform Xyce what output variables the external simulator wants to see
  virtual void requestedOutputs(std::vector<std::string> & outputVars) {};

  /// Inform external simulator about status of output string parsing
  ///
  /// Xyce will use this method after attempting to decode the strings
  /// provided to it by requestedOutputs.  It will pass a vector of bools
  /// indicating the parse status of each string.  A "true" in this vector
  /// means that the string at that index was parseable into the low-level
  /// data structure that Xyce uses for print variables.  A "false" indicates
  /// that some syntax problem was detected, and the string at that position
  /// was ignored.
  ///
  /// Note that at this stage, we don't know if the requested output is
  /// actually valid, only that we were able to parse it.  It could be the
  /// case that a string that is parseable at this stage is still invalid,
  /// but we don't actually know that until later, when we try to create
  /// access ops.  
  virtual void reportParseStatus(std::vector<bool> & statusVec) {};

  
  /// Inform external simulator of names of fields that will actually be output
  ///
  /// This method is intentionally pure-virtual, to assure that every
  /// derived class reimplement it.
  ///
  /// After all output set-up is complete, but before any simulation output is
  /// sent, Xyce will call this function to give the external code the
  /// actual list of output fields that it will generate.  These output
  /// names will be in the same order as the data that will be delivered.
  ///
  /// Ideally, if all requested fields were valid, there will be as many
  /// names here as requested.  Unparseable names will have been ignored,
  /// though, and it is the external code's responsibility to check.
  ///
  /// Furthermore, Xyce may prepend an independent variable to the
  /// output list.  This method MUST be implemented.
  ///
  /// For transient, this will be TIME, for AC, FREQUENCY,
  /// etc.  For DC simulations, no variable will be prepended, just as for
  /// .print lines --- the user must explicitly request DC sweep variables
  /// to be output.  
  ///
  /// It is therefore essential that external codes using this interface
  /// check the names that come back, so they can properly relate them to
  /// the original request.
  ///
  /// These names are taken from the actual ops generated to process the
  /// output, and will generally have been upcased.  Further, any
  /// brace-delimited expressions will have names that are constructed by
  /// traversing the expression parse tree, and therefore will not only be
  /// upcased, but may have had spaces removed.
  ///
  /// These are the same names that would normally appear in the header
  /// of a standard Xyce output file with the same .print arguments.
  virtual void outputFieldNames(std::vector<std::string> & outputNames) = 0;

  /// Deliver a vector of real outputs to the external code
  ///
  /// This is the output routine that will be called when
  /// doing transient or DC simulations.  
  virtual void outputReal(std::vector<double> & outputData) {};

  /// Deliver a vector of complex outputs to the external code
  ///
  /// This would be usable to deliver AC or NOISE output
  virtual void outputComplex(std::vector<complex> & outputData) {};


  // Other output functions may be desirable here for specific analysis types
  // like HB, which might need to deliver both time- and frequency-domain
  // data in a single call.

  // These functions are used to signal special events

  /// Function to signal that Xyce has completed all output
  virtual void finishOutput() {};

  /// Function to signal that Xyce is doing .STEP processing, and has
  /// begun stepNumber out of maxStep steps.
  virtual void newStepOutput(int stepNumber, int maxStep) {};

  /// Function to signal that Xyce has completed all steps of a .step
  /// analysis.
  virtual void finishedStepping() {};

};
}
}

#endif
