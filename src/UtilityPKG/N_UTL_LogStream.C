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
// Purpose        : Describe the purpose of the contents of the file. If the
//                  contains the header file of a class, provide a clear
//                  description of the nature of the class.
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : David Baur
//
// Creation Date  : 3/28/2013
//
//
//
//
//-------------------------------------------------------------------------

/**
 * @file   N_UTL_LogStream.C
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 06:38:04 2013
 *
 * @brief  Logging and diagnostic stream management
 *
 * There are three output stream to help manage logging for a multi-processor run.  These
 * streams provide basic logging, diagnostic logging and message forwarding.  The logging
 * stream provide a few features for massively parallel processing using Tee and
 * Indentation stream buffers.
 *
 * lout() - Regular logging output.
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), writes to a tee stream buffer with std::cout
 *     added on rank 0 and nothing added to other rank processors.
 *
 *   - If openLogFile(path, false), writes to path on rank 0 and nothing on the
 *     other rank processors.
 *
 *   - If outputLogFile(path, true), writes to path on rank 0 and to path.<rank>.<size>
 *     on other rank processors.
 *
 *
 * dout() - Diagnostic output
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), write to indent stream buffer which then
 *     writes to lout().
 *
 *
 * pout() - Log messages to be forwarded to the rank 0 processor.
 *
 *   - Before initializeLogStream(rank, size), writes to std::cout.
 *
 *   - After initializeLogStream(rank, size), writes to stream buffer with lout() added on
 *     all processors.
 *
 *   - When pout(comm) is called, all messages for each processor are forwarded to the
 *     rank 0 processor and written to the log file, clearing the backlog.
 *
 *
 * The end result is logging happens on all processors.  The rank 0 processor is the
 * regular log output.  If per processor logging is enabled, then each processor gets log
 * output.  Diagnostic messages are written to the log file for each processor.  Processor
 * unique messages are written to pout() and logged when communication is possible and
 * makes sense.
 *
 */

// NEVER, leave this out!
#include <Xyce_config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_LogStream.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_TeeStreamBuf.h>

#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

namespace Xyce {

const char *section_divider  = "------------------------------------------------------------";
const char *subsection_divider = "----------------------------------------";

namespace {

int                             s_rank = 0;                                     ///< Processor rank
int                             s_size = 1;                                     ///< Number of processors
Util::tee_streambuf             s_loutStreambuf;                                ///< Tee'd stream buffer for logging output
Util::indent_streambuf          s_doutStreambuf(std::cout.rdbuf());             ///< Indentation stream buffer for diagnostic output
Util::tee_streambuf             s_poutStreambuf;                                ///< Processor zero forwarder
std::ofstream *                 s_logFileStream;                                ///< Log file output stream
std::ofstream *                 s_diagnosticFileStream;                         ///< Diagnostic log file output stream

std::ostream                    s_lout(std::cout.rdbuf());                      ///< Logging output stream
std::ostream                    s_dout(std::cout.rdbuf());                      ///< Diagnostic output stream
std::ostream                    s_pout(std::cout.rdbuf());                      ///< Processor zero forwarded output stream
std::ostringstream              s_poutBacklog;                                  ///< Stream holding messages to be forwarded

const char *section_divider  = "------------------------------------------------------------";
const char *subsection_divider = "----------------------------------------";


/**
 * Initializes output streams and stream buffers.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 06:37:28 2013
 */
struct INIT
{

  /**
   * Connects the output stream via tee and indentation stream buffers.
   *
   * lout() to std::cout, dout() to lout() via indentation stream buffer, and pout() to
   * lout().  pout() is not initialized to write to the backlog.
   * initializeLogStream(rank, size) handles that since that is the indication that
   * parallel processing is happening and that message forwarding will be done.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Thu Oct  3 07:05:15 2013
   */
  INIT()
  {
    s_loutStreambuf.add(&std::cout);
    s_doutStreambuf.redirect(&s_loutStreambuf);
    s_poutStreambuf.add(&s_lout);

    s_lout.rdbuf(&s_loutStreambuf);
    s_dout.rdbuf(&s_doutStreambuf);
    s_pout.rdbuf(&s_poutStreambuf);
  }

  /**
   * Disconnect the output streams from the tee and indentation stream buffers.
   *
   * lout(), dout() and pout() once again write to std::cout.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Thu Oct  3 07:07:49 2013
   */
  ~INIT()
  {
    s_dout.flush();
    s_pout.flush();
    s_lout.flush();

    s_loutStreambuf.clear();
    s_poutStreambuf.clear();

    delete s_logFileStream;

    s_lout.rdbuf(std::cout.rdbuf());
    s_dout.rdbuf(std::cout.rdbuf());
    s_pout.rdbuf(std::cout.rdbuf());

    s_logFileStream = 0;
  }
};

/**
 * Let the C++ abi handle the initialization and cleanup.
 *
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:10:01 2013
 */
inline void init()
{
  static INIT s_init;
}

} // namespace <unnamed>


/**
 * Initialize the log streams based on processor rank and size.
 *
 * On rank 0 processor, lout() write to std::cout, nowhere on other ranked processors.
 * pout() writes to the forward message backlog.
 *
 * @param rank This processor's rank
 * @param size Number of processor's in parallel run
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:10:38 2013
 */
void initializeLogStream(int rank, int size)
{
  init();

  s_rank = rank;
  s_size = size;

  if (rank != 0)
  {
    s_poutStreambuf.add(&s_poutBacklog);
    s_loutStreambuf.remove(&std::cout);
  }
}

/**
 * Logging output stream.
 *
 * @return reference to the logging output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &lout()
{
  init();

  return s_lout;
}

/**
 * Diagnostic output stream.
 *
 * @return reference to the diagnostic output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &dout()
{
  init();

  return s_dout;
}

/**
 * Processor zero forwarding output stream.
 *
 * @return reference to the processor zero forwarding output stream.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:15:07 2013
 */
std::ostream &pout()
{
  init();

  return s_pout;
}

/**
 * Open output file for logging.
 *
 * Open the specified path on processor rank 0 for the lout() stream.  If per_processor is
 * specified, then the other ranked processors lout() writes to the file path.r.s where r
 * is rank and s is size from initializeLogStream(rank, size).
 *
 * @param path Path of the logging file to open
 * @param per_processor true if per processor logging is to occur
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:16:28 2013
 */
bool openLogFile(const std::string &path, bool per_processor)
{
  // Allocate the output stream
  if (path == "cout")
  {
    s_loutStreambuf.remove(&std::cout);
    if (per_processor || s_rank == 0)
      s_loutStreambuf.add(&std::cout);

    return true;
  }

  std::ofstream *ofs = 0;

  if (per_processor && s_rank != 0)
  {
    std::ostringstream perprocess_path;
    perprocess_path << path << "." << s_rank << "." << s_size;

    ofs = new std::ofstream();
    ofs->open(perprocess_path.str().c_str());
    if ( !ofs->good() )
    {
      Report::UserFatal() << "Could not open log file " << perprocess_path.str() << " for output";
    } 
  }
  else if (s_rank == 0)
  {
    ofs = new std::ofstream();
    ofs->open(path.c_str());
    if ( !ofs->good() )
    {
      Report::UserFatal() << "Could not open log file " << path << " for output";
    } 
  }

  if (ofs)
  {
    s_logFileStream = ofs;
    s_loutStreambuf.add(s_logFileStream);
    s_loutStreambuf.remove(&std::cout);
  }

  return ofs != 0;
}

/**
 * Open output file for diagnostic logging.
 *
 * Open the specified path on processor rank 0 for the lout() stream.  If per_processor is
 * specified, then the other ranked processors lout() writes to the file path.r.s where r
 * is rank and s is size from initializeDiagnosticStream(rank, size).
 *
 * @param path Path of the diagnostic file to open
 * @param per_processor true if per processor diagnostic file is to occur
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:16:28 2013
 */
bool openDiagnosticFile(const std::string &path, bool per_processor)
{
  std::ofstream *ofs = 0;

  if (per_processor && s_rank != 0)
  {
    std::ostringstream perprocess_path;
    perprocess_path << path << "." << s_rank << "." << s_size;

    ofs = new std::ofstream();
    ofs->open(perprocess_path.str().c_str());
  }
  else if (s_rank == 0)
  {
    ofs = new std::ofstream();
    ofs->open(path.c_str());
  }

  // If opening the stream was unsuccessful, silently ignore diagnostic logging to the file.
  if (ofs && ofs->fail())
  {
    delete ofs;
    ofs = 0;
  }

  if (ofs)
  {
    s_diagnosticFileStream = ofs;
    s_doutStreambuf.redirect(s_diagnosticFileStream->rdbuf());
    s_loutStreambuf.add(s_diagnosticFileStream);
  }

  return ofs != 0;
}

/**
 * Close output log file.
 *
 * Closes the output log file on all processors.  On processor rank 0, the output is
 * written to cout.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Thu Oct  3 07:23:25 2013
 */
void closeLogFile()
{
  pout().flush();
  dout().flush();
  lout().flush();

  s_loutStreambuf.remove(s_logFileStream);

  if (s_diagnosticFileStream) {
    s_loutStreambuf.remove(s_diagnosticFileStream);

    delete s_diagnosticFileStream;

    s_diagnosticFileStream = 0;
    s_doutStreambuf.redirect(&s_loutStreambuf);
  }

  delete s_logFileStream;

  s_logFileStream = 0;

  if (s_rank != 0)
  {
    s_loutStreambuf.remove(&std::cout);
  }
}

/** 
 * Gather the stored per-processor output to processor rank zero.
 *
 * 
 * @param comm  Communicator to use to gather the messages
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Dec  4 10:15:26 2013
 */
void pout(Parallel::Machine comm)
{
  Parallel::AllWriteString(comm, lout(), s_poutBacklog.str());

  // Clear the backlog
  s_poutBacklog.str("");
  s_poutBacklog.clear();
}

} // namespace Xyce
