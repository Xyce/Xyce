//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 1/7/2014
//
//
//
//
//-------------------------------------------------------------------------


#ifndef Xyce_N_ERH_Message_h
#define Xyce_N_ERH_Message_h

#include <string>
#include <sstream>
#include <cstddef>

#include <N_ERH_fwd.h>

#include <N_UTL_NetlistLocation.h>

namespace Xyce {
namespace Report {

static const bool PARSABLE_LINE = false;        ///< Message useful for Emacs compilation

///
/// @brief Typedef <b>MessageId</b> defines a message identifier. 
//
/// Message identifiers must be consist from reference to reference,
/// unique for each instance, and yet consist within each instance
/// across multiple processors.  To meet these criteria, the message
/// identifier is implemented as a static memory location.  It must be
/// declared by the application developer as static.  This results in
/// the linker selecting an address for each instance, which never
/// changes from reference to reference, is unique for each instance
/// and the same regardless of executable (assuming the executable is
/// mapped into the same memory location for each process).  In order
/// to remove the pointer-ness of the static memory location, the
/// address is cast to a pointer difference type which is an integral
/// type by subtracting the zero pointer from it.
///
typedef std::ptrdiff_t MessageId;

///
/// @brief Enumeration <b>ThrottleGroup</b> lists defined throttling groups.
///
/// When messages are throttled, the throttle count may be reset at
/// varior points during an application run.  Some throttles defined
/// for the application, while other may be reset at each time step or
/// other interval.  This allows warnings to be repeated at each time
/// step rather than cut off.
///
enum ThrottleGroup
{
  MSG_APPLICATION       = 0,
  MSG_TIME_STEP         = 1,
  MSG_SOLVER            = 2
};

///
/// @brief Class <b>Throttle</b> describes the cutoff limits for a message throttle.
///
///
struct Throttle
{
  ///
  ///Creates a new <b>Throttle</b> instance.
  ///
  /// @param cutoff a <b>size_t</b> value to display before the message is no longer displayed.
  ///
  /// @param group an <b>int</b> value to identify the throttle group that this message belongs to.
  ///
  ///
  Throttle(size_t cutoff, int group)
    : m_cutoff(cutoff),
      m_group(group),
      m_count(0)
  {}

  size_t        m_cutoff;                       ///< Maximum number to display
  int           m_group;                        ///< Throttle group of message
  size_t        m_count;                        ///< Number which have been displayed
};

///
/// @brief Class <b>MessageCode</b> declares a message identifier and
/// throttle characteristics for a message.  THESE MUST BE DECLARED
/// STATIC.
///
/// All messages have an associated message code.  This message code
/// is used to identify a message for throttling and aggregation.
///
/// Message identifiers must be consist from reference to reference,
/// unique for each instance, and yet consist within each instance
/// across multiple processors.  To meet these criteria, the message
/// identifier is implemented as a static memory location.  It must be
/// declared by the application developer as static.  This results in
/// the linker selecting an address for each instance, which never
/// changes from reference to reference, is unique for each instance
/// and the same regardless of executable (assuming the executable is
/// mapped into the same memory location for each process).  In order
/// to remove the pointer-ness of the static memory location, the
/// address is cast to a pointer difference type which is an integral
/// type by subtracting the zero pointer from it.
///
struct MessageCode
{
  ///
  /// Creates a new <b>MessageCode</b> instance.
  ///
  /// @param throttle_cutoff	a <b>size_t</b> value to display before the
  ///                           message is no longer displayed.
  ///
  /// @param throttle_group	an <b>int</b> value to identify the throttle
  ///                           group that this message belongs to.
  ///
  MessageCode(size_t throttle_cutoff = 5, int throttle_group = MSG_APPLICATION)
    : m_id(&m_id - (MessageId *) 0),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  ///
  /// Creates a new <b>MessageCode</b> instance.  Be particularly
  /// careful when using this constructor.  The message_id value must
  /// be the same on all processors for deferred message reporting to
  /// work properly.
  ///
  /// @param message_id a <b>MessageId</b> value of the message id.
  ///                   This value must be the same for each message across all
  ///                   processors.
  ///
  /// @param throttle_cutoff a <b>size_t</b> value to display before
  ///                            the message is no longer displayed.
  ///
  /// @param throttle_group an <b>int</b> value to identify the
  ///                            throttle group that this message
  ///                            belongs to.
  ///
  MessageCode(MessageId message_id, size_t throttle_cutoff, int throttle_group)
    : m_id(message_id),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  static MessageCode    s_defaultMessageCode;   ///< Default message code

  MessageId             m_id;                   ///< Message identifier
  Throttle              m_throttle;             ///< Throttle characteristics
};

typedef std::ostream &(*OStreamFunctionPtr)(std::ostream &);
typedef std::ios_base &(*IOSBaseFunctionPtr)(std::ios_base &);

///
/// @brief Enumeration <b>MessageType</b> declares the global message types.
///
/// Currently warning and error message types are defined.  Additional
/// types may be added after MSG_FATAL.  The MSG_SYMMETRIC bit
/// indicates that the message was generated identically across all
/// processors.
///
/// Note that while the MSG_SYMMETRIC bit is intended to signal that a 
/// given error was generated on all processors, it is the developer's
/// responsibility to assure that the bit is only set when this is in fact
/// true.
///
/// For the FATAL errors, the symmetry or asymmetry of an error is a critical
/// piece of information used by the report handler that causes the code
/// to exit after printing such a message.  In parallel, MPI_Abort must
/// be called, and incorrectly setting a fatal error message to "MSG_SYMMETRIC"
/// when it is really only emitted on one processor, or FAILING to set
/// MSG_SYMMETRIC when the message is emitted on all processors can lead to
/// a hang at exit time.  This happens either because MPI_Abort is never called
/// (the first case) or MPI_Abort being called from ALL processors (the second).
/// In the case of MPI_Abort being called from all processors, the result
/// is generally a random hang, because it can trip a race condition in mpirun.
///
enum MessageType
{
  MSG_TYPE_MASK         = 0x000000FF,           ///< Mask of levels
  MSG_WARNING           =          0,           ///< Message is a warning
  MSG_ERROR             =          1,           ///< Message is an error, but processing can proceed
  MSG_FATAL             =          2,           ///< Message is a fatal error, seg fault imminent
  MSG_EXCEPTION         =          3,           ///< Message is an exception
  MSG_INFORMATION       =          4,           ///< Message is informational
  MSG_DEBUG             =          5,           ///< Message is debug

  MSG_USER              = 0x00000100,
  MSG_DEVEL             = 0x00000200,

#ifdef Xyce_PARALLEL_MPI
  MSG_PARALLEL          = 0x80000000,                   ///< Compiled in a parallel environment
#else
  MSG_PARALLEL          = 0x00000000,                   ///< Compiled in a parallel environment
#endif
  MSG_SYMMETRIC         = MSG_PARALLEL | 0x40000000,    ///< Message is symmetrical, same on all processors
  MSG_DEFERRED          = 0x20000000,                   ///< Message is deferred, forward to 0, aggregated
  MSG_UNUSED0           = 0x10000000,                   ///< Unused

  MSG_TERMINATE         = 0x00010000,

  // user error types:
  USR_FATAL        = MSG_USER | MSG_FATAL | MSG_TERMINATE,
  USR_ERROR        = MSG_USER | MSG_ERROR,
  USR_WARNING      = MSG_USER | MSG_WARNING,
  USR_FATAL_0      = USR_FATAL | MSG_SYMMETRIC,
  USR_ERROR_0      = USR_ERROR | MSG_SYMMETRIC,
  USR_WARNING_0    = USR_WARNING | MSG_SYMMETRIC,

  // developer error types:
  DEV_FATAL        = MSG_DEVEL | MSG_FATAL | MSG_TERMINATE,
  DEV_WARNING      = MSG_DEVEL | MSG_WARNING,
  DEV_FATAL_0      = DEV_FATAL | MSG_SYMMETRIC,
  DEV_WARNING_0    = DEV_WARNING | MSG_SYMMETRIC,
};

///
/// @brief Function <b>operator<<</b> writes the message type name to
/// the output stream.  If the symmetric bit is set, "parallel" is
/// prefixed to the name.
///
/// @param os		a <b>std::ostream</b> reference to the output stream.
///
/// @param message_type a <b>MessageType</b> const reference to write
///                     the name of.
///
/// @return		a <b>std::ostream</b> reference to os
///
std::ostream &operator<<(std::ostream &os, const MessageType &message_type);

///
/// Message creates an error message to be submitted to the error manager
/// on object destruction.
///
/// Operator<< (put-to) can be used to populate a std::ostringstream
/// contained in the message class to make assembly of a meaningful error
/// meessage painless.
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
/// @date   Tue Sep 17 09:28:54 2013
///
class Message
{
public:
  /// Constructor
  Message(MessageType message_type, MessageCode message_code = MessageCode::s_defaultMessageCode);

  ///
  /// @brief Destructor
  ///
  /// On destruction, the text of the message is actually output through the Xyce::Report::report_message function.
  ///
  virtual ~Message();

private:
  /// Disabled copy constructor
  Message(const Message &right);
  /// Disabled assignment operator
  Message &operator=(const Message &);

public:
  /// Add netlist location information to error message
  Message &at(const NetlistLocation &netlist_location) 
  {
    netlistLocation_ = netlist_location;

    return *this;
  }

  /// Add file/line information to error message
  Message &at(const std::string &path, int line_number) 
  {
    netlistLocation_ = NetlistLocation(path, line_number);

    return *this;
  }

  /// Convert message to immediate termination (Fatal) message
  Message &die()
  {
    messageType_ |= MSG_TERMINATE;
    return *this;
  }

  /// Get internal stringstream from message
  std::ostringstream &os()
  {
    return oss_;
  }

  /// General "output to message" operator
  template<class T>
  Message &operator<<(const T &t)
  {   
    os() << t;
    return *this;
  }

  ///  "output to message" operator with ostream pointer arg
  Message &operator<<(OStreamFunctionPtr f)
  {
    f(os());
    return *this;
  }

  ///  "output to message" operator with iosbase pointer arg
  Message &operator<<(IOSBaseFunctionPtr f)
  {
    f(os());
    return *this;
  }

private:
  unsigned            messageType_;
  MessageCode         messageCode_;
  std::ostringstream  oss_;
  NetlistLocation     netlistLocation_;

protected:
  const char *        functionName_;
};

// Special-case message types:
// User messages:  Messages that are the result of user input.
// Developer messages: These are the "this should never happen to a user"
//                     that generally indicate a bug.


///
/// @brief An asymmetric (one-processor) fatal error that must cause
///             an immediate exit.
///
/// Use this version to report a fatal error due to invalid user input
/// resulting in a single processor being unable to proceed.
///
struct UserFatal : public Message
{
  UserFatal(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_FATAL, message_code)
  {}
};

///
/// A symmetric (ALL-processor) fatal error that must cause an
/// immediate exit.  Never call this unless all processors are
/// encountering the same error, generally at a barrier.
///
/// Use this version to report a fatal error due to invalid user input
/// resulting in ALL processors detecting the same error and being unable to
/// proceed.
///
struct UserFatal0 : public Message
{
  UserFatal0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_FATAL_0, message_code)
  {}
};

///
/// An error due to user input that must ultimately abort
/// the run, but which does not require an immediate exit.
///
/// This asymmetric version will emit the error message on any processor
/// that encounters it.
///
/// It is the caller's responsibility to assure that the code exits after
/// the error is emitted.  This is most often done by using safeBarrier.
///
struct UserError : public Message
{
  UserError(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_ERROR, message_code)
  {}
};

///
/// An error due to user input that must ultimately abort the run, but
/// which does not require an immediate exit.
///
/// This symmetric version will emit the error message only from proc 0
/// in parallel, no matter how many processors call the function with the
/// same message.
///
/// It is the caller's responsibility to assure that the code exits after
/// the error is emitted.  This is most often done by using safeBarrier.
///
 struct UserError0 : public Message
{
  UserError0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_ERROR_0, message_code)
  {}
};

///
/// A non-fatal warning due to user input
///
/// This asymmetric version will emit the error message on any
/// processor that encounters it.
///
struct UserWarning : public Message
{
  UserWarning(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_WARNING, message_code)
  {}
};

///
/// A non-fatal warning due to user input
///
/// This symmetric version will emit the warning message only from proc 0,
/// no matter how many processors emit the same message.
///
struct UserWarning0 : public Message
{
  UserWarning0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_WARNING_0, message_code)
  {}
};

///
/// An asymmetric (one-processor) fatal error that must cause an
/// immediate exit.
///
/// Use this function to report a fatal error due to internal errors
/// resulting in a single processor being unable to proceed.  This is used
/// to report "developer errors" not user input errors.
///
struct DevelFatal : public Message
{
  DevelFatal(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_FATAL, message_code)
  {}

  /// Add function name to error message
  Message &in(const char *function_name)
  {
    functionName_ = function_name;

    return *this;
  }
};

///
/// A symmetric (ALL-processor) fatal error that must cause an
/// immediate exit.  Never call this unless all processors are
/// encountering the same error on the same data at the same place.
///
/// Use this function to report a fatal error due to internal errors
/// resulting in ALL processors detecting the same error and being
/// unable to proceed.  This is used to report a BUG, not invalid user
/// input.
///
struct DevelFatal0 : public Message
{
  DevelFatal0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_FATAL_0, message_code)
  {}

  /// Add function name to error message
  Message &in(const char *function_name)
  {
    functionName_ = function_name;

    return *this;
  }
};

// The devel warnings are not used much, and are meant to signal a
// warning condition indicating a possible bug.  It is primarily a debugging
// tool.
struct DevelWarning : public Message
{
  DevelWarning(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_WARNING, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = function_name;

    return *this;
  }
};

struct DevelWarning0 : public Message
{
  DevelWarning0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_WARNING_0, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = function_name;

    return *this;
  }
};

} // namespace Report
} // namespace Xyce

#endif // Xyce_N_ERH_Message_h
