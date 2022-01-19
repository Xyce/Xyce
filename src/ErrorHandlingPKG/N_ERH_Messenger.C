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
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 1/7/2014
//
//-------------------------------------------------------------------------


#include <Xyce_config.h>

#include <list>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>

#include <N_ERH_Messenger.h>
#include <N_ERH_Message.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_ReportHandler.h>
#include <N_PDS_ParallelMachine.h>

namespace Xyce {
namespace Report {

namespace {

void bootstrap()
{
  register_message_type(MSG_WARNING, 10000000, "warning");
  register_message_type(MSG_ERROR, 10000000, "error");
  register_message_type(MSG_FATAL, 10000000, "fatal");
  register_message_type(MSG_EXCEPTION, 1000000, "exception");
  register_message_type(MSG_INFORMATION, 1000000, "informational");
}

typedef std::pair<MessageId, std::string> MessageKey;

typedef std::map<MessageKey, Throttle> MessageIdMap;

MessageIdMap s_messageIdMap;

MessageIdMap s_deferredMessageIdMap;

struct DeferredMessage
{
  DeferredMessage()
  {}

  DeferredMessage(
    size_t              type,
    MessageId           message_id,
    size_t              throttle_cutoff,
    int                 throttle_group,
    const std::string & header,
    const std::string & aggregate)
    : m_type(type),
      m_messageId(message_id),
      m_rank(0),
      m_throttleCutoff(throttle_cutoff),
      m_throttleGroup(throttle_group),
      m_header(header),
      m_aggregate(aggregate)
  {}

  size_t                m_type;
  MessageId             m_messageId;
  int                   m_rank;
  size_t                m_throttleCutoff;
  int                   m_throttleGroup;
  std::string           m_header;
  std::string           m_aggregate;
};

typedef std::vector<DeferredMessage> DeferredMessageVector;

struct DeferredMessageLess : public std::binary_function<DeferredMessage, DeferredMessage, bool>
{
  bool operator()(const DeferredMessage &key_1, const DeferredMessage &key_2) const
  {
    return (key_1.m_type < key_2.m_type)
      || (!(key_2.m_type < key_1.m_type) && key_1.m_messageId < key_2.m_messageId)
      || (!(key_2.m_type < key_1.m_type) && !(key_2.m_messageId < key_1.m_messageId) && key_1.m_header < key_2.m_header);
  }
};

DeferredMessageVector s_deferredMessageVector;

struct MessageTypeInfo
{
  MessageTypeInfo()
    : m_count(0),
      m_maxCount(10000000),
      m_name("unknown")
  {}

  unsigned              m_count;
  unsigned              m_maxCount;
  std::string           m_name;
};

typedef std::map<unsigned, MessageTypeInfo> MessageTypeInfoMap;

MessageTypeInfoMap s_messageTypeInfo;

MessageTypeInfo &
get_message_type_info(
  unsigned              type)
{
  if (s_messageTypeInfo.empty())
    bootstrap();
  
  MessageTypeInfoMap::iterator it = s_messageTypeInfo.find(type & MSG_TYPE_MASK);
  if (it != s_messageTypeInfo.end())
    return (*it).second;
  else
    return s_messageTypeInfo[type & MSG_TYPE_MASK];
}

std::ostream &
operator<<(
  std::ostream &        os,
  const MessageType &   message_type)
{
  os << get_message_type_info(message_type).m_name;

  return os;
}

enum CutoffStatus
  {
    MSG_DISPLAY           = 0,
    MSG_CUTOFF            = 1,
    MSG_CUTOFF_EXCEEDED   = 2
  };


CutoffStatus
count_message(
  MessageId             message_id,
  const char *          message,
  const Throttle &      throttle)
{
  std::pair<MessageIdMap::iterator, bool> res = s_messageIdMap.insert(MessageIdMap::value_type(MessageIdMap::key_type(message_id, message), throttle));
  size_t count = ++(*res.first).second.m_count;

  if (count < (*res.first).second.m_cutoff)
    return MSG_DISPLAY;
  else if (count == (*res.first).second.m_cutoff)
    return MSG_CUTOFF;
  else
    return MSG_CUTOFF_EXCEEDED;
}

Util::Marshal &operator<<(Util::Marshal &mout, const DeferredMessage &s)
{
  mout << s.m_type << s.m_messageId << s.m_rank << s.m_throttleGroup << s.m_throttleCutoff << s.m_header << s.m_aggregate;
  return mout;
}

Util::Marshal &operator>>(Util::Marshal &min, DeferredMessage &s)
{
  min >> s.m_type >> s.m_messageId >> s.m_rank >> s.m_throttleGroup >> s.m_throttleCutoff >> s.m_header >> s.m_aggregate;
  return min;
}

} // namespace <empty>


void
register_message_type(
  unsigned              message_type,
  unsigned              max_count,
  const char *          name)
{
  MessageTypeInfo &message_info = s_messageTypeInfo[message_type];

  message_info.m_maxCount = max_count;
  message_info.m_name = name;
}


unsigned
get_message_count(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_count;
}


unsigned
increment_message_count(
  unsigned              message_type)
{
  return ++get_message_type_info(message_type).m_count;
}


void
reset_message_count(
  unsigned              message_type)
{
  get_message_type_info(message_type).m_count = 0;
}

void
reset_message_counts()

{
  reset_message_count(MSG_WARNING);
  reset_message_count(MSG_ERROR);
  reset_message_count(MSG_FATAL);
  reset_message_count(MSG_EXCEPTION);
  reset_message_count(MSG_INFORMATION);
}

const std::string &
get_message_name(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_name;
}


void
set_max_message_count(
  unsigned              message_type,
  unsigned              max_count)
{
  get_message_type_info(message_type).m_maxCount = max_count;
}


unsigned
get_max_message_count(
  unsigned              message_type)
{
  return get_message_type_info(message_type).m_maxCount;
}


void
report_message(
  const char *          message,
  unsigned              message_type,
  const MessageCode &   message_code)
{
  if (message_type & MSG_DEFERRED)
    Xyce::report(message, message_type);

  else
  {
    unsigned count = increment_message_count(message_type);
    unsigned max_count = get_max_message_count(message_type);

    if (count == max_count)
    {
      Xyce::report(message, message_type);

      std::ostringstream s;
      s << "Maximum " << get_message_name(message_type) << " count has been exceeded and will no longer be displayed";
      Xyce::report(s.str().c_str(), MSG_WARNING | MSG_SYMMETRIC);
    }

    else if (count < max_count)
    {
      CutoffStatus cutoff = count_message(message_code.m_id, "", message_code.m_throttle);

      if (cutoff == MSG_CUTOFF)
      {
        Xyce::report(message, message_type);

        std::ostringstream s;
        s << "Maximum count for this " << get_message_name(message_type) << " has been exceeded and will no longer be displayed";
        Xyce::report(s.str().c_str(), MSG_WARNING | MSG_SYMMETRIC);
      }

      else if (cutoff == MSG_DISPLAY)
        Xyce::report(message, message_type);
    }
  }
}


void
reset_throttle_group(
  int                   throttle_group)
{
  for (MessageIdMap::iterator it = s_messageIdMap.begin(); it != s_messageIdMap.end(); ++it)
    if ((*it).second.m_group == throttle_group)
      (*it).second.m_count = 0;
}


void
add_deferred_message(
  int                   message_type,
  MessageId             message_id,
  size_t                throttle_cutoff,
  int                   throttle_group,
  const char *          header,
  const char *          aggegrate)
{
  std::ostringstream s;
  s << header << " " << aggegrate;

  Xyce::report(s.str().c_str(), message_type | MSG_DEFERRED);

  std::pair<MessageIdMap::iterator, bool> res = s_deferredMessageIdMap.insert(MessageIdMap::value_type(MessageIdMap::key_type(message_id, header), Throttle(throttle_cutoff, throttle_group)));
  size_t count = ++(*res.first).second.m_count;

  if (count <= throttle_cutoff)
    s_deferredMessageVector.push_back(DeferredMessage(message_type, message_id, throttle_cutoff, throttle_group, header, aggegrate));
}


/// @todo REFACTOR Should the message counts be broadcast to the other processors?

void
report_deferred_messages(
  Parallel::Machine       comm)
{
#ifdef Xyce_PARALLEL_MPI
  const int p_root = 0 ;
  int p_size = Parallel::size(comm);
  int p_rank = Parallel::rank(comm);

  for (DeferredMessageVector::iterator it = s_deferredMessageVector.begin(); it != s_deferredMessageVector.end(); ++it)
    (*it).m_rank = p_rank;

  Util::Marshal mout;
  mout << s_deferredMessageVector;

  DeferredMessageVector deferred_message_vector;

  // Gather the send counts on root processor
  std::string send_string(mout.stream.str());
  int send_count = send_string.size();
  std::vector<int> recv_count(p_size, 0);
  int * const recv_count_ptr = &recv_count[0] ;

  int result = MPI_Gather(&send_count, 1, MPI_INT,
                          recv_count_ptr, 1, MPI_INT,
                          p_root, comm);
  if (MPI_SUCCESS != result)
  {
    std::ostringstream message ;
    message << "stk::report_deferred_messages FAILED: MPI_Gather = " << result ;
    throw std::runtime_error(message.str());
  }

  // Receive counts are only non-zero on the root processor:
  std::vector<int> recv_displ(p_size + 1, 0);

  for (int i = 0 ; i < p_size ; ++i)
  {
    recv_displ[i + 1] = recv_displ[i] + recv_count[i] ;
  }

  const int recv_size = recv_displ[p_size] ;

  std::vector<char> buffer(recv_size);

  {
    const char * const send_ptr = send_string.data();
    char * const recv_ptr = recv_size ? & buffer[0] : (char *) NULL ;
    int * const recv_displ_ptr = & recv_displ[0] ;

    result = MPI_Gatherv((void *) send_ptr, send_count, MPI_CHAR,
                         recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                         p_root, comm);
    if (MPI_SUCCESS != result)
    {
      std::ostringstream message ;
      message << "stk::report_deferred_messages FAILED: MPI_Gatherv = " << result ;
      throw std::runtime_error(message.str());
    }


    if (p_rank == p_root)
    {
      for (int i = 0; i < p_size; ++i)
      {
        Util::Marshal min(std::string(recv_ptr + recv_displ[i], recv_ptr + recv_displ[i + 1]));
        min >> deferred_message_vector;
      }

      std::stable_sort(deferred_message_vector.begin(), deferred_message_vector.end(), DeferredMessageLess());

      DeferredMessageVector::const_iterator current_message_it = deferred_message_vector.begin();
      while (current_message_it != deferred_message_vector.end())
      {
        const DeferredMessage &current_message = (*current_message_it);

        DeferredMessageVector::const_iterator end = current_message_it + 1;
        while (end != deferred_message_vector.end()
               && current_message.m_messageId == (*end).m_messageId
               && current_message.m_header == (*end).m_header)
          ++end;

        std::ostringstream s;

        s << current_message.m_header << current_message.m_aggregate;

        for (DeferredMessageVector::const_iterator it1 = current_message_it + 1; it1 != end; ++it1)
        {
          bool print = true;
          for (DeferredMessageVector::const_iterator it2 = current_message_it; it2 != it1; ++it2)
            if ((*it1).m_aggregate == (*it2).m_aggregate)
            {
              print = false;
              break;
            }
          if (print)
          {
            if (!(*it1).m_aggregate.find('\n'))
              s << ", ";
            s << (*it1).m_aggregate;
          }
        }

        report_message(s.str().c_str(), current_message.m_type | MSG_SYMMETRIC, MessageCode(current_message.m_messageId, current_message.m_throttleCutoff, current_message.m_throttleGroup));

        current_message_it = end;
      }
    }
  }

  s_deferredMessageIdMap.clear();
  s_deferredMessageVector.clear();
#endif
}


void
aggregate_messages(
  Parallel::Machine       comm,
  std::ostringstream &  os,
  const char *          separator)
{
#ifdef Xyce_PARALLEL_MPI
  std::string message = os.str();
  os.str("");

  const int p_root = 0 ;
  int p_size = Parallel::size(comm);
  int p_rank = Parallel::rank(comm);

  int result ;

  // Gather the send counts on root processor

  int send_count = message.size();

  std::vector<int> recv_count(p_size, 0);

  int * const recv_count_ptr = & recv_count[0] ;

  result = MPI_Gather(& send_count, 1, MPI_INT,
                      recv_count_ptr, 1, MPI_INT,
                      p_root, comm);

  if (MPI_SUCCESS != result)
  {
    std::ostringstream s;
    s << "stk::all_write FAILED: MPI_Gather = " << result ;
    throw std::runtime_error(s.str());
  }

  // Receive counts are only non-zero on the root processor:
  std::vector<int> recv_displ(p_size + 1, 0);

  for (int i = 0 ; i < p_size ; ++i)
  {
    recv_displ[i + 1] = recv_displ[i] + recv_count[i] ;
  }

  const int recv_size = recv_displ[ p_size ] ;

  std::vector<char> buffer(recv_size);

  {
    const char * const send_ptr = message.c_str();
    char * const recv_ptr = recv_size ? & buffer[0] : (char *) NULL ;
    int * const recv_displ_ptr = & recv_displ[0] ;

    result = MPI_Gatherv((void*) send_ptr, send_count, MPI_CHAR,
                         recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                         p_root, comm);
  }

  if (MPI_SUCCESS != result)
  {
    std::ostringstream s ;
    s << "stk::all_write FAILED: MPI_Gatherv = " << result ;
    throw std::runtime_error(s.str());
  }

  if (p_root == (int) p_rank)
  {
    bool first = true;
    for (int i = 0 ; i < p_size ; ++i)
    {
      if (recv_count[i])
      {
        if (!first)
          os << separator;
        first = false;
        char * const ptr = & buffer[ recv_displ[i] ];
        os.write(ptr, recv_count[i]);
      }
    }
    os.flush();
  }
  else
    os << message;
#endif
}

} // namespace Report
} // namespace Xyce
