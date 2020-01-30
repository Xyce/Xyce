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
// Purpose        : Stream buffer that performs indentationx
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_TeeStreambuf_h
#define Xyce_N_UTL_TeeStreambuf_h

#include <string>
#include <streambuf>
#include <ostream>
#include <set>
#include <map>

namespace Xyce{
namespace Util {

/**
 * @brief Maintains a list of destination output streams to send written characters to.
 *
 * Many destination output streams may be added.  For each character written to this stream buffer, the same character
 * is written to the stream buffer of each destination stream.
 *
 */
template<class Ch, class Tr = std::char_traits<Ch> >
class basic_tee_streambuf : public std::basic_streambuf<Ch, Tr>
{
  typedef std::set<std::ostream *> StreamSet;
  typedef std::map<std::ostream *, int> StreamErrorMap;

public:
  /**
   * Creates a new <b>basic_tee_streambuf</b> instance.
   *
   */
  basic_tee_streambuf()
  {}

  /**
   * Creates a new <b>basic_tee_streambuf</b> instance and adds the specified destination output
   * stream.
   *
   */
  explicit basic_tee_streambuf(std::basic_ostream<Ch, Tr> *os) 
  {
    add(os);
  }

  /**
   * Destroys a <b>basic_tee_streambuf</b> instance.
   *
   */
  virtual ~basic_tee_streambuf()
  {}

  /**
   * @brief Member function <b>eof</b> returns the current end-of-file status.
   *
   * @return			an <b>int</b> value of the current end-of-file status.
   */
  int eof() 
  {
    return std::basic_streambuf<Ch, Tr>::traits_type::eof();
  }

  /**
   * @brief Member function <b>add</b> adds the specified destination output stream.
   *
   * @param sb			a <b>std::stream</b> pointer to the output stream to add.
   *
   */
  void add(std::ostream *os) 
  {
    m_destinations.insert(os);
  }

  /**
   * @brief Member function <b>remove</b> removes the specified destination output stream.
   *
   * @param sb			a <b>std::stream</b> pointer to the output stream buffer to
   *                            remove.
   *
   */
  void remove(std::ostream *os) 
  {
    m_destinations.erase(os);
  }

  /**
   * @brief Member function <b>clear</b> removes all destination output streams.
   *
   */
  void clear() 
  {
    m_destinations.clear();
  }

private:
  /**
   * @brief Member function <b>sync</b> syncs the destination stream buffers of each output stream.
   *
   * @return			an <b>int</b> value of 1 if successful.
   */
  virtual int sync() 
  {
    if (m_destinations.empty())
      return 1;

    StreamErrorMap return_code;

    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) 
    {
      if ((*it)->rdbuf() != this) 
      {
        int ret = (*it)->rdbuf()->pubsync();
        return_code[*it] = ret;
      }
    }

    // Remove streambufs with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it)
      if (return_code[*it] == eof())
        m_destinations.erase(it);

    if (m_destinations.empty())
      return 1;

    return 1;
  }

  /**
   * @brief Member function <b>overflow</b> writes the specified character to all the destination
   * output stream buffers.
   *
   * @param c			an <b>int</b> const value of the character to write.
   *
   * @return			an <b>int</b> value of the character written.
   */
  virtual typename std::basic_streambuf<Ch, Tr>::int_type overflow(const int c) 
  {
    if (m_destinations.empty())
      return 1;

    StreamErrorMap return_code;

    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) 
    {
      int ret = (*it)->rdbuf()->sputc(c);
      return_code[*it] = ret;
    }

    // Remove streambufs with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it)
      if (return_code[*it] == eof())
        m_destinations.erase(it);

    if (m_destinations.empty())
      return 1;

    return 1;
  }

  /**
   * @brief Member function <b>xsputn</b> writes the specified characters to all the destination
   *
   * @param buffer		a <b>char</b> const pointer to the character string to write.
   *
   * @param n			a <b>std::streamsize</b> value of the number of characters to write.
   *
   * @return			a <b>std::streamsize</b> value of the number of characters written.
   */
  virtual std::streamsize xsputn(char const *buffer, std::streamsize n) 
  {
    if (m_destinations.empty())
      return n;

    StreamErrorMap return_code;

    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) 
    {
      std::ostream *os = (*it);
      int ret = os->rdbuf()->sputn(buffer, n);
      return_code[*it] = ret;
    }

    // Remove ostreams with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) 
    {
      if (return_code[*it] < 0) 
      {
        m_destinations.erase(it);
      }
    }

    if (m_destinations.empty())
      return n;

    return n;
  }

private:
  StreamSet             m_destinations;    ///< Destination output streams to write to
};

typedef Xyce::Util::basic_tee_streambuf<char> tee_streambuf;

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_TeeStreambuf_h
