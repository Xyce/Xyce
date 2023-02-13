//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Listener_h
#define Xyce_N_UTL_Listener_h

#include <vector>
#include <algorithm>

namespace Xyce {
namespace Util {

template <class T>
class Notifier;

//-----------------------------------------------------------------------------
// Class         : Listener
// Purpose       :
// Special Notes :
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:07:15 2014
//-----------------------------------------------------------------------------
/// Template class <b>Listener</b> describes an interface to be called
/// by a <b>Notifier</b> to notify the <b>Listener</b> that event
/// <b>T</b> has occurred. 
///
template <class T>
class Listener
{
public:
  Listener() {}
  virtual ~Listener() {}

  virtual void notify(const T &event) = 0;

private:
  Listener(const Listener &);
  Listener &operator=(const Listener &);  
};

//-----------------------------------------------------------------------------
// Class         : ListenerAutoSubscribe
// Purpose       :
// Special Notes :
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:09:31 2014
//-----------------------------------------------------------------------------
///
/// Template class <b>ListenerAutoSubscribe</b> describes an interface
/// to be called by a <b>Notifier</b> to notify the <b>Listener</b> that
/// event <b>T</b> has occurred. 
///
template <class T>
class ListenerAutoSubscribe : public Listener<T>
{
public:
  explicit ListenerAutoSubscribe(Notifier<T> &notifier)
    : m_notifier(notifier)
  {
    m_notifier.subscribe(*this);
  }

  explicit ListenerAutoSubscribe(Notifier<T> *notifier)
    : m_notifier(*notifier)
  {
    m_notifier.subscribe(*this);
  }

  virtual ~ListenerAutoSubscribe() {
    m_notifier.unsubscribe(*this);
  }

private:
  ListenerAutoSubscribe(const ListenerAutoSubscribe &);
  ListenerAutoSubscribe &operator=(const ListenerAutoSubscribe &);
  
private:
  Notifier<T> &		m_notifier;                     ///< Notifier subscribed to
};

//-----------------------------------------------------------------------------
// Class         : Notifier
// Purpose       :
// Special Notes :
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:09:31 2014
//-----------------------------------------------------------------------------
///
/// Maintains a list for subscribers for notification events of type T.
/// When the an event is published, the notify(event) function is called
/// for all listeners.  Unsubscribing results in the zeroing of the
/// pointer so that iterators that may be active are not invalidated.
///
/// @invariant Listener list iterators are not invalidated during
/// notification even if unsubcribe is called during notification
///
template <class T>
class Notifier
{
public:
  typedef std::vector<Listener<T> *> ListenerList;	///< Registered listener list type

  Notifier() {}
  virtual ~Notifier() {}

  void subscribe(Listener<T> &listener) { m_listenerList.push_back(&listener); }

  void unsubscribe(Listener<T> &listener) 
  {
    for (typename ListenerList::iterator it = m_listenerList.begin(); it != m_listenerList.end(); ++it)
      if (*it == &listener)
	(*it) = typename ListenerList::value_type(0);
  }

  void publish(const T &event = T()) 
  {
    for (typename ListenerList::iterator it = m_listenerList.begin(); it != m_listenerList.end(); ++it)
      if (*it)
	(*it)->notify(event);

    m_listenerList.erase(std::remove(m_listenerList.begin(), m_listenerList.end(), typename ListenerList::value_type(0)), m_listenerList.end());
  }

private:
  ListenerList          m_listenerList;	///< List of <b>Listeners</b> to notify
};

//-----------------------------------------------------------------------------
// Function      : subscribe
// Purpose       : 
// Special Notes :
// Scope         : 
// Creator       : Dave Baur
// Creation Date : 2014
//-----------------------------------------------------------------------------
template<class T>
void subscribe(Notifier<T> &notifier, Listener<T> &listener) {
  notifier.subscribe(listener);
}

//-----------------------------------------------------------------------------
// Function      : publish
// Purpose       : 
// Special Notes :
// Scope         : 
// Creator       : Dave Baur
// Creation Date : 2014
//-----------------------------------------------------------------------------
template<class T>
void publish(Notifier<T> &notifier, const T &event) {
  notifier.publish(event);
}

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Listener_h
