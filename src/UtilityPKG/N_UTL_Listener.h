//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

///
/// @addtogroup ListenerNotifierDetail
/// @{
///

template <class T>
class Notifier;

//-----------------------------------------------------------------------------
// Function      : Listener
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:07:15 2014
//-----------------------------------------------------------------------------
///
/// Template class <b>Listener</b> describes an interface to be called
/// by a <b>Notifier</b> to notify the <b>Listener</b> that event
/// <b>T</b> has occurred. 
///
template <class T>
class Listener
{
public:
//-----------------------------------------------------------------------------
// Function      : Listener
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:07:52 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>Listener</b> instance.
///
  Listener()
  {}

//-----------------------------------------------------------------------------
// Function      : ~Listener
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:08:08 2014
//-----------------------------------------------------------------------------
///
/// Destroys a <b>Listener</b> instance.
///
/// @invariant
///
///
/// @return 
///
///
  virtual ~Listener()
  {}

//-----------------------------------------------------------------------------
// Function      : notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:08:28 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>update</b> is the interface member function which
/// the <b>Notifier</b> uses to notify the the Listener that event
/// <b>T</b> has occurred. 
///
/// @invariant
///
/// @param event      a <b>T</b> reference to the event.
///
///
virtual void notify(const T &event) = 0;

private:
  Listener(const Listener &);
  Listener &operator=(const Listener &);  
};

//-----------------------------------------------------------------------------
// Function      : ListenerAutoSubscribe
// Purpose       :
// Special Notes :
// Scope         : public
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
//-----------------------------------------------------------------------------
// Function      : ListenerAutoSubscribe
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:10:06 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>Listener</b> instance.
///
/// @invariant
///
/// @param notifier     Notifier to subscribe to
///
  explicit ListenerAutoSubscribe(Notifier<T> &notifier)
    : m_notifier(notifier)
  {
    m_notifier.subscribe(*this);
  }

//-----------------------------------------------------------------------------
// Function      : ListenerAutoSubscribe
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:10:06 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>Listener</b> instance.
///
/// @invariant
///
/// @param notifier     Notifier to subscribe to
///
  explicit ListenerAutoSubscribe(Notifier<T> *notifier)
    : m_notifier(*notifier)
  {
    m_notifier.subscribe(*this);
  }

//-----------------------------------------------------------------------------
// Function      : ~ListenerAutoSubscribe
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:11:23 2014
//-----------------------------------------------------------------------------
///
/// Destroys a <b>Listener</b> instance.
///
/// @invariant Unsubscribes from the subscribed notifier
///
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
// Class         : ListenerProxy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:12:41 2014
//-----------------------------------------------------------------------------
///
/// Template class <b>ListenerProxy</b> describes an adapter interface
/// to be called by a <b>Notifier</b> to notify the object <b>U</b>
/// contained within <b>ListenerProxy</b> that event <b>T</b> has
/// occurred. 
///
template <class T, class U>
class ListenerProxy : public Listener<T>
{
public:
  typedef void (U::*Function)(const T &event);		///< Function signature

//-----------------------------------------------------------------------------
// Function      : ListenerProxy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:13:14 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>Listener</b> instance.
///
/// @invariant
///
/// @param object       Object to call member function on notication
/// @param function     Member function to call on notication
///
  ListenerProxy(U &object, Function function)
    : m_object(object),
      m_function(function)
  {}

//-----------------------------------------------------------------------------
// Function      : ~ListenerProxy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:14:04 2014
//-----------------------------------------------------------------------------
///
/// Destroys a <b>Listener</b> instance.
///
  virtual ~ListenerProxy()
  {}

//-----------------------------------------------------------------------------
// Function      : notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:14:26 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>update</b> is the interface member function which
/// the <b>Notifier</b> uses to notify the the Listener that event
/// <b>T</b> has occurred. 
///
/// @param event a <b>T</b> reference to the event.
///
///
  virtual void notify(const T &event) {
    (m_object.*m_function)(event);
  }

private:
  U &		m_object;                       ///< Object 
  Function	m_function;                     ///< Function to call
};


//-----------------------------------------------------------------------------
// Function      : void
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:16:35 2014
//-----------------------------------------------------------------------------
///
/// Template class <b>ListenerAdapter</b> describes an adapter interface
/// to be called by a <b>Notifier</b> to notify the object <b>U</b>
/// contained within <b>ListenerAdapter</b> that event <b>T</b> has
/// occurred. 
///
template <class T, class U, typename F = void (U::*)()>
class ListenerAdapter : public Listener<T>
{
public:
  typedef F Function;				///< Function signature

//-----------------------------------------------------------------------------
// Function      : ListenerAdapter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:16:59 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>ListenerAdapter</b> instance.
///
/// @param object       Object to call member function on notication
/// @param function     Member function to call on notication
///
  ListenerAdapter(U &object, Function function)
    : m_object(object),
      m_function(function)
  {}

//-----------------------------------------------------------------------------
// Function      : ~ListenerAdapter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:18:00 2014
//-----------------------------------------------------------------------------
///
/// Destroys a <b>ListenerAdapter</b> instance.
///
/// @invariant
///
///
/// @return 
///
///
  virtual ~ListenerAdapter()
  {}

//-----------------------------------------------------------------------------
// Function      : notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:18:22 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>update</b> is the interface member function which
/// the <b>Notifier</b> uses to notify the the Listener that event
/// <b>T</b> has* occurred. 
///
/// @param event a <b>T</b> reference to the event.
///
///
  virtual void notify(const T &event) {
    (m_object.*m_function)();
  }

private:
  U &		m_object;                       ///< Object 
  Function	m_function;                     ///< Function to call
};


//-----------------------------------------------------------------------------
// Function      : Notifier
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:19:25 2014
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

//-----------------------------------------------------------------------------
// Function      : Notifier
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:23:19 2014
//-----------------------------------------------------------------------------
///
/// Creates a new <b>Notifier</b> instance.
///
  Notifier()
  {}

//-----------------------------------------------------------------------------
// Function      : ~Notifier
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:23:39 2014
//-----------------------------------------------------------------------------
///
/// Destroys a <b>Notifier</b> instance.
///
  virtual ~Notifier()
  {}

  /**
   * @brief Member function <b>subscribe</b> registers the <b>Listener</b> with
   * the <b>Notifier</b> so that it may receive notification messages.
   *
   * @param listener	an <b>Listener</b> reference which is requesting
   *			notification from the Notifier when event <b>T</b> occurs.
   *
   */
//-----------------------------------------------------------------------------
// Function      : subscribe
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:23:58 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>subscribe</b> registers the <b>Listener</b> with
/// the <b>Notifier</b> so that it may receive notification messages. 
///
/// @param listener 	an <b>Listener</b> reference which is requesting
///                     notification from the Notifier when event <b>T</b> occurs.
///
///
  void subscribe(Listener<T> &listener) {
    m_listenerList.push_back(&listener);
  }

//-----------------------------------------------------------------------------
// Function      : unsubscribe
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:24:50 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>unsubscribe</b> unregisters the <b>Listener</b>
/// from the <b>Notifier</b> so that is will not longer recieve
/// messages. 
///
/// @invariant
///
/// @param listener 	an <b>Listener</b> reference which no longer wishes to
///			receive notification when event <b>T</b> occurs.
///
///
  void unsubscribe(Listener<T> &listener) {
    for (typename ListenerList::iterator it = m_listenerList.begin(); it != m_listenerList.end(); ++it)
      if (*it == &listener)
	(*it) = typename ListenerList::value_type(0);
  }

//-----------------------------------------------------------------------------
// Function      : publish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Jul 15 13:25:20 2014
//-----------------------------------------------------------------------------
///
/// Member function <b>publish</b> sends the event to all of the
/// registered listeners. 
///
/// @invariant
///
/// @param event 	a <b>T</b> const reference to the event describing the
///			notification.
///
///
  void publish(const T &event = T()) {
    for (typename ListenerList::iterator it = m_listenerList.begin(); it != m_listenerList.end(); ++it)
      if (*it)
	(*it)->notify(event);

    m_listenerList.erase(std::remove(m_listenerList.begin(), m_listenerList.end(), typename ListenerList::value_type(0)), m_listenerList.end());
  }

private:
  ListenerList          m_listenerList;	///< List of <b>Listeners</b> to notify
};

template<class T>
void subscribe(Notifier<T> &notifier, Listener<T> &listener) {
  notifier.subscribe(listener);
}

template<class T>
void publish(Notifier<T> &notifier, const T &event) {
  notifier.publish(event);
}

///
/// @}
///

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Listener_h
