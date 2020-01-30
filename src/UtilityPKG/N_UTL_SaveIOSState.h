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
//                  
//
// Special Notes  : 
//                  
//
// Creator        : David Baur
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------

//  Boost io/ios_state.hpp header file  --------------------------------------

//  Copyright 2002, 2005 Daryle Walker.  Use, modification, and distribution
//  are subject to the Boost Software License, Version 1.0.  (See accompanying
//  file LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/io/> for the library's home page.


#ifndef Xyce_IOState_h
#define Xyce_IOState_h

#include <ios>        // for std::ios_base, std::basic_ios, etc.
#include <locale>     // for std::locale
#include <ostream>    // for std::basic_ostream
#include <streambuf>  // for std::basic_streambuf
#include <string>     // for std::char_traits

namespace Xyce {

class ios_flags_saver;

template < typename Ch, class Tr = ::std::char_traits<Ch> >
    class basic_ios_all_saver;

typedef basic_ios_all_saver<char>            ios_all_saver;
typedef basic_ios_all_saver<wchar_t>        wios_all_saver;

//  Basic stream state saver class declarations  -----------------------------//

class ios_flags_saver
{
public:
    typedef ::std::ios_base            state_type;
    typedef ::std::ios_base::fmtflags  aspect_type;

    explicit  ios_flags_saver( state_type &s )
	: s_save_( s ), a_save_( s.flags() )
	{}
    ios_flags_saver( state_type &s, aspect_type const &a )
	: s_save_( s ), a_save_( s.flags(a) )
	{}
    ~ios_flags_saver()
	{ this->restore(); }

    void  restore()
	{ s_save_.flags( a_save_ ); }

private:
    state_type &       s_save_;
    aspect_type const  a_save_;
};


//  Combined stream state saver class (template) declarations  ---------------//


template < typename Ch, class Tr >
class basic_ios_all_saver
{
public:
    typedef ::std::basic_ios<Ch, Tr>  state_type;

    explicit  basic_ios_all_saver( state_type &s )
	: s_save_( s ), a1_save_( s.flags() ), a2_save_( s.precision() )
	, a3_save_( s.width() ), a4_save_( s.rdstate() )
	, a5_save_( s.exceptions() ), a6_save_( s.tie() )
	, a7_save_( s.rdbuf() ), a8_save_( s.fill() )
	, a9_save_( s.getloc() )
	{}

    ~basic_ios_all_saver()
	{ this->restore(); }

    void  restore()
    {
	s_save_.imbue( a9_save_ );
	s_save_.fill( a8_save_ );
	s_save_.rdbuf( a7_save_ );
	s_save_.tie( a6_save_ );
	s_save_.exceptions( a5_save_ );
	s_save_.clear( a4_save_ );
	s_save_.width( a3_save_ );
	s_save_.precision( a2_save_ );
	s_save_.flags( a1_save_ );
    }

private:
    state_type &                            s_save_;
    typename state_type::fmtflags const     a1_save_;
    ::std::streamsize const                 a2_save_;
    ::std::streamsize const                 a3_save_;
    typename state_type::iostate const      a4_save_;
    typename state_type::iostate const      a5_save_;
    ::std::basic_ostream<Ch, Tr> * const    a6_save_;
    ::std::basic_streambuf<Ch, Tr> * const  a7_save_;
    typename state_type::char_type const    a8_save_;
    ::std::locale const                     a9_save_;
};

} // namespace Xyce

#endif  // Xyce_IOState_h

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
