//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  N_UTL_PARAM_H
#define  N_UTL_PARAM_H

#include <stdexcept>
#include <iosfwd>
#include <string>
#include <vector>
#include <typeinfo>
#include <iostream>
#include <unordered_set>
using std::unordered_set;

#include <N_UTL_fwd.h>
#include <N_UTL_Pack.h>
#include <N_UTL_NoCase.h>

namespace Xyce {
namespace Util {

class Expression;


/** 
 * Parameter type enumeration.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 15:02:32 2013
 */

enum { STR, DBLE, INT, LNG, EXPR, BOOL, STR_VEC, INT_VEC, DBLE_VEC, DBLE_VEC_IND, COMPOSITE };

}

/** 
 * DataTypeTraits casts a C++ data type to the corresponding enumerated value.
 *
 * 
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 14:55:02 2013
 */
template<>
struct DataTypeTrait<std::string>
{
  enum {type = Util::STR};
};

template<>
struct DataTypeTrait<double>
{
  enum {type = Util::DBLE};
};

template<>
struct DataTypeTrait<int>
{
  enum {type = Util::INT};
};

template<>
struct DataTypeTrait<long>
{
  enum {type = Util::LNG};
};

template<>
struct DataTypeTrait<bool>
{
  enum {type = Util::BOOL};
};

template<>
struct DataTypeTrait<std::vector<std::string> >
{
  enum {type = Util::STR_VEC};
};

template<>
struct DataTypeTrait<std::vector<int> >
{
  enum {type = Util::INT_VEC};
};

template<>
struct DataTypeTrait<std::vector<double> >
{
  enum {type = Util::DBLE_VEC};
};


template<>
struct DataTypeTrait<Util::Expression>
{
  enum {type = Util::EXPR};
};


template<>
struct DataTypeTrait<Util::Param>
{
  enum {type = -1};
};

namespace Util {

/** 
 * Print the parameters value to the output stream.
 *
 * @param os output stream
 * @param value value to print
 *
 * @return output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 14:56:12 2013
 */
template<class T>
inline std::ostream &printValue(std::ostream &os, const T &value) 
{
  return os << value;
}

/** 
 * Print the vector parameters value to the output stream as a comma separated list.
 *
 * @param os output stream
 * @param vector vector of values to print
 *
 * @return output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 14:56:12 2013
 */
template<class T>
inline std::ostream &printValue(std::ostream &os, const std::vector<T> &vector) 
{
  for (typename std::vector<T>::const_iterator it = vector.begin(); it != vector.end(); ++it) 
  {
    if (it != vector.begin())
      os << ", ";
    os << *it;
  }
  
  return os;
}

/** 
 * Print an expressions value (current does nothing) to the output stream.
 *
 * 
 * @param os output stream
 * @param expression expression to print
 *
 * @return 
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 14:58:01 2013
 */
inline std::ostream &printValue(std::ostream &os, const Expression &expression)
{
  return os;
}

template<class T>
class ParamData;

/** 
 * ParamData<void> is the base class for storing a parameters value.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 14:58:58 2013
 */template<>
 class ParamData<void>
   {
   public:
     ParamData()
     {}

   private:
     ParamData(const ParamData &param_data);
     ParamData &operator=(const ParamData &param_data);

   public:
     virtual ~ParamData()
     {}

     /** 
      * Returns the typeid of the parameter's value.
      *
      * 
      *
      * @return 
      *
      * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
      * @date   Tue Dec 17 15:00:24 2013
      */
     virtual const std::type_info &type() const = 0;

     /** 
      * Returns the type enumeration value of the parameter's value.
      *
      * @return type enumeration value
      *
      * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
      * @date   Tue Dec 17 15:01:08 2013
      */
     virtual int enumType() const = 0;

     /** 
      * Return the parameter's value converted to a string.
      *
      * @return string of the parameter's value
      *
      * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
      * @date   Tue Dec 17 15:02:01 2013
      */
     virtual std::string stringValue() = 0;

     /** 
      * Deep copy of the parameter's value.
      *
      * @return copy of the parameter's value
      *
      * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
      * @date   Tue Dec 17 15:03:44 2013
      */
     virtual ParamData *clone() = 0;

     /** 
      * Returns true if the parameter's value is of the specifed type.
      *
      * @return true if the parameter's value is of the specifed type.
      *
      * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
      * @date   Tue Dec 17 15:04:29 2013
      */
     template <class U>
     bool isType()
     {
       return typeid(U) == type();
     }
 };


/** 
 * Stores a parameter value of any type.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Dec 17 15:05:08 2013
 */
template <class T>
class ParamData : public ParamData<void>
{
public:
  ParamData()
    : value_(),
      enumType_(-1)
  {}

  ParamData(const T &t)
    : enumType_(DataTypeTrait<T>::type),
      value_(t)
  {}

  virtual ~ParamData()
  {}

private:
  ParamData(const ParamData &param_data);
  ParamData &operator=(const ParamData &param_data);

public:
  virtual const std::type_info &type() const
  {
    return typeid(T);
  }

  virtual int enumType() const
  {
    return enumType_;
  }

  T &getValue()
  {
    return value_;
  }

  const T &getValue() const
  {
    return value_;
  }

  virtual ParamData *clone()
  {
    return new ParamData(value_);
  }

  virtual std::string stringValue()
  {
    std::ostringstream oss;
    printValue(oss, value_);
    return oss.str();
  }

private:
  int         enumType_;                      ///< Parameter's type enumeration value
  T           value_;                         ///< Parameter's value
};

//-----------------------------------------------------------------------------
// Class         : Param
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
class Param
{
  friend class Pack<Param>;
  friend std::ostream &operator<<(std::ostream & os, const Param & p);

public:
  Param()
    : tag_(),
      data_(0)
  {}

  template <class T>
  Param( const std::string & t, const T & v )
    : tag_(t),
      data_(new ParamData<T>(v))
  {}

  Param(const std::string & t, const char * v)
    : tag_(t),
      data_(new ParamData<std::string>(std::string(v)))
  {}

  Param(const Param &rhsParam)
    : tag_(rhsParam.tag_),
      data_(rhsParam.data_ ? rhsParam.data_->clone() : 0)
  {}


  Param &operator=(const Param & rhsParam)
  {
    if (this != &rhsParam)
    {
      tag_ = rhsParam.tag_;
      delete data_;
      data_ = rhsParam.data_ ? rhsParam.data_->clone() : 0;
    }

    return *this;
  }

  virtual ~Param()
  {
    delete data_;
  }

  void setTag(const std::string & tag) 
  {
    tag_ = tag;
  }

  void setVal(const Param &param)
  {
    delete data_;
    data_ = param.data_->clone();
  }

  template <class T>
  void setVal(const T &t) 
  {
    delete data_;
    data_ = new ParamData<T>(t);
  }

  void setVal(const char *t) 
  {
    delete data_;
    data_ = new ParamData<std::string>(std::string(t));
  }

  void setVal(const ExtendedString &t);

  // Methods to set value.
  template <class T>
  Param &set(const std::string & tag, const T &val) 
  {
    tag_ = tag;
    setVal(val);
    return *this;
  }

  template<class T>
  T &getValue() 
  {
    if (data_->isType<T>())
      return static_cast<ParamData<T> &>(*data_).getValue();
    throw std::runtime_error("Wrong type");
  }

  template<class T>
  const T &getValue() const 
  {
    if (data_->isType<T>())
      return static_cast<ParamData<T> &>(*data_).getValue();
    throw std::runtime_error("Wrong type");
  }

  template<class T>
  T getImmutableValue() const;

  // Methods to get tag.
  const std::string &tag() const 
  {
    return tag_;
  }

  std::string uTag() const;

  // Methods to get the "val" in the desired format.
  std::string stringValue() const;

  // Special accessor methods.
  std::string usVal() const;
  std::string lsVal() const;

  int getType() const 
  {
    return data_->enumType();
  }

  // Method for checking if the parameter is expression valued.
  bool hasExpressionValue() const;

  // Method to check if the parameter value is enclosed in double quotes.
  bool isQuoted();
  
  // Method to check if the parameter value starts with tablefile" and ends with ".
  bool isTableFileTypeQuoted();
  
  // Method to check if the parameter value starts with string" and ends with ".
  bool isStringTypeQuoted();

  // Method to check whether a parameter has a legal real or integer numeric
  // value.
  bool isNumeric() const;
  bool isInteger() const;
  bool isBool() const;

  // Methods for working with time dependency of parameters.
  void setTimeDependent( bool timeDependent );
  bool isTimeDependent() const;

  bool operator==(const Param &right) const
  {
    return (compare_nocase(tag_.c_str(), right.tag_.c_str()) == 0);
  }

  bool operator!=(const Param &right) const
  {
    return !(compare_nocase(tag_.c_str(), right.tag_.c_str()) == 0);
  }

private:
  std::string         tag_;
  ParamData<void> *   data_;
};

struct LessParam : public std::binary_function<Param, std::string, bool>
{
  bool operator()(const Param &param, const std::string &tag) const {
    return compare_nocase(param.tag().c_str(), tag.c_str()) < 0;
  }
};

struct EqualParam : public std::binary_function<Param, std::string, bool>
{
  EqualParam(const std::string &tag = "")
    : tag_(tag)
  {}
  
  bool operator()(const Param &param, const std::string &tag) const {
    return compare_nocase(param.tag().c_str(), tag.c_str()) == 0;
  }

  bool operator()(const Param &param) const {
    return compare_nocase(param.tag().c_str(), tag_.c_str()) == 0;
  }

  const std::string tag_;
};

inline void setParamValue(Param &param, const Param &from_param) 
{
  param.setVal(from_param);
}

inline void setParam(Param &param, const std::string &tag, const Param &from_param) 
{
  param.set(tag, from_param);
}

template<class T>
inline void setParamValue(Param &param, const T &t) 
{
  param.setVal(t);
}

template<class T>
inline void setParam(Param &param, const std::string &tag, const T &t) 
{
  param.set(tag, t);
}

//-----------------------------------------------------------------------------
// Function      : hasExpressionTag
// Purpose       : Determine if the Param tag is an expression.
// Special Notes : This checks the tag.
// Scope         :
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 08/15/04
//-----------------------------------------------------------------------------
inline bool hasExpressionTag(const std::string &tag)
{
  return !tag.empty() && *tag.begin() == '{' && *(tag.end() - 1) == '}';
}

inline bool hasExpressionTag(const Param &param)
{
  return hasExpressionTag(param.tag());
}

template <class It>
typename It::value_type *findParameter(It begin, It end, const std::string &name)
{
  It paramIter = std::find_if(begin, end, EqualParam(name));
  if ( paramIter != end)
    return (typename It::value_type *) &(*paramIter);
  else
    return NULL;
}

// template <class T, class F>
// bool setValue(const Param &param, const char *name, T &value) {
//   if (compare_nocase(param.tag().c_str(), name) == 0) {
//     value = param.getImmutableValue<T>();
//     return true;
//   }
//   return false;
// }

// deepCompare -- compare TAG and Value if TAGS are the same
bool deepCompare(const Param &s0, const Param &s1);

typedef unordered_set<Param> UParamList;
typedef std::list<Param> ParamList;
typedef unordered_map<std::string, Param, HashNoCase, EqualNoCase> ParamMap;
typedef unordered_map<std::string, ParamMap, HashNoCase, EqualNoCase> OptionsMetadataMap;

bool isVectorParam(const Xyce::Util::Param &param, std::string &name, int &index);

template <class T>
bool setValue(const Param &param, const std::string &tag, T &t) {
  if (equal_nocase(tag, param.tag())) {
    t = param.getImmutableValue<T>();
    return true;
  }

  return false;
}

template <class T>
bool setValue(const Param &param, const std::string &tag, T &t, bool &given) {
  if (equal_nocase(tag, param.tag())) {
    t = param.getImmutableValue<T>();
    given = true;
    return true;
  }

  return false;
}

template <class T>
bool setValue(const Param &param, const std::string &tag, void (&f)(T)) {
  if (equal_nocase(tag, param.tag())) {
    f(param.getImmutableValue<T>());
    return true;
  }

  return false;
}


template<class T>
bool setValue(const Xyce::Util::Param &param, const char *tag, std::vector<T> &value) {
  std::string name;
  int index;

  if (isVectorParam(param, name, index)) {
    if (Xyce::equal_nocase(tag, name.c_str())) {
      value.push_back(param.getImmutableValue<T>());
      return true;
    }
  }

  return false;
}

template<class T>
bool setValue(const Xyce::Util::Param &param, const char *tag, std::vector<T> &value, bool &given) {
  std::string name;
  int index;

  if (isVectorParam(param, name, index)) {
    if (Xyce::equal_nocase(tag, name.c_str())) {
      value.push_back(param.getImmutableValue<T>());
      given = true;
      
      return true;
    }
  }

  return false;
}

} // namespace Util
} // namespace Xyce


namespace std
{
template<>
struct hash<Xyce::Util::Param>
{
  std::size_t operator()(const Xyce::Util::Param& p) const
  {
    std::size_t res = std::hash<std::string>{}(p.uTag());
    return res;
  }
};
}

#endif
