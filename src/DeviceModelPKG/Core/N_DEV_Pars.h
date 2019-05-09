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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Netlist device and model parameter management
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Pars_h
#define Xyce_N_DEV_Pars_h

#include <Sacado.hpp>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_Units.h>

///
/// ParameterType::ExprAccess is enumeration of parameter usage types and masks
///
namespace ParameterType {
enum ExprAccess
  {
    NO_DEP    = 0x0,              ///< Parameter can only be set to a constant from netlist
    TIME_DEP  = 0x1,              ///< Parameter may be specified as time dependent expression from netlist
    SOLN_DEP  = 0x2,              ///< Parameter may be specified as a solution dependent expression from netlist
    LOG_T_DEP = 0x8,              ///< Parameter uses temperature interpolation based on log of value
    MIN_RES   = 0x10,             ///< Parameter is subject to being set to minimum lead resistance
    MIN_CAP   = 0x20,             ///< Parameter is subject to being set to minimum junction capacitance
    NO_DOC    = 0x40              ///< Parameter is not to be documented
  };
};

using namespace ParameterType;


namespace Xyce {
namespace Device {

///
/// @addtogroup xyce_device_parameters_detail
///
/// @brief
///
/// Device models, device instances and composite parameters are
/// collectively called an <b>entity</b>, which are all derived from
/// the ParameterBase tag class, and manage the mapping from a
/// parameter's string name to a member variable in each created
/// entity object.  The mapping is refered to as parameter binding and
/// the object's bound member variable is refered to as the parameter
/// member variable.  The object may have an additional member
/// variable for each binding that indicates if the value was
/// specified in the netlist device and is known as the given member
/// variable.
///
/// Each class of entity has a ParametricData<T> class, where T is the
/// entity class, and a singleton object of that class.  This
/// singleton manages the parameter binding.  A parameter descriptor,
/// Descriptor, has the parameter member variable pointer for entity
/// object of type T that sets and retrieves the current value of an
/// entity's member variable using the parameters string name.
///
/// The Descriptor maintains other information about the parameter.
/// It has a serial number for each parameter and a boolean member
/// variable pointer that determines if the parameter has been given
/// in the netlist.  The given value map, GivenValueMap, is a mapping
/// from entity object pointer and serial number to a boolean given
/// member varaible indicating if the parameter was provided by the
/// netlist.
///
/// Any entity may also wish to restore a parameter member variable
/// value to the value originally set during initialization.  The
/// Descriptor maintains an index into an original value map to store
/// and retrieve the parameters original value.  The original value
/// map, OriginalValueMap, is a mapping from entity object pointer and
/// original value index to the parameter's origin value.  If the
/// original value index is -1, then no original value information is
/// maintained.
///
/// A parameter may also be vectorized, allowing it to maintain a
/// vector values.  The Descriptor maintains that vector length.
///
/// A parameter has specific usages as described by the
/// ParameterType::ExprAccess enumeration.
///
/// A parameter may also serves as an aggregation of parametric values
/// which is stored as a composite parameter, CompositeParam, entity.
///
/// And a parameter has units, a documentation category and a
/// documentation description.
///
/// Entity - DeviceModel, DeviceInstance or CompositeParam
///
/// Given - flag and possibly and entity member variable that
/// indicates if a parameter was given in the netlist.
///
/// OriginalValue - parameters value after initialization
///

typedef std::map<int, double> OriginalValueMap;
typedef std::set<int> GivenValueSet;

///
/// Base sensitivity functor
///
class baseSensitivity
{
public:
  baseSensitivity() {}
  virtual ~baseSensitivity() {}

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &param,
    std::vector<double> & dfdp,
    std::vector<double> & dqdp,
    std::vector<double> & dbdp,
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const =0;
};

///
/// Base AC sensitivity functor 
///
class baseACSensitivity
{
public:
  baseACSensitivity() {}
  virtual ~baseACSensitivity() {}

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &param,
    std::vector< std::complex<double> > & dbdp,
    std::vector<int> & Bindices
    ) const =0;
};

///
/// Base matrix sensitivity functor
///
class baseMatrixSensitivity
{
public:
  baseMatrixSensitivity() {}
  virtual ~baseMatrixSensitivity() {}

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &param,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLids,
    std::vector< std::vector<int> > & Q_jacLids
    ) const =0;
};

///
/// Base class for all parameters
///
class ParameterBase
{
public:
  ParameterBase()
    : originalValueMap_(),
      givenValueSet_()
  {}

  virtual ~ParameterBase()
  {}

private:
  ParameterBase(const ParameterBase &);
  ParameterBase &operator=(const ParameterBase &);

public:
  double getOriginalValue(int serial_number)
  {
    return originalValueMap_[serial_number];
  }

  void setOriginalValue(int serial_number, double value)
  {
    originalValueMap_[serial_number] = value;
  }

  bool wasValueGiven(int serial_number) const
  {
    return givenValueSet_.find(serial_number) != givenValueSet_.end();
  }

  void setValueGiven(int serial_number, bool value)
  {
    if (value)
      givenValueSet_.insert(serial_number);
    else
      givenValueSet_.erase(serial_number);
  }

private:
  OriginalValueMap            originalValueMap_;      ///< Map from device entity and original value index to original value
  GivenValueSet               givenValueSet_;         ///< Map from device entity and serial number to value given flag
};

template<class C>
class ParametricData;

template<class T>
class Entry;

template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<T> &entry);

template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::vector<T> > &entry);

template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::string> &entry);

template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<bool> &entry);

inline std::ostream &printEntry(std::ostream &os, const Entry<CompositeMap> &entry);

void typeMismatch(const std::type_info &from_type, const std::type_info &to_type);
void nonexistentParameter(const std::string &name, const std::type_info &entity_type);

///
/// Class Entry<void> defines the parameter binding value entry interface.
///
/// This defines the interface to check the data type of an Entry<T>
/// object and to print the value.  Type specific Entry classes
/// inherit from this class.  The entry_cast<T>() function is used to
/// cast an object of this type to the derived Entry<T> class safely.
///
template<>
class Entry<void>
{
protected:
  ///
  /// Constructs the Entry base class
  ///
  /// @note The construct is protected so it may only be constructed by Entry<T> classes.
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 08:32:36 2013
  ///
  Entry()
  {}

public:
  ///
  /// Destroys Entry
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 08:34:16 2013
  virtual ~Entry()
  {}

private:
  Entry(const Entry &);
  Entry &operator=(const Entry &);

public:
  ///
  /// Returns the type_info of the data type being stored in the entry.
  ///
  /// @return const reference to the type_info of the data type being stored in the entry.
  ///
  virtual const std::type_info &type() const = 0;

  ///
  /// Prints the value of the entry to the output stream
  ///
  /// @param os output stream to write to
  ///
  /// @return reference to the output stream
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:28:02 2013
  std::ostream &print(std::ostream &os) const
  {
    doPrint(os);

    return os;
  }

private:
  ///
  /// Prints the value of the entry to the output stream
  ///
  /// @param os output stream to write to
  ///
  /// @return reference to the output stream
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:29:29 2013
  virtual std::ostream &doPrint(std::ostream &os) const = 0;
};

///
/// Class Entry<T> defines the parameter member variable access for
/// parameter member variable of type T
///
/// The pointer to the parameter member variable and the default value
/// to set are contained in this class.
///
template<class T>
class Entry : public Entry<void>
{
public:
  ///
  /// Constructs an Entry.
  ///
  /// Initializes the parameter member variable pointer of type T.
  /// Note the member variable pointer are actually offsets into an
  /// object where the data for the member exists.  So, by storing
  /// this pointer the value of any object of type can be retrieved.
  ///
  /// The default value ot type T is also contained here.
  ///
  /// @param member_value_ptr parameter member variable pointer
  /// @param default_value Default value of parameter
  ///
  Entry(T ParameterBase::*member_value_ptr, const T &default_value = T())
    : defaultValue_(default_value),
      memberValuePtr_(member_value_ptr)
  {}

  ///
  /// Destroys the entry
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:30:29 2013
  ///
  virtual ~Entry()
  {}

private:
  Entry(const Entry &);               ///< No copying
  Entry &operator=(const Entry &);    ///< No assignment

public:
  ///
  /// Returns type_info of this entry.
  ///
  /// @return const reference to type_info of this entry
  ///
  virtual const std::type_info &type() const
  {
    return typeid(T);
  }

private:
  ///
  /// Prints the value of the entry to the output stream
  ///
  /// @param os output stream to write to
  ///
  /// @return reference to the output stream
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:29:29 2013
  ///
  virtual std::ostream &doPrint(std::ostream &os) const
  {
    return printEntry(os, *this);
  }

public:
  ///
  /// Return the member pointer to the data member variable that holds
  /// the value associated with this parameter
  ///
  /// @return member pointer to the data member variable
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:31:21 2013
  ///
  T ParameterBase::*getMemberPtr() const
  {
    return memberValuePtr_;
  }

  ///
  /// Return the value of the entity's parameter
  ///
  /// @param entity device class or device instance
  ///
  /// @return const reference to the value of the entity's parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:32:58 2013
  ///
  const T &getValue(const ParameterBase &entity) const
  {
    return entity.*memberValuePtr_;
  }

  ///
  /// Return the value of the entity's parameter
  ///
  /// @param entity device class or device instance
  ///
  /// @return reference to the value of the entity's parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:32:58 2013
  ///
  T &getValue(ParameterBase &entity) const
  {
    return entity.*memberValuePtr_;
  }

  ///
  /// Sets the value of the entity's parameter
  ///
  /// @param entity device class or device instance
  /// @param value value to set the entity's parameter to
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:32:58 2013
  ///
  void setValue(ParameterBase &entity, const T &value) const
  {
    entity.*memberValuePtr_ = value;
  }

  ///
  /// Return the default value of the parameter
  ///
  /// All parameters provide a default value when created.  The
  /// parameter is set to this value and the given flag is cleared on
  /// construction.  If the value is provided by the net list, the
  /// parameter's value is set accordingly and given flag is set to
  /// true.
  ///
  /// @return const reference to the parameter's default value
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:35:18 2013
  ///
  const T &getDefaultValue() const
  {
    return defaultValue_;
  }

  ///
  /// Sets the parameter's default value
  ///
  /// All parameters provide a default value when created.  The
  /// parameter is set to this value and the given flag is cleared on
  /// construction.  If the value is provided by the net list, the
  /// parameter's value is set accordingly and given flag is set to
  /// true.
  ///
  /// @param value default value of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:37:50 2013
  ///
  void setDefaultValue(const T &value)
  {
    defaultValue_ = value;
  }

private:
  T                           defaultValue_;        ///< Default value of parameter
  T ParameterBase::*          memberValuePtr_;      ///< Member pointer containing value
};

///
/// Casts the entry to type T
///
/// If the entry if not of the specified type, typeMismatch() is
/// called to emit a fatal error.
///
/// @tparam T type to cast to
/// @param entry entry to cast
///
/// @return const reference to Entry<T>
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Thu Jul 25 13:21:09 2013
///
template<class T>
const Entry<T> &entry_cast(const Entry<void> &entry)
{
  if (entry.type() != typeid(T))
    typeMismatch(entry.type(), typeid(T));

  return static_cast<const Entry<T> &>(entry);
}

///
/// Casts the entry to type T
///
/// If the entry if not of the specified type, typeMismatch() is called to emit a fatal error.
///
/// @tparam T type to cast to
/// @param entry entry to cast
///
/// @return reference to Entry<T>
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Thu Jul 25 13:21:09 2013
///
template<class T>
Entry<T> &entry_cast(Entry<void> &entry)
{
  if (entry.type() != typeid(T))
    typeMismatch(entry.type(), typeid(T));

  return static_cast<Entry<T> &>(entry);
}

///
/// Class Descriptor describes the parameters stored in the
/// ParametricData parameter map.
///
/// The descriptor contains the Entry for the parameter member
/// variable.  Each parameter is assigned a serialNumber_ when it is
/// created and is used by manage the given value map.  The parameter
/// may have a given member variable as well.  The parameter can be
/// marked to store its original initialized value and the
/// serialNumber_ is used to store this value in the originalValueMap.
///
/// If the parameter is vectorized, the index to this parameter is in
/// vec_.
///
/// The expressionAccess_ is used to indicate the parameter's usage.
///
/// The Descriptor also contains the units, catagory and description
/// for documentation generation.  serial number, its original value
/// index if the initiazes value needs to be restored, its usage
/// (ExprAccess), units,
///
class Descriptor
{
public:
  ///
  /// Constructs Descriptor
  ///
  /// @param entry the parameter member variable pointer, type and default value
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:38:43 2013
  ///
  Descriptor(Entry<void> *const entry)
    : serialNumber_(0),
      originalValueFlag_(false),
      vec_(0),
      expressionAccess_(ParameterType::NO_DEP),
      entry_(entry),
      unit_(U_NONE),
      category_(CAT_NONE),
      description_(""),
      compositeParametricData_(0),
      given_(static_cast<bool ParameterBase::*>(0)),
      sensitivityAvailable_(false),
      acSensitivityAvailable_(false),
      matrixSensitivityAvailable_(false),
      sensPtr_(0),
      acSensPtr_(0),
      matSensPtr_(0),
      autoConvertTemp_(true)
  {}

  ///
  /// Destroy Descriptor
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 15:42:55 2013
  ///
  virtual ~Descriptor()
  {
    delete entry_;
  }

private:
  Descriptor(const Descriptor &descriptor);
  Descriptor &operator=(const Descriptor &descriptor);

public:
  ///
  /// Tests entry data type
  ///
  /// @return true if entry is of type T
  ///
  /// @date   Wed Aug  7 10:31:04 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  template <class T>
  bool isType() const
  {
    if (entry_)
      return entry_->type() == typeid(T);

    return typeid(T) == typeid(int);
  }

  bool hasCompositeData() const
  {
    return compositeParametricData_ != 0;
  }

  ///
  /// Sets a boolean marking an original value having been stored.
  ///
  /// @param original_value_flag
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:37:03 2013
  ///
  Descriptor &setOriginalValueStored(bool original_value_flag)
  {
    originalValueFlag_ = original_value_flag;
    return *this;
  }

  ///
  /// Returns whether an original value has been stored
  ///
  /// @return True if original value has been stored.  False otherwise.
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:19:59 2013
  ///
  bool hasOriginalValueStored() const
  {
    return originalValueFlag_;
  }

  ///
  /// Sets the expression access which describe the usage of the parameter
  ///
  /// @param expression_access usage of the paraemeter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:41:56 2013
  ///
  Descriptor &setExpressionAccess(ExprAccess expression_access)
  {
    expressionAccess_ = expression_access;
    return *this;
  }

  ///
  /// Gets the expression access which describes the usage of the paramter
  ///
  /// @return usage of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:41:26 2013
  ///
  ExprAccess getExpressionAccess() const
  {
    return expressionAccess_;
  }

  ///
  /// Sets the units of the parameter, only used to document the parameter
  ///
  /// @param unit units of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:43:32 2013
  ///
  Descriptor &setUnit(ParameterUnit unit)
  {
    unit_ = unit;
    return *this;
  }

  ///
  /// Gets the units of the parameter, only used to document the parameter
  ///
  /// @return units of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:44:16 2013
  ///
  ParameterUnit getUnit() const
  {
    return unit_;
  }

  ///
  /// Sets the category of the parameter, only used to document the parameter
  ///
  /// @param category category of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:45:09 2013
  Descriptor &setCategory(ParameterCategory category)
  {
    category_ = category;
    return *this;
  }

  ///
  /// Gets the category of the parameter, only used to document the parameter
  ///
  /// @return categort of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:45:40 2013
  ///
  ParameterCategory getCategory() const
  {
    return category_;
  }

  ///
  /// Sets the description of the parameter, only used to document the parameter
  ///
  /// @param description description of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:46:21 2013
  ///
  Descriptor &setDescription(const std::string &description)
  {
    description_ = description;
    return *this;
  }

  ///
  /// Gets the description of the parameter, only used to document the parameter
  ///
  /// @return description of the parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:46:46 2013
  ///
  const std::string &getDescription() const
  {
    return description_;
  }

  ///
  /// Sets the composite parametric
  ///
  /// A composite parameter is a named aggregation of parameters
  ///
  /// @param composite_parametric_data composite parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:47:19 2013
  ///
  Descriptor &setCompositeParametricData(const ParametricData<void> *composite_parametric_data)
  {
    compositeParametricData_ = composite_parametric_data;
    return *this;
  }

  ///
  /// Return the composite parameter
  ///
  /// A composite parameter is a named aggregation of parameters
  ///
  /// @return const pointer to the composite parameter
  ///
  /// @note I think the template parameter is really only ever CompositeParam, though sometimes its
  /// currently void.
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:09:01 2013
  ///
  template <class U>
  const ParametricData<U> *getCompositeParametricData() const
  {
    return static_cast<const ParametricData<U> *>(compositeParametricData_);
  }

  ///
  /// sets the vector length of an vectorized parameter
  ///
  /// @param index
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:50:05 2013
  ///
  Descriptor &setVec(int index)
  {
    vec_ = index;
    return *this;
  }

  ///
  /// Gets the vector length of a vectorized parameter
  ///
  /// @return length of the vectorized parameter
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:51:10 2013
  ///
  int getVec() const
  {
    return vec_;
  }

  ///
  /// Gets the entry object of the parameter
  ///
  /// @return const reference to the Entry<void>
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 15:29:16 2013
  ///
  const Entry<void> &getEntry() const
  {
    return *entry_;
  }

  ///
  /// Gets the entry object of the parameter
  ///
  /// @return reference to the Entry<void>
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 15:29:16 2013
  ///
  Entry<void> &getEntry()
  {
    return *entry_;
  }

  ///
  /// Sets the serial number used to store and retrieve given boolean from the GivenValueMap
  ///
  /// @param serial_number serial number of the parameter in the GivenValueMap
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:38:18 2013
  ///
  Descriptor &setSerialNumber(int serial_number)
  {
    serialNumber_ = serial_number;
    return *this;
  }

  ///
  /// Gets the serial number used to store and retireve given boolean fromt he GivenValueMap
  ///
  /// @return parameter serial number
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 14:40:30 2013
  ///
  int getSerialNumber() const
  {
    return serialNumber_;
  }

public:
  ///
  /// Returns the value of the parameter for the entity
  ///
  /// @param entity device class or device instance
  ///
  /// @return reference to the value of the parameter for the entity
  ///
  /// @date   Wed Aug  7 10:36:10 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  template <class T>
  const T &value(const ParameterBase &entity) const
  {
    const Entry<T> &entry = entry_cast<T>(*entry_);

    return entry.getValue(entity);
  }

  ///
  /// Returns the value of the parameter for the entity
  ///
  /// @param entity device class or device instance
  ///
  /// @return reference to the value of the parameter for the entity
  ///
  /// @date   Wed Aug  7 10:36:10 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  template <class T>
  T &value(ParameterBase &entity) const
  {
    const Entry<T> &entry = entry_cast<T>(*entry_);

    return entry.getValue(entity);
  }

  ///
  /// Returns the parameter member variable pointer of the enrtry
  ///
  /// @return parameter member pointer of type T
  ///
  /// @date   Wed Aug  7 10:54:58 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  template <class T>
  T ParameterBase::*getMemberPtr() const
  {
    const Entry<T> &entry = entry_cast<T>(*entry_);

    return entry.getMemberPtr();
  }

  ///
  /// Tests if parameter has a given data member
  ///
  /// Parameters may provide a boolean member variable that is set true if the netlist provides
  /// the value.
  ///
  /// @return true if the parameter has a boolean member variable to set if it is provided by the netlist
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:03:55 2013
  ///
  bool hasGivenMember() const
  {
    return given_ != static_cast<bool ParameterBase::*>(0);
  }

  ///
  /// Sets the boolean member variable to set if the netlist provides the value.
  ///
  /// Parameters may provide a boolean member variable that is set true if the netlist provides
  /// the value.
  ///
  /// @param given boolean member variable to be set
  ///
  /// @date   Wed Aug  7 11:02:38 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  template<class U>
  Descriptor &setGivenMember(bool U::*given)
  {
    given_ = static_cast<bool ParameterBase::*>(given);
    return *this;
  }

  ///
  /// Tests if the parameter has been given by the netlist.
  ///
  /// @return true if the parameter was given by the netline
  ///
  /// @date   Wed Aug  7 11:01:42 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  bool getGiven(ParameterBase &entity) const
  {
    if (hasGivenMember())
    {
      return entity.*given_;
    }

    return false;
  }

  ///
  /// Sets the given state of the parameter to value
  ///
  /// The parameter's boolean member variable that is set true if exists and the netlist provides
  /// the value.
  ///
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:09:00 2013
  ///
  void setGiven(ParameterBase &entity, bool value) const
  {
    if (hasGivenMember())
    {
      entity.*given_ = value;
    }
  }

  ///
  /// Sets a boolean to indicate if analytic sensitivities are available w.r.t. this parameter.
  //
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author Eric R. Keiter, Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:09:00 2013
  ///
  Descriptor &setAnalyticSensitivityAvailable(bool value)
  {
    sensitivityAvailable_ = value;
    return *this;
  }

  Descriptor &setAnalyticACSensitivityAvailable(bool value)
  {
    acSensitivityAvailable_ = value;
    return *this;
  }

  Descriptor &setAnalyticMatrixSensitivityAvailable(bool value)
  {
    matrixSensitivityAvailable_ = value;
    return *this;
  }

  Descriptor &setSensitivityFunctor(const baseSensitivity * sensPtr)
  {
    sensPtr_ = sensPtr;
    return *this;
  }

  Descriptor &setACSensitivityFunctor(const baseACSensitivity * sensPtr)
  {
    acSensPtr_ = sensPtr;
    return *this;
  }
  
  Descriptor &setMatrixSensitivityFunctor(const baseMatrixSensitivity * matSensPtr)
  {
    matSensPtr_ = matSensPtr;
    return *this;
  }

  ///
  /// returns a boolean to indicate if analytic sensitivities are available w.r.t. this parameter.
  //
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author Eric R. Keiter, Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:09:00 2013
  ///
  bool getAnalyticSensitivityAvailable() const
  {
    return sensitivityAvailable_;
  }

  bool getAnalyticACSensitivityAvailable() const
  {
    return acSensitivityAvailable_;
  }

  bool getAnalyticMatrixSensitivityAvailable() const
  {
    return matrixSensitivityAvailable_;
  }

  ///
  /// returns analytic sensitivity, evaluated at the most recent solve
  /// point.
  //
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author Eric R. Keiter, Sandia National Laboratories 1355
  /// @date   Wed Aug  7 11:09:00 2013
  ///
  bool getAnalyticSensitivity(
      const ParameterBase &entity1,
      const std::string & name,
      std::vector<double> & dfdp, std::vector<double> & dqdp, std::vector<double> & dbdp,
      std::vector<int> & Findices, std::vector<int> & Qindices, std::vector<int> & Bindices) const
  {
    if (sensPtr_)
    {
      sensPtr_->operator()(entity1, name, dfdp, dqdp, dbdp, Findices, Qindices, Bindices);
      return true;
    }
    else
    {
      return false;
    }
  }

  bool getAnalyticMatrixSensitivity(
      const ParameterBase &entity1,
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLids,
    std::vector< std::vector<int> > & Q_jacLids
      ) const
  {
    if (matSensPtr_)
    {
      matSensPtr_->operator()(entity1, name, d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLids, Q_jacLids);
      return true;
    }
    else
    {
      return false;
    }
  }

  bool getAnalyticBSensVectorsforAC(
      const ParameterBase &entity1,
      const std::string & name,
      std::vector< std::complex<double> > & dbdp,
      std::vector<int> &        Bindices) const
  {
    if (acSensPtr_)
    {
      acSensPtr_->operator()(entity1, name, dbdp, Bindices);
      return true;
    }
    else
    {
      return false;
    }
  }

  ///
  /// Sets a boolean to indicate if this parameter may be autoconverted
  /// from Celsius to Kelvin (or vice-versa) if it's a temperature param.
  //
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author Tom Russo, Sandia National Laboratories 1355
  /// @date   Wed Mar  30
  ///
  Descriptor &setAutoConvertTemperature(bool value)
  {
    autoConvertTemp_ = value;
    return *this;
  }

  ///
  /// returns a boolean to indicate if analytic sensitivities are available w.r.t. this parameter.
  //
  /// @param entity device class or device instance
  /// @param value true if provided by netlist
  ///
  /// @author Tom Russo, Sandia National Laboratories 1355
  /// @date   Wed Mar 30
  ///
  bool getAutoConvertTemperature() const
  {
    return autoConvertTemp_;
  }

  
private:
  int                                 serialNumber_;          ///< Unique identifier of descriptor
  bool                                originalValueFlag_;     ///< Flag indicating original value was stored
  int                                 vec_;                   ///< If > 0 specifies a vector of params.(eg: if = 3 then IC becomes IC1, IC2, IC3)
  ExprAccess                          expressionAccess_;      ///< Flags for parameter attributes, such as whether can be input by user, may depend on time, etc.
  Entry<void> * const                 entry_;                 ///< Pointer to entry which contains the value
  ParameterUnit                       unit_;                  ///< Unit designator for documentation
  ParameterCategory                   category_;              ///< Category designator for documentation
  std::string                         description_;           ///< Description of parameter for documentation
  const ParametricData<void> *        compositeParametricData_; ///< If a composite, then this points to the ParameterData of the composite
  bool ParameterBase::*               given_;                 ///< Pointer to given bool, usually 0
  bool                                sensitivityAvailable_;  ///< Flag indicating that analytic sensitivity is available for this param.
  bool                                acSensitivityAvailable_;  ///< Flag indicating that analytic AC sensitivity is available for this param.
  bool                                matrixSensitivityAvailable_;  ///< Flag indicating that analytic matrix sensitivity is available for this param.
  const baseSensitivity *             sensPtr_;               ///< pointer to the sensitivity functor
  const baseACSensitivity *           acSensPtr_;             ///< pointer to the AC sensitivity functor
  const baseMatrixSensitivity *       matSensPtr_;            ///< pointer to the matrix sensitivity functor
  bool                                autoConvertTemp_;  ///< Flag indicating that temperature should be converted to/from Kelvin
};

///
/// Prints the entry default value to the output stream
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:32:26 2013
///
template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<T> &entry)
{
  os << entry.getDefaultValue();

  return os;
}

///
/// Prints the entry default values of a vectorized parameter
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:33:24 2013
///
template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::vector<T> > &entry)
{
  for (typename std::vector<T>::const_iterator it = entry.getDefaultValue().begin(); it != entry.getDefaultValue().end(); ++it)
    os << (*it) << std::endl;

  return os;
}

///
/// Prints the entry default values of a map parameter
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Thu Feb  6 16:37:16 2014
///
template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::map<std::string, T> > &entry)
{
  for (typename std::map<std::string, T>::const_iterator it = entry.getDefaultValue().begin(); it != entry.getDefaultValue().end(); ++it)
    os << (*it).first << std::endl;

  return os;
}

///
/// Prints the entry default string value, within single quotes
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:34:18 2013
///
template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::string> &entry)
{
  os << "'" << entry.getDefaultValue() << "'";

  return os;
}

///
/// Prints the entry default boolean value, printed as true or false
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:34:18 2013
///
template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<bool> &entry)
{
  os << (entry.getDefaultValue() ? "true" : "false");

  return os;
}

///
/// Prints the entry composite value as newline terminated list of colon separated name, value pairs
///
/// @param os output stream
/// @param entry entry to print
///
/// @return reference to the output stream
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:34:18 2013
///
inline std::ostream &printEntry(std::ostream &os, const Entry<CompositeMap> &entry)
{
  for (CompositeMap::const_iterator it = entry.getDefaultValue().begin(); it != entry.getDefaultValue().end(); ++it)
    os << (*it).first << ": " << (*it).second << std::endl;

  return os;
}

///
/// Gets the default value of the parameter
///
/// @param descriptor descriptor of the parameter
///
/// @return const reference to the default value of the entry of the descriptor
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:37:33 2013
///
template <class T>
inline const T &getDefaultValue(const Descriptor &descriptor)
{
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getDefaultValue();
}

///
/// Sets the default value of the parameter
///
/// @param descriptor descriptor of the parameter
/// @param t default value
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:38:49 2013
///
template <class T>
inline void setDefaultValue(Descriptor &descriptor, const T &t)
{
  Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  entry.setDefaultValue(t);
}

///
/// Returns the value of the parameter for the entity
///
/// @param entity device class or device instance
/// @param descriptor descriptor of the parameter
///
/// @return const reference to the value of the entry of the descriptor
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:39:40 2013
///
template <class T>
inline const T &value(const ParameterBase &entity, const Descriptor &descriptor)
{
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}

///
/// Returns the value of the parameter for the entity
///
/// @param entity device class or device instance
/// @param descriptor descriptor of the parameter
///
/// @return reference to the value of the entry of the descriptor
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:39:40 2013
///
template <class T>
inline T &value(ParameterBase &entity, const Descriptor &descriptor)
{
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}


///
/// Gets the value of the parameter for the entity
///
/// @param entity device class or device instance
/// @param descriptor descriptor of the parameter
///
/// @return const reference to the value of the entry of the descriptor
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:39:40 2013
///
template <class T, class U>
inline const T &getValue(const ParameterBase &entity, const Descriptor &descriptor)
{
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}

///
/// Sets the value of the parameter
///
/// @param entity device class or device instance
/// @param descriptor descriptor of the parameter
/// @param value value
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Fri Aug  9 15:42:13 2013
///
template <class T, class U>
inline void setValue(ParameterBase &entity, const Descriptor &descriptor, const T &value)
{
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  entry.setValue(entity, value);
}

///
/// Class ParametricData<void> manages the configuration information and the parameter binding map
///
/// Parametric data associated with a device instance, device model or composite parameter
///
/// The Parametric data class manages the mapping of parameter string names to descriptors and the
/// general configuration information associated with a device model.
///
/// To restore original values during perturbation, the originalValueCount_ and serialNumber_ members
/// maintain counts of original values to be stored and of parameters declared.
///
/// @date   Tue Aug  6 13:10:21 2013
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
///
template<>
class ParametricData<void>
{
public:
  ///
  /// Constructs a ParametricData object
  ///
  /// @date   Tue Aug  6 13:10:21 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  ParametricData()
    : map_()
  {}

  ///
  /// Destroys a ParametricData object
  ///
  /// @date   Tue Aug  6 13:13:47 2013
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  ///
  virtual ~ParametricData()
  {
    for (ParameterMap::iterator it = map_.begin(); it != map_.end(); ++it)
      delete (*it).second;
  }

private:
  ParametricData(const ParametricData &parametric_data);              ///< No copying
  ParametricData &operator=(const ParametricData &parametric_data);   ///< No assignment

public:
  ///
  /// Gets the parameter binding map map
  ///
  /// @return reference to the parameter binding map
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 15:50:06 2013
  ///
  ParameterMap &getMap()
  {
    return map_;
  }

  ///
  /// Returns the parameter binding map
  ///
  /// @return const reference to the parameter binding map
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 15:50:06 2013
  ///
  const ParameterMap &getMap() const
  {
    return map_;
  }

protected:
  ///
  /// Adds the parameter to the parameter binding map
  ///
  /// @param name parameter name
  /// @param descriptor descriptor created for the parameter
  /// @param parameter_data_class typeinfo to get the class name for diagnostics
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Fri Aug  9 16:23:45 2013
  ///
  void addDescriptor(const std::string &name, Descriptor *descriptor, const std::type_info &parameter_data_class);

protected:
  ParameterMap  map_;                 ///< Mapping from parameter name to descriptor
};

void checkExprAccess(const std::string &name, ParameterType::ExprAccess &expr_access, const std::type_info &parameter_data_class);

///
/// Manages parameter binding for class C.
///
///
template<class C>
class ParametricData : public ParametricData<void>
{
public:
  ///
  /// Constructs the parameter data map
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:39:06 2014
  ///
  ParametricData()
    : ParametricData<void>()
  {}

  ///
  /// Destroys the parameter data map
  ///
  ///
  ///
  /// @return
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:39:34 2014
  ///
  virtual ~ParametricData()
  {}

private:
  ParametricData(const ParametricData &parametric_data);              ///< No copying
  ParametricData &operator=(const ParametricData &parametric_data);   ///< No assignment

public:
  ///
  /// Adds the parameter description to the parameter map
  ///
  /// @tparam T                data type of the parameter
  /// @tparam U                class containing the member data storing the parameter's value
  /// @param parName           const pointer to the parameter name
  /// @param default_value     default value
  /// @param varPtr            member pointer to the parameter's value
  ///
  /// @return
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:33:48 2014
  ///
  template<class T, class U>
  Descriptor &addPar(const char *parName, T default_value, T U::*varPtr)
  {
    Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), default_value));

    addDescriptor(parName, descriptor, typeid(C));

    return *descriptor;
  }

  ///
  /// Adds the parameter description to the parameter map
  ///
  /// This is specialization to allow const char * (quoted strings) to create a std::string type parameter.
  ///
  /// @tparam U                class containing the member data storing the parameter's value
  /// @param parName           const pointer to the parameter name
  /// @param default_value     default value
  /// @param varPtr            member pointer to the parameter's value
  ///
  /// @return
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:33:48 2014
  ///
  template<class U>
  Descriptor &addPar(const char *parName, const char *default_value, std::string U::*varPtr)
  {
    Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), default_value));

    addDescriptor(parName, descriptor, typeid(C));

    return *descriptor;
  }

  ///
  /// Adds a composite parameter to the parameter map
  ///
  /// THIS IS SUPER DANGEROUS U * MAY NOT ACTUALLY POINT TO CompositeParam *.  IT COULD BE SLICED
  ///
  /// @param comp_name
  /// @param composite_pars
  /// @param composite_map
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:47:59 2014
  ///
  template<class U, class V>
  void addComposite(const char *comp_name, const ParametricData<U> &composite_pars, std::map<std::string, U *> V::*composite_map)
  {
    Descriptor *descriptor = new Descriptor(new Entry<CompositeMap>(reinterpret_cast<CompositeMap ParameterBase::*>(composite_map)));

    descriptor->setUnit(U_INVALID);
    descriptor->setCategory(CAT_NONE);
    descriptor->setCompositeParametricData(&composite_pars);

    addDescriptor(comp_name, descriptor, typeid(C));
  }

  ///
  /// Adds a composite vector parameter to the parameter map
  ///
  /// THIS IS SUPER DANGEROUS U * MAY NOT ACTUALLY POINT TO CompositeParam *.  IT COULD BE SLICED
  ///
  /// @param comp_name
  /// @param composite_pars
  /// @param composite_vector
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:47:59 2014
  ///
  template<class U, class V>
  void addComposite(const char *comp_name, const ParametricData<U> &composite_pars, std::vector<U *> V::*composite_vector)
  {
    Descriptor *descriptor = new Descriptor(new Entry<CompositeVector>(reinterpret_cast<CompositeVector ParameterBase::*>(composite_vector)));

    descriptor->setUnit(U_INVALID);
    descriptor->setCategory(CAT_NONE);
    descriptor->setCompositeParametricData(&composite_pars);
    addDescriptor(comp_name, descriptor, typeid(C));
  }

  ///
  /// Allows the parameter to be specified as a vector
  ///
  ///
  /// @param cname
  /// @param len
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   Thu Feb  6 16:50:29 2014
  ///
  void makeVector(const std::string &cname, int len)
  {
    for (int i = 1; i <= len; ++i)
    {
      std::ostringstream oss;
      oss << cname << i;
      std::string param = oss.str();

      ParameterMap::iterator it = map_.find(param);
      if (it == map_.end())
        nonexistentParameter(param, typeid(C));

      Descriptor &descriptor = *(*it).second;
      descriptor.setVec(i);
    }
  }
};

///
/// Set the default values for the parameter.
///
/// Note that any values provided in the objects constructor that are defined as parameters have their initialized value
/// replaced with the default value given in the addPar() call
///
/// @param parameter_base
/// @param begin
/// @param end
/// @param device_options
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Thu Feb  6 16:51:00 2014
///
void setDefaultParameters(ParameterBase &parameter_base, ParameterMap::const_iterator begin, ParameterMap::const_iterator end, const DeviceOptions &device_options);

///
/// Retrieve a parameter's original value
///
/// @param parameter_base        device entity holding parameter
/// @param serial_number         entity's serial number
///
/// @return original value of the parameter
///
/// @date   Tue Aug  6 13:51:16 2013
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
///
inline double getOriginalValue(ParameterBase &parameter_base, int serial_number)
{
  return parameter_base.getOriginalValue(serial_number);
}

///
/// Set a parameter's original value
///
/// @param parameter_base        device entity or composite holding parameter
/// @param serial_number         entity's serial number
/// @param value                 value to be stored
///
/// @date   Tue Aug  6 13:53:29 2013
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
///
inline void setOriginalValue(ParameterBase &parameter_base, int serial_number, double value)
{
  parameter_base.setOriginalValue(serial_number, value);
}


///
/// Return true if a value was provided for the device
///
/// @param parameter_base        device entity or composite holding parameter
/// @param serial_number         serial number of parameter
///
/// @return true if a value was provided for the device
///
/// @date   Tue Aug  6 13:54:25 2013
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
///
inline bool wasValueGiven(const ParameterBase &parameter_base, int serial_number)
{
  return parameter_base.wasValueGiven(serial_number);
}


///
/// Set the given value state of a parameter
///
/// @param parameter_base        device entity or composite holding parameter
/// @param serial_number         serial number of parameter
/// @param value                 true if the value was given
///
/// @date   Tue Aug  6 13:55:34 2013
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
///
inline void setValueGiven(ParameterBase &parameter_base, int serial_number, bool value)
{
  parameter_base.setValueGiven(serial_number, value);
}

///
/// Returns true if the name is TNOM or TEMP
///
/// @param name  parameter name
///
/// @return true if the name is TNOM or TEMP
///
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
/// @date   Thu Feb  6 16:52:26 2014
///
inline bool isTempParam(const std::string &name)
{
  return equal_nocase(name, "TNOM") || equal_nocase(name, "TEMP");
}

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Pars_h
