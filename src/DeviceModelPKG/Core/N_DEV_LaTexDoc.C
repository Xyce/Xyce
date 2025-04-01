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
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 3/20/2013
//
//
//
//
//-------------------------------------------------------------------------

// All Xyce source files must include this!
#include <Xyce_config.h>

#include <string>
#include <map>
#include <vector>
#include <fstream>

#include <cctype>

#include <N_DEV_LaTexDoc.h>

#include <N_DEV_fwd.h>
#include <N_DEV_Units.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace Device {

namespace {

typedef std::map<std::string, Descriptor *, LessNoCase> OrderedParameterMap;

std::ostream &laTexComposite(std::ostream &os, const std::string &composite_name, const std::string &composite_description, const ParametricData<void> &parameters);

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::escape
// Purpose       : Add escape backslash before underscore for LaTeX
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/16/06
//-----------------------------------------------------------------------------
struct laTexProtect {
  laTexProtect(const std::string &t)
    : text(t)
  {}

  const std::string &text;
};


std::ostream &operator<<(std::ostream &os, const laTexProtect &laTex_protect)
{
  for (std::string::const_iterator cit = laTex_protect.text.begin(); cit != laTex_protect.text.end(); ++cit) {
    if ((*cit) == '_')
      os << "\\_";
    else
      os << (*cit);
  }

  return os;
}


struct laTexClean {
  laTexClean(const std::string &t)
    : text(t)
  {}

  const std::string &text;
};


std::ostream &operator<<(std::ostream &os, const laTexClean &laTex_clean)
{
  for (std::string::const_iterator cit = laTex_clean.text.begin(); cit != laTex_clean.text.end(); ++cit) {
    if ((*cit) != '_')
      os << (*cit);
  }

  return os;
}


const std::string &categoryName(int category) {
  static std::vector<std::string> s_categoryList;

  if (s_categoryList.empty()) {
    s_categoryList.resize(CAT_MAX);
    s_categoryList[CAT_NONE] = "";
    s_categoryList[CAT_AC] = "AC Parameters";
    s_categoryList[CAT_BASIC] = "Basic Parameters";
    s_categoryList[CAT_BIN] = "Bin Parameters";
    s_categoryList[CAT_CAP] = "Capacitance Parameters";
    s_categoryList[CAT_CONTROL] = "Control Parameters";
    s_categoryList[CAT_CURRENT] = "Current Parameters";
    s_categoryList[CAT_DC] = "DC Parameters";
    s_categoryList[CAT_DEPENDENCY] = "Dependency Parameters";
    s_categoryList[CAT_DOPING] = "Doping Parameters";
    s_categoryList[CAT_FLICKER] = "Flicker and Thermal Noise Parameters";
    s_categoryList[CAT_GEOMETRY] = "Geometry Parameters";
    s_categoryList[CAT_INITIAL] = "Initial Condition Parameters";
    s_categoryList[CAT_MATERIAL] = "Material Parameters";
    s_categoryList[CAT_NQS] = "NQS Parameters";
    s_categoryList[CAT_RADP] = "Radiation Pulse Parameters";
    s_categoryList[CAT_RES] = "Resistance Parameters";
    s_categoryList[CAT_PROCESS] = "Process Parameters";
    s_categoryList[CAT_RF] = "RF Parameters";
    s_categoryList[CAT_RAD] = "Radiation Parameters";
    s_categoryList[CAT_TEMP] = "Temperature Parameters";
    s_categoryList[CAT_TUNNEL] = "Tunnelling Parameters";
    s_categoryList[CAT_VBI] = "Built-in Potential Lowering Parameters";
    s_categoryList[CAT_VOLT] = "Voltage Parameters";
    s_categoryList[CAT_ASYMRDS] = "Asymmetric and Bias-Dependent $R_{ds}$ Parameters";
    s_categoryList[CAT_ASYMDDS] = "Asymmetric Source/Drain Junction Diode Parameters";
    s_categoryList[CAT_IMPACT] = "Impact Ionization Current Parameters";
    s_categoryList[CAT_GDLEAKAGE] = "Gate-induced Drain Leakage Model Parameters";
    s_categoryList[CAT_STRESS] = "Stress Effect Model Parameters";
    s_categoryList[CAT_WELL] = "Well-Proximity Effect Model Parameters";

    s_categoryList[CAT_STATIC] = "Static Model Parameters";
    s_categoryList[CAT_DYNAMIC] = "Dynamic Model Parameters";
    s_categoryList[CAT_CARRIER] = "Carrier Model Parameters";
    s_categoryList[CAT_OUTPUT] = "Model Output Parameters";
    s_categoryList[CAT_PULSE] = "Pulse Parameters";

    s_categoryList[CAT_SCALING] = "Scaling Parameters";
    s_categoryList[CAT_BOUNDARYCONDITIONS] = "Boundary Condition Parameters";
  }

  return s_categoryList.size() ? s_categoryList[category] : s_categoryList[CAT_UNKNOWN];
}

void unitDescription(
  const std::string &           name,
  const OrderedParameterMap &   parameter_map,
  const Descriptor &            descriptor,
  ParameterUnit &               unit,
  ParameterCategory &           category,
  std::string &                 description,
  std::string &                 versionRange)
{
  unit = descriptor.getUnit();
  if (unit == STANDARD)
  {
    for (const StdDescription *it = Units::descriptionTable; it != Units::descriptionTable + Units::descriptionTableSize; ++it) {
      if (name == (*it).Name)
      {
        description = (*it).Description;
        category = (*it).Category;
        unit = (*it).Unit;
        break;
      }
    }

    if (description.empty())
    {
      Report::UserWarning0() << "No description given for " << name << ".  Setting category and unit as UNKNOWN.";

      description = "Unspecified";
      category = CAT_UNKNOWN;
      unit = U_UNKNOWN;
    }
  }
  else
  {
    category = descriptor.getCategory();
    if (category == CAT_DEPENDENCY)
    {
      int num;
      std::string base_name = name.substr(1);
      if (name[0] == 'L')
      {
        description = "Length";
        num = 1;
      }
      else if (name[0] == 'W')
      {
        description = "Width";
        num = 1;
      }
      else if (name[0] == 'P')
      {
        description = "Cross-term";
        num = 2;
      }
      description += " dependence of ";
      description += base_name;

      OrderedParameterMap::const_iterator it = parameter_map.find(base_name);
      if (it == parameter_map.end())
      {
        Report::DevelFatal().in("DeviceEntity::getUnitDescription")
          << "Base parameter not found for " << name << " in " <<  base_name;
      }

      else {
        const Descriptor &base_descriptor = *(*it).second;

        std::string base_description;
        ParameterCategory base_category = CAT_UNKNOWN;
        std::string versionRange;

        unitDescription(base_name, parameter_map, base_descriptor, unit, base_category, base_description,versionRange);
        if (unit != U_INVALID && unit != U_UNKNOWN)
        {
          for (int i = 0; i < num; ++i)
          {
            for (const UnitInfo *it = Units::unitTable; it != Units::unitTable + Units::unitTableSize; ++it)
            {
              if ((*it).Unit == unit)
              {
                unit = (*it).UnitM;
                if (unit == U_INVALID)
                {
                  Report::UserWarning0() << "Need unit for " << (*it).description << " times m";
                }
                break;
              }
            }
          }
        }
      }
    }
    else
    {
      description = descriptor.getDescription();
    }
  }
  {
    std::ostringstream verRange;
    versionRange="";
    if (descriptor.isMinVersionSet() && descriptor.isMaxVersionSet())
    {
      verRange << "[Only between versions " << descriptor.getMinimumVersion() << " and " << descriptor.getMaximumVersion() << "]" ;
      versionRange = verRange.str();
    }
    else if (descriptor.isMinVersionSet())
    {
      verRange << "[Only for versions starting with " << descriptor.getMinimumVersion() << "]" ;
      versionRange = verRange.str();
    }
    else if (descriptor.isMaxVersionSet())
    {
      verRange << "[Only for versions up to " << descriptor.getMaximumVersion() << "]" ;
      versionRange = verRange.str();
    }
  }

}

const UnitInfo &findUnit(const int unit)
{
  for (const UnitInfo *it = Units::unitTable; it != Units::unitTable + Units::unitTableSize; ++it)
    if ((*it).Unit == unit)
      return *it;
  return *(Units::unitTable);
}

std::ostream &documentParameter(std::ostream &os, const std::string &name, const int parameter_unit, const std::string &description, const std::string &versionRange, const Descriptor &descriptor)
{
  if (descriptor.getCompositeParametricData<void>()) {
    os << laTexProtect(name) << " & " << laTexProtect(description) << " & ";

    const UnitInfo &unit_info = findUnit(parameter_unit);
    os << "\\multicolumn{2}{c}{See Table~\\ref{" << name << "_Composite_Params}} ";
  }
  else {
    os << laTexProtect(name);
    if (versionRange != "")
    {
      os << "\\newline" << "{\\normalfont " << versionRange << "}";
    }
    os << " & " << laTexProtect(description) << " & ";

    const UnitInfo &unit_info = findUnit(parameter_unit);
    os << unit_info.doc;
    os << " & ";

    if ((descriptor.isType<double>() && name == "TNOM") || name == "TEMP") 
    {
      double default_value = getDefaultValue<double>(descriptor);
      // os << default_value - CONSTCtoK;
      // os << CONSTREFTEMP - CONSTCtoK;
      if (default_value == 0.0)
      {
        os << "Ambient Temperature";
      }
      else if (default_value - CONSTCtoK < 0.0)
      {
        os << default_value;
      }
      else
      {
        os << default_value - CONSTCtoK;
      }
    }
    else
    {
      descriptor.getEntry().print(os);
    }
  }

  if (descriptor.getCompositeParametricData<void>()) {
    const ParametricData<void> &composite_parametric_data = *descriptor.getCompositeParametricData<void>();
    const OrderedParameterMap composite_parametric_map(composite_parametric_data.getMap().begin(), composite_parametric_data.getMap().end());

    std::ostringstream path;
    path << name << "_Composite_Params.tex";

    std::ofstream os(path.str().c_str(), std::ios_base::out);

    os << "% This table was generated by Xyce" << std::endl;

    laTexComposite(os, name, name, composite_parametric_data);
  }

  os << " \\\\ \\hline" << std::endl;

  return os;
}

std::ostream &
laTexComposite(
  std::ostream &                os,
  const std::string &           composite_name,
  const std::string &           composite_description,
  const ParametricData<void> &  parameters)
{
  // Place caption on table and add index entry
  std::string composite_description_lc = composite_description;
  std::transform(composite_description_lc.begin(), composite_description_lc.end(), composite_description_lc.begin(), (int (*)(int)) std::tolower);

  os << "\\index{" << laTexClean(composite_description_lc) << "!" << "composite parameters}" << std::endl
     << "\\begin{CompositeParamTableGenerated}{" << laTexProtect(composite_description) << " " << "Composite Parameters}"
     << "{" << composite_name << "_Composite_Params}" << std::endl;

  const OrderedParameterMap parameter_map(parameters.getMap().begin(), parameters.getMap().end());
  for (OrderedParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
    std::string parameter_name = (*it).first;
    const Descriptor &descriptor = *(*it).second;

    if ((descriptor.getExpressionAccess() & ParameterType::NO_DOC) == 0) {
      std::string parameter_description;
      ParameterUnit parameter_unit = U_INVALID;
      ParameterCategory parameter_category = CAT_UNKNOWN;
      std::string versionRange;
      unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description,versionRange);

      documentParameter(os, parameter_name, parameter_unit, parameter_description, versionRange, descriptor);
    }
  }

  os << "\\end{CompositeParamTableGenerated}" << std::endl;

  return os;
}

} // namespace <unnamed>

std::ostream &
laTexDevice(
  std::ostream &                os,
  const std::string &           device_name,
  const int                     device_level,
  const int                     type,
  const std::string &           device_description,
  const ParametricData<void> &  parameters,
  OutputMode::Mode              format)
{
  // Place caption on table and add index entry
  std::string device_description_lc = device_description;
  std::transform(device_description_lc.begin(), device_description_lc.end(), device_description_lc.begin(), (int (*)(int)) std::tolower);

  os << "\\index{" << laTexClean(device_description_lc) << "!" << (type ? "device model parameters" : "device instance parameters") << "}" << std::endl
     << "\\begin{DeviceParamTableGenerated}{" << laTexProtect(device_description) << " " << (type ? "Device Model Parameters" : "Device Instance Parameters") << "}"
     << "{" << device_name << "_" << device_level << (type ? "_Device_Model_Params" : "_Device_Instance_Params") << "}" << std::endl;

  // If catagorical listing, group by category
  if (format == OutputMode::DOC_CAT) {
    const OrderedParameterMap parameter_map(parameters.getMap().begin(), parameters.getMap().end());
    for (int category = CAT_NONE; category != CAT_MAX; ++category) {
      const std::string &header = categoryName(category);
      bool header_printed = false;

      for (OrderedParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
        std::string parameter_name = (*it).first;
        const Descriptor &descriptor = *(*it).second;

        if ((descriptor.getExpressionAccess() & ParameterType::NO_DOC) == 0) {
          std::string parameter_description;
          ParameterUnit parameter_unit = U_INVALID;
          ParameterCategory parameter_category = CAT_UNKNOWN;
          std::string versionRange;
          unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description,versionRange);

          if (parameter_category == category) {
            if (!header.empty() && !header_printed) {
              header_printed = true;

              os << std::endl
                 << "\\category{" << header << "}" << "\\\\ \\hline" << std::endl;
            }
            documentParameter(os, parameter_name, parameter_unit, parameter_description, versionRange, descriptor);
          }
        }
      }
    }
  }
  else {
    const OrderedParameterMap parameter_map(parameters.getMap().begin(), parameters.getMap().end());
    for (OrderedParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
      std::string parameter_name = (*it).first;
      const Descriptor &descriptor = *(*it).second;

      if ((descriptor.getExpressionAccess() & ParameterType::NO_DOC) == 0) {
        std::string parameter_description;
        ParameterUnit parameter_unit = U_INVALID;
        ParameterCategory parameter_category = CAT_UNKNOWN;
        std::string versionRange;
        unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description,versionRange);

        documentParameter(os, parameter_name, parameter_unit, parameter_description, versionRange, descriptor);
      }
    }
  }

  os << "\\end{DeviceParamTableGenerated}" << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce
