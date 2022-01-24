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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iomanip>
#include <set>

#include <unordered_map>
using std::unordered_map;

#include <unordered_set>
using std::unordered_set;

#include <N_DEV_Const.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Param.h>
#include <N_DEV_Message.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_NumericalJacobian.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_MachDepParams.h>

#include <N_ERH_ErrorMgr.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceModel::DeviceModel
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceModel::DeviceModel(
  const ModelBlock &            model_block,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceEntity(parametric_data, factory_block.solverState_, factory_block.deviceOptions_, model_block.getNetlistLocation()),
    name_(model_block.getName()),
    type_(model_block.getType()),
    level_(model_block.getLevel()),
    temperatureModel(""),
    doseModel(""),
    iModel(TEMP),
    iMethod(QUAD),
    base_temp(CONSTREFTEMP)
{}

std::ostream &
DeviceModel::printName(std::ostream &os) const
{
  return os << "model " << name_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::setModParams
// Purpose       : Set up parameter fits from model line
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/26/05
//-----------------------------------------------------------------------------
void DeviceModel::setModParams(const std::vector<Param> &params)
{
  std::vector<int> m_start;
  unordered_set<std::string> pname;
  std::vector<std::string> ptype_expr;

  int param_index = 0;
  m_start.push_back(param_index);
  for (std::vector<Param>::const_iterator mp = params.begin(); mp != params.end(); ++mp)
  {
    ++param_index;
    if ((*mp).tag() == "INDEPENDENT;PARAM")
    {
      m_start.push_back(param_index);
      pname.clear();
    }
    else
    {
      std::pair<unordered_set<std::string>::iterator,bool> ret = pname.insert((*mp).tag());
      if (ret.second==false)
      {
        UserError(*this) << "Duplicate specification of parameter " << (*mp).tag();
        return;
      }
      if (m_start.size() == 1)
      {
        // Collect all parameters that are expressions.
        if ((*mp).getType() == Util::EXPR)
        {
          ptype_expr.push_back( (*mp).tag() );
        }
      }
    }
  }
  
  // Make sure list of expression model parameters is unique and sorted.
  std::sort( ptype_expr.begin(), ptype_expr.end() );
  ptype_expr.erase(std::unique(ptype_expr.begin(), ptype_expr.end() ), ptype_expr.end() );

  if (m_start.size() == 1)
  {
    setParams(params);
  }
  else
  {
    Xyce::lout() << "\nThe temperature interpolation capability, invoked by TEMPMODEL=QUADRATIC in a .MODEL statement, is considered deprecated.  It will be removed from a future version of Xyce\n" << std::endl;;

    m_start.push_back(param_index + 1);

    // An interpolation method is present, first figure out what it
    // is and make sure that all models agree on the method.
    std::string tmod("");
    std::string dmod("");
    std::string modName("");

    for (int i = 0; i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        if (params[j].tag() == "TEMPMODEL" && params[j].stringValue() != "NONE")
        {
          if (i == 0)
          {
            tmod = params[j].stringValue();
          }
          else
          {
            if (tmod != params[j].stringValue())
            {
              UserError(*this) << "Inconsistent or missing TEMPMODEL parameter, " << params[j].stringValue() << " specified here, " << tmod << " specified previously";
              return;
            }
          }
        }
        if (params[j].tag() == "DOSEMODEL" && params[j].stringValue() != "NONE")
        {
          if (i == 0)
          {
            dmod = params[j].stringValue();
          }
          else
          {
            if (dmod != params[j].stringValue())
            {
              UserError(*this) << "Inconsistent or missing DOSEMODEL parameter, " << params[j].stringValue() << " specified here, " << dmod << " specified previously";
              return;
            }
          }
        }
      }
    }

    if (tmod == "" && dmod == "")
    {
      UserError(*this) << "Duplicate model specification implies parameter interpolation, TEMPMODEL or DOSEMODEL parameters must specified";
      return;
    }
    else if (tmod != "" && dmod != "")
    {
      UserError(*this) << "Only one of TEMPMODEL or DOSEMODEL parameters may be specified";
      return;
    }
    else if (tmod != "")
    {
      modName = tmod;
      iModel = TEMP;
    }
    else if (dmod != "")
    {
      modName = dmod;
      iModel = DOSE;
    }

    if (modName == "QUADRATIC")
    {
      iMethod = QUAD;
      if (m_start.size() != 4)
      {
        UserError(*this) << "Three model specifications required for QUADRATIC fit";
        return;
      }
      fit.resize(3);
    }
    else if (modName == "PWL")
    {
      iMethod = PWL;
      fit.resize(m_start.size() - 1);
    }
    else
    {
      UserError(*this) << "Only QUADRATIC or PWL interpolation method is supported";
      return;
    }

    // First find the params that vary between the specified models:
    unordered_map<std::string,double> basePars;
    std::vector<double> t;
    std::string par;

    for (int i = 0 ;i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        par = params[j].tag();
        ParameterMap::const_iterator parIt = getParameterMap().find(par);
        if (parIt == getParameterMap().end())
          DevelFatal(*this).in("DeviceModel::setModParams") << "Parameter " << par << " not found";

        const Descriptor &nPar = *(*parIt).second;
        if (params[j].given() && nPar.isType<double>())
        {
          if (iModel == TEMP)
          {
            if (par == "TNOM")
            {
              t.push_back(params[j].getImmutableValue<double>());
            }
          }
          else if (iModel == DOSE)
          {
            if (par == "DOSE")
            {
              t.push_back(params[j].getImmutableValue<double>());
            }
          }
          if (i == 0)
          {
            if (params[j].getType() != Util::EXPR)
            {
              basePars[par] = params[j].getImmutableValue<double>();
            }
          }
          else
          {
            if (std::binary_search( ptype_expr.begin(), ptype_expr.end(), par))
            {
              UserWarning(*this) << "Non-constant expression for parameter " << par << ", it will not interpolated";
            }
            else
            {
              if (basePars.find(par) == basePars.end())
              {
                UserError(*this) << "Unknown parameter " << params[j].tag() <<  " in temperature compensation .MODEL statement";
                return;
              }
              if (basePars[par] != params[j].getImmutableValue<double>() && params[j].given())
              {
                if (fitMap.find(par) == fitMap.end())
                {
                  fitMap[par] = fit[0].size();
                  fitParams.push_back(nPar.getMemberPtr<double>());
                  fit[0].push_back(basePars[par]);
                }
              }
            }
          }
        }
      }
    }

    if (t.size() != m_start.size()-1)
    {
      UserError(*this) << (iModel == TEMP ? "TNOM" : "DOSE") << " not specified in all .MODEL statements";
      return;
    }

    for (int i = 1; i < t.size(); ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        if (t[i] == t[j])
        {
          UserError(*this) << "Identical " << (iModel == TEMP ? "TNOM" : "DOSE") << " values in .MODEL statements";
          return;
        }
      }
    }

    // Now, collect the values to use for the fits:
    int nFit = fitMap.size();
    //int nSet = t.size();
    std::vector<std::vector<double> > temp(nFit);
    std::vector<std::vector<double> > val(nFit);
    oldParams.resize(nFit);

    for (int i = 1; i < fit.size(); ++i)
    {
      fit[i].resize(nFit);
    }
    min_par.resize(nFit);
    max_par.resize(nFit);
    parType.resize(nFit);
    std::map<std::string, int>::iterator fm = fitMap.begin();
    std::map<std::string, int>::iterator fm_end = fitMap.end();

    for ( ; fm != fm_end ; ++fm)
    {
      par = (*fm).first;
      ParameterMap::const_iterator parIt = getParameterMap().find(par);
      if (parIt == getParameterMap().end())
        throw std::runtime_error(std::string("Parameter ") + par + " not found");
      {
        const Descriptor &nPar = *(*parIt).second;

        if (nPar.getExpressionAccess() & ParameterType::LOG_T_DEP)
        {
          parType[fitMap[par]] = LOG_FIT;
        }
        else
        {
          parType[fitMap[par]] = LINEAR_FIT;
        }
      }
    }

    base_temp = t[0];
    if (iModel == TEMP)
    {
      base_temp += CONSTCtoK;
    }

    for (int i = 0; i < nFit; ++i)
    {
      temp[i].push_back(t[0]);
      val[i].push_back(fit[0][i]);
      if (parType[i] == LOG_FIT)
      {
        if (val[i][0] <= 0)
        {
          UserError message(*this);
          message << "Negative parameter value for log interpolation based parameter ";
          fm = fitMap.begin();
          for ( ; fm != fm_end ; ++fm)
          {
            if ((*fm).second == i)
            {
              par = (*fm).first;
              break;
            }
          }
          message << par;
          return;
        }
        val[i][0] = log(val[i][0]);
      }
      min_par[i] = val[i][0];
      max_par[i] = val[i][0];
    }

    for (int i = 1; i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        if (params[j].getType() == Util::DBLE && params[j].given())
        {
          par = params[j].tag();
          std::map<std::string, int>::iterator fm1 = fitMap.find(par);
          if (fm1 != fitMap.end())
          {
            int k = fm1->second;
            temp[k].push_back(t[i]);
            double p = params[j].getImmutableValue<double>();
            if (parType[k] == LOG_FIT)
            {
              if (p <= 0)
              {
                UserError message(*this);
                message << "Negative parameter value for log interpolation based parameter ";
                fm = fitMap.begin();
                for ( ; fm != fm_end ; ++fm)
                {
                  if ((*fm).second == k)
                  {
                    par = (*fm).first;
                    break;
                  }
                }
                message << par;
                return;
              }
              p = log(p);
            }
            val[k].push_back(p);
            if (p > max_par[k])
              max_par[k] = p;
            if (p < min_par[k])
              min_par[k] = p;
          }
        }
      }
    }

    // Finally, do the actual fits:
    if (fitMap.size() > 0)
    {
      if (iMethod == QUAD)
      {
        std::map<std::string, int>::iterator f;
        for (f=fitMap.begin() ; f!=fitMap.end() ; ++f)
        {
          int i = (*f).second;
          if (temp[i].size() == 2)
          {
            fit[0][i] = val[i][0];
            fit[1][i] = (val[i][1] - val[i][0])/(temp[i][1] - temp[i][0]);
            fit[2][i] = 0;
          }
          else if (temp[i].size() == 3)
          {
            fit[0][i] = val[i][0];
            double x1,x2,y1,y2;
            x1 = temp[i][1] - temp[i][0];
            y1 = val[i][1];
            x2 = temp[i][2] - temp[i][0];
            y2 = val[i][2];
            fit[2][i] = (x2*y1-x1*y2-fit[0][i]*(x2-x1))/(x2*x1*x1-x1*x2*x2);
            fit[1][i] = (y1-fit[2][i]*x1*x1-fit[0][i])/x1;
          }
          else
          {
            DevelFatal(*this).in("setModParams") << "Internal error in DeviceModel, illegal number of fit points for parameter " << (*f).first;
          }
          if ((*f).first == "TNOM")
          {
            fit[0][i] += CONSTCtoK;
          }
        }
      }
      else if (iMethod == PWL)
      {
        int nT = fit.size();
        std::map<double,int> tOrd;
        for (int i = 0; i < nT; ++i)
        {
          tOrd[t[i]] = 0;
        }
        int i = 0;
        std::map<double,int>::iterator tOrd_i = tOrd.begin();
        std::map<double,int>::iterator tOrd_end = tOrd.end();
        for ( ; tOrd_i != tOrd_end; ++tOrd_i)
        {
          if (iModel == TEMP)
          {
            base.push_back((*tOrd_i).first+CONSTCtoK);
          }
          else
          {
            base.push_back((*tOrd_i).first);
          }
          (*tOrd_i).second = i++;
        }
        std::map<std::string, int>::iterator f;
        std::map<std::string, int>::iterator f_end;
        std::vector<bool> p(nT,false);
        f=fitMap.begin();
        f_end=fitMap.end();
        for ( ; f!=f_end; ++f)
        {
          i = (*f).second;
          for (int j = 0 ;j < nT; ++j)
          {
            p[j] = false;
          }
          for (int j = 0; j < temp[i].size() ; ++j)
          {
            int ind = tOrd[temp[i][j]];
            p[ind] = true;
            fit[ind][i] = val[i][j];
          }
          for (int j = 0; j < nT ; ++j)
          {
            if (!p[j])
            {
              int k_lo = j;
              int k_hi = j;
              while (k_lo >= 0 && !p[k_lo])
              {
                k_lo--;
              }
              while (k_hi <nT && !p[k_hi])
              {
                ++k_hi;
              }
              if (k_lo == -1)
              {
                if (k_hi < nT)
                {
                  fit[j][i] = fit[k_hi][i];
                }
                else
                {
                  DevelFatal(*this).in("DeviceModel::setModParams") <<"DeviceModel::setModParams: Internal error forming PWL interpolation";
                }
              }
              else
              {
                if (k_hi < nT)
                {
                  double frac = (base[j]-base[k_lo])/(base[k_hi]-base[k_lo]);
                  fit[j][i] = fit[k_hi][i]*frac+fit[k_lo][i]*(1-frac);
                }
                else
                {
                  fit[j][i] = fit[k_lo][i];
                }
              }
            }
            if ((*f).first == "TNOM")
            {
              fit[j][i] += CONSTCtoK;
            }
          }
        }
      }
    }

    // params.resize(m_start[1]-1);
    // setParams(params);
    std::vector<Param> remaining_params(&params[0], &params[m_start[1] - 1]);
    setParams(remaining_params);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::~DeviceModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceModel::~DeviceModel()
{
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::saveParams
// Purpose       : save existing param values before fitting params to temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
void DeviceModel::saveParams()
{
  int nFit = fitMap.size();
  int i;

  if (nFit == 0)
  {
    return;
  }

  for (i=0 ; i<nFit ; ++i)
  {
    oldParams[i] = this->*(fitParams[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolateTNOM
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolateTNOM(double t)
{
  if (iModel != TEMP)
  {
    return false;
  }

  return interpolate(t);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolateDOSE
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/11/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolateDOSE(double d)
{
  if (iModel != DOSE)
  {
    return false;
  }

  return interpolate(d);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolated
// Purpose       : returns true if an interpolation is in effect
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
bool DeviceModel::interpolated()
{
  return (fitMap.size() > 0);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolate
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/11/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolate(double t)
{
  int nFit = fitMap.size();
  int i, j, k_hi, k_lo;
  double del;
  double frac, p;

  if (nFit == 0)
  {
    return false;
  }

  if (iMethod == QUAD)
  {
    del = t - base_temp;
    //    for (i=0 ; i<nFit ; ++i)
    std::map<std::string,int>::iterator fp;
    std::map<std::string,int>::iterator fm_begin=fitMap.begin();
    std::map<std::string,int>::iterator fm_end=fitMap.end();
    for (fp=fm_begin; fp != fm_end; fp++)
    {
      i=fp->second;
      p = (fit[2][i]*del + fit[1][i])*del + fit[0][i];

      if (p > max_par[i] && i>0)
      {
        this->*(fitParams[i]) = max_par[i];
      }
      else if (p < min_par[i] && i>0)
      {
        this->*(fitParams[i]) = min_par[i];
      }
      else
      {
        this->*(fitParams[i]) = p;
      }
    }
  }
  else if (iMethod == PWL)
  {
    del = t;
    k_hi = 0;
    for (j=0 ; j<fit.size() ; ++j)
    {
      if (base[j] >= del)
      {
        break;
      }
      k_hi = j+1;
    }
    if (k_hi == 0)
    {
      frac = 0;
    }
    else if (k_hi == fit.size())
    {
      k_hi = fit.size()-1;
      frac = 1;
    }
    else
    {
      k_lo = k_hi-1;
      frac = (del-base[k_lo])/(base[k_hi]-base[k_lo]);
    }
    if (frac == 1)
    {
      for (i=0 ; i<nFit ; ++i)
      {
        this->*(fitParams[i]) = fit[k_hi][i];
      }
    }
    else
    {
      for (i=0 ; i<nFit ; ++i)
        this->*(fitParams[i]) = fit[k_hi][i]*frac + fit[k_lo][i]*(1-frac);
    }
  }
  for (i=0 ; i<nFit ; ++i)
  {
    if (parType[i] == LOG_FIT)
    {
      this->*(fitParams[i]) = exp(this->*(fitParams[i]));
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::restoreParams
// Purpose       : restore previously saved param values
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
void DeviceModel::restoreParams()
{
  int nFit = fitMap.size();
  int i;

  if (nFit == 0)
  {
    return;
  }

  for (i=0 ; i<nFit ; ++i)
  {
    this->*(fitParams[i]) = oldParams[i];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::setupBaseInstanceContainer
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/25/2018
//-----------------------------------------------------------------------------
void DeviceModel::setupBaseInstanceContainer()
{
  UserError(*this) << "DeviceModel::setupBaseInstanceContainer: not implemented for device " << getName() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::getNumericalSensitivity
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/24/2018
//-----------------------------------------------------------------------------
bool DeviceModel::getNumericalSensitivity ( const std::string & name,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec )
{
  std::string paramName = Util::paramNameFromFullParamName(name);
  double origParamValue = 0.0;
  bool found = getParam(paramName, origParamValue); // const?
  if (found)
  {
    if(baseInstanceContainer.empty() ) { setupBaseInstanceContainer(); }

    if(baseInstanceContainer.empty() ) 
    { 
      UserError(*this) << "DeviceModel::getNumericalSensitivity: Failed to setup baseInstance container!" << std::endl;
    }

    std::vector<Xyce::Device::DeviceInstance *>::iterator begin = baseInstanceContainer.begin();
    std::vector<Xyce::Device::DeviceInstance *>::iterator end   = baseInstanceContainer.end();
    std::vector<Xyce::Device::DeviceInstance *>::iterator iterInst;

    int numInst = baseInstanceContainer.size();

    int lidSizeMax = -1;
    int lidSizeMin = 1e9;

    int lidStateSizeMax = -1;
    int lidStateSizeMin = 1e9;

    for (iterInst=begin; iterInst != end; ++iterInst)
    {
      Xyce::Device::DeviceInstance * instPtr = *(iterInst);
      instPtr->consolidateDevLIDs();

      int tmpSize = instPtr->getDevLIDs().size();
      if(tmpSize > lidSizeMax) { lidSizeMax = tmpSize; }
      if(tmpSize < lidSizeMin) { lidSizeMin = tmpSize; }

      int tmpStateSize = instPtr->getStaLIDVec().size();
      if(tmpStateSize > lidStateSizeMax) { lidStateSizeMax = tmpStateSize; }
      if(tmpStateSize < lidStateSizeMin) { lidStateSizeMin = tmpStateSize; }
    }

    if (lidSizeMax != lidSizeMin)
    {
      UserError(*this) << "Unable to compute numerical sensitivity for model " << getName() << ".  Number of instances = " << numInst << std::endl;
    }

    int lidSize = lidSizeMax;
    int lidStateSize = lidStateSizeMax;

    if (lidSize>0)
    {
      MatrixLoadData & mlData = (*begin)->getMatrixLoadData();
      const ExternData & extData = (*begin)->getExternData();

      if (!(FindicesVec.empty())) {  FindicesVec.clear(); }
      if (!(QindicesVec.empty())) {  QindicesVec.clear(); }
      if (!(BindicesVec.empty())) {  BindicesVec.clear(); }

      std::vector<int> stateIndicesVec;

      unordered_map<int,int> lidCheckMap;
      unordered_map<int,int> stateLidCheckMap;
      for (iterInst=begin; iterInst != end; ++iterInst)
      {
        Xyce::Device::DeviceInstance * instPtr = *(iterInst);
        const IdVector & devLIDs = instPtr->getDevLIDs();

        for (int i=0;i<devLIDs.size();++i)
        {
          int newID = devLIDs[i];
          if (lidCheckMap.find(newID) == lidCheckMap.end())
          {
            FindicesVec.push_back(newID);
            lidCheckMap[newID] = i;
          }
        }

        const IdVector & stateDevLIDs = instPtr->getStaLIDVec();
        for (int i=0;i<stateDevLIDs.size();++i)
        {
          int newID = stateDevLIDs[i];
          if (stateLidCheckMap.find(newID) == stateLidCheckMap.end())
          {
            stateIndicesVec.push_back(newID);
            stateLidCheckMap[newID] = i;
          }
        } 
      }

      // this sort is a convenience for debug output.
      std::sort(FindicesVec.begin(),FindicesVec.end());

      QindicesVec.insert(QindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 
      BindicesVec.insert(BindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 

      int solSize=FindicesVec.size();
      int stateSize=stateIndicesVec.size();

      Linear::Vector & Fvec          = (*extData.daeFVectorPtr);
      Linear::Vector & Qvec          = (*extData.daeQVectorPtr);
      Linear::Vector & Bvec          = (*extData.daeBVectorPtr);

      Linear::Vector & currSol      = (*extData.currSolVectorPtr);
      Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

      Linear::Vector & lastSta      = (*extData.lastStaVectorPtr);
      Linear::Vector & currSta      = (*extData.currStaVectorPtr);
      Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
      Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

      mlData.resizeSolnSizedVectors (solSize);
      mlData.resizeStateSizedVectors (stateSize);

      std::vector<double> & saveF = mlData.saveF;
      std::vector<double> & pertF = mlData.pertF;
      std::vector<double> & origF = mlData.origF;

      std::vector<double> & saveQ = mlData.saveQ;
      std::vector<double> & pertQ = mlData.pertQ;
      std::vector<double> & origQ = mlData.origQ;

      std::vector<double> & saveB = mlData.saveB;
      std::vector<double> & pertB = mlData.pertB;
      std::vector<double> & origB = mlData.origB;

      std::vector<double> & saveLastState = mlData.saveLastState;
      std::vector<double> & saveCurrState = mlData.saveCurrState;
      std::vector<double> & saveNextState = mlData.saveNextState;
      std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

      saveF.assign(saveF.size(),0.0);
      pertF.assign(pertF.size(),0.0);

      saveQ.assign(saveQ.size(),0.0);
      pertQ.assign(pertQ.size(),0.0);

      saveB.assign(saveB.size(),0.0);
      pertB.assign(pertB.size(),0.0);

      dfdpVec.resize(saveF.size(),0.0);
      dqdpVec.resize(saveQ.size(),0.0);
      dbdpVec.resize(saveB.size(),0.0);

      // save RHS information
      std::vector<int> origFlagVec(numInst,-1);
      int i=0;
      for (iterInst=begin; iterInst != end; ++iterInst, ++i)
      {
        Xyce::Device::DeviceInstance * instPtr = *(iterInst);
        origFlagVec[i] = static_cast<int>(instPtr->getOrigFlag());
      }

      for (i=0; i<saveF.size(); ++i)
      {
        saveF[i]      = Fvec[FindicesVec[i]];
        saveQ[i]      = Qvec[QindicesVec[i]];
        saveB[i]      = Bvec[BindicesVec[i]];
      }
      
      for (i=0; i<saveLastState.size(); ++i)
      {
        saveLastState[i] = lastSta[stateIndicesVec[i]];
        saveCurrState[i] = currSta[stateIndicesVec[i]];
        saveNextState[i] = nextSta[stateIndicesVec[i]];
        saveStateDerivs[i] = nextStaDeriv[stateIndicesVec[i]];
      }

      // zero out the RHS, and re-load, so that we have only the
      // elements from one device present.
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = 0.0;
        Qvec[QindicesVec[i]] = 0.0;
        Bvec[BindicesVec[i]] = 0.0;
      }
      // load new values for F,Q,B to get the "orig" values
      for (iterInst=begin; iterInst != end; ++iterInst, ++i)
      {
        Xyce::Device::DeviceInstance * instPtr = *(iterInst);
        instPtr->numJacPtr->loadLocalDAEVectorsIncludingB(*instPtr);
      } 

      // Save RHS for just this instance:
      for (int i=0 ; i<origF.size(); ++i)
      {
        origF[i] = Fvec[FindicesVec[i]];
        origQ[i] = Qvec[QindicesVec[i]];
        origB[i] = Bvec[BindicesVec[i]];
      }

      // perturb the parameter
      double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
      double sqrtEta= std::sqrt(epsilon);
      double dP = sqrtEta * fabs( origParamValue );
      double minDouble = Util::MachineDependentParams::DoubleMin();
      if (dP < minDouble)
      {
        dP = sqrtEta;
      }
  
      double newParamValue = origParamValue + dP;

      if (paramName == "")
      {
        setDefaultParam(newParamValue, false);
      }
      else
      {
        setParam(paramName, newParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();
      
      // zero out F,Q,B
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = 0.0;
        Qvec[QindicesVec[i]] = 0.0;
        Bvec[BindicesVec[i]] = 0.0;
      }

      // load new values for F,Q,B
      for (iterInst=begin; iterInst != end; ++iterInst, ++i)
      {
        Xyce::Device::DeviceInstance * instPtr = *(iterInst);
        instPtr->numJacPtr->loadLocalDAEVectorsIncludingB(*instPtr);
      }

      // save the perturbed DAE vectors
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        pertF[i] = Fvec[FindicesVec[i]];
        pertQ[i] = Qvec[BindicesVec[i]];
        pertB[i] = Bvec[QindicesVec[i]];
      }

      // compute the derivatives
      for (int i=0; i<pertF.size(); ++i)
      {
        dfdpVec[i] = (pertF[i] - origF[i])/dP;
        dqdpVec[i] = (pertQ[i] - origQ[i])/dP;
        dbdpVec[i] = (pertB[i] - origB[i])/dP;
      }

      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() << paramName << ": ";
        Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "deviceSens_dp = " << dP << std::endl;

        for (int k1=0; k1<pertF.size(); ++k1)
        {
          int ii=FindicesVec[k1];
          Xyce::dout() 
            <<"pertF["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertF[k1]
            <<" origF["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origF[k1]
            <<" dfdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dfdpVec[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl;
        for (int k1=0; k1<pertQ.size(); ++k1)
        {
          int ii=QindicesVec[k1];
          Xyce::dout() 
            <<"pertQ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertQ[k1]
            <<" origQ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origQ[k1]
            <<" dqdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dqdpVec[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl;
        for (int k1=0; k1<pertB.size(); ++k1)
        {
          int ii=BindicesVec[k1];
          Xyce::dout() 
            <<"pertB["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertB[k1]
            <<" origB["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origB[k1]
            <<" dbdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dbdpVec[k1]
            <<std::endl;
        }
      }

      // restore everything:
      if (paramName == "")
      {
        setDefaultParam(origParamValue, false);
      }
      else
      {
        setParam(paramName, origParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();      
      
      for (i=0, iterInst=begin; iterInst != end; ++iterInst, ++i)
      {
        Xyce::Device::DeviceInstance * instPtr = *(iterInst);
        bool origFlag = static_cast<bool> (origFlagVec[i]);
        instPtr->setOrigFlag(origFlag);
      }

      for (i=0; i<saveF.size(); ++i)
      {
        Fvec[FindicesVec[i]] = saveF[i];
        Qvec[QindicesVec[i]] = saveQ[i];
        Bvec[BindicesVec[i]] = saveB[i];
      }

      for (i=0; i<saveLastState.size(); ++i)
      {
        lastSta[stateIndicesVec[i]] = saveLastState[i];
        currSta[stateIndicesVec[i]] = saveCurrState[i];
        nextSta[stateIndicesVec[i]] = saveNextState[i];
        nextStaDeriv[stateIndicesVec[i]] = saveStateDerivs[i];
      }
    }
  }

  return found;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::getNumericalMatrixSensitivity
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/19/2019
//-----------------------------------------------------------------------------
bool DeviceModel::getNumericalMatrixSensitivity ( 
    const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & FindicesVec,
    std::vector<int> & QindicesVec,
    std::vector< std::vector<int> > & FjacLIDs,
    std::vector< std::vector<int> > & QjacLIDs )
{
  double origParamValue = 0.0;
  bool found = getParam(paramName, origParamValue); // const?

  if (found)
  {
    if(baseInstanceContainer.empty() ) { setupBaseInstanceContainer(); }
    if(baseInstanceContainer.empty() ) 
    { 
      UserError(*this) << "DeviceModel::getNumericalMatrixSensitivity: Failed to setup baseInstance container!" << std::endl;
    }

    std::vector<Xyce::Device::DeviceInstance *>::iterator begin = baseInstanceContainer.begin();
    std::vector<Xyce::Device::DeviceInstance *>::iterator end   = baseInstanceContainer.end();
    std::vector<Xyce::Device::DeviceInstance *>::iterator iterInst;

    for (iterInst=begin; iterInst != end; ++iterInst)
    {
      Xyce::Device::DeviceInstance * instPtr = *(iterInst);
      instPtr->consolidateDevLIDs();
    }

    MatrixLoadData & mlData = (*begin)->getMatrixLoadData();

    int numRows=0, numState=0;
    // get the TOTAL size for the numerical jacobian from ALL instances in this model
    for (iterInst=begin;iterInst!=end;++iterInst)
    {
      numRows += (*iterInst)->jacobianStamp().size();
      numState += (*iterInst)->getStaLIDVec().size();
    }

    // check if numRows = zero.  if so, return, as there is nothing to do
    if (numRows==0) { return true; }
    if (numRows > mlData.saveF.size())
    {
      mlData.resizeTestJacSolData(numRows);
      mlData.resizeTestJacQData(numRows);
    }

    if (numState > mlData.saveCurrState.size())
    {
      mlData.resizeTestJacStateData(numState);
    }

    if (!(FindicesVec.empty())) {  FindicesVec.clear(); }
    if (!(QindicesVec.empty())) {  QindicesVec.clear(); }
    if (!(FjacLIDs.empty())) {  FjacLIDs.clear(); }
    if (!(QjacLIDs.empty())) {  QjacLIDs.clear(); }

    // NOTE:  the use of the lidCheckMap was an attempt to make sure that 
    // I didn't double count anything.  I eventually gave up on that, however,
    // as I screwed up the logic and that caused problems.
    //
    // Now, I go ahead and compute things more than once, if the same 
    // element appears multiple times in the LIDs.  That should be OK,
    // as the AC class that uses this information will just overwrite 
    // any duplicates.
    unordered_map<int,int> lidCheckMap;
    unordered_map<int,int> stateLidCheckMap;

    std::vector<int> stateIndicesVec;

    for (iterInst=begin;iterInst!=end;++iterInst)
    {
      const IdVector & devLIDs = (*iterInst)->getDevLIDs();
      const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();
      for (int i=0;i<devLIDs.size();++i)
      {
        int newID = devLIDs[i];
        const std::vector<int> & jacRow = devJacLIDs[i]; 

        //if (lidCheckMap.find(newID) == lidCheckMap.end())
        {
          FindicesVec.push_back(newID);
          lidCheckMap[newID] = i;
          FjacLIDs.push_back(jacRow); 
        }
      }

      const IdVector & stateDevLIDs = (*iterInst)->getStaLIDVec();
      for (int i=0;i<stateDevLIDs.size();++i)
      {
        int newID = stateDevLIDs[i];
        //if (stateLidCheckMap.find(newID) == stateLidCheckMap.end())
        {
          stateIndicesVec.push_back(newID);
          stateLidCheckMap[newID] = i;
        }
      }
    }
    QindicesVec.insert(QindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 

    if (FjacLIDs.size() != FindicesVec.size() ) // Houston we have a problem
    {
      Report::DevelFatal().in("DeviceInstance::getNumericalMatrixSensitivity") << "Internal Error 1";
    }

    QjacLIDs.resize(FjacLIDs.size());
    for (int i=0 ; i<FjacLIDs.size(); ++i)
    {
      int jCol = FjacLIDs[i].size(); QjacLIDs[i].resize(jCol);
      for (int j=0 ; j<jCol ; ++j) { QjacLIDs[i][j] = FjacLIDs[i][j]; }
    }

    // grab some references
    const ExternData & extData = (*begin)->getExternData();
    Linear::Vector & Fvec          = (*extData.daeFVectorPtr);
    Linear::Vector & Qvec          = (*extData.daeQVectorPtr);
    Linear::Vector & Bvec          = (*extData.daeBVectorPtr);
    Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

    Linear::Vector & lastSta      = (*extData.lastStaVectorPtr);
    Linear::Vector & currSta      = (*extData.currStaVectorPtr);
    Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
    Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

    Linear::Matrix & dQdxMat      = (*extData.dQdxMatrixPtr);
    Linear::Matrix & dFdxMat      = (*extData.dFdxMatrixPtr);

    int solSize=FindicesVec.size();
    int stateSize=stateIndicesVec.size();
    mlData.resizeSolnSizedVectors (solSize);
    mlData.resizeStateSizedVectors (stateSize);

    std::vector< std::vector<double> > & numJacF = mlData.numJac;
    std::vector< std::vector<double> > & saveJacF = mlData.saveJac;
    std::vector< std::vector<double> > & devJacF = mlData.devJac;

    std::vector< std::vector<double> > & numJacQ = mlData.numJacQ;
    std::vector< std::vector<double> > & saveJacQ = mlData.saveJacQ;
    std::vector< std::vector<double> > & devJacQ = mlData.devJacQ;

    std::vector< std::vector<int> > & stencil = mlData.stencil;

    std::vector<double> & saveF = mlData.saveF;
    std::vector<double> & saveQ = mlData.saveQ;
    std::vector<double> & saveB = mlData.saveB;

    std::vector<double> & saveLastState = mlData.saveLastState;
    std::vector<double> & saveCurrState = mlData.saveCurrState;
    std::vector<double> & saveNextState = mlData.saveNextState;
    std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

    // initialize/ zero things out
    saveF.assign(saveF.size(),0.0);
    saveQ.assign(saveQ.size(),0.0);
    saveB.assign(saveB.size(),0.0);

    // zero out the matrix stuff
    for (int i=0;i<numJacF.size();++i)
    {
      numJacF[i].assign(numJacF[i].size(),0.0);
      saveJacF[i].assign(saveJacF[i].size(),0.0);
      devJacF[i].assign(devJacF[i].size(),0.0);
      stencil[i].assign(stencil[i].size(),0);

      numJacQ[i].assign(numJacQ[i].size(),0.0);
      saveJacQ[i].assign(saveJacQ[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
    }

    d_dfdx_dp.clear();
    d_dqdx_dp.clear();
    d_dfdx_dp.resize(FindicesVec.size());
    d_dqdx_dp.resize(FindicesVec.size());

    for (int i=0;i<d_dfdx_dp.size();++i)
    {
      d_dfdx_dp[i].resize(d_dfdx_dp.size(),0.0);
      d_dqdx_dp[i].resize(d_dqdx_dp.size(),0.0);
    }

    // save RHS information
    std::vector<int> origFlags( baseInstanceContainer.size() );
    int ii=0;
    for (iterInst=begin;iterInst!=end;++iterInst,++ii)
    {
      origFlags[ii] = static_cast<int> (  (*iterInst)->getOrigFlag() );
    }

    for (int i=0; i<saveF.size(); ++i)
    {
      saveF[i]      = Fvec[FindicesVec[i]];
      saveQ[i]      = Qvec[FindicesVec[i]];
      saveB[i]      = Bvec[FindicesVec[i]];
    }
    //const std::vector<int> & devStateLIDs = getStaLIDVec();
    for (int i=0; i<saveLastState.size(); ++i)
    {
      saveLastState[i] = lastSta[stateIndicesVec[i]];
      saveCurrState[i] = currSta[stateIndicesVec[i]];
      saveNextState[i] = nextSta[stateIndicesVec[i]];
      saveStateDerivs[i] = nextStaDeriv[stateIndicesVec[i]];
    }

    // re-load for just this instance:
    for (iterInst=begin;iterInst!=end;++iterInst,++ii)
    {
      (*iterInst)->numJacPtr->loadLocalDAEVectorsIncludingB(*(*iterInst));
    }

    // save the Jacobian matrix information. 
    //
    // the number of rows in saveF is equal to a union of all the rows from all the instances.
    // So this number (-1) is what "ii" should end up at.
    //
    // Or it should be the sum of numRows from each instance.  If done this way, then some rows might
    // be represented more than once.  The code below is written with this logic.
    // Essentially this will be a vector of jacobian stamps, push into the same STL object.
    {
      int ii=0; // master row index
      for (iterInst=begin;iterInst!=end;++iterInst)
      {
        const IdVector & devLIDs = (*iterInst)->getDevLIDs();
        const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();

        for (int i=0 ; i<devLIDs.size(); ++i, ++ii)
        {
          int jCol = devJacLIDs[i].size();
          for (int j=0 ; j<jCol ; ++j)
          {
            double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
            saveJacF[ii][j] = valF;
            double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
            saveJacQ[ii][j] = valQ;
          }
        }
      }
    }

    // Zeroing needs to be done after all saved values are
    // recorded because there can be multiple references
    // to the same element
    {
      for (iterInst=begin;iterInst!=end;++iterInst)
      {
        const IdVector & devLIDs = (*iterInst)->getDevLIDs();
        const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();

        for (int i=0 ; i<devLIDs.size(); ++i)
        {
          int jCol = devJacLIDs[i].size();
          for (int j=0 ; j<jCol ; ++j)
          {
            dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
            dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
          }
        }
      }
    }

    // Now that the original load has been zeroed out, re-load the
    // analytic contributions, to get the contributions from *just* this
    // device.
    {
      for (iterInst=begin;iterInst!=end;++iterInst)
      {
        const IdVector & devLIDs = (*iterInst)->getDevLIDs();
        const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();

        (*iterInst)->loadDAEdQdx ();
        (*iterInst)->loadDAEdFdx ();
      }
    }

    for (int i=0 ; i<devJacF.size() ; ++i)
    {
      devJacF[i].assign(devJacF[i].size(),0.0);
      devJacQ[i].assign(devJacQ[i].size(),0.0);
      stencil[i].assign(stencil[i].size(),0);
    }

    for (iterInst=begin;iterInst!=end;++iterInst)
    {
      const IdVector & devLIDs = (*iterInst)->getDevLIDs();
      const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();
      const std::vector< std::vector<int> > & jacStamp   = (*iterInst)->jacobianStamp();

      int ii=0; // master row index
      for (int i=0 ; i<devJacLIDs.size() ; ++i,++ii)
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
          devJacF[ii][jacStamp[i][j]] = valF;
          double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
          devJacQ[ii][jacStamp[i][j]] = valQ;
          stencil[ii][jacStamp[i][j]] = 1;
        }
      }
    }

    // zero again
    {
      for (iterInst=begin;iterInst!=end;++iterInst)
      {
        const IdVector & devLIDs = (*iterInst)->getDevLIDs();
        const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();

        for (int i=0 ; i<devLIDs.size(); ++i)
        {
          int jCol = devJacLIDs[i].size();
          for (int j=0 ; j<jCol ; ++j)
          {
            dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
            dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
          }
        }
      }
    }

    // perturb the parameter
    double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
    double sqrtEta= std::sqrt(epsilon);
    double dP = sqrtEta * fabs( origParamValue );
    double minDouble = Util::MachineDependentParams::DoubleMin();
    if (dP < minDouble)
    {
      dP = sqrtEta;
    }

    double newParamValue = origParamValue + dP;

    if (paramName == "")
    {
      setDefaultParam(newParamValue, false);
    }
    else
    {
      setParam(paramName, newParamValue, false);
    }

    processParams (); // if this "entity" is a model, then need to
    // also do a  "processParams" on the related instances.
    processInstanceParams();

    // load new values for F,Q,B  (this is done here b/c Jacobian elements computed here)
    for (iterInst=begin; iterInst != end; ++iterInst)
    {
      Xyce::Device::DeviceInstance * instPtr = *(iterInst);
      instPtr->numJacPtr->loadLocalDAEVectorsIncludingB(*instPtr);
      instPtr->loadDAEdQdx ();
      instPtr->loadDAEdFdx ();
    }

    // Save the new matrix, compute the derivatives
    for (iterInst=begin;iterInst!=end;++iterInst)
    {
      const IdVector & devLIDs = (*iterInst)->getDevLIDs();
      const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();
      const std::vector< std::vector<int> > & jacStamp   = (*iterInst)->jacobianStamp();

      int ii=0; // will be the master row index

      for (int i=0 ; i<devJacLIDs.size() ; ++i,++ii)      
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          double pertValF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
          numJacF[ii][j] = pertValF;
          double pertValQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
          numJacQ[ii][j] = pertValQ;

          d_dfdx_dp[i][j] = (pertValF - devJacF[ii][j])/dP;
          d_dqdx_dp[i][j] = (pertValQ - devJacQ[ii][j])/dP;
        }
      }
    }

    // restore everything:
    if (paramName == "")
    {
      setDefaultParam(origParamValue, false);
    }
    else
    {
      setParam(paramName, origParamValue, false);
    }

    processParams (); // if this "entity" is a model, then need to
    // also do a  "processParams" on the related instances.
    processInstanceParams();

    {
      int ii=0;
      for (iterInst=begin;iterInst!=end;++iterInst,++ii)
      {
        (*iterInst)->setOrigFlag( static_cast<bool>(origFlags[ii])  );
      }
    } 
    
    for (int i=0; i<FindicesVec.size(); ++i)
    {
      Fvec[FindicesVec[i]] = saveF[i];
      Qvec[FindicesVec[i]] = saveQ[i];
      Bvec[FindicesVec[i]] = saveB[i];
    }

    for (int i=0; i<saveLastState.size(); ++i)
    {
      lastSta[stateIndicesVec[i]] = saveLastState[i];
      currSta[stateIndicesVec[i]] = saveCurrState[i];
      nextSta[stateIndicesVec[i]] = saveNextState[i];
      nextStaDeriv[stateIndicesVec[i]] = saveStateDerivs[i];
    }

    // restore the Jacobian matrix information
    {
      int ii=0; // master row index
      for (iterInst=begin;iterInst!=end;++iterInst)
      {
        const IdVector & devLIDs = (*iterInst)->getDevLIDs();
        const std::vector<IdVector > & devJacLIDs = (*iterInst)->getDevJacLIDs();

        for (int i=0 ; i<devLIDs.size(); ++i, ++ii)
        {
          int jCol = devJacLIDs[i].size();
          for (int j=0 ; j<jCol ; ++j)
          {
            dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacF[ii][j];
            dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = saveJacQ[ii][j];
          }
        }
      }
    }


  } // found

  return found;
}

} // namespace Device
} // namespace Xyce
