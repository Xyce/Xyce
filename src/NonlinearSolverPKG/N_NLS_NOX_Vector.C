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

//-------------------------------------------------------------------------
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_Vector.h"
#include "N_NLS_NOX.h"
#include "N_LAS_Vector.h"
#include "N_LAS_System.h"
#include "N_LAS_Builder.h"
#include "N_ERH_ErrorMgr.h"

// ----------   Namespaces   ----------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

// ----------   Code   ----------

Vector::Vector(Linear::Vector& vector, Linear::System& lasSys) :
  vectorPtr_(&vector),
  lasSys_(lasSys),
  doDelete_(false)
{
  
}


Vector::Vector(const Vector& source, NOX::CopyType type) :
  vectorPtr_(0),
  lasSys_(source.lasSys_),
  doDelete_(true)
{
  vectorPtr_ = lasSys_.builder().createVector();
  if (vectorPtr_ == 0) {
    error("Vector Copy Constructor - unable to create vector");
  }

  if (type == NOX::DeepCopy)
    *vectorPtr_ = *(source.vectorPtr_);
}

Vector::~Vector()
{
  if (doDelete_)
    delete vectorPtr_;
}

NOX::size_type Vector::length() const
{
  return vectorPtr_->globalLength();
}

NOX::Abstract::Vector& Vector::init(double value)
{
  vectorPtr_->putScalar(value);
  return *this;
}

NOX::Abstract::Vector& Vector::abs(const NOX::Abstract::Vector& source)
{
  return abs(dynamic_cast<const Vector&>(source));
}

NOX::Abstract::Vector& Vector::abs(const Vector& source)
{
  vectorPtr_->absValue(source.getNativeVectorRef_());
  return *this;
}

NOX::Abstract::Vector& Vector::operator=(const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

NOX::Abstract::Vector& Vector::operator=(const Vector& source)
{
  *vectorPtr_ = *(source.vectorPtr_);
  return *this;
}

NOX::Abstract::Vector& Vector::reciprocal(const NOX::Abstract::Vector& source)
{
  return reciprocal(dynamic_cast<const Vector&>(source));
}

NOX::Abstract::Vector& Vector::reciprocal(const Vector& source)
{
  vectorPtr_->reciprocal(source.getNativeVectorRef_());
  return *this;
}

NOX::Abstract::Vector& Vector::scale(double gamma)
{
  vectorPtr_->scale(gamma);
  return *this;
}

NOX::Abstract::Vector& Vector::scale(const NOX::Abstract::Vector& y)
{
  return scale(dynamic_cast<const Vector&>(y));
}

NOX::Abstract::Vector& Vector::scale(const Vector& y)
{
  vectorPtr_->multiply(y.getNativeVectorRef_());
  return *this;
}

NOX::Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a,
				      double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}


NOX::Abstract::Vector& Vector::update(double alpha, const Vector& a,
				      double gamma)
{
  vectorPtr_->update(alpha, a.getNativeVectorRef_(), gamma);
  return *this;
}

NOX::Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a,
				      double beta, const NOX::Abstract::Vector& b,
				      double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a),
		beta, dynamic_cast<const Vector&>(b),
		gamma);
}



NOX::Abstract::Vector& Vector::update(double alpha, const Vector& a,
				      double beta, const Vector& b,
				      double gamma)
{
  vectorPtr_->update(alpha, a.getNativeVectorRef_(),
                     beta, b.getNativeVectorRef_(), gamma);
  return *this;
}

NOX::Abstract::Vector& Vector::random(bool useSeed, int seed) 
{
  vectorPtr_->random();
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector> 
Vector::clone(NOX::CopyType type) const
{
  Teuchos::RCP<Vector> ptr = Teuchos::rcp(new Vector(*this, type));
  return ptr;
}

double Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  double tmp = 0.0;
  switch(type)
  {
  case NOX::Abstract::Vector::TwoNorm:
    vectorPtr_->lpNorm(2, &tmp);
    break;
  case NOX::Abstract::Vector::OneNorm:
    vectorPtr_->lpNorm(1, &tmp);
    break;
  case NOX::Abstract::Vector::MaxNorm:
    vectorPtr_->infNorm(&tmp);
    break;
  default:
    error("Vector::norm - invalid norm type");
  }
  return tmp;
}

double Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double Vector::norm(const Vector& weights) const
{
  error("N_NLS::NOX::Vector::norm with weights is not supported");
  return 0.0;
}

double Vector::innerProduct(const Vector& y) const
{
  return vectorPtr_->dotProduct(y.getNativeVectorRef_());
}

double Vector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const Vector&>(y));
}

void Vector::print(std::ostream &os) const
{
  vectorPtr_->print(os);
}

}}}
