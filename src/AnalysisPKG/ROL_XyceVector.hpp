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


#ifndef ROL_XYCEVECTOR_H
#define ROL_XYCEVECTOR_H

// Xyce includes
#include "Xyce_config.h"
#include "N_UTL_fwd.h"
#include "N_ERH_ErrorMgr.h"
#include "N_LAS_SystemHelpers.h"
#include "N_LAS_MultiVector.h"
#include "N_LAS_Vector.h"
#include "N_PDS_Comm.h"
#include "N_PDS_ParMap.h"
#include "N_UTL_FeatureTest.h"

/*
// Epetra includes
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Comm.h"
*/

#include "ROL_Vector.hpp"

#include <cmath>
#include <fstream>

/** \class Xyce::Linear::XyceVector
    \brief Implements the ROL::Vector interface for a Xyce::Linear::MultiVector.
*/

namespace Xyce  {
namespace Linear {

template <class Real>
class ROL_XyceVector : public ::ROL::Vector<Real> {
private:
  // a vector of pointers to Xyce vectors
  //Teuchos::RCP< std::vector< MultiVector *> >  xyce_multi_vec_;
  Teuchos::RCP< std::vector<Teuchos::RCP<Vector> > >  xyce_multi_vec_;
  int size_;

public:
  virtual ~ROL_XyceVector() {}

  // constructor
  ROL_XyceVector( const int size, const MultiVector & xyce_multi_vec ){
    size_ = size;
    xyce_multi_vec_ = Teuchos::rcp( new std::vector< Teuchos::RCP<Vector> >(size_));
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i] = Teuchos::rcp( createVector( *xyce_multi_vec.pmap(), *xyce_multi_vec.omap() ));
    }
  }

  // copy constructor
  ROL_XyceVector(const Teuchos::RCP<std::vector<Teuchos::RCP<Vector> > > & xyce_multi_vec) : xyce_multi_vec_(xyce_multi_vec) {
    size_ = xyce_multi_vec_->size();
  }

  // copy constructor accepting vectors of raw pointers
  ROL_XyceVector(std::vector<Vector *> & xyce_multi_vec){
    size_ = xyce_multi_vec.size();
    xyce_multi_vec_ = Teuchos::rcp( new std::vector< Teuchos::RCP<Vector> >(size_));
    //std::vector<Vector *> & xmv = const_cast<std::vector<Vector *> &>(xyce_multi_vec);
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i] = Teuchos::rcpFromRef(*xyce_multi_vec[i]);
    }
  }

  // access functions
  Teuchos::RCP<const std::vector<Teuchos::RCP<Vector> > > getVector() const {
    return this->xyce_multi_vec_;
  }


  Teuchos::RCP<std::vector<Teuchos::RCP<Vector> > > getVector() {
    return this->xyce_multi_vec_;
  }


  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const ::ROL::Vector<Real> &x ) {
    const ROL_XyceVector &ex = Teuchos::dyn_cast<const ROL_XyceVector>(x);

    // assuming x and this are the same size
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->update( (double)1.0, *( (*ex.getVector())[i] ) );
    }
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->scale( (double)alpha );
    }
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const ::ROL::Vector<Real> &x ) const {
    double val=0.0;
    const ROL_XyceVector &ex = Teuchos::dyn_cast<const ROL_XyceVector>(x);
    //xyce_multi_vec_->dotProduct( *ex.getVector() );
    for (int i=0;i<size_;i++){
      val += (*xyce_multi_vec_)[i]->dotProduct( *(*ex.getVector())[i] );
    }
    return (Real)val;
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    std::vector<double> vals(size_,0.0);
    double res = 0.0;
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->lpNorm(2,&vals[i]);      
      res += vals[i]*vals[i];
    }
    return (Real) sqrt(res);
  } 

  /** \brief Clone to make a new (uninitialized) vector.
  */
  Teuchos::RCP< ::ROL::Vector<Real> > clone() const{
    // return Teuchos::rcp(new ROL_XyceVector( Teuchos::rcp( MultiVector( *(xyce_multi_vec_->pmap()), xyce_multi_vec_->numVectors() ) )) );
    //return Teuchos::rcp(new ROL_XyceVector( Teuchos::rcp( new std::vector< MultiVector *>(MultiVector(*((*xyce_multi_vec_)[0]) ) ) ) ) );

    Teuchos::RCP< std::vector<Teuchos::RCP<Vector> > > x = Teuchos::rcp( new std::vector<Teuchos::RCP<Vector> >(size_));

    //Teuchos::RCP< std::vector< MultiVector * > > x = Teuchos::rcp( new std::vector< MultiVector * >(size_, &*Teuchos::rcp(new EpetraVector( *((*xyce_multi_vec_)[0]->pmap()), *((*xyce_multi_vec_)[0]->omap()) ) )  ));

    //Teuchos::RCP< std::vector< MultiVector * > > x = Teuchos::rcp( new std::vector< MultiVector * >(size_, dynamic_cast< MultiVector* >(&*newMV(0)) ));
    
    //Teuchos::RCP< Vector > temp; 
    for (int i=0;i<size_;i++){
      (*x)[i] = Teuchos::rcp( (*xyce_multi_vec_)[i]->cloneVector() );

      // Using copy constructor (segfault)
      //*((*x)[i]) = MultiVector(*((*xyce_multi_vec_)[i])); 
    }
    //if (is_null(x)) std::cout << "x is null!" << std::endl;
    return Teuchos::rcp( new ROL_XyceVector( x ));
    
  }

  void randomize() {
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->random();
    }
  }

  void print(std::ostream & outStream=std::cout) {
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->print(outStream);
    }
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void axpy( const Real alpha, const ::ROL::Vector<Real> &x ) {
    const ROL_XyceVector &ex = Teuchos::dyn_cast<const ROL_XyceVector>(x);
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->update( (double)alpha, *( (*ex.getVector())[i] ) );
    }
  }

  /**  \brief Set to zero vector.
  */
  virtual void zero() {
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->putScalar( 0.0 );
    }
  }

  virtual void putScalar(Real alpha) {
    for (int i=0;i<size_;i++){
      (*xyce_multi_vec_)[i]->putScalar( (double)alpha );
    }
  }

  // /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  // */
  // virtual void set( const Vector<Real> &x ) {
  //   const EpetraMultiVector &ex = Teuchos::dyn_cast<const EpetraMultiVector>(x);
  //   epetra_vec_->Scale(1.0,*ex.getVector());
  // }

  // Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  //   Teuchos::RCP<EpetraMultiVector> e = Teuchos::rcp( new EpetraMultiVector( Teuchos::rcp(new Epetra_MultiVector(epetra_vec_->Map(),epetra_vec_->NumVectors(),true)) ));
  //   const Epetra_BlockMap & domainMap = e->getVector()->Map();

  //   Epetra_Map linearMap(domainMap.NumGlobalElements(), domainMap.NumMyElements(), 0, domainMap.Comm());
  //   int lid = linearMap.LID(i);
  //   if(lid >=0)
  //     (*e->getVector())[0][lid]= 1.0;

  //   return e;

  //   /*


  //   // Build IntVector of GIDs on all processors.
  //   const Epetra_Comm & comm = domainMap.Comm();
  //   int numMyElements = domainMap.NumMyElements();
  //   Epetra_BlockMap allGidsMap(-1, numMyElements, 1, 0, comm);
  //   Epetra_IntVector allGids(allGidsMap);
  //   for (int j=0; j<numMyElements; j++) {allGids[j] = domainMap.GID(j);}

  //   // Import my GIDs into an all-inclusive map. 
  //   int numGlobalElements = domainMap.NumGlobalElements();
  //   Epetra_LocalMap allGidsOnRootMap(numGlobalElements, 0, comm);
  //   Epetra_Import importer(allGidsOnRootMap, allGidsMap);
  //   Epetra_IntVector allGidsOnRoot(allGidsOnRootMap);
  //   allGidsOnRoot.Import(allGids, importer, Insert);
  //   Epetra_Map rootDomainMap(-1, allGidsOnRoot.MyLength(), allGidsOnRoot.Values(), domainMap.IndexBase(), comm);

  //   for (int j = 0; j < this->dimension(); j++) {
  //     // Put 1's in slots
  //     int curGlobalCol = rootDomainMap.GID(i); // Should return same value on all processors
  //     if (domainMap.MyGID(curGlobalCol)){
  //       int curCol = domainMap.LID(curGlobalCol);
  //       (*e->getVector())[0][curCol]= 1.0;
  //     }
  //   }

  //   return e;

  //   */
  // }

  // int dimension() const {return epetra_vec_->GlobalLength();}


}; // class XyceVector

} // namespace Linear
} // namespace Xyce


#endif

