///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_DiagMatrixArithFunc.h"
#include "TMV_DiagMatrix.h"
#include "TMV_Vector.h"
#include "TMV_VectorArithFunc.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_Matrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> ConstVectorView<T> DiagMatrixComposite<T>::cdiag() const
  {
    if (!inst.get()) inst.reset(new DiagMatrix<T>(*this));
    return inst->diag();
  }

  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x 
    // yi (+)= alpha * Ai * xi
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
#ifdef XDEBUG
    //cout<<"MultMV: \n";
    //cout<<"alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"x = "<<Type(x)<<"  "<<x<<endl;
    //cout<<"y = "<<Type(y)<<"  "<<y<<endl;
    Vector<T> y0 = y;
    Vector<Tx> x0 = x;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = alpha*A0*x0;
    if (add) y2 += y0;
#endif

    if (y.size() > 0) {
      if (alpha == T(0)) {
	if (!add) y.Zero();
      } 
      else if (!add) {
	if (SameStorage(A.diag(),y)) {
	  if (y.SameAs(A.diag())) ElementProd(alpha,x,y);
	  else if (SameStorage(x,y)) {
	    if (y.SameAs(x)) ElementProd(alpha,A.diag(),y);
	    else {
	      Vector<T> yy = x;
	      ElementProd(alpha,A.diag(),yy.View());
	      y = yy;
	    }
	  } 
	  else {
	    y = A.diag();
	    ElementProd(alpha,x,y);
	  }
	}
	else {
	  y = x;
	  ElementProd(alpha,A.diag(),y);
	}
      }
      else AddElementProd(alpha,A.diag(),x,y);
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultDV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


