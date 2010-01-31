///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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


//#define XDEBUG


#include "tmv/TMV_DiagMatrixArithFunc.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArithFunc.h"

#ifdef XDEBUG
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    ConstVectorView<T> DiagMatrixComposite<T>::cdiag() const
    {
        if (!inst.get()) inst.reset(new DiagMatrix<T>(*this));
        return inst->diag();
    }

    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
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
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
        //cout<<"x = "<<TMV_Text(x)<<"  "<<x<<endl;
        //cout<<"y = "<<TMV_Text(y)<<"  "<<y<<endl;
        Vector<T> y0 = y;
        Vector<Tx> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = alpha*A0*x0;
        if (add) y2 += y0;
#endif

        if (y.size() > 0) {
            if (alpha == T(0)) {
                if (!add) y.setZero();
            } 
            else if (!add) {
                if (SameStorage(A.diag(),y)) {
                    if (y.isSameAs(A.diag())) ElementProd(alpha,x,y);
                    else if (SameStorage(x,y)) {
                        if (y.isSameAs(x)) ElementProd(alpha,A.diag(),y);
                        else {
                            Vector<T> yy = x;
                            ElementProd(alpha,A.diag(),yy.view());
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
        if (Norm(y-y2) > 0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                                (add?Norm(y0):TMV_RealType(T)(0)))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"-> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultDV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


