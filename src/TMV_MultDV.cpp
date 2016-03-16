///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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
        VectorView<T> y)
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

        ElemMultVV<add>(alpha,A.diag(),x,y); 

#ifdef XDEBUG
        if (!(Norm(y-y2) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
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


