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



#include "tmv/TMV_SymHouseholder.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

namespace tmv {

    template <class T1, class T2> 
    void HouseholderLRMult(
        const GenVector<T1>& v, T1 beta, SymMatrixView<T2> m)
    {
        // The input vector, v, is taken to be the vector for a  
        // Householder matrix, H.  This routine takes m <- H m Ht
        // if m is Hermitian or H m HT if m is symmetric.
        TMVAssert(m.size() > 0);
        TMVAssert(v.size() == m.size()-1);

        // If m is Hermitian:
        //
        // H m Ht = (I - beta v vt) m (I - beta* v vt)
        //        = m - beta v vt m - beta* m v vt + |beta|^2 v vt m v vt
        //        = m - beta v (mv)t - beta* (mv) vt + (|beta|^2 vt (mv)) v vt
        //
        // If m is symmetric:
        //
        // H m HT = (I - beta v vt) m (I - beta v* vT)
        //        = m - beta v vt m - beta m v* vT + beta^2 v vt m v* vT
        //        = m - beta v (mv*)T - beta (mv*) vT + (beta^2 vt (mv*)) v vT
        //
        ptrdiff_t N = m.size();
        if (N > 0 && beta != T1(0)) {
            // Normally, I take the unit first element of v to be implicit, 
            // but this calculation is more complicated, so for now I just 
            // copy the vector to a temporary Vector (vv).  
            // Someday maybe I'll change this to not need this temporary.
            Vector<T1> vv(N);
            vv(0) = T1(1);
            vv.subVector(1,N) = v;

            T2 betasqvtmv;
            Vector<T2> mv(N);
            if (m.isherm()) {
                mv = m*vv;
                betasqvtmv = TMV_NORM(beta)*vv.conjugate()*mv;
            } else {
                mv = m*vv.conjugate();
                betasqvtmv = beta*beta*vv.conjugate()*mv;
            }
            Rank2Update<true>(T2(-beta),vv,mv,m);
            if (m.issym()) {
                m += betasqvtmv * (vv^vv);
            } else {
                // imag part is 0 - make it exact.
                m += TMV_REAL(betasqvtmv) * (vv^vv.conjugate());
            }
        }
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymHouseholder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


