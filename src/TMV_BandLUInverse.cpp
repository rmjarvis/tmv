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


#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Matrix.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif


#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#else
#define TRI_DIV_BLOCKSIZE 64
#endif

namespace tmv {

    template <class T, class T1> 
    void LU_Inverse(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> minv)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(minv.isSquare());
        TMVAssert(minv.colsize() == LUx.colsize());
#ifdef XDEBUG
        LowerTriMatrix<T,UnitDiag> L0(LUx.colsize());
        LU_PackedPL_Unpack(LUx,p,L0.view());
        UpperTriMatrix<T> U0 = BandMatrixViewOf(LUx,0,LUx.nhi());
        Matrix<T> PLU = L0 * U0;
        if (LUx.nlo() > 0) PLU.reversePermuteRows(p);
        Matrix<T> minv2 = PLU.inverse();
#endif

        if (minv.colsize() > 0) {
            if ( !(minv.iscm() || minv.isrm())) {
                Matrix<T,ColMajor> temp(minv.colsize(),minv.colsize());
                LU_Inverse(LUx,p,temp.view());
                minv = temp;
            } else {
                minv.setZero();
                UpperTriMatrixView<T> U = minv.upperTri();
                U = BandMatrixViewOf(LUx,0,LUx.nhi());
                TriInverse(U,LUx.nhi());
                LU_PackedPL_RDivEq(LUx,p,minv);
            }
        }

#ifdef XDEBUG
        TMV_RealType(T) normdiff = Norm(PLU*minv - T(1));
        TMV_RealType(T) kappa = Norm(PLU)*Norm(minv);
        if (normdiff > 0.001*kappa*minv.colsize()) {
            cerr<<"LUInverse:\n";
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"p = ";
            for(ptrdiff_t i=0;i<LUx.colsize();i++) cerr<<p[i]<<" ";
            cerr<<endl;
            cerr<<"PLU = "<<PLU<<endl;
            cerr<<"minv = "<<minv<<endl;
            cerr<<"minv2 = "<<minv2<<endl;
            cerr<<"m*minv = "<<PLU*minv<<endl;
            cerr<<"minv*m = "<<minv*PLU<<endl;
            cerr<<"Norm(m*minv - 1) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
    }

    template <bool unit, class T> 
    static void RecursiveInverse(UpperTriMatrixView<T> U, ptrdiff_t nhi)
    {
        TMVAssert(U.iscm() || U.isrm());
        TMVAssert(unit == U.isunit());

        const ptrdiff_t N = U.size();
        const ptrdiff_t nb = TRI_DIV_BLOCKSIZE;

        if (N == 1) {
            if (!unit) {
                T*const Uptr = U.ptr();
                if (*Uptr == T(0))  {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<T>(U);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(Uptr >= U._first);
                TMVAssert(Uptr < U._last);
#endif
                *Uptr = TMV_RealType(T)(1) / (*Uptr);
            }
        } else {
            ptrdiff_t k = N/2;
            if (k > nb) k = k/nb*nb;

            UpperTriMatrixView<T> U00 = U.subTriMatrix(0,k);
            MatrixView<T> U01 = U.subMatrix(0,k,k,N);
            UpperTriMatrixView<T> U11 = U.subTriMatrix(k,N);

            // U00 U01' + U01 U11' = 0
            // U00 U01' = -U01 U11'
            // U01' = -U00' U01 U11'

            RecursiveInverse<unit>(U00,nhi);
            RecursiveInverse<unit>(U11,nhi);

            if (nhi < N-k)
                U01.colRange(0,nhi) = -U00 * U01.colRange(0,nhi);
            else
                U01 = -U00 * U01;
            U01 *= U11;
        }
    }

    template <class T> 
    static inline void DoInverse(UpperTriMatrixView<T> U, ptrdiff_t nhi)
    {
#ifndef NOTHROW
        try {
#endif
            if (U.isunit()) RecursiveInverse<true>(U,nhi);
            else RecursiveInverse<false>(U,nhi);
#ifndef NOTHROW
        } catch (Singular) {
            throw SingularUpperTriMatrix<T>(U);
        }
#endif
    }

    template <class T> 
    void TriInverse(UpperTriMatrixView<T> U, ptrdiff_t nhi)
    {
#ifdef XDEBUG
        Matrix<T> U0(U);
#endif

        if (U.size() > 0) {
            if (nhi == 0) DiagMatrixViewOf(U.diag()).invertSelf();
            else if (U.iscm() || U.isrm()) DoInverse(U,nhi);
            else {
                UpperTriMatrix<T> temp = U;
                DoInverse(temp.view(),nhi);
                U = temp;
            }
        }
#ifdef XDEBUG
        Matrix<T> eye = U*U0;
        if (Norm(eye-T(1)) > 0.0001*(Norm(U0)+Norm(U))) {
            cerr<<"UpperTriMatrix Inverse:\n";
            cerr<<"U = "<<TMV_Text(U)<<"  "<<U0<<endl;
            cerr<<"Uinv = "<<U<<endl;
            cerr<<"Uinv*U = "<<U*U0<<endl;
            cerr<<"U*Uinv = "<<U0*U<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandLUInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


