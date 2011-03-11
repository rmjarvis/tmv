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

//#define PRINTALGO_UV

#include "TMV_Blas.h"
#include "tmv/TMV_MultUV.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_ProdXV.h"

namespace tmv {

    // To save on compile size, I only compile the MultEq function:
    // v *= U, v *= L
    // This requires a temporary for add=true, but those are actually 
    // fairly rare, and it's not much of a hit on the execution time.
    // If people really want it without a temporary, they should 
    // call InlineMultMV directly, rather than use the Inst version.
    //
    // In fact, I wonder if it would be worth not even having an Inst
    // version for the Add option.  Let that be done inline all the time.
    // I'll have the think about that...
    
    template <class M1, class T2>
    static void NonBlasMultEq(const M1& m1, VectorView<T2> v2)
    {
        TMVAssert(m1.iscm() || m1.isrm());
        if (v2.step() == 1) {
            const Scaling<1,typename Traits<T2>::real_type> one;
            VectorView<T2,1> v2u = v2.unitView();
            if (m1.iscm()) 
                InlineMultMV<false>(one,m1.cmView(),v2u.constView(),v2u);
            else 
                InlineMultMV<false>(one,m1.rmView(),v2u.constView(),v2u);
        } else {
            Vector<T2> v2c = v2;
            DoMultEq(m1,v2c.xView());
            InstCopy(v2c.constView().xView(),v2);
        }
    }

#ifdef BLAS
    template <class M1, class T2, class T1> 
    static inline void BlasMultEq(const M1& A, VectorView<T2> x, T1)
    { NonBlasMultEq(A,x); }
#ifdef TMV_INST_DOUBLE
    template <class M1>
    static void BlasMultEq(const M1& A, VectorView<double> x, double)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        double* xp = x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        BLASNAME(dtrmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, VectorView<std::complex<double> > x, std::complex<double>)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<double>* xp = x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        if (A.iscm() && A.isconj()) {
            x.conjugateSelf();
            BLASNAME(ztrmv) (
                BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
        } else {
            BLASNAME(ztrmv) (
                BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M1>
    static void BlasMultEq( 
        const M1& A, VectorView<std::complex<double> > x, double)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        double* xp = (double*) x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        BLASNAME(dtrmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif
#ifdef TMV_INST_FLOAT
    template <class M1>
    static void BlasMultEq(const M1& A, VectorView<float> x, float)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        float* xp = x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        BLASNAME(strmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <class M1>
    static void BlasMultEq(
        const M1& A, VectorView<std::complex<float> > x, std::complex<float>)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<float>* xp = x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        if (A.iscm() && A.isconj()) {
            x.conjugateSelf();
            BLASNAME(ctrmv) (
                BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
        } else {
            BLASNAME(ctrmv) (
                BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M1>
    static void BlasMultEq( 
        const M1& A, VectorView<std::complex<float> > x, float)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        float* xp = (float*) x.ptr();
        if (xs < 0) xp += (n-1)*xs;
        BLASNAME(strmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(strmv) (
            BLASCM A.iscm()==M1::mupper?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif
#endif // BLAS

    template <class M1, class T2>
    static void DoMultEq(const M1& m1, VectorView<T2> v2)
    {
#ifdef BLAS
        const typename M1::value_type t1(0);
        if ( (m1.isrm() && m1.stepi()>0) || (m1.iscm() && m1.stepj()>0) ) {
            BlasMultEq(m1,v2,t1);
        } else {
            if (m1.isunit()) {
                BlasMultEq(
                    m1.copy().viewAsUnitDiag().constView().xdView(),v2,t1);
            } else {
                BlasMultEq(m1.copy().constView().xdView(),v2,t1);
            }
        }
#else
        if (m1.isrm() || m1.iscm()) {
            NonBlasMultEq(m1,v2);
        } else {
            if (m1.isunit()) {
                NonBlasMultEq(
                    m1.copy().viewAsUnitDiag().constView().xdView(),v2);
            } else {
                NonBlasMultEq(m1.copy().constView().xdView(),v2);
            }
        }
#endif
    }

    template <class T, class M1, class V2>
    void DoInstMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T> v3)
    {
        InstMultXV(x,v2,v3);
        DoMultEq(m1,v3);
    }

    template <class T, class M1, class V2>
    void DoInstAddMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T> v3)
    {
        Vector<T> v2c = v2;
        DoMultEq(m1,v2c.xView());
        InstAddMultXV(x,v2c.constView().xView(),v3);
    }

    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { DoInstMultMV(x,m1,v2,v3); }
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { DoInstAddMultMV(x,m1,v2,v3); }

    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { DoInstMultMV(x,m1,v2,v3); }
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { DoInstAddMultMV(x,m1,v2,v3); }

#define InstFile "TMV_MultUV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


