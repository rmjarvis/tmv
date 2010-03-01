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

#include "TMV_Blas.h"
#include "tmv/TMV_MultUV.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_ProdXV.h"

// The CBLAS trick of using RowMajor with ConjTrans when we have a 
// case of A.conjugate() * x doesn't seem to be working with MKL 10.2.2.
// I haven't been able to figure out why.  (e.g. Is it a bug in the MKL
// code, or am I doing something wrong?)  So for now, just disable it.
#ifdef CBLAS
#undef CBLAS
#endif


namespace tmv {

    // To save on compile size, I only compile the MultEq function:
    // v *= U, v *= L
    // This requires a temporary for add=true, but those are actually 
    // fairly rare, and it's not much of a hit on the execution time.
    // If people really want it without a temporary, they should 
    // call InlineMultMV directly, rather than use the Inst version.
    
    template <class T, class M1>
    static void DoMultEqUV(const M1& m1, VectorView<T> v3)
    {
        if (v3.step() == 1) {
            const Scaling<1,typename Traits<T>::real_type> one;
            ConstVectorView<T,1> v2u = v3.unitView();
            VectorView<T,1> v3u = v3.unitView();
            if (m1.isrm()) InlineMultMV<false>(one,m1.rmView(),v2u,v3u);
            else InlineMultMV<false>(one,m1.cmView(),v2u,v3u);
        } else {
            Vector<T> v3c = v3;
            DoMultEqUV(m1,v3c.xView());
            InstCopy(v3c.constView().xView(),v3);
        }
    }

#ifdef BLAS
#ifdef INST_DOUBLE
    template <DiagType D>
    static void DoMultEqUV( 
        const GenUpperTriMatrix<double,D>& A, const VectorView<double>& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        double* xp = x.ptr();
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV( 
        const GenLowerTriMatrix<double,D>& A, const VectorView<double>& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        double* xp = x.ptr();
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV(
        const GenUpperTriMatrix<std::complex<double>,D>& A,
        const VectorView<std::complex<double> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<double>* xp = x.ptr();
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztrmv) (
                BLASRM BLASCH_LO, BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ztrmv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            BLASNAME(ztrmv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <DiagType D>
    static void DoMultEqUV(
        const GenLowerTriMatrix<std::complex<double>,D>& A,
        const VectorView<std::complex<double> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<double>* xp = x.ptr();
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztrmv) (
                BLASRM BLASCH_UP, BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ztrmv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            BLASNAME(ztrmv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <DiagType D>
    static void DoMultEqUV( 
        const GenUpperTriMatrix<double,D>& A,
        const VectorView<std::complex<double> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        double* xp = (double*) x.ptr();
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV( 
        const GenLowerTriMatrix<double,D>& A,
        const VectorView<std::complex<double> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        double* xp = (double*) x.ptr();
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif
#endif
#ifdef INST_FLOAT
    template <DiagType D>
    static void DoMultEqUV( 
        const GenUpperTriMatrix<float,D>& A, const VectorView<float>& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        float* xp = x.ptr();
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV( 
        const GenLowerTriMatrix<float,D>& A, const VectorView<float>& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        float* xp = x.ptr();
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV(
        const GenUpperTriMatrix<std::complex<float>,D>& A,
        const VectorView<std::complex<float> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<float>* xp = x.ptr();
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ctrmv) (
                BLASRM BLASCH_LO, BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ctrmv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            BLASNAME(ctrmv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
    template <DiagType D>
    static void DoMultEqUV(
        const GenLowerTriMatrix<std::complex<float>,D>& A,
        const VectorView<std::complex<float> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=x.step();
        std::complex<float>* xp = x.ptr();
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ctrmv) (
                BLASRM BLASCH_UP, BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
            x.conjugateSelf();
            BLASNAME(ctrmv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
            x.conjugateSelf();
#endif
        } else {
            BLASNAME(ctrmv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),
                BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        }
    }
#ifdef TMV_INST_MIX
    template <DiagType D>
    static void DoMultEqUV( 
        const GenUpperTriMatrix<float,D>& A,
        const VectorView<std::complex<float> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        float* xp = (float*) x.ptr();
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
    template <DiagType D>
    static void DoMultEqUV( 
        const GenLowerTriMatrix<float,D>& A,
        const VectorView<std::complex<float> >& x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(x.ct() == NonConj);

        int n=A.size();
        int lda=A.isrm()?A.stepi():A.stepj();
        int xs=2*x.step();
        float* xp = (float*) x.ptr();
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs) BLAS1 BLAS1 BLAS1);
        BLASNAME(strmv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp+1),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
#endif
#endif // FLOAT
#endif // BLAS

    template <class T, class M1, class V2>
    void GenInstMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T> v3)
    {
        if (m1.isrm() || m1.iscm()) {
            InstMultXV(x,v2,v3);
            DoMultEqUV(m1,v3);
        } else {
            GenInstMultMV(x,m1.copy().constView().xView(),v2,v3);
        }
    }

    template <class T, class M1, class V2>
    void GenInstAddMultMV(
        const T x, const M1& m1, const V2& v2, VectorView<T> v3)
    {
        if (m1.isrm() || m1.iscm()) {
            Vector<T> v2c = v2;
            DoMultEqUV(m1,v2c.xView());
            InstAddMultXV(x,v2c.constView().xView(),v3);
        } else {
            GenInstAddMultMV(x,m1.copy().constView().xView(),v2,v3);
        }
    }

    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { GenInstMultMV(x,m1,v2,v3); }
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { GenInstAddMultMV(x,m1,v2,v3); }

    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { GenInstMultMV(x,m1,v2,v3); }
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3)
    { GenInstAddMultMV(x,m1,v2,v3); }

#define InstFile "TMV_MultUV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


