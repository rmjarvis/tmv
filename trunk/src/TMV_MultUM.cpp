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
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ProdXM.h"

namespace tmv {

    template <class M1, class M3>
    void DoMultEqUM(const M1& m1, M3& m3)
    {
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m3.isrm() || m3.iscm());
        typedef typename M3::real_type RT;
        Scaling<1,RT> one;
        if (m3.iscm()) {
            typename M3::cmview_type m3cm = m3.cmView();
            if (m1.iscm())
                InlineMultMM<false>(one,m1.cmView(),m3cm,m3cm);
            else
                InlineMultMM<false>(one,m1.rmView(),m3cm,m3cm);
        } else {
            typename M3::rmview_type m3rm = m3.rmView();
            if (m1.iscm())
                InlineMultMM<false>(one,m1.cmView(),m3rm,m3rm);
            else
                InlineMultMM<false>(one,m1.rmView(),m3rm,m3rm);
        }
    }

#ifdef BLAS
#ifdef INST_DOUBLE
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<double>& A, 
        const MatrixView<double>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const double alpha(1);

        BLASNAME(dtrmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<double>& A, 
        const MatrixView<double>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const double alpha(1);

        BLASNAME(dtrmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_LO:BLASCH_UP, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<std::complex<double> >& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const std::complex<double> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<std::complex<double> >& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const std::complex<double> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ztrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<double>& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = 2*B.stepi();
            double alpha(1);
            BLASNAME(dtrmm) (
                BLASCM BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.isrm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((double*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<double,ColMajor> B1 = B.realPart();
            DoMultEqMM(A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            DoMultEqMM(A,B1.view());
            B.imagPart() = B1;
        }
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<double>& A,
        const MatrixView<std::complex<double> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m=2*B.rowsize();
            int n=B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.stepi();
            double alpha(1);
            BLASNAME(dtrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((double*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<double,ColMajor> B1 = B.realPart();
            DoMultEqMM(1.,A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            DoMultEqMM(1.,A,B1.view());
            B.imagPart() = B1;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<float>& A, 
        const MatrixView<float>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        float alpha(1);

        BLASNAME(strmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<float>& A, 
        const MatrixView<float>& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        float alpha(1);

        BLASNAME(strmm) (
            BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
            A.iscm()?BLASCH_LO:BLASCH_UP, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
            BLASV(alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
            BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<std::complex<float> >& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const std::complex<float> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<std::complex<float> >& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        int m=B.iscm()?B.colsize():B.rowsize();
        int n=B.iscm()?B.rowsize():B.colsize();
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        const std::complex<float> alpha(1);

        if (A.iscm()==B.iscm() && A.isconj()) {
            B.conjugateSelf();
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
            B.conjugateSelf();
        } else {
            BLASNAME(ctrmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASP(&alpha),BLASP(A.cptr()),BLASV(lda), BLASP(B.ptr()), 
                BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void DoMultEqMM(
        const GenUpperTriMatrix<float>& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m = 2*B.rowsize();
            int n = B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = 2*B.stepi();
            float alpha(1);
            BLASNAME(strmm) (
                BLASCM BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.isrm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda), 
                BLASP((float*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<float,ColMajor> B1 = B.realPart();
            DoMultEqMM(1.F,A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            DoMultEqMM(1.F,A,B1.view());
            B.imagPart() = B1;
        }
    }
    template <> 
    void DoMultEqMM(
        const GenLowerTriMatrix<float>& A,
        const MatrixView<std::complex<float> >& B)
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(B.rowsize() > 0);
        TMVAssert(B.colsize() > 0);
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) {
            int m=2*B.rowsize();
            int n=B.colsize();
            int lda = A.isrm()?A.stepi():A.stepj();
            int ldb = B.stepi();
            float alpha(1);
            BLASNAME(strmm) (
                BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
                A.iscm()?BLASCH_LO:BLASCH_UP, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU, BLASV(m),BLASV(n),
                BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
                BLASP((float*)B.ptr()), BLASV(ldb)
                BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<float,ColMajor> B1 = B.realPart();
            DoMultEqMM(1.F,A,B1.view());
            B.realPart() = B1;
            B1 = B.imagPart();
            DoMultEqMM(1.F,A,B1.view());
            B.imagPart() = B1;
        }
    }
#endif
#endif // BLAS

    template <class T, class M1, class M2>
    void GenInstMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        if (!(m1.isrm() || m1.iscm())) 
            GenInstMultMM(x,m1.copy().constView().xdView(),m2,m3);
        else if (!(m3.iscm() || m3.isrm())) {
            Matrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
            InstMultXM(x,m2,m3c.xView());
            DoMultEqUM(m1,m3c);
            InstCopy(m3c.constView().xView(),m3);
        } else {
            InstMultXM(x,m2,m3);
            DoMultEqUM(m1,m3);
        }
    }

    template <class T, class M1, class M2>
    void GenInstAddMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        if (!(m1.isrm() || m1.iscm()))
            GenInstAddMultMM(x,m1.copy().constView().xdView(),m2,m3);
        else {
            Matrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
            InstMultXM(x,m2,m3c.xView());
            DoMultEqUM(m1,m3c);
            InstAddMultXM(T(1),m3c.xView().constView(),m3);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }


#define InstFile "TMV_MultUM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


