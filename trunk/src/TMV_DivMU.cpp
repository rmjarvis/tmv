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

//#define PRINTALGO_DIVU

#include "TMV_Blas.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_MultMM.h"

namespace tmv {

    template <class M1, class M2>
    static void NonBlasLDivEq(M1& m1, const M2& m2)
    {
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m2.isrm() || m2.iscm());
        if (m1.iscm()) {
            typename M1::cmview_type m1cm = m1.cmView();
            if (m2.iscm())
                InlineLDivEq(m1cm,m2.cmView());
            else
                InlineLDivEq(m1cm,m2.rmView());
        } else {
            typename M1::rmview_type m1rm = m1.rmView();
            if (m2.iscm())
                InlineLDivEq(m1rm,m2.cmView());
            else
                InlineLDivEq(m1rm,m2.rmView());
        }
    }

#ifdef BLAS
    template <class T1, class M2, class T2> 
    static inline void BlasLDivEq(MatrixView<T1> A, const M2& B, T2)
    { NonBlasLDivEq(A,B); }
#ifdef TMV_INST_DOUBLE
    template <class M2> 
    static void BlasLDivEq(MatrixView<double> A, const M2& B, double )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        double alpha(1);
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(dtrsm) (
            BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
            B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            B.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
            BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M2>
    static void BlasLDivEq(
        MatrixView<std::complex<double> > A, const M2& B, std::complex<double> )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        std::complex<double> alpha(1);
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        if (A.iscm()==B.iscm() && B.isconj()) {
            A.conjugateSelf();
            BLASNAME(ztrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
            A.conjugateSelf();
        } else {
            BLASNAME(ztrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2>
    static void BlasLDivEq(
        MatrixView<std::complex<double> > A, const M2& B, double )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        if (A.isrm()) {
            int m=2*A.rowsize();
            int n=A.colsize();
            double alpha(1);
            int lda = 2*A.stepi();
            int ldb = B.isrm()?B.stepi():B.stepj();

            BLASNAME(dtrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP((double*)A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<double,ColMajor> A1 = A.realPart();
            BlasLDivEq(A1.xView(),B,double(0));
            A.realPart() = A1;
            A1 = A.imagPart();
            BlasLDivEq(A1.xView(),B,double(0));
            A.imagPart() = A1;
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    template <class M2> 
    static void BlasLDivEq(MatrixView<float> A, const M2& B, float )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        float alpha(1);
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();

        BLASNAME(strsm) (
            BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
            B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
            A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
            B.isunit()?BLASCH_U:BLASCH_NU, 
            BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
            BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
    }
    template <class M2>
    static void BlasLDivEq(
        MatrixView<std::complex<float> > A, const M2& B, std::complex<float> )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        int m=A.iscm()?A.colsize():A.rowsize();
        int n=A.iscm()?A.rowsize():A.colsize();
        std::complex<float> alpha(1);
        int lda = A.isrm()?A.stepi():A.stepj();
        int ldb = B.isrm()?B.stepi():B.stepj();
        if (A.iscm()==B.iscm() && B.isconj()) {
            A.conjugateSelf();
            BLASNAME(ctrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, BLASCH_NT,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
            A.conjugateSelf();
        } else {
            BLASNAME(ctrsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASP(&alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        }
    }
    template <class M2>
    static void BlasLDivEq(
        MatrixView<std::complex<float> > A, const M2& B, float )
    {
        TMVAssert(B.size() == A.colsize());
        TMVAssert(A.colsize()>0);
        TMVAssert(A.rowsize()>0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(B.isrm() || B.iscm());

        if (A.isrm()) {
            int m=2*A.rowsize();
            int n=A.colsize();
            float alpha(1);
            int lda = 2*A.stepi();
            int ldb = B.isrm()?B.stepi():B.stepj();

            BLASNAME(strsm) (
                BLASCM A.iscm()?BLASCH_L:BLASCH_R, 
                B.iscm()==M2::mupper?BLASCH_UP:BLASCH_LO, 
                A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
                B.isunit()?BLASCH_U:BLASCH_NU, 
                BLASV(m),BLASV(n),BLASV(alpha),BLASP(B.cptr()),BLASV(ldb),
                BLASP((float*)A.ptr()),BLASV(lda) BLAS1 BLAS1 BLAS1 BLAS1);
        } else {
            Matrix<float,ColMajor> A1 = A.realPart();
            BlasLDivEq(A1.xView(),B,float(0));
            A.realPart() = A1;
            A1 = A.imagPart();
            BlasLDivEq(A1.xView(),B,float(0));
            A.imagPart() = A1;
        }
    }
#endif
#endif // BLAS

    template <class T1, class M2>
    static void DoLDivEq(MatrixView<T1> m1, const M2& m2)
    {
#ifdef BLAS
        const typename M2::value_type t2(0);
        if ((m1.isrm() && m1.stepi()>0) || (m1.iscm() && m1.stepj()>0) ) {
            if ((m2.isrm() && m2.stepi()>0) || (m2.iscm() && m2.stepj()>0)) {
                BlasLDivEq(m1,m2,t2);
            } else {
                if (m2.isunit()) 
                    BlasLDivEq(
                        m1,m2.copy().viewAsUnitDiag().constView().xdView(),t2);
                else 
                    BlasLDivEq(m1,m2.copy().constView().xdView(),t2);
            }
        } else {
            Matrix<T1,ColMajor> m1c(m1);
            DoLDivEq(m1c.xView(),m2);
            InstCopy(m1c.constView().xView(),m1);
        }
#else
        if (m1.iscm() || m1.isrm()) {
            if (m2.iscm() || m2.isrm()) {
                NonBlasLDivEq(m1,m2);
            } else {
                if (m2.isunit()) 
                    NonBlasLDivEq(
                        m1,m2.copy().viewAsUnitDiag().constView().xdView());
                else 
                    NonBlasLDivEq(m1,m2.copy().constView().xdView());
            }
        } else {
            Matrix<T1,ColMajor> m1c(m1);
            DoLDivEq(m1c.xView(),m2);
            InstCopy(m1c.constView().xView(),m1);
        }
#endif
    }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        MatrixView<T1> m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoLDivEq(m1,m2); }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        MatrixView<T1> m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoLDivEq(m1,m2); }

#define InstFile "TMV_DivMU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


