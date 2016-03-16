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


#include "TMV_Blas.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <bool ua, class T, class Ta> 
    static void RowUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(B.isrm());
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(A.nlo() == 0);
        TMVAssert(B.ct() == NonConj);

        ptrdiff_t N = B.colsize();

        ptrdiff_t k = A.nhi();
        if (ua) {
            for(ptrdiff_t i=N-1; i>=0; --i) {
                B.row(i) -= A.row(i,i+1,N) * B.rowRange(i+1,N);
                if (k > 0) --k; else --N;
            }
        } else {
            const ptrdiff_t ds = A.diagstep();
            const Ta* Aii = A.cptr() + (N-1)*ds;
            for(ptrdiff_t i=N-1; i>=0; --i,Aii-=ds) {
                B.row(i) -= A.row(i,i+1,N) * B.rowRange(i+1,N);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
                B.row(i) /= (A.isconj() ? TMV_CONJ(*Aii) : *Aii);
                if (k > 0) --k; else --N;
            }
        } 
    }

    template <bool ua, class T, class Ta> 
    static void ColUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(B.isrm());
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(A.nlo() == 0);
        TMVAssert(A.colsize()>A.nhi());
        TMVAssert(B.ct() == NonConj);

        ptrdiff_t N = A.colsize();

        ptrdiff_t i1 = N-1-A.nhi();
        if (ua) {
            for(ptrdiff_t j=N-1; j>0; --j) {
                B.rowRange(i1,j) -= A.col(j,i1,j) ^ B.row(j);
                if (i1 > 0) --i1;
            }
        } else {
            const ptrdiff_t ds = A.diagstep();
            const Ta* Ajj = A.cptr() + (N-1)*ds;
            for(ptrdiff_t j=N-1; j>=0; --j,Ajj-=ds) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
                B.row(j) /= (A.isconj() ? TMV_CONJ(*Ajj) : *Ajj);
                B.rowRange(i1,j) -= A.col(j,i1,j) ^ B.row(j);
                if (i1 > 0) --i1;
            }
        } 
    }

    template <bool ua, class T, class Ta> 
    static void RowLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(B.isrm());
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(A.nhi() == 0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t N = B.colsize();

        ptrdiff_t i1=0;
        ptrdiff_t k=A.nlo();
        if (ua) {
            for(ptrdiff_t i=0; i<N; ++i) {
                B.row(i) -= A.row(i,i1,i) * B.rowRange(i1,i);
                if (k>0) --k; else ++i1;
            }
        } else {
            const ptrdiff_t ds = A.diagstep();
            const Ta* Aii = A.cptr();
            for(ptrdiff_t i=0; i<N; ++i,Aii+=ds) {
                B.row(i) -= A.row(i,i1,i) * B.rowRange(i1,i);
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
                B.row(i) /= (A.isconj() ? TMV_CONJ(*Aii) : *Aii);
                if (k>0) --k; else ++i1;
            }
        }
    }

    template <bool ua, class T, class Ta> 
    static void ColLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(B.isrm());
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 0);
        TMVAssert(A.nhi() == 0);
        TMVAssert(B.ct() == NonConj);

        const ptrdiff_t N = B.colsize();

        ptrdiff_t i2=A.nlo()+1;
        if (ua) {
            for(ptrdiff_t j=0; j<N; ++j) {
                B.rowRange(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
                if (i2 < N) ++i2;
            }
        } else {
            const ptrdiff_t ds = A.diagstep();
            const Ta* Ajj = A.cptr();
            for(ptrdiff_t j=0; j<N; ++j,Ajj+=ds) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
                B.row(j) /= (A.isconj() ? TMV_CONJ(*Ajj) : *Ajj);
                B.rowRange(j+1,i2) -= A.col(j,j+1,i2) ^ B.row(j);
                if (i2 < N) ++i2;
            }
        }
    }

    template <bool ua, class T, class Ta> 
    static void NonLapTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B)
    {
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 1);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        TMVAssert(B.ct() == NonConj);

        if (B.isrm()) 
            if (A.nlo()==0) 
                if (A.isrm()) RowUpperTriLDivEq<ua>(A,B);
                else if (A.iscm()) ColUpperTriLDivEq<ua>(A,B);
                else RowUpperTriLDivEq<ua>(A,B);
            else 
                if (A.isrm()) RowLowerTriLDivEq<ua>(A,B);
                else if (A.iscm()) ColLowerTriLDivEq<ua>(A,B);
                else RowLowerTriLDivEq<ua>(A,B);
        else {
            const ptrdiff_t N = B.rowsize();
            for(ptrdiff_t j=0;j<N;++j) {
                TriLDivEq(A,B.col(j),ua?UnitDiag:NonUnitDiag);
            }
        }
    }

#ifdef LAP
    template <class T, class Ta> 
    static inline void LapTriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B, DiagType dt)
    { 
        if (dt == UnitDiag) NonLapTriLDivEq<true>(A,B);
        else NonLapTriLDivEq<false>(A,B);
    }
#ifdef INST_DOUBLE
    template <> void LapTriLDivEq(
        const GenBandMatrix<double>& A, MatrixView<double> B, DiagType dt)
    {
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 1);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(B.iscm());
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int nrhs = B.rowsize();
        int aoffset = 
            (A.nlo()==0 && A.iscm()) ? A.nhi() :
            (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
        int lda = A.diagstep();
        int ldb = B.stepj();
        int Lap_info=0;
        LAPNAME(dtbtrs) (
            LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
            A.iscm()?LAPCH_NT:LAPCH_T, dt==UnitDiag?LAPCH_U:LAPCH_NU,
            LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
            LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
        LAP_Results(Lap_info,"dtbtrs");
    }
    template <> void LapTriLDivEq(
        const GenBandMatrix<std::complex<double> >& A,
        MatrixView<std::complex<double> > B, DiagType dt)
    {
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 1);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        TMVAssert(B.ct() == NonConj);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(B.iscm());
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int nrhs = B.rowsize();
        int aoffset = 
            (A.nlo()==0 && A.iscm()) ? A.nhi() :
            (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
        int lda = A.diagstep();
        int ldb = B.stepj();
        int Lap_info=0;
        if (A.iscm() && A.isconj()) {
            B.conjugateSelf();
            LAPNAME(ztbtrs) (
                LAPCM A.nlo()==0?LAPCH_UP:LAPCH_LO,
                LAPCH_NT, dt==UnitDiag?LAPCH_U:LAPCH_NU,
                LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
                LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
            B.conjugateSelf();
        } else {
            LAPNAME(ztbtrs) (
                LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
                A.iscm()?LAPCH_NT:A.isconj()?LAPCH_CT:LAPCH_T, 
                dt==UnitDiag?LAPCH_U:LAPCH_NU,
                LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
                LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
        }
        LAP_Results(Lap_info,"ztbtrs");
    }
#endif
#ifdef INST_FLOAT
    template <> void LapTriLDivEq(
        const GenBandMatrix<float>& A, MatrixView<float> B, DiagType dt)
    {
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 1);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(B.ct() == NonConj);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(B.iscm());
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int nrhs = B.rowsize();
        int aoffset = 
            (A.nlo()==0 && A.iscm()) ? A.nhi() :
            (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
        int lda = A.diagstep();
        int ldb = B.stepj();
        int Lap_info=0;
        LAPNAME(stbtrs) (
            LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
            A.iscm()?LAPCH_NT:LAPCH_T, dt==UnitDiag?LAPCH_U:LAPCH_NU,
            LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
            LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
        LAP_Results(Lap_info,"stbtrs");
    }
    template <> void LapTriLDivEq(
        const GenBandMatrix<std::complex<float> >& A,
        MatrixView<std::complex<float> > B, DiagType dt)
    {
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        TMVAssert(B.rowsize() > 1);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        TMVAssert(B.ct() == NonConj);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(B.iscm());
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int nrhs = B.rowsize();
        int aoffset = 
            (A.nlo()==0 && A.iscm()) ? A.nhi() :
            (A.nhi()==0 && A.isrm()) ? A.nlo() : 0;
        int lda = A.diagstep();
        int ldb = B.stepj();
        int Lap_info=0;
        if (A.iscm() && A.isconj()) {
            B.conjugateSelf();
            LAPNAME(ctbtrs) (
                LAPCM A.nlo()==0?LAPCH_UP:LAPCH_LO,
                LAPCH_NT, dt==UnitDiag?LAPCH_U:LAPCH_NU,
                LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
                LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
            B.conjugateSelf();
        } else {
            LAPNAME(ctbtrs) (
                LAPCM (A.iscm()?A.nlo():A.nhi())==0?LAPCH_UP:LAPCH_LO,
                A.iscm()?LAPCH_NT:A.isconj()?LAPCH_CT:LAPCH_T, 
                dt==UnitDiag?LAPCH_U:LAPCH_NU,
                LAPV(n),LAPV(kd),LAPV(nrhs), LAPP(A.cptr()-aoffset),LAPV(lda),
                LAPP(B.ptr()),LAPV(ldb) LAPINFO LAP1 LAP1 LAP1);
        }
        LAP_Results(Lap_info,"ctbtrs");
    }
#endif
#endif // LAP

    template <class T, class Ta> void TriLDivEq(
        const GenBandMatrix<Ta>& A, MatrixView<T> B, DiagType dt)
    {
#ifdef XDEBUG
        Matrix<T> B0 = B;
        Matrix<T> A0 = A;
        if (dt == UnitDiag) A0.diag().setAllTo(T(1));
#endif
        TMVAssert(A.isSquare());
        TMVAssert(B.colsize() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.colsize() > 0);
        if (B.rowsize() == 0) return;
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());

        if (B.rowsize() == 1) TriLDivEq(A,B.col(0),dt);
        else if (B.isconj()) 
            TriLDivEq(A.conjugate(),B.conjugate(),dt);
        else
#ifdef LAP
            if (!(B.iscm() && B.stepj()>0)) {
                Matrix<T,ColMajor> BB=B;
                TriLDivEq(A,BB.view(),dt);
                B = BB;
            } else if ( !((A.iscm() && A.stepj()>0) || 
                          (A.isrm() && A.stepi()>0))) {
                BandMatrix<T,ColMajor> AA = A;
                LapTriLDivEq(AA,B,dt);
            } else 
                LapTriLDivEq(A,B,dt);
#else
            if (dt == UnitDiag) NonLapTriLDivEq<true>(A,B);
            else NonLapTriLDivEq<false>(A,B);
#endif

#ifdef XDEBUG
            Matrix<T> BB = A0*B;
            if (Norm(BB-B0) > 0.001*Norm(B0)*Norm(A0)*Norm(A0.inverse())) {
                cerr<<"TriLDivEq Matrix:\n";
                cerr<<"A = "<<TMV_Text(A)<<A0<<endl;
                cerr<<"B = "<<TMV_Text(B)<<B0<<endl;
                cerr<<"--> B = "<<B<<endl;
                cerr<<"A*B = "<<BB<<endl;
                abort();
            }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandTriDiv_M.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


