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
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"

// Using CblasRowMajor in all the other files isn't working correcly with
// MKL 10.2.2.  The usage in here hasn't given me any problems, but 
// just to be safe, I'm disabling it here as well.  My guess is that I 
// haven't tested the V/B code as thoroughly as the others, so it might
// turn up as a problem in more complicated code.
#ifdef CBLAS
#undef CBLAS
#endif

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // TriLDivEq
    //

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.step()==1);
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(A.nlo() == 0);
        TMVAssert(A.nhi() > 0);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = b.size();

        const ptrdiff_t sj = (rm?1:A.stepj());
        const ptrdiff_t ds = A.diagstep();
        const Ta* Aii = A.cptr() + (ua ? N-2 : N-1)*ds;
        T* bi = b.ptr() + (ua ? N-2 : N-1);

        if (!ua) {
            if (*Aii==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularBandLU<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bi >= b._first);
            TMVAssert(bi < b._last);
#endif
            *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aii) : *Aii));
            Aii -= ds;
            --bi;
        }
        if (N==1) return;

        ptrdiff_t k = A.nhi()-1;
        for(ptrdiff_t i=N-1,len=1;i>0;--i,Aii-=ds,--bi) {
            // Actual row being done is i-1, not i

            // *bi -= A.row(i,i+1,j2) * b.subVector(i+1,j2);
            const T* bj = bi+1;
            const Ta* Aij = Aii+sj;
            for(ptrdiff_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (!ua) {
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aii) : *Aii));
            }
            if (k > 0) { --k; ++len; }
        } 
    }

    template <bool rm, bool ua, class T, class Ta> 
    static inline void RowUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            DoRowUpperTriLDivEq<rm,true,ua>(A,b);
        else
            DoRowUpperTriLDivEq<rm,false,ua>(A,b);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.step()==1);
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(A.nlo() == 0);
        TMVAssert(A.nhi() > 0);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = b.size();

        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t ds = A.diagstep();
        const ptrdiff_t hi = A.nhi();

        ptrdiff_t i1 = N-1;
        if (i1 > hi) i1 -= hi;
        else i1 = 0;

        const Ta* Ai1j = A.cptr()+(N-1)*sj+i1*si;
        const Ta* Ajj = (ua ? 0 : A.cptr()+(N-1)*ds); // if unit, this isn't used
        T* bi1 = b.ptr()+i1;
        T* bj = b.ptr()+N-1;

        for(ptrdiff_t len=N-1-i1;len>0;--bj) {
            if (*bj != T(0)) {
                if (!ua) {
                    if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                        std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                        exit(1); 
#else
                        throw SingularBandLU<Ta>(A);
#endif
                    }
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= b._first);
                    TMVAssert(bj < b._last);
#endif
                    *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
                    Ajj -= ds;
                }

                // b.subVector(i1,j) -= (*bj) * A.col(j,i1,j);
                T* bi = bi1;
                const Ta* Aij = Ai1j;
                for(ptrdiff_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(bi >= b._first);
                    TMVAssert(bi < b._last);
#endif
                    *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
            else if (!ua) Ajj -= ds;

            if (i1 > 0) { --i1; --bi1; Ai1j-=ds; } 
            else { --len; Ai1j-=sj; }
        } 
        if (!ua && *bj != T(0)) {
            if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular BandUpperTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularBandLU<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bj >= b._first);
            TMVAssert(bj < b._last);
#endif
            *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
        }
    }

    template <bool cm, bool ua, class T, class Ta> 
    static inline void ColUpperTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            DoColUpperTriLDivEq<cm,true,ua>(A,b);
        else
            DoColUpperTriLDivEq<cm,false,ua>(A,b);
    }

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.step()==1);
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(A.nhi() == 0);
        TMVAssert(A.nlo() > 0);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = A.colsize();

        const ptrdiff_t sj = (rm ? 1 : A.stepj());
        const ptrdiff_t si = A.stepi();
        const ptrdiff_t ds = A.diagstep();

        const Ta* Aij1 = A.cptr();
        const T* bj1 = b.cptr();
        T* bi = b.ptr();

        if (!ua) {
            if (*Aij1==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularBandLU<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bi >= b._first);
            TMVAssert(bi < b._last);
#endif
            *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aij1) : *Aij1));
        }

        ++bi;
        Aij1 += si;
        ptrdiff_t k=A.nlo()-1;

        for(ptrdiff_t i=1,len=1;i<N;++i,++bi) {
            // *bi -= A.row(i,j1,i) * b.subVector(j1,i);
            const Ta* Aij = Aij1;
            const T* bj = bj1;
            for(ptrdiff_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(bj >= b._first);
                TMVAssert(bj < b._last);
#endif
                *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
            if (!ua) {
                // Aij is Aii after the above for loop
                if (*Aij == Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularBandLU<Ta>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aij) : *Aij));
            }
            if (k>0) { --k; ++len; Aij1+=A.stepi(); } 
            else { ++bj1; Aij1+=ds; }
        }
    }

    template <bool rm, bool ua, class T, class Ta> 
    static inline void RowLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            DoRowLowerTriLDivEq<rm,true,ua>(A,b);
        else
            DoRowLowerTriLDivEq<rm,false,ua>(A,b);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.step()==1);
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(A.nhi() == 0);
        TMVAssert(A.nlo() > 0);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());

        const ptrdiff_t N = A.colsize();

        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t ds = A.diagstep();

        const Ta* Ajj= A.cptr();
        T* bj = b.ptr();

        ptrdiff_t i2=TMV_MIN(A.nlo()+1,A.colsize());

        for(ptrdiff_t len=i2-1;len>0;++bj,Ajj+=ds) {
            if (*bj != T(0)) {
                if (!ua) {
                    if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                        std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                        exit(1); 
#else
                        throw SingularBandLU<Ta>(A);
#endif
                    }
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= b._first);
                    TMVAssert(bj < b._last);
#endif
                    *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
                }
                // b.subVector(j+1,i2) -= *bj * A.col(j,j+1,i2)
                T* bi = bj+1;
                const Ta* Aij = Ajj+si;
                for(ptrdiff_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(bi >= b._first);
                    TMVAssert(bi < b._last);
#endif
                    *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
            if (i2 < N) ++i2; else --len;
        }
        if (!ua && *bj != T(0)) {
            if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular BandLowerTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularBandLU<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bj >= b._first);
            TMVAssert(bj < b._last);
#endif
            *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
        }
    }

    template <bool rm, bool ua, class T, class Ta> 
    static inline void ColLowerTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            DoColLowerTriLDivEq<rm,true,ua>(A,b);
        else
            DoColLowerTriLDivEq<rm,false,ua>(A,b);
    }

    template <class T, class Ta> 
    static void UpperTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b, DiagType dt)
    {
        TMVAssert(A.nlo() == 0);
        if (A.nhi() == 0) {
            if (dt == NonUnitDiag) b /= DiagMatrixViewOf(A.diag());
        } else {
            if (dt == UnitDiag)
                if (A.isrm()) RowUpperTriLDivEq<true,true>(A,b);
                else if (A.iscm()) ColUpperTriLDivEq<true,true>(A,b);
                else RowUpperTriLDivEq<false,true>(A,b);
            else
                if (A.isrm()) RowUpperTriLDivEq<true,false>(A,b);
                else if (A.iscm()) ColUpperTriLDivEq<true,false>(A,b);
                else RowUpperTriLDivEq<false,false>(A,b);
        }
    }

    template <class T, class Ta> 
    static void LowerTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b, DiagType dt)
    {
        TMVAssert(A.nhi() == 0);
        if (A.nlo() == 0) {
            if (dt == NonUnitDiag) b /= DiagMatrixViewOf(A.diag());
        } else {
            if (dt == UnitDiag)
                if (A.isrm()) RowLowerTriLDivEq<true,true>(A,b);
                else if (A.iscm()) ColLowerTriLDivEq<true,true>(A,b);
                else RowLowerTriLDivEq<false,true>(A,b);
            else
                if (A.isrm()) RowLowerTriLDivEq<true,false>(A,b);
                else if (A.iscm()) ColLowerTriLDivEq<true,false>(A,b);
                else RowLowerTriLDivEq<false,false>(A,b);
        }
    }

    template <class T, class Ta> 
    static void NonBlasTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b, DiagType dt)
    {
        //cout<<"Start NonBlasTriLDvEq\n";
        // Solve A x = y  where A is an upperor lower band triangle matrix
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(b.size() > 0);
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(b.ct() == NonConj);

        const ptrdiff_t N = b.size();
        if (b.step() == 1) {
            if (A.nlo() == 0) {
                ptrdiff_t i2 = N;
                for(const T* b2 = b.cptr()+i2-1; i2>0 && *b2==T(0); --i2,--b2);
                if (i2==0) return;
                else if (i2 == N) 
                    UpperTriLDivEq(A,b,dt);
                else if (A.nhi() < i2)
                    UpperTriLDivEq(A.subBandMatrix(0,i2,0,i2,0,A.nhi()),
                                   b.subVector(0,i2),dt);
                else
                    UpperTriLDivEq(A.subBandMatrix(0,i2,0,i2,0,i2-1),
                                   b.subVector(0,i2),dt);
            } else {
                ptrdiff_t i1 = 0;
                for(const T* b1 = b.ptr(); i1<N && *b1==T(0); ++i1,++b1);
                if (i1==N) return;
                else if (i1 == 0) 
                    LowerTriLDivEq(A,b,dt);
                else if (A.nlo() < (N-i1))
                    LowerTriLDivEq(A.subBandMatrix(i1,N,i1,N,A.nlo(),0),
                                   b.subVector(i1,N),dt);
                else
                    LowerTriLDivEq(A.subBandMatrix(i1,N,i1,N,N-i1-1,0),
                                   b.subVector(i1,N),dt);
            }
        } else {
            Vector<T> bb = b;
            NonBlasTriLDivEq(A,bb.view(),dt);
            b = bb;
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasTriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b, DiagType dt)
    { NonBlasTriLDivEq(A,b,dt); }
#ifdef INST_DOUBLE
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<double>& A, VectorView<double> b,
        DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = b.step();
        double* bp = b.ptr();
        if (s < 0) bp += (n-1)*s;
        BLASNAME(dtbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > b, DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = b.step();
        std::complex<double>* bp = b.ptr();
        if (s < 0) bp += (n-1)*s;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztbsv) (
                BLASRM A.nlo()==0 ? BLASCH_LO : BLASCH_UP, BLASCH_CT,
                dt==UnitDiag ? BLASCH_U : BLASCH_NU, 
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ztbsv) (
                BLASCM A.nlo()==0?BLASCH_UP:BLASCH_LO, BLASCH_NT, 
                dt==UnitDiag ? BLASCH_U : BLASCH_NU,
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ztbsv) (
                BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
                A.isrm() ? A.isconj() ? BLASCH_CT : BLASCH_T : BLASCH_NT, 
                dt==UnitDiag ? BLASCH_U : BLASCH_NU,
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<double>& A, 
        VectorView<std::complex<double> > b, DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = 2*b.step();
        double* bp = (double*) b.ptr();
        if (s < 0) bp += (n-1)*s;
        BLASNAME(dtbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
        BLASNAME(dtbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp+1), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<float>& A, VectorView<float> b, DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = b.step();
        float* bp = b.ptr();
        if (s < 0) bp += (n-1)*s;
        BLASNAME(stbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > b, DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = b.step();
        std::complex<float>* bp = b.ptr();
        if (s < 0) bp += (n-1)*s;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ctbsv) (
                BLASRM A.nlo()==0 ? BLASCH_LO : BLASCH_UP, BLASCH_CT,
                dt==UnitDiag ? BLASCH_U : BLASCH_NU, 
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ctbsv) (
                BLASCM A.nlo()==0?BLASCH_UP:BLASCH_LO, BLASCH_NT, 
                dt==UnitDiag ? BLASCH_U : BLASCH_NU,
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ctbsv) (
                BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
                A.isrm() ? A.isconj() ? BLASCH_CT : BLASCH_T : BLASCH_NT, 
                dt==UnitDiag ? BLASCH_U : BLASCH_NU,
                BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
                BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenBandMatrix<float>& A, 
        VectorView<std::complex<float> > b, DiagType dt)
    {
        int n = A.colsize();
        int kd = A.nlo()==0 ? A.nhi() : A.nlo();
        int aoffset = A.isrm() ? A.nlo() : A.nhi();
        int ds = A.diagstep();
        int s = 2*b.step();
        float* bp = (float*) b.ptr();
        if (s < 0) bp += (n-1)*s;
        BLASNAME(stbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp), BLASV(s) BLAS1 BLAS1 BLAS1);
        BLASNAME(stbsv) (
            BLASCM ((A.nlo()==0) == A.isrm()) ?BLASCH_LO:BLASCH_UP, 
            A.isrm() ? BLASCH_T : BLASCH_NT, 
            dt==UnitDiag ? BLASCH_U : BLASCH_NU,
            BLASV(n), BLASV(kd), BLASP(A.cptr()-aoffset), BLASV(ds),
            BLASP(bp+1), BLASV(s) BLAS1 BLAS1 BLAS1);
    }
#endif 
#endif // BLAS

    template <class T, class Ta> 
    void TriLDivEq(
        const GenBandMatrix<Ta>& A, VectorView<T> b, DiagType dt)
    {
#ifdef XDEBUG
        Vector<T> b0 = b;
        Matrix<T> A0 = A;
        if (dt == UnitDiag) A0.diag().setAllTo(T(1));
#endif
        TMVAssert(A.isSquare());
        TMVAssert(b.size() == A.colsize());
        TMVAssert(A.nlo() == 0 || A.nhi() == 0);
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.nhi() < A.colsize());
        TMVAssert(A.nlo() < A.colsize());
        if (b.isconj())
            TriLDivEq(A.conjugate(),b.conjugate(),dt);
        else {
#ifdef BLAS
            if (A.isrm() || A.iscm())
                BlasTriLDivEq(A,b,dt);
            else {
                BandMatrix<Ta,ColMajor> AA = A;
                BlasTriLDivEq(AA,b,dt);
            }
#else
            NonBlasTriLDivEq(A,b,dt);
#endif
        }

#ifdef XDEBUG
        Vector<T> bb = A0*b;
        if (Norm(bb-b0) > 0.001*Norm(b0)*Norm(A0)*Norm(A0.inverse())) {
            cerr<<"TriLDivEq Vector:\n";
            cerr<<"A = "<<TMV_Text(A)<<A0<<endl;
            cerr<<"b = "<<TMV_Text(b)<<b0<<endl;
            cerr<<"--> b = "<<b<<endl;
            cerr<<"A*b = "<<bb<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandTriDiv_V.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


