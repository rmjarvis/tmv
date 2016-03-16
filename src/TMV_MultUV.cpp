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
#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "TMV_MultMV.h"

// CBLAS trick of using RowMajor with ConjTrans when we have a 
// case of A.conjugate() * x doesn't seem to be working with MKL 10.2.2.
// I haven't been able to figure out why.  (e.g. Is it a bug in the MKL
// code, or am I doing something wrong?)  So for now, just disable it.
#ifdef CBLAS
#undef CBLAS
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    const T* UpperTriMatrixComposite<T>::cptr() const
    {
        if (!itsm.get()) {
            itsm.resize(this->size()*this->size());
            UpperTriMatrixView<T>(
                itsm.get(),this->size(),stepi(),stepj(),this->dt(),NonConj) =
                *this;
        }
        return itsm.get();
    }

    template <class T> 
    ptrdiff_t UpperTriMatrixComposite<T>::stepi() const 
    { return 1; }

    template <class T> 
    ptrdiff_t UpperTriMatrixComposite<T>::stepj() const 
    { return this->size(); }

    template <class T> 
    const T* LowerTriMatrixComposite<T>::cptr() const
    {
        if (!itsm.get()) {
            itsm.resize(this->size()*this->size());
            LowerTriMatrixView<T>(
                itsm.get(),this->size(),stepi(),stepj(),this->dt(),NonConj) =
                *this;
        }
        return itsm.get();
    }

    template <class T> 
    ptrdiff_t LowerTriMatrixComposite<T>::stepi() const 
    { return 1; }

    template <class T> 
    ptrdiff_t LowerTriMatrixComposite<T>::stepj() const 
    { return this->size(); }

    // 
    // MultEqMV
    //

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    {
        //cout<<"RowMultEqMV Upper\n";
        TMVAssert(x.step()==1);
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = x.size();

        const ptrdiff_t sj = rm ? 1 : A.stepj();
        const ptrdiff_t ds = A.stepi()+sj;
        T* xi = x.ptr();
        const Ta* Aii = A.cptr();
        ptrdiff_t len = N-1;

        for(; len>0; --len,++xi,Aii+=ds) {
            // i = 0..N-2
            // x(i) = A.row(i,i,N) * x.subVector(i,N);
            if (!ua) *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
            const T* xj = xi+1;
            const Ta* Aij = Aii+sj;
            for(ptrdiff_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi += (*xj) * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
        }
#ifdef TMVFLDEBUG
        TMVAssert(xi >= x._first);
        TMVAssert(xi < x._last);
#endif
        if (!ua) *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
    }

    template <bool rm, class T, class Ta> 
    static void RowMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            if (A.isunit())
                DoRowMultEqMV<rm,true,true>(A,x);
            else
                DoRowMultEqMV<rm,true,false>(A,x);
        else
            if (A.isunit())
                DoRowMultEqMV<rm,false,true>(A,x);
            else
                DoRowMultEqMV<rm,false,false>(A,x);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    {
        //cout<<"ColMultEqMV Upper\n";
        TMVAssert(x.step()==1);
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = x.size();

        T* x0 = x.ptr();
        const T* xj = x0+1;
        const ptrdiff_t si = cm ? 1 : A.stepi();
        const Ta* A0j = A.cptr();

#ifdef TMVFLDEBUG
        TMVAssert(x0 >= x._first);
        TMVAssert(x0 < x._last);
#endif
        if (!ua) *x0 *= (ca ? TMV_CONJ(*A0j) : *A0j);
        A0j += A.stepj();

        for(ptrdiff_t len=1; len<N; ++len,++xj,A0j+=A.stepj()) if (*xj != T(0)) {
            // j = 1..N-1
            // x.subVector(0,j) += *xj * A.col(j,0,j);
            const Ta* Aij = A0j;
            T* xi = x0;
            for(ptrdiff_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
            // Now Aij == Ajj, xi == xj
            // so this next statement is really *xj *= *Ajj
#ifdef TMVFLDEBUG
            TMVAssert(xi >= x._first);
            TMVAssert(xi < x._last);
#endif
            if (!ua) *xi *= (ca ? TMV_CONJ(*Aij) : *Aij);
        }
    }

    template <bool cm, class T, class Ta> 
    static void ColMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            if (A.isunit())
                DoColMultEqMV<cm,true,true>(A,x);
            else
                DoColMultEqMV<cm,true,false>(A,x);
        else
            if (A.isunit())
                DoColMultEqMV<cm,false,true>(A,x);
            else
                DoColMultEqMV<cm,false,false>(A,x);
    }

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    {
        //cout<<"RowMultEqMV Lower\n";
        TMVAssert(x.step()==1);
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = x.size();
        const ptrdiff_t si = A.stepi();
        const ptrdiff_t sj = rm ? 1 : A.stepj();
        const ptrdiff_t ds = si+sj;

        const T* x0 = x.cptr();
        T* xi = x.ptr() + N-1;
        const Ta* Ai0 = A.cptr()+(N-1)*si;
        const Ta* Aii = Ai0 + (N-1)*sj;

        for(ptrdiff_t len=N-1; len>0; --len,--xi,Ai0-=si,Aii-=ds) {
            // i = N-1..1
            // x(i) = A.row(i,0,i+1) * x.subVector(0,i+1);
            T xx = *xi;
            if (!ua) xx *= (ca ? TMV_CONJ(*Aii) : *Aii);
            const Ta* Aij = Ai0;
            const T* xj = x0;
            for(ptrdiff_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
                xx += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
#ifdef TMVFLDEBUG
            TMVAssert(xi >= x._first);
            TMVAssert(xi < x._last);
#endif
            *xi = xx;
        }
#ifdef TMVFLDEBUG
        TMVAssert(xi >= x._first);
        TMVAssert(xi < x._last);
#endif
        if (!ua) *xi *= (ca ? TMV_CONJ(*Aii) : *Aii);
    }

    template <bool rm, class T, class Ta> 
    static void RowMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            if (A.isunit())
                DoRowMultEqMV<rm,true,true>(A,x);
            else
                DoRowMultEqMV<rm,true,false>(A,x);
        else
            if (A.isunit())
                DoRowMultEqMV<rm,false,true>(A,x);
            else
                DoRowMultEqMV<rm,false,false>(A,x);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    {
        //cout<<"ColMultEqMV Lower\n";
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.ct() == NonConj);
        TMVAssert(x.step() == 1);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = x.size();

        T* xj = x.ptr() + N-2;
        const ptrdiff_t si = cm ? 1 : A.stepi();
        const ptrdiff_t ds = A.stepj()+si;
        const Ta* Ajj = A.cptr()+(N-2)*ds;

#ifdef TMVFLDEBUG
        TMVAssert(xj+1 >= x._first);
        TMVAssert(xj+1 < x._last);
#endif
        if (!ua) *(xj+1) *= (ca ? TMV_CONJ(*(Ajj+ds)) : *(Ajj+ds));
        for(ptrdiff_t jj=N-1,len=1;jj>0;--jj,++len,--xj,Ajj-=ds) if (*xj!=T(0)) {
            // j = N-2..0, jj = j+1
            // x.subVector(j+1,N) += *xj * A.col(j,j+1,N);
            T* xi = xj+1;
            const Ta* Aij = Ajj+si;
            for (ptrdiff_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                TMVAssert(xi >= x._first);
                TMVAssert(xi < x._last);
#endif
                *xi += *xj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
#ifdef TMVFLDEBUG
            TMVAssert(xj >= x._first);
            TMVAssert(xj < x._last);
#endif
            if (!ua) *xj *= (ca ? TMV_CONJ(*Ajj) : *Ajj);
        }
    }

    template <bool cm, class T, class Ta> 
    static void ColMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    { 
        if (A.isconj())
            if (A.isunit())
                DoColMultEqMV<cm,true,true>(A,x);
            else
                DoColMultEqMV<cm,true,false>(A,x);
        else
            if (A.isunit())
                DoColMultEqMV<cm,false,true>(A,x);
            else
                DoColMultEqMV<cm,false,false>(A,x);
    }

    template <class T, class Ta> 
    static inline void DoMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    // x = A * x
    {
        if (A.isrm()) RowMultEqMV<true>(A,x);
        else if (A.iscm()) ColMultEqMV<true>(A,x);
        else RowMultEqMV<false>(A,x);
    }

    template <class T, class Ta> 
    static void NonBlasMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        //     [ A11 A12 A13 ] [ 0  ]   [ A12 x2 ]
        // x = [  0  A22 A23 ] [ x2 ] = [ A22 x2 ]
        //     [  0   0  A33 ] [ 0  ]   [   0    ]

        const ptrdiff_t N = x.size(); // = A.size()
        ptrdiff_t j2 = N;
        for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
        if (j2 == 0) return;
        ptrdiff_t j1 = 0;
        for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) {
            DoMultEqMV(A,x);
        } else {
            TMVAssert(j1 < j2);
            ConstUpperTriMatrixView<Ta> A22 = A.subTriMatrix(j1,j2);
            VectorView<T> x2 = x.subVector(j1,j2);

            if (j1 != 0) {
                ConstMatrixView<Ta> A12 = A.subMatrix(0,j1,j1,j2);
                VectorView<T> x1 = x.subVector(0,j1);
                UnitAMultMV1<true,false>(A12,x2,x1);
            }
            DoMultEqMV(A22,x2);
        }
    }

    template <class T, class Ta> 
    static inline void DoMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    {
        if (A.isrm()) RowMultEqMV<true>(A,x);
        else if (A.iscm() && !SameStorage(A,x))
            ColMultEqMV<true>(A,x);
        else RowMultEqMV<false>(A,x);
    }

    template <class T, class Ta> 
    static void NonBlasMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    // x = A * x
    {
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);
        TMVAssert(x.ct() == NonConj);

        //     [ A11  0   0  ] [ 0  ]   [   0    ]
        // x = [ A21 A22  0  ] [ x2 ] = [ A22 x2 ]
        //     [ A31 A32 A33 ] [ 0  ]   [ A32 x2 ]

        const ptrdiff_t N = x.size(); // = A.size()
        ptrdiff_t j2 = N;
        for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
        if (j2 == 0) return;
        ptrdiff_t j1 = 0;
        for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
        if (j1 == 0 && j2 == N) {
            DoMultEqMV(A,x);
        } else {
            TMVAssert(j1 < j2);
            ConstLowerTriMatrixView<Ta> A22 = A.subTriMatrix(j1,j2);
            VectorView<T> x2 = x.subVector(j1,j2);

            if (j2 != N) {
                ConstMatrixView<Ta> A32 = A.subMatrix(j2,N,j1,j2);
                VectorView<T> x3 = x.subVector(j2,N);
                UnitAMultMV1<true,false>(A32,x2,x3);
            }
            DoMultEqMV(A22,x2);
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasMultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    { NonBlasMultEqMV(A,x); }
    template <class T, class Ta> 
    static inline void BlasMultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    { NonBlasMultEqMV(A,x); }
#ifdef INST_DOUBLE
    template <> 
    void BlasMultEqMV( 
        const GenUpperTriMatrix<double>& A, VectorView<double> x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenLowerTriMatrix<double>& A, VectorView<double> x)
    {
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
    template <> 
    void BlasMultEqMV(
        const GenUpperTriMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > x)
    {
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
    template <> 
    void BlasMultEqMV(
        const GenLowerTriMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenUpperTriMatrix<double>& A,
        VectorView<std::complex<double> > x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenLowerTriMatrix<double>& A,
        VectorView<std::complex<double> > x)
    {
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
#ifdef INST_FLOAT
    template <> 
    void BlasMultEqMV( 
        const GenUpperTriMatrix<float>& A, VectorView<float> x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenLowerTriMatrix<float>& A, VectorView<float> x)
    {
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
    template <> 
    void BlasMultEqMV(
        const GenUpperTriMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > x)
    {
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
    template <> 
    void BlasMultEqMV(
        const GenLowerTriMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenUpperTriMatrix<float>& A,
        VectorView<std::complex<float> > x)
    {
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
    template <> 
    void BlasMultEqMV( 
        const GenLowerTriMatrix<float>& A,
        VectorView<std::complex<float> > x)
    {
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
#endif // FLOAT
#endif // BLAS

    template <class T, class Ta> 
    static void MultEqMV(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> x)
    {
#ifdef XDEBUG
        Vector<T> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> x2 = A0 * x0;
        //cout<<"MultEqMV: A = "<<A<<"x = "<<x<<endl;
#endif
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);

        if (x.isconj()) MultEqMV(A.conjugate(),x.conjugate());
        else {
#ifdef BLAS
            if ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) 
                BlasMultEqMV(A,x);
            else {
                if (A.isunit()) {
                    UpperTriMatrix<Ta,UnitDiag|RowMajor> A2(A);
                    BlasMultEqMV(A2,x);
                } else {
                    UpperTriMatrix<Ta,NonUnitDiag|RowMajor> A2(A);
                    BlasMultEqMV(A2,x);
                }
            }
#else
            NonBlasMultEqMV(A,x);
#endif
        }
#ifdef XDEBUG
        //cout<<"-> x = "<<x<<endl<<"x2 = "<<x2<<endl;
        if (!(Norm(x-x2) <= 0.001*(Norm(A0)*Norm(x0)))) {
            cerr<<"MultEqMV: \n";
            cerr<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<x.cptr()<<"  "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"-> x = "<<x<<endl;
            cerr<<"x2 = "<<x2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    static void MultEqMV(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> x)
    {
#ifdef XDEBUG
        Vector<T> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> x2 = A0 * x0;
#endif
        TMVAssert(A.size() == x.size());
        TMVAssert(x.size() > 0);
        TMVAssert(x.step() == 1);

        if (x.isconj()) MultEqMV(A.conjugate(),x.conjugate());
        else {
#ifdef BLAS
            if ( (A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0) )
                BlasMultEqMV(A,x);
            else {
                if (A.isunit()) {
                    LowerTriMatrix<Ta,UnitDiag|RowMajor> A2(A);
                    BlasMultEqMV(A2,x);
                } else {
                    LowerTriMatrix<Ta,NonUnitDiag|RowMajor> A2(A);
                    BlasMultEqMV(A2,x);
                }
            }
#else
            NonBlasMultEqMV(A,x);
#endif
        }

#ifdef XDEBUG
        if (!(Norm(x-x2) <= 0.001*(Norm(A0)*Norm(x0)))) {
            cerr<<"MultEqMV: \n";
            cerr<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<x.cptr()<<"  "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"-> x = "<<x<<endl;
            cerr<<"x2 = "<<x2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y)
    // y (+)= alpha * A * x
    { 
#ifdef XDEBUG
        Vector<Tx> x0 = x;
        Vector<T> y0 = y;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = alpha*A0*x0;
        if (add) y2 += y0;
#endif
        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());

        if (y.size() > 0) {
            if (alpha==T(0)) {
                if (!add) y.setZero();
            } else if (!add && y.step() == 1) {
                y = x;
                MultEqMV(A,y);
                y *= alpha;
            } else {
                Vector<T> xx = alpha*x;
                MultEqMV(A,xx.view());
                if (add) y += xx;
                else y = xx;
            }
        } 
#ifdef XDEBUG
        if (!(Norm(y-y2) <=
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"-> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            abort();
        }
#endif
    }

    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y)
    // y (+)= alpha * A * x
    { 
#ifdef XDEBUG
        Vector<T> y0 = y;
        Vector<Tx> x0 = x;
        Matrix<Ta> A0 = A;
        Vector<T> y2 = alpha*A0*x0;
        if (add) y2 += y0;
#endif

        TMVAssert(A.size() == x.size());
        TMVAssert(A.size() == y.size());

        if (y.size() > 0) {
            if (alpha==T(0)) {
                if (!add) y.setZero();
            } else if (!add && y.step() == 1) {
                y = x;
                MultEqMV(A,y);
                if (alpha != T(1)) y *= alpha;
            } else {
                Vector<T> xx = alpha*x;
                MultEqMV(A,xx.view());
                if (add) y += xx;
                else y = xx;
            }
        }
#ifdef XDEBUG
        if (!(Norm(y-y2) <= 
              0.001*(TMV_ABS(alpha)*Norm(A0)*Norm(x0)+
                     (add?Norm(y0):TMV_RealType(T)(0))))) {
            cerr<<"MultMV: alpha = "<<alpha<<endl;
            cerr<<"add = "<<add<<endl;
            cerr<<"A = "<<A.cptr()<<"  "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"x = "<<x.cptr()<<"  "<<TMV_Text(x)<<" step "<<x.step()<<"  "<<x0<<endl;
            cerr<<"y = "<<y.cptr()<<"  "<<TMV_Text(y)<<" step "<<y.step()<<"  "<<y0<<endl;
            cerr<<"-> y = "<<y<<endl;
            cerr<<"y2 = "<<y2<<endl;
            abort();
        }
#endif
    }

#define InstFile "TMV_MultUV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


