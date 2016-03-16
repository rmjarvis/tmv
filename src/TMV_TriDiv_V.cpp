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
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"

// Using CblasRowMajor in all the other files isn't working correcly with
// MKL 10.2.2.  The usage in here hasn't given me any problems, but 
// just to be safe, I'm disabling it here as well.  My guess is that I 
// haven't tested the V/U code as thoroughly as the others, so it might
// turn up as a problem in more complicated code.
#ifdef CBLAS
#undef CBLAS
#endif

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // TriLDivEq V
    //

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        // Solve A x = y  where A is an upper triangle matrix
        //cout<<"Row Upper\n";
        TMVAssert(b.step()==1);
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = A.size();

        const ptrdiff_t sj = (rm?1:A.stepj());
        const ptrdiff_t ds = A.stepi()+sj;
        const Ta* Aii = A.cptr() + (ua ? N-2 : N-1)*ds;
        T* bi = b.ptr() + (ua ? N-2 : N-1);

        if (!ua) {
            if (*Aii==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular UpperTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularUpperTriMatrix<Ta>(A);
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

        for(ptrdiff_t i=N-1,len=1; i>0; --i,++len,Aii-=ds,--bi) {
            // Actual row being done is i-1, not i

            // *bi -= A.row(i,i+1,N) * b.subVector(i+1,N);
            const T* bj = bi+1;
            const Ta* Aij = Aii + sj;
            for(ptrdiff_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi -= (*bj) * (ca ? TMV_CONJ(*Aij) : *Aij);
            }

            if (!ua) {
                if (*Aii==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular UpperTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularUpperTriMatrix<Ta>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aii) : *Aii));
            }
        }
    }

    template <bool rm, class T, class Ta> 
    static inline void RowTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            if (A.isunit())
                DoRowTriLDivEq<rm,true,true>(A,b);
            else
                DoRowTriLDivEq<rm,true,false>(A,b);
        else
            if (A.isunit())
                DoRowTriLDivEq<rm,false,true>(A,b);
            else
                DoRowTriLDivEq<rm,false,false>(A,b);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        //cout<<"colmajor upper\n";
        // Solve A x = y  where A is an upper triangle matrix
        TMVAssert(b.step()==1);
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = A.size();

        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t sj = A.stepj();
        const ptrdiff_t ds = si+sj;
        const Ta* A0j = A.cptr()+(N-1)*sj;
        const Ta* Ajj = (ua ? 0 : A0j+(N-1)*si); // if unit, this isn't used.
        T*const b0 = b.ptr();
        T* bj = b0 + N-1;

        for(ptrdiff_t j=N-1; j>0; --j,--bj,A0j-=sj) {
            if (*bj != T(0)) {
                if (!ua) {
                    if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                        std::cerr<<"Singular UpperTriMatrix found\n"; 
                        exit(1); 
#else
                        throw SingularUpperTriMatrix<Ta>(A);
#endif
                    }
#ifdef TMVFLDEBUG
                    TMVAssert(bj >= b._first);
                    TMVAssert(bj < b._last);
#endif
                    *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
                    Ajj-=ds;
                }

                // b.subVector(0,j) -= *bj * A.col(j,0,j);
                T* bi = b0;
                const Ta* Aij = A0j;
                for(ptrdiff_t i=j;i>0;--i,++bi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
                    TMVAssert(bi >= b._first);
                    TMVAssert(bi < b._last);
#endif
                    *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
                }
            }
            else if (!ua) Ajj -= ds;
        }
        if (!ua && *bj != T(0)) {
            if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular UpperTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularUpperTriMatrix<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bj >= b._first);
            TMVAssert(bj < b._last);
#endif
            *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
        } 
    }

    template <bool cm, class T, class Ta> 
    static inline void ColTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            if (A.isunit())
                DoColTriLDivEq<cm,true,true>(A,b);
            else
                DoColTriLDivEq<cm,true,false>(A,b);
        else
            if (A.isunit())
                DoColTriLDivEq<cm,false,true>(A,b);
            else
                DoColTriLDivEq<cm,false,false>(A,b);
    }

    template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        // Solve A x = y  where A is a lower triangle matrix
        TMVAssert(b.step()==1);
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(rm == A.isrm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = A.size();

        const ptrdiff_t sj = (rm ? 1 : A.stepj());
        const ptrdiff_t si = A.stepi();

        const Ta* Ai0 = A.cptr();
        T* b0 = b.ptr();

        if (!ua) {
            if (*Ai0==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular LowerTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularLowerTriMatrix<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(b0 >= b._first);
            TMVAssert(b0 < b._last);
#endif
            *b0 = TMV_Divide(*b0,(ca ? TMV_CONJ(*Ai0) : *Ai0));
        }

        T* bi = b0+1;
        Ai0 += si;
        for(ptrdiff_t i=1,len=1;i<N;++i,++len,++bi,Ai0+=si) {
            // *bi -= A.row(i,0,i) * b.subVector(0,i);
            const Ta* Aij = Ai0;
            const T* bj = b0;
            for(ptrdiff_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)){
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi -= (*bj) * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
            if (!ua) {
                // Aij is Aii after the above for loop
                if (*Aij==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular LowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularLowerTriMatrix<Ta>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi = TMV_Divide(*bi,(ca ? TMV_CONJ(*Aij) : *Aij));
            }
        }
    }

    template <bool rm, class T, class Ta> 
    static inline void RowTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            if (A.isunit())
                DoRowTriLDivEq<rm,true,true>(A,b);
            else
                DoRowTriLDivEq<rm,true,false>(A,b);
        else
            if (A.isunit())
                DoRowTriLDivEq<rm,false,true>(A,b);
            else
                DoRowTriLDivEq<rm,false,false>(A,b);
    }

    template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        // Solve A x = y  where A is a lower triangle matrix
        TMVAssert(b.step()==1);
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size() > 0);
        TMVAssert(b.ct() == NonConj);
        TMVAssert(cm == A.iscm());
        TMVAssert(ca == A.isconj());
        TMVAssert(ua == A.isunit());

        const ptrdiff_t N = A.size();

        const ptrdiff_t si = (cm ? 1 : A.stepi());
        const ptrdiff_t ds = A.stepj()+si;
        const Ta* Ajj = A.cptr();
        T* bj = b.ptr();

        for(ptrdiff_t j=0,len=N-1;len>0;++j,--len,++bj,Ajj+=ds) if (*bj != T(0)) {
            if (!ua) {
                if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular LowerTriMatrix found\n"; 
                    exit(1); 
#else
                    throw SingularLowerTriMatrix<Ta>(A);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(bj >= b._first);
                TMVAssert(bj < b._last);
#endif
                *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
            }
            // b.subVector(j+1,N) -= *bj * A.col(j,j+1,N);
            T* bi = bj+1;
            const Ta* Aij = Ajj+si;
            for(ptrdiff_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si)){
#ifdef TMVFLDEBUG
                TMVAssert(bi >= b._first);
                TMVAssert(bi < b._last);
#endif
                *bi -= *bj * (ca ? TMV_CONJ(*Aij) : *Aij);
            }
        }
        if (!ua && *bj != T(0)) {
            if (*Ajj==Ta(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular LowerTriMatrix found\n"; 
                exit(1); 
#else
                throw SingularLowerTriMatrix<Ta>(A);
#endif
            }
#ifdef TMVFLDEBUG
            TMVAssert(bj >= b._first);
            TMVAssert(bj < b._last);
#endif
            *bj = TMV_Divide(*bj,(ca ? TMV_CONJ(*Ajj) : *Ajj));
        } 
    }

    template <bool cm, class T, class Ta> 
    static inline void ColTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isconj())
            if (A.isunit())
                DoColTriLDivEq<cm,true,true>(A,b);
            else
                DoColTriLDivEq<cm,true,false>(A,b);
        else
            if (A.isunit())
                DoColTriLDivEq<cm,false,true>(A,b);
            else
                DoColTriLDivEq<cm,false,false>(A,b);
    }

    template <class T, class Ta> 
    static inline void DoTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isrm()) RowTriLDivEq<true>(A,b);
        else if (A.iscm()) ColTriLDivEq<true>(A,b);
        else RowTriLDivEq<false>(A,b); 
    }

    template <class T, class Ta> 
    static inline void DoTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        if (A.isrm()) RowTriLDivEq<true>(A,b);
        else if (A.iscm()) ColTriLDivEq<true>(A,b);
        else RowTriLDivEq<false>(A,b); 
    }

    template <class T, class Ta> 
    static void NonBlasTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        //cout<<"Upper LDivEq vect\n";
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size()>0);
        TMVAssert(b.ct() == NonConj);

        if (b.step() == 1) {
            ptrdiff_t i2 = b.size();
            for(const T* b2 = b.cptr()+i2-1; i2>0 && *b2==T(0); --i2,--b2);
            if (i2==0) return;
            else DoTriLDivEq(A.subTriMatrix(0,i2),b.subVector(0,i2));
        } else {
            Vector<T> bb = b;
            NonBlasTriLDivEq(A,bb.view());
            b = bb;
        }
    }

    template <class T, class Ta> 
    static void NonBlasTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        //cout<<"Lower LDivEq vect\n";
        TMVAssert(A.size() == b.size());
        TMVAssert(b.size()>0);
        TMVAssert(b.ct() == NonConj);

        if (b.step() == 1) {
            const ptrdiff_t N = b.size();
            ptrdiff_t i1 = 0;
            for(const T* b1 = b.cptr(); i1<N && *b1==T(0); ++i1,++b1);
            if (i1==N) return;
            else
                DoTriLDivEq(A.subTriMatrix(i1,N),b.subVector(i1,N));
        } else {
            Vector<T> bb = b;
            NonBlasTriLDivEq(A,bb.view());
            b = bb;
        }
    }

#ifdef BLAS
    template <class T, class Ta> 
    static inline void BlasTriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    { NonBlasTriLDivEq(A,b); }
    template <class T, class Ta> 
    static inline void BlasTriLDivEq(
        const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    { NonBlasTriLDivEq(A,b); }
#ifdef INST_DOUBLE
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<double>& A, VectorView<double> b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        double* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<double>& A, VectorView<double> b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        double* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        std::complex<double>* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztrsv) (
                BLASRM BLASCH_LO,BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ztrsv) (
                BLASCM BLASCH_UP,BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ztrsv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<std::complex<double> >& A,
        VectorView<std::complex<double> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        std::complex<double>* bp = b.ptr();
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ztrsv) (
                BLASRM BLASCH_UP,BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ztrsv) (
                BLASCM BLASCH_LO,BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ztrsv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<double>& A,
        VectorView<std::complex<double> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = 2*b.step();
        double* bp = (double*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<double>& A, 
        VectorView<std::complex<double> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = 2*b.step();
        double* bp = (double*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(dtrsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<float>& A, VectorView<float> b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        float* bp = (b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<float>& A, VectorView<float> b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        float* bp = (b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        std::complex<float>* bp = (b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ctrsv) (
                BLASRM BLASCH_LO,BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ctrsv) (
                BLASCM BLASCH_UP,BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ctrsv) (
                BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<std::complex<float> >& A,
        VectorView<std::complex<float> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = b.step();
        std::complex<float>* bp = (b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        if (A.iscm() && A.isconj()) {
#ifdef CBLAS
            BLASNAME(ctrsv) (
                BLASRM BLASCH_UP,BLASCH_CT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
#else
            b.conjugateSelf();
            BLASNAME(ctrsv) (
                BLASCM BLASCH_LO,BLASCH_NT,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
            b.conjugateSelf();
#endif
        } else {
            BLASNAME(ctrsv) (
                BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
                A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
                A.isunit()?BLASCH_U:BLASCH_NU,
                BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
                BLAS1 BLAS1 BLAS1);
        }
    }
    template <> 
    void BlasTriLDivEq(
        const GenUpperTriMatrix<float>& A,
        VectorView<std::complex<float> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = 2*b.step();
        float* bp = (float*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
    template <> 
    void BlasTriLDivEq(
        const GenLowerTriMatrix<float>& A, 
        VectorView<std::complex<float> > b)
    {
        int n=A.size();
        int lda = A.isrm()?A.stepi():A.stepj();
        int bs = 2*b.step();
        float* bp = (float*)(b.ptr());
        if (bs < 0) bp += (n-1)*bs;
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
        BLASNAME(strsv) (
            BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
            A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
            BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(bp+1),BLASV(bs)
            BLAS1 BLAS1 BLAS1);
    }
#endif // FLOAT
#endif // BLAS

    template <class T, class Ta> 
    void TriLDivEq(
        const GenUpperTriMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.size() == A.size());
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Vector<T> b0(b);
#endif
        if (b.size() > 0) {
            if (b.isconj()) TriLDivEq(A.conjugate(),b.conjugate());
            else 
#ifdef BLAS
                if (isComplex(T()) && isReal(Ta()))
                    BlasTriLDivEq(A,b);
                else if ( !((A.isrm() && A.stepi()>0) || 
                            (A.iscm() && A.stepj()>0)) ) {
                    if (A.isunit()) {
                        UpperTriMatrix<Ta,UnitDiag|ColMajor> AA = A;
                        BlasTriLDivEq(AA,b);
                    } else {
                        UpperTriMatrix<Ta,NonUnitDiag|ColMajor> AA = A;
                        BlasTriLDivEq(AA,b);
                    }
                }
                else BlasTriLDivEq(A,b);
#else
            NonBlasTriLDivEq(A,b);
#endif
        }
#ifdef XDEBUG
        Vector<T> b2 = A0*b;
        if (Norm(b2-b0) > 0.001*Norm(A0)*Norm(b0)) {
            cerr<<"TriLDivEq: v/Upper\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"b = "<<TMV_Text(b)<<"  "<<b0<<endl;
            cerr<<"Done: b = "<<b<<endl;
            cerr<<"A*b = "<<b2<<endl;
            abort();
        }
#endif
    }

    template <class T, class Ta> 
    void TriLDivEq(const GenLowerTriMatrix<Ta>& A, VectorView<T> b)
    {
        TMVAssert(b.size() == A.size());
#ifdef XDEBUG
        Matrix<Ta> A0(A);
        Vector<T> b0(b);
#endif
        if (b.size() > 0) {
            if (b.isconj()) TriLDivEq(A.conjugate(),b.conjugate());
            else {
#ifdef BLAS
                if (isComplex(T()) && isReal(Ta()))
                    BlasTriLDivEq(A,b);
                else if ( !((A.isrm() && A.stepi()>0) ||
                            (A.iscm() && A.stepj()>0)) ) {
                    if (A.isunit()) {
                        LowerTriMatrix<Ta,UnitDiag|ColMajor> AA = A;
                        BlasTriLDivEq(AA,b);
                    } else {
                        LowerTriMatrix<Ta,NonUnitDiag|ColMajor> AA = A;
                        BlasTriLDivEq(AA,b);
                    }
                } else {
                    BlasTriLDivEq(A,b);
                }
#else
                NonBlasTriLDivEq(A,b);
#endif
            }
        }
#ifdef XDEBUG
        Vector<T> b2 = A0*b;
        if (Norm(b2-b0) > 0.001*Norm(A0)*Norm(b0)) {
            cerr<<"TriLDivEq: v/Lower\n";
            cerr<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
            cerr<<"b = "<<TMV_Text(b)<<"  "<<b0<<endl;
            cerr<<"Done: b = "<<b<<endl;
            cerr<<"A*b = "<<b2<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_TriDiv_V.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


