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


//---------------------------------------------------------------------------
//
// This file defines the TMV SmallMatrix class.
//
// Constructors:
//
//    SmallMatrix<T,M,N,A>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,A>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,A>(const GenMatrix<T>& m)
//        Make a SmallMatrix which copies the elements of m.
//
// SmallMatrix doesn't have views like a regular Matrix.
// All the normal viewing kinds of routines just return a regular MatrixView.
// It is mostly useful for fast element access and simple tasks
// like multiplication and addition.  For most division routines,
// SmallMatrix just sends the task to a regular matrix.  The exception
// is 2x2, which has specialization for det, inverse, etc.


#ifndef TMV_SmallMatrix_H
#define TMV_SmallMatrix_H

#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_Array.h"
#include "tmv/TMV_IOStyle.h"

namespace tmv {

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline void Copy(
        const SmallMatrix<T1,M,N,A1>& m1, SmallMatrix<T2,M,N,A2>& m2);

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline void Copy(
        const SmallMatrix<std::complex<T>,M,N,A1>& , SmallMatrix<T,M,N,A2>& )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotXm_1;

    template <typename T, ptrdiff_t M, ptrdiff_t N>
    class SmallMatrixComposite;

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline T DoDet(const SmallMatrix<T,M,N,A>& m);

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A2>
    inline void DoInverse(
        const SmallMatrix<T,M,N,A>& m, SmallMatrix<T2,N,M,A2>& minv);

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A2>
    inline void DoInverseATA(
        const SmallMatrix<T,M,N,A>& m, SmallMatrix<T,N,N,A2>& ata);

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class SmallMatrix
    {
    public:

        enum { S = A & AllStorageType };
        enum { I = A & FortranStyle };
        enum { Si = (S == int(RowMajor) ? N : 1) };
        enum { Sj = (S == int(RowMajor) ? 1 : M) };

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef SmallMatrix<T,M,N,A> type;
        typedef type copy_type;
        typedef ConstMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef ConstMatrixView<RT,I> const_realpart_type;
        typedef MatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef MatrixView<RT,I> realpart_type;
        typedef T& reference;
        typedef typename MatrixIterHelper<S,type>::rowmajor_iterator
            rowmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::const_rowmajor_iterator
            const_rowmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::colmajor_iterator
            colmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::const_colmajor_iterator
            const_colmajor_iterator;


        //
        // Constructors
        //

        inline SmallMatrix()
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        explicit inline SmallMatrix(const T& x)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            if (x == T(0)) setZero();
            else setAllTo(x);
        }

        inline SmallMatrix(const type& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            for(ptrdiff_t i=0;i<M*N;++i) itsm[i] = m2.itsm[i];
        }

        template <typename T2, int A2>
        inline SmallMatrix(const SmallMatrix<T2,M,N,A2>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
        }

        inline SmallMatrix(const GenMatrix<T>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        template <typename T2>
        inline SmallMatrix(const GenMatrix<T2>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const AssignableToMatrix<RT>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const AssignableToMatrix<CT>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const SmallMatrixComposite<RT,M,N>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            m2.assignTom(*this);
        }

        inline SmallMatrix(const SmallMatrixComposite<CT,M,N>& m2)
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(Attrib<A>::matrixok);
            m2.assignTom(*this);
        }

        inline ~SmallMatrix()
        {
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            if (&m2 != this) Copy(m2,*this);
            return *this;
        }

        template <typename T2, int A2>
        inline type& operator=(const SmallMatrix<T2,M,N,A2>& m2)
        {
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
            return *this;
        }

        inline type& operator=(const T& x)
        { return setToIdentity(x); }

        template <typename T2>
        inline type& operator=(const GenMatrix<T2>& m2)
        {
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            TMVAssert(isComplex(T()) || isReal(T2()));
            view() = m2;
            return *this;
        }

        inline type& operator=(const AssignableToMatrix<RT>& m2)
        {
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
            return *this;
        }

        inline type& operator=(const AssignableToMatrix<CT>& m2)
        {
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            TMVAssert(isComplex(T()));
            view() = m2;
            return *this;
        }

        inline type& operator=(const SmallMatrixComposite<RT,M,N>& m2)
        {
            m2.assignTom(*this);
            return *this;
        }

        inline type& operator=(const SmallMatrixComposite<CT,M,N>& m2)
        {
            m2.assignTom(*this);
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(rowmajor_begin(),M*N,x); }

        //
        // Access
        //

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            if (I == int(CStyle)) {
                TMVAssert(i>=0 && i<M);
                TMVAssert(j>=0 && j<N);
                return cref(i,j);
            } else {
                TMVAssert(i>=1 && i<=M);
                TMVAssert(j>=1 && j<=N);
                return cref(i-1,j-1);
            }
        }

        inline T& operator()(ptrdiff_t i,ptrdiff_t j)
        {
            if (I == int(CStyle)) {
                TMVAssert(i>=0 && i<M);
                TMVAssert(j>=0 && j<N);
                return ref(i,j);
            } else {
                TMVAssert(i>=1 && i<=M);
                TMVAssert(j>=1 && j<=N);
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(ptrdiff_t i) const
        { return view().row(i); }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return view().row(i,j1,j2); }

        inline const_vec_type operator[](ptrdiff_t i) const
        { return view().row(i); }

        inline const_vec_type col(ptrdiff_t j) const
        { return view().col(j); }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        { return view().col(j,i1,i2); }

        inline const_vec_type diag() const
        { return view().diag(); }

        inline const_vec_type diag(ptrdiff_t i) const
        { return view().diag(i); }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return view().diag(i,j1,j2); }

        inline vec_type row(ptrdiff_t i)
        { return view().row(i); }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        { return view().row(i,j1,j2); }

        inline vec_type operator[](ptrdiff_t i)
        { return view().row(i); }

        inline vec_type col(ptrdiff_t j)
        { return view().col(j); }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        { return view().col(j,i1,i2); }

        inline vec_type diag()
        { return view().diag(); }

        inline vec_type diag(ptrdiff_t i)
        { return view().diag(i); }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        { return view().diag(i,j1,j2); }

        //
        // Functions of Matrix
        //

        inline T trace() const
        {
            TMVAssert(M == N);
            T sum(0);
            for(ptrdiff_t i=0; i<M*N; i+=M+1) sum += itsm[i];
            return sum;
        }

        inline T sumElements() const
        {
            T sum(0);
            for(ptrdiff_t i=0;i<M*N; ++i) sum += itsm[i];
            return sum;
        }

        inline RT sumAbsElements() const
        {
            RT sum(0);
            for(ptrdiff_t i=0;i<M*N; ++i) sum += TMV_ABS(itsm[i]);
            return sum;
        }

        inline RT sumAbs2Elements() const
        {
            RT sum(0);
            for(ptrdiff_t i=0;i<M*N; ++i) sum += TMV_ABS2(itsm[i]);
            return sum;
        }

        inline RT norm() const
        { return normF(); }

        inline RT normF() const
        { return TMV_SQRT(normSq()); }

        // normF()^2
        inline RT normSq(RT scale=RT(1)) const
        {
            RT sum(0);
            if (scale == RT(1))
                for(ptrdiff_t i=0;i<M*N; ++i) sum += TMV_NORM(itsm[i]);
            else
                for(ptrdiff_t i=0;i<M*N; ++i) sum += TMV_NORM(itsm[i]*scale);
            return sum;
        }

        // 1-norm = max_j (sum_i |a_ij|)
        inline RT norm1() const
        {
            RT max(0);
            for(ptrdiff_t j=0;j<N;++j) {
                RT temp(0);
                for(ptrdiff_t i=0;i<M;++i) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // inf-norm = max_i (sum_j |a_ij|)
        inline RT normInf() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<M;++i) {
                RT temp(0);
                for(ptrdiff_t j=0;j<N;++j) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|a_ij|)
        inline RT maxAbsElement() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<M*N; ++i) {
                RT temp = TMV_ABS(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        inline RT maxAbs2Element() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<M*N; ++i) {
                RT temp = TMV_ABS2(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        inline T det() const
        {
            TMVAssert(M == N);
            return DoDet(*this);
        }
        inline RT logDet(T* sign=0) const
        {
            TMVAssert(M==N);
            if (M <= 3) {
                T d = det();
                RT absd = TMV_ABS(d);
                if (sign) {
                    if (isReal(T())) *sign = TMV_REAL(d) > 0 ?
                        RT(1) : RT(-1);
                    else *sign = d / absd;
                }
                return TMV_LOG(absd);
            }
            else return view().logDet(sign);
        }
        inline RT norm2() const
        { return view().doNorm2(); }
        inline RT doNorm2() const
        { return view().doNorm2(); }
        inline bool isSingular() const
        { return view().isSingular(); }
        inline RT condition() const
        { return view().doCondition(); }
        inline RT doCondition() const
        { return view().doCondition(); }

        //
        // Division Control
        //

        inline QuotXm_1<T,T,M,N,A> inverse() const
        { return QuotXm_1<T,T,M,N,A>(T(1),*this); }

        inline void makeInverse(MatrixView<T> minv) const
        { view().makeInverse(minv); }

        template <typename T1>
        inline void makeInverse(MatrixView<T1> minv) const
        { view().makeInverse(minv); }

        template <typename T2, int A2>
        inline void makeInverse(Matrix<T2,A2>& minv) const
        { view().makeInverse(minv); }

        inline void makeInverseATA(MatrixView<T> ata) const
        { view().makeInverseATA(ata); }

        template <int A2>
        inline void makeInverseATA(Matrix<T,A2>& ata) const
        { view().makeInverseATA(ata); }

        template <typename T2, int A2>
        inline void makeInverse(SmallMatrix<T2,N,M,A2>& minv) const
        { DoInverse(*this,minv); }

        template <typename T2, int A2>
        inline void makeInverseATA(SmallMatrix<T2,N,N,A2>& ata) const
        { DoInverseATA(*this,ata); }

        //
        // Modifying Functions
        //

        inline type& setZero()
        { for(ptrdiff_t i=0;i<M*N;++i) itsm[i] = T(0); return *this; }

        inline type& clip(RT thresh)
        {
            for(ptrdiff_t i=0; i<M*N; ++i)
                if (TMV_ABS(itsm[i]) < thresh) itsm[i] = T(0);
            return *this;
        }

        inline type& setAllTo(const T& x)
        {
            for(ptrdiff_t i=0; i<M*N; ++i) itsm[i] = x;
            return *this;
        }

        inline type& addToAll(const T& x)
        {
            for(ptrdiff_t i=0; i<M*N; ++i) itsm[i] += x;
            return *this;
        }

        inline type& transposeSelf()
        {
            TMVAssert(M == N);
            for(ptrdiff_t i=1; i<M; ++i)
                for(ptrdiff_t j=0; j<i; ++j) TMV_SWAP(ref(i,j),ref(j,i));
            return *this;
        }

        inline type& conjugateSelf()
        {
            if (isComplex(T())) {
                RT* itsmi = reinterpret_cast<RT*>(ptr())+1;
                for(ptrdiff_t i=0;i<2*M*N;i+=2) itsmi[i] = -itsmi[i];
            }
            return *this;
        }

        inline type& setToIdentity(const T& x=T(1))
        {
            TMVAssert(M == N);
            setZero();
            for(ptrdiff_t i=0; i<M*N; i+=M+1) itsm[i] = x;
            return *this;
        }

        inline type& swapRows(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(CStyle)) {
                TMVAssert(i1 >= 0 && i1 < M);
                TMVAssert(i2 >= 0 && i2 < M);
                if (i1 != i2)
                    for(ptrdiff_t j=0; j<N; ++j) TMV_SWAP(ref(i1,j),ref(i2,j));
            } else {
                TMVAssert(i1 >= 1 && i1 <= M);
                TMVAssert(i2 >= 1 && i2 <= M);
                if (i1 != i2)
                    for(ptrdiff_t j=0; j<N; ++j) TMV_SWAP(ref(i1-1,j),ref(i2-1,j));
            }
            return *this;
        }

        inline type& swapCols(ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(CStyle)) {
                TMVAssert(j1 >= 0 && j1 < N);
                TMVAssert(j2 >= 0 && j2 < N);
                if (j1 != j2)
                    for(ptrdiff_t i=0; i<M; ++i) TMV_SWAP(ref(i,j1),ref(i,j2));
            } else {
                TMVAssert(j1 >= 1 && j1 <= N);
                TMVAssert(j2 >= 1 && j2 <= N);
                if (j1 != j2)
                    for(ptrdiff_t i=0; i<M; ++i) TMV_SWAP(ref(i,j1-1),ref(i,j2-1));
            }
            return *this;
        }

        inline type& permuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(CStyle)) {
                TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
                for(ptrdiff_t i=i1;i<i2;++i) swapRows(i,p[i]);
            } else {
                TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
                for(ptrdiff_t i=i1-1;i<i2;++i) swapRows(i,p[i]);
            }
            return *this;
        }

        inline type& permuteRows(const ptrdiff_t* p)
        { permuteRows(p,I==int(CStyle)?0:1,M); return *this; }

        inline type& permuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(CStyle)) {
                TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
                for(ptrdiff_t j=j1;j<j2;++j) swapCols(j,p[j]);
            } else {
                TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
                for(ptrdiff_t j=j1-1;j<j2;++j) swapCols(j,p[j]);
            }
            return *this;
        }

        inline type& permuteCols(const ptrdiff_t* p)
        { permuteCols(p,I==int(CStyle)?0:1,N); return *this; }

        inline type& reversePermuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(CStyle)) {
                TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
                for(ptrdiff_t i=i2-1;i>=i1;--i) swapRows(i,p[i]);
            } else {
                TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
                for(ptrdiff_t i=i2-1;i>=i1-1;--i) swapRows(i,p[i]);
            }
            return *this;
        }

        inline type& reversePermuteRows(const ptrdiff_t* p)
        { reversePermuteRows(p,I==int(CStyle)?0:1,M); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(CStyle)) {
                TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
                for(ptrdiff_t j=j2-1;j>=j1;--j) swapCols(j,p[j]);
            } else {
                TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
                for(ptrdiff_t j=j2-1;j>=j1-1;--j) swapCols(j,p[j]);
            }
            return *this;
        }

        inline type& reversePermuteCols(const ptrdiff_t* p)
        { reversePermuteCols(p,I==int(CStyle)?0:1,N); return *this; }

        //
        // SubMatrix
        //

        inline const_view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return view().cSubMatrix(i1,i2,j1,j2); }

        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return view().subMatrix(i1,i2,j1,j2); }

        inline const_view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return view().cSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return view().subMatrix(i1,i2,j1,j2,istep,jstep); }

        inline const_vec_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        { return view().cSubVector(i,j,istep,jstep,s); }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        { return view().subVector(i,j,istep,jstep,s); }

        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        { return view().colPair(j1,j2); }

        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        { return view().rowPair(i1,i2); }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        { return view().colRange(j1,j2); }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        { return view().rowRange(i1,i2); }

        inline const_realpart_type realPart() const
        { return view().realPart(); }

        inline const_realpart_type imagPart() const
        { return view().imagPart(); }

        inline view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        { return view().cSubMatrix(i1,i2,j1,j2); }

        inline view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        { return view().subMatrix(i1,i2,j1,j2); }

        inline view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        { return view().cSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        { return view().subMatrix(i1,i2,j1,j2,istep,jstep); }

        inline vec_type cSubVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s)
        { return view().cSubVector(i,j,istep,jstep,s); }

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s)
        { return view().subVector(i,j,istep,jstep,s); }

        inline view_type colPair(ptrdiff_t j1, ptrdiff_t j2)
        { return view().colPair(j1,j2); }

        inline view_type rowPair(ptrdiff_t i1, ptrdiff_t i2)
        { return view().rowPair(i1,i2); }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        { return view().colRange(j1,j2); }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        { return view().rowRange(i1,i2); }

        inline realpart_type realPart()
        { return view().realPart(); }

        inline realpart_type imagPart()
        { return view().imagPart(); }

        //
        // Views
        //

        inline const_view_type view() const
        { return const_view_type(cptr(),M,N,Si,Sj,NonConj,M*N); }

        inline const_view_type transpose() const
        {
            return const_view_type(cptr(),N,M,Sj,Si,NonConj,M*N);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),M,N,Si,Sj,isReal(T())?NonConj:Conj,M*N);
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),N,M,Sj,Si,isReal(T())?NonConj:Conj,M*N);
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return const_uppertri_type(cptr(),N,Si,Sj,dt,NonConj); }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return const_lowertri_type(cptr(),M,Si,Sj,dt,NonConj); }

        inline const_vec_type constLinearView() const
        { return const_vec_type(cptr(),M*N,1,NonConj); }

        inline view_type view()
        { return view_type(ptr(),M,N,Si,Sj,NonConj,M*N); }

        inline view_type transpose()
        { return view_type(ptr(),N,M,Sj,Si,NonConj,M*N); }

        inline view_type conjugate()
        { return view_type(ptr(),M,N,Si,Sj,isReal(T())?NonConj:Conj,M*N); }

        inline view_type adjoint()
        { return view_type(ptr(),N,M,Sj,Si,isReal(T())?NonConj:Conj,M*N); }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag)
        { return uppertri_type(ptr(),N,Si,Sj,dt,NonConj); }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        { return lowertri_type(ptr(),M,Si,Sj,dt,NonConj); }

        inline vec_type linearView()
        { return vec_type(ptr(),M*N,1,NonConj); }


        inline ptrdiff_t colsize() const { return M; }
        inline ptrdiff_t rowsize() const { return N; }
        inline bool isSquare() const { return M == N; }
        inline ptrdiff_t stepi() const { return Si; }
        inline ptrdiff_t stepj() const { return Sj; }
        inline bool isrm() const { return S == int(RowMajor); }
        inline bool iscm() const { return S == int(ColMajor); }
        inline bool isconj() const { return false; }
        inline ConjType ct() const { return NonConj; }
        inline ptrdiff_t ls() const { return M*N; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(ptrdiff_t i, ptrdiff_t j) const
        { return S == int(RowMajor) ? itsm[i*N+j] : itsm[j*M+i]; }

        inline T& ref(ptrdiff_t i, ptrdiff_t j)
        { return S == int(RowMajor) ? itsm[i*N+j] : itsm[j*M+i]; }

        inline ptrdiff_t rowstart(ptrdiff_t ) const { return 0; }
        inline ptrdiff_t rowend(ptrdiff_t ) const { return N; }
        inline ptrdiff_t colstart(ptrdiff_t ) const { return 0; }
        inline ptrdiff_t colend(ptrdiff_t ) const { return M; }

        inline rowmajor_iterator rowmajor_begin()
        { return MatrixIterHelper<S,type>::rowmajor_begin(this); }
        inline rowmajor_iterator rowmajor_end()
        { return MatrixIterHelper<S,type>::rowmajor_end(this); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return MatrixIterHelper<S,type>::rowmajor_begin(this); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return MatrixIterHelper<S,type>::rowmajor_end(this); }

        inline colmajor_iterator colmajor_begin()
        { return MatrixIterHelper<S,type>::colmajor_begin(this); }
        inline colmajor_iterator colmajor_end()
        { return MatrixIterHelper<S,type>::colmajor_end(this); }

        inline const_colmajor_iterator colmajor_begin() const
        { return MatrixIterHelper<S,type>::colmajor_begin(this); }
        inline const_colmajor_iterator colmajor_end() const
        { return MatrixIterHelper<S,type>::colmajor_end(this); }


    protected :

        StackArray<T,M*N> itsm;

    }; // SmallMatrix

    //-------------------------------------------------------------------------

    //
    // Copy Matrices
    //

    template <ptrdiff_t M, ptrdiff_t N, StorageType S1, StorageType S2, typename T1, typename T2>
    struct DoCopym {};
    template <ptrdiff_t M, ptrdiff_t N, StorageType S, typename T1, typename T2>
    struct DoCopym<M,N,S,S,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(ptrdiff_t i=0;i<M*N;++i) m2[i] = m1[i]; }
    };
    template <ptrdiff_t M, ptrdiff_t N, typename T1, typename T2>
    struct DoCopym<M,N,ColMajor,RowMajor,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(ptrdiff_t i=0;i<M;++i) for(ptrdiff_t j=0;j<N;++j) m2[i*N+j] = m1[j*M+i]; }
    };
    template <ptrdiff_t M, ptrdiff_t N, typename T1, typename T2>
    struct DoCopym<M,N,RowMajor,ColMajor,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(ptrdiff_t i=0;i<M;++i) for(ptrdiff_t j=0;j<N;++j) m2[j*M+i] = m1[i*N+j]; }
    };

    template <ptrdiff_t M, ptrdiff_t N, StorageType S, typename T>
    struct DoCopym<M,N,S,S,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };
    template <ptrdiff_t M, ptrdiff_t N, typename T>
    struct DoCopym<M,N,ColMajor,RowMajor,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };
    template <ptrdiff_t M, ptrdiff_t N, typename T>
    struct DoCopym<M,N,RowMajor,ColMajor,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };

    template <ptrdiff_t M, ptrdiff_t N, StorageType S1, StorageType S2, typename T1, typename T2>
    inline void SmallMatrixCopy(const T1* m1, T2* m2)
    { DoCopym<M,N,S1,S2,T1,T2>(m1,m2); }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline void Copy(
        const SmallMatrix<T1,M,N,A1>& m1, SmallMatrix<T2,M,N,A2>& m2)
    {
        TMVAssert(isComplex(T2()) || isReal(T1()));
        const StorageType S1 = static_cast<StorageType>(A1 & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
        SmallMatrixCopy<M,N,S1,S2>(m1.cptr(),m2.ptr());
    }

    //
    // Swap Matrices
    //

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline void Swap(
        SmallMatrix<T,M,N,A1>& m1, SmallMatrix<T,M,N,A2>& m2)
    {
        if ((A1&AllStorageType)==(A2&AllStorageType))
            for(ptrdiff_t i=0;i<M*N;++i) TMV_SWAP(m1.ptr()[i],m2.ptr()[i]);
        else
            for(ptrdiff_t i=0;i<M;++i) for(ptrdiff_t j=0;j<N;++j)
                TMV_SWAP(m1.ref(i,j),m2.ref(i,j));
    }

    //
    // Functions of Matrices:
    //

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline T Det(const SmallMatrix<T,M,N,A>& m)
    { return m.det(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) LogDet(const SmallMatrix<T,M,N,A>& m)
    { return m.logDet(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline T Trace(const SmallMatrix<T,M,N,A>& m)
    { return m.trace(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) Norm(const SmallMatrix<T,M,N,A>& m)
    { return m.norm(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) NormSq(const SmallMatrix<T,M,N,A>& m)
    { return m.normSq(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) NormF(const SmallMatrix<T,M,N,A>& m)
    { return m.normF(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) Norm1(const SmallMatrix<T,M,N,A>& m)
    { return m.norm1(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) Norm2(const SmallMatrix<T,M,N,A>& m)
    { return m.norm2(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) NormInf(const SmallMatrix<T,M,N,A>& m)
    { return m.normInf(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) MaxAbsElement(const SmallMatrix<T,M,N,A>& m)
    { return m.maxAbsElement(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline TMV_RealType(T) MaxAbs2Element(const SmallMatrix<T,M,N,A>& m)
    { return m.maxAbs2Element(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline ConstMatrixView<T,A&FortranStyle> Transpose(
        const SmallMatrix<T,M,N,A>& m)
    { return m.transpose(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline MatrixView<T,A&FortranStyle> Transpose(SmallMatrix<T,M,N,A>& m)
    { return m.transpose(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline ConstMatrixView<T,A&FortranStyle> Conjugate(
        const SmallMatrix<T,M,N,A>& m)
    { return m.conjugate(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline MatrixView<T,A&FortranStyle> Conjugate(SmallMatrix<T,M,N,A>& m)
    { return m.conjugate(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline ConstMatrixView<T,A&FortranStyle> Adjoint(
        const SmallMatrix<T,M,N,A>& m)
    { return m.adjoint(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline MatrixView<T,A&FortranStyle> Adjoint(SmallMatrix<T,M,N,A>& m)
    { return m.adjoint(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline QuotXm_1<T,T,M,N,A> Inverse(const SmallMatrix<T,M,N,A>& m)
    { return m.inverse(); }

    //
    // Matrix ==, != Matrix
    //

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline bool operator==(
        const SmallMatrix<T1,M,N,A1>& m1, const SmallMatrix<T2,M,N,A2>& m2)
    {
        if ((A1&AllStorageType)==(A2&AllStorageType))
            for(ptrdiff_t i=0;i<M*N;++i) {
                if (m1.cptr()[i] != m2.cptr()[i]) return false;
            }
        else
            for(ptrdiff_t i=0;i<M;++i) for(ptrdiff_t j=0;j<N;++j) {
                if (m1(i,j) != m2(i,j)) return false;
            }
        return true;
    }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline bool operator!=(
        const SmallMatrix<T1,M,N,A1>& m1, const SmallMatrix<T2,M,N,A2>& m2)
    { return !(m1 == m2); }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline bool operator==(
        const GenMatrix<T1>& m1, const SmallMatrix<T2,M,N,A>& m2)
    { return m1 == m2.view(); }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline bool operator==(
        const SmallMatrix<T1,M,N,A>& m1, const GenMatrix<T2>& m2)
    { return m1.view() == m2; }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline bool operator!=(
        const GenMatrix<T1>& m1, const SmallMatrix<T2,M,N,A>& m2)
    { return m1 != m2.view(); }

    template <typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline bool operator!=(
        const SmallMatrix<T1,M,N,A>& m1, const GenMatrix<T2>& m2)
    { return m1.view() != m2; }



    //
    // I/O
    //

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const SmallMatrix<T,M,N,A>& m)
    { m.view().write(writer); return writer.getos(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline std::ostream& operator<<(
        std::ostream& os, const SmallMatrix<T,M,N,A>& m)
    { return os << IOStyle() << m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, SmallMatrix<T,M,N,A>& m)
    { m.view().read(reader); return reader.getis(); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline std::istream& operator>>(
        std::istream& is, SmallMatrix<T,M,N,A>& m)
    { return is >> IOStyle() >> m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline std::string TMV_Text(const SmallMatrix<T,M,N,A>& )
    {
        std::ostringstream s;
        s << std::string("SmallMatrix<")<<TMV_Text(T())<<','<<M<<','<<N;
        s <<','<<Attrib<A>::text()<<'>';
        return s.str();
    }

} // namespace tmv

#endif
