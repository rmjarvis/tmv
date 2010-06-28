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


//---------------------------------------------------------------------------
//
// This file defines the TMV SmallMatrix class.
//
// Constructors:
//
//    SmallMatrix<T,M,N,stor,I>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,stor,I>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,stor,I>(const vector<vector<T> >& m)
//        Makes a SmallMatrix with a_ij = m[i][j]
//
//    SmallMatrix<T,M,N,stor,I>(const T* m)
//    SmallMatrix<T,M,N,stor,I>(const vector<T>& m)
//    SmallMatrix<T,M,N,stor,I>(const GenMatrix<T>& m)
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

namespace tmv {

    template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Copy(
        const SmallMatrix<T1,M,N,S1,I1>& m1,
        SmallMatrix<T2,M,N,S2,I2>& m2);

    template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Copy(
        const SmallMatrix<std::complex<T>,M,N,S1,I1>& ,
        SmallMatrix<T,M,N,S2,I2>& )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Tm, int M, int N, StorageType S, IndexStyle I> 
    class QuotXm_1;

    template <class T, int M, int N> 
    class SmallMatrixComposite;

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline T DoDet(const SmallMatrix<T,M,N,S,I>& m);

    template <class T, class T2, int M, int N, StorageType S, StorageType S2, IndexStyle I, IndexStyle I2> 
    inline void DoInverse(
        const SmallMatrix<T,M,N,S,I>& m, SmallMatrix<T2,N,M,S2,I2>& minv);

    template <class T, int M, int N, StorageType S, StorageType S2, IndexStyle I, IndexStyle I2> 
    inline void DoInverseATA(
        const SmallMatrix<T,M,N,S,I>& m, SmallMatrix<T,N,N,S2,I2>& ata);


    // defaults S=ColMajor and I=CStyle are set in TMV_BaseMatrix.h
    template <class T, int M, int N, StorageType S, IndexStyle I> 
    class SmallMatrix 
    {
    public:

        enum { Si = (S == RowMajor ? N : 1) };
        enum { Sj = (S == RowMajor ? 1 : M) };

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef SmallMatrix<T,M,N,S,I> type;
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

        //
        // Constructors
        //

        inline SmallMatrix() 
        {
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        explicit inline SmallMatrix(const T& x) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            if (x == T(0)) setZero();
            else setAllTo(x);
        }

        inline SmallMatrix(const T* vv) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            for(int i=0;i<M*N;++i) itsm[i] = vv[i];
        }

        inline SmallMatrix(const std::vector<T>& vv) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == M*N);
            for(int i=0;i<M*N;++i) itsm[i] = vv[i];
        }

        explicit inline SmallMatrix(const std::vector<std::vector<T> >& vv) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            for(int i=0;i<M;++i) {
                TMVAssert(vv[i].size() == N);
                for(int j=0;j<N;++j) ref(i,j) = vv[i][j];
            }
        }

        inline SmallMatrix(const type& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            for(int i=0;i<M*N;++i) itsm[i] = m2.itsm[i];
        }

        template <IndexStyle I2> 
        inline SmallMatrix(const SmallMatrix<T,M,N,S,I2>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            for(int i=0;i<M*N;++i) itsm[i] = m2.cptr()[i];
        }

        template <class T2, StorageType S2, IndexStyle I2> 
        inline SmallMatrix(const SmallMatrix<T2,M,N,S2,I2>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
        }

        inline SmallMatrix(const GenMatrix<T>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        template <class T2> 
        inline SmallMatrix(const GenMatrix<T2>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const AssignableToMatrix<RT>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const AssignableToMatrix<CT>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            view() = m2;
        }

        inline SmallMatrix(const SmallMatrixComposite<RT,M,N>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignTom(*this);
        }

        inline SmallMatrix(const SmallMatrixComposite<CT,M,N>& m2) 
        { 
            TMVAssert(M>0);
            TMVAssert(N>0);
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignTom(*this);
        }

        inline ~SmallMatrix()
        {
#ifdef TMVDEBUG
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

        template <IndexStyle I2> 
        inline type& operator=(const SmallMatrix<T,M,N,S,I2>& m2)
        { 
            Copy(m2,*this);
            return *this; 
        }

        template <class T2, StorageType S2, IndexStyle I2> 
        inline type& operator=(const SmallMatrix<T2,M,N,S2,I2>& m2)
        { 
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
            return *this; 
        }

        inline type& operator=(const T& x) 
        { return setToIdentity(x); }

        template <class T2> 
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

        typedef ListAssigner<T,VIt<T,Unit,NonConj> > MyListAssigner;
        TMV_DEPRECATED(inline MyListAssigner operator=(ListInitClass))
        { return MyListAssigner(VIt<T,Unit,NonConj>(ptr(),1),M*N); }

        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(VIt<T,Unit,NonConj>(ptr(),1),M*N,x); }


        //
        // Access
        //

        inline T operator()(int i,int j) const
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<M);
                TMVAssert(j>=0 && j<N);
                return cref(i,j); 
            } else {
                TMVAssert(i>=1 && i<=M);
                TMVAssert(j>=1 && j<=N);
                return cref(i-1,j-1); 
            }
        }

        inline T& operator()(int i,int j) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<M);
                TMVAssert(j>=0 && j<N);
                return ref(i,j);
            } else {
                TMVAssert(i>=1 && i<=M);
                TMVAssert(j>=1 && j<=N);
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(int i) const 
        { return view().row(i); }

        inline const_vec_type row(int i, int j1, int j2) const 
        { return view().row(i,j1,j2); }

        inline const_vec_type operator[](int i) const
        { return view().row(i); }

        inline const_vec_type col(int j) const
        { return view().col(j); }

        inline const_vec_type col(int j, int i1, int i2) const
        { return view().col(j,i1,i2); }

        inline const_vec_type diag() const
        { return view().diag(); }

        inline const_vec_type diag(int i) const
        { return view().diag(i); }

        inline const_vec_type diag(int i, int j1, int j2) const
        { return view().diag(i,j1,j2); }

        inline vec_type row(int i)
        { return view().row(i); }

        inline vec_type row(int i, int j1, int j2)
        { return view().row(i,j1,j2); }

        inline vec_type operator[](int i)
        { return view().row(i); }

        inline vec_type col(int j)
        { return view().col(j); }

        inline vec_type col(int j, int i1, int i2)
        { return view().col(j,i1,i2); }

        inline vec_type diag()
        { return view().diag(); }

        inline vec_type diag(int i)
        { return view().diag(i); }

        inline vec_type diag(int i, int j1, int j2) 
        { return view().diag(i,j1,j2); }

        //
        // Functions of Matrix
        //

        inline T trace() const
        { 
            TMVAssert(M == N);
            T sum(0);
            for(int i=0; i<M*N; i+=M+1) sum += itsm[i];
            return sum;
        }

        inline T sumElements() const
        {
            T sum(0);
            for(int i=0;i<M*N; ++i) sum += itsm[i];
            return sum;
        }

        inline RT sumAbsElements() const
        {
            RT sum(0);
            for(int i=0;i<M*N; ++i) sum += TMV_ABS(itsm[i]);
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
                for(int i=0;i<M*N; ++i) sum += TMV_NORM(itsm[i]);
            else
                for(int i=0;i<M*N; ++i) sum += TMV_NORM(itsm[i]*scale);
            return sum;
        }

        // 1-norm = max_j (sum_i |a_ij|)
        inline RT norm1() const
        {
            RT max(0);
            for(int j=0;j<N;++j) {
                RT temp(0);
                for(int i=0;i<M;++i) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // inf-norm = max_i (sum_j |a_ij|)
        inline RT normInf() const
        {
            RT max(0);
            for(int i=0;i<M;++i) {
                RT temp(0);
                for(int j=0;j<N;++j) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|a_ij|)
        inline RT maxAbsElement() const
        {
            RT max(0);
            for(int i=0;i<M*N; ++i) {
                RT temp = TMV_ABS(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        inline RT maxAbs2Element() const
        {
            RT max(0);
            for(int i=0;i<M*N; ++i) {
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

        TMV_DEPRECATED(T Det() const)
        { return det(); }
        TMV_DEPRECATED(RT LogDet(T* sign=0) const)
        { return logDet(sign); }
        TMV_DEPRECATED(T Trace() const)
        { return trace(); }
        TMV_DEPRECATED(RT Norm() const )
        { return norm(); }
        TMV_DEPRECATED(RT NormF() const)
        { return normF(); }
        TMV_DEPRECATED(RT NormSq(const RT scale = RT(1)) const)
        { return normSq(); }
        TMV_DEPRECATED(RT Norm1() const)
        { return norm1(); }
        TMV_DEPRECATED(RT NormInf() const)
        { return normInf(); }
        TMV_DEPRECATED(RT MaxAbsElement() const)
        { return maxAbsElement(); }
        TMV_DEPRECATED(bool Singular() const)
        { return isSingular(); }
        TMV_DEPRECATED(RT Norm2() const)
        { return norm2(); }
        TMV_DEPRECATED(RT DoNorm2() const)
        { return doNorm2(); }
        TMV_DEPRECATED(RT Condition() const)
        { return condition(); }
        TMV_DEPRECATED(RT DoCondition() const)
        { return doCondition(); }

        // 
        // Division Control
        //

        inline QuotXm_1<T,T,M,N,S,I> inverse() const
        { return QuotXm_1<T,T,M,N,S,I>(T(1),*this); }

        inline void makeInverse(const MatrixView<T>& minv) const
        { view().makeInverse(minv); }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        { view().makeInverse(minv); }

        template <class T2, StorageType S2, IndexStyle I2> 
        inline void makeInverse(Matrix<T2,S2,I2>& minv) const
        { view().makeInverse(minv); }

        inline void makeInverseATA(const MatrixView<T>& ata) const
        { view().makeInverseATA(ata); }

        template <StorageType S2, IndexStyle I2> 
        inline void makeInverseATA(Matrix<T,S2,I2>& ata) const
        { view().makeInverseATA(ata); }

        template <class T2, StorageType S2, IndexStyle I2> 
        inline void makeInverse(SmallMatrix<T2,N,M,S2,I2>& minv) const
        { DoInverse(*this,minv); }

        template <class T2, StorageType S2, IndexStyle I2> 
        inline void makeInverseATA(SmallMatrix<T2,N,N,S2,I2>& ata) const
        { DoInverseATA(*this,ata); }

        TMV_DEPRECATED(void Inverse(const MatrixView<T>& minv) const)
        { makeInverse(minv); }
        TMV_DEPRECATED(void InverseATA(const MatrixView<T>& ata) const)
        { makeInverseATA(ata); }
        typedef QuotXm_1<T,T,M,N,S,I> MyQuotXm;
        TMV_DEPRECATED(MyQuotXm Inverse() const)
        { return inverse(); }
        template <class T2>
        TMV_DEPRECATED(void Inverse(const MatrixView<T2>& minv) const);
        template <class T2, StorageType S2, IndexStyle I2>
        TMV_DEPRECATED(void Inverse(Matrix<T2,S2,I2>& minv) const);
        template <StorageType S2, IndexStyle I2>
        TMV_DEPRECATED(void InverseATA(Matrix<T,S2,I2>& ata) const);
        template <class T2, StorageType S2, IndexStyle I2> 
        TMV_DEPRECATED(void Inverse(SmallMatrix<T2,N,M,S2,I2>& minv) const);
        template <class T2, StorageType S2, IndexStyle I2> 
        TMV_DEPRECATED(void InverseATA(SmallMatrix<T2,N,N,S2,I2>& ata) const);


        //
        // Modifying Functions
        //

        inline type& setZero() 
        { for(int i=0;i<M*N;++i) itsm[i] = T(0); return *this; }

        inline type& clip(RT thresh)
        { 
            for(int i=0; i<M*N; ++i)
                if (TMV_ABS(itsm[i]) < thresh) itsm[i] = T(0);
            return *this;
        }

        inline type& setAllTo(const T& x) 
        {
            for(int i=0; i<M*N; ++i) itsm[i] = x;
            return *this;
        }

        inline type& transposeSelf() 
        {
            TMVAssert(M == N);
            for(int i=1; i<M; ++i) 
                for(int j=0; j<i; ++j) TMV_SWAP(ref(i,j),ref(j,i));
            return *this;
        }

        inline type& conjugateSelf() 
        {
            if (isComplex(T())) {
                RT* itsmi = reinterpret_cast<RT*>(ptr())+1;
                for(int i=0;i<2*M*N;i+=2) itsmi[i] = -itsmi[i];
            }
            return *this;
        }

        inline type& setToIdentity(const T& x=T(1)) 
        { 
            TMVAssert(M == N);
            setZero();
            for(int i=0; i<M*N; i+=M+1) itsm[i] = x;
            return *this;
        }

        inline type& swapRows(int i1, int i2)
        {
            if (I == CStyle) {
                TMVAssert(i1 >= 0 && i1 < M);
                TMVAssert(i2 >= 0 && i2 < M);
                if (i1 != i2)
                    for(int j=0; j<N; ++j) TMV_SWAP(ref(i1,j),ref(i2,j));
            } else {
                TMVAssert(i1 >= 1 && i1 <= M);
                TMVAssert(i2 >= 1 && i2 <= M);
                if (i1 != i2)
                    for(int j=0; j<N; ++j) TMV_SWAP(ref(i1-1,j),ref(i2-1,j));
            }
            return *this;
        }

        inline type& swapCols(int j1, int j2)
        {
            if (I == CStyle) {
                TMVAssert(j1 >= 0 && j1 < N);
                TMVAssert(j2 >= 0 && j2 < N);
                if (j1 != j2)
                    for(int i=0; i<M; ++i) TMV_SWAP(ref(i,j1),ref(i,j2));
            } else {
                TMVAssert(j1 >= 1 && j1 <= N);
                TMVAssert(j2 >= 1 && j2 <= N);
                if (j1 != j2)
                    for(int i=0; i<M; ++i) TMV_SWAP(ref(i,j1-1),ref(i,j2-1));
            }
            return *this;
        }

        inline type& permuteRows(const int* p, int i1, int i2)
        {
            if (I == CStyle) {
                TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
                for(int i=i1;i<i2;++i) swapRows(i,p[i]);
            } else {
                TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
                for(int i=i1-1;i<i2;++i) swapRows(i,p[i]);
            }
            return *this;
        }

        inline type& permuteRows(const int* p)
        { permuteRows(p,I==CStyle?0:1,M); return *this; }

        inline type& permuteCols(const int* p, int j1, int j2)
        {
            if (I == CStyle) {
                TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
                for(int j=j1;j<j2;++j) swapCols(j,p[j]);
            } else {
                TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
                for(int j=j1-1;j<j2;++j) swapCols(j,p[j]);
            }
            return *this;
        }

        inline type& permuteCols(const int* p)
        { permuteCols(p,I==CStyle?0:1,N); return *this; }

        inline type& reversePermuteRows(const int* p, int i1, int i2)
        {
            if (I == CStyle) {
                TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
                for(int i=i2-1;i>=i1;--i) swapRows(i,p[i]);
            } else {
                TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
                for(int i=i2-1;i>=i1-1;--i) swapRows(i,p[i]);
            }
            return *this;
        }

        inline type& reversePermuteRows(const int* p)
        { reversePermuteRows(p,I==CStyle?0:1,M); return *this; }

        inline type& reversePermuteCols(const int* p, int j1, int j2)
        {
            if (I == CStyle) {
                TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
                for(int j=j2-1;j>=j1;--j) swapCols(j,p[j]);
            } else {
                TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
                for(int j=j2-1;j>=j1-1;--j) swapCols(j,p[j]);
            }
            return *this;
        }

        inline type& reversePermuteCols(const int* p)
        { reversePermuteCols(p,I==CStyle?0:1,N); return *this; }

        TMV_DEPRECATED(type& Zero())
        { return setZero(); }
        TMV_DEPRECATED(type& SetAllTo(const T& x))
        { return setAllTo(x); }
        TMV_DEPRECATED(type& Clip(RT thresh))
        { return clip(thresh); }
        TMV_DEPRECATED(type& TransposeSelf())
        { return transposeSelf(); }
        TMV_DEPRECATED(type& ConjugateSelf())
        { return conjugateSelf(); }
        TMV_DEPRECATED(type& SetToIdentity(const T& x=T(1)))
        { return setToIdentity(x); }
        TMV_DEPRECATED(type& SwapRows(int i1, int i2))
        { return swapRows(i1,i2); }
        TMV_DEPRECATED(type& SwapCols(int j1, int j2))
        { return swapCols(j1,j2); }
        TMV_DEPRECATED(type& PermuteRows(const int* p, int i1, int i2))
        { return permuteRows(p,i1,i2); }
        TMV_DEPRECATED(type& PermuteRows(const int* p))
        { return permuteRows(p); }
        TMV_DEPRECATED(type& PermuteCols(const int* p, int j1, int j2))
        { return permuteCols(p,j1,j2); }
        TMV_DEPRECATED(type& PermuteCols(const int* p))
        { return permuteCols(p); }
        TMV_DEPRECATED(type& ReversePermuteRows(const int* p, int i1, int i2))
        { return reversePermuteRows(p,i1,i2); }
        TMV_DEPRECATED(type& ReversePermuteRows(const int* p))
        { return reversePermuteRows(p); }
        TMV_DEPRECATED(type& ReversePermuteCols(const int* p, int j1, int j2))
        { return reversePermuteCols(p,j1,j2); }
        TMV_DEPRECATED(type& ReversePermuteCols(const int* p))
        { return reversePermuteCols(p); }


        //
        // SubMatrix
        //

        inline const_view_type cSubMatrix(int i1, int i2, int j1, int j2) const
        { return view().cSubMatrix(i1,i2,j1,j2); }

        inline const_view_type subMatrix(int i1, int i2, int j1, int j2) const
        { return view().subMatrix(i1,i2,j1,j2); }

        inline const_view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return view().cSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline const_view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return view().subMatrix(i1,i2,j1,j2,istep,jstep); }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return view().cSubVector(i,j,istep,jstep,s); }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        { return view().subVector(i,j,istep,jstep,s); }

        inline const_view_type colPair(int j1, int j2) const
        { return view().colPair(j1,j2); }

        inline const_view_type rowPair(int i1, int i2) const
        { return view().rowPair(i1,i2); }

        inline const_view_type colRange(int j1, int j2) const
        { return view().colRange(j1,j2); }

        inline const_view_type rowRange(int i1, int i2) const
        { return view().rowRange(i1,i2); }

        inline const_realpart_type realPart() const
        { return view().realPart(); }

        inline const_realpart_type imagPart() const
        { return view().imagPart(); }

        inline view_type cSubMatrix(int i1, int i2, int j1, int j2)
        { return view().cSubMatrix(i1,i2,j1,j2); }

        inline view_type subMatrix(int i1, int i2, int j1, int j2)
        { return view().subMatrix(i1,i2,j1,j2); }

        inline view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        { return view().cSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        { return view().subMatrix(i1,i2,j1,j2,istep,jstep); }

        inline vec_type cSubVector(int i, int j, int istep, int jstep, int s) 
        { return view().cSubVector(i,j,istep,jstep,s); }

        inline vec_type subVector(int i, int j, int istep, int jstep, int s) 
        { return view().subVector(i,j,istep,jstep,s); }

        inline view_type colPair(int j1, int j2) 
        { return view().colPair(j1,j2); }

        inline view_type rowPair(int i1, int i2) 
        { return view().rowPair(i1,i2); }

        inline view_type colRange(int j1, int j2) 
        { return view().colRange(j1,j2); }

        inline view_type rowRange(int i1, int i2) 
        { return view().rowRange(i1,i2); }

        inline realpart_type realPart() 
        { return view().realPart(); }

        inline realpart_type imagPart() 
        { return view().imagPart(); }

        TMV_DEPRECATED(const_view_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(const_view_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(const_vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(const_uppertri_type UpperTri(
                DiagType dt=NonUnitDiag) const)
        { return upperTri(dt); }
        TMV_DEPRECATED(const_lowertri_type LowerTri(
                DiagType dt=NonUnitDiag) const)
        { return lowerTri(dt); }
        TMV_DEPRECATED(const_view_type ColPair(int j1, int j2) const)
        { return colPair(j1,j2); }
        TMV_DEPRECATED(const_view_type RowPair(int i1, int i2) const)
        { return rowPair(i1,i2); }
        TMV_DEPRECATED(const_view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(const_view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(view_type SubMatrix(int i1, int i2, int j1, int j2))
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(view_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep))
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s))
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(uppertri_type UpperTri(DiagType dt=NonUnitDiag))
        { return upperTri(); }
        TMV_DEPRECATED(lowertri_type LowerTri(DiagType dt=NonUnitDiag))
        { return lowerTri(); }
        TMV_DEPRECATED(view_type ColPair(int j1, int j2))
        { return colPair(j1,j2); }
        TMV_DEPRECATED(view_type RowPair(int i1, int i2))
        { return rowPair(i1,i2); }
        TMV_DEPRECATED(view_type Cols(int j1, int j2))
        { return colRange(j1,j2); }
        TMV_DEPRECATED(view_type Rows(int i1, int i2))
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(realpart_type Real())
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag())
        { return imagPart(); }


        //
        // Views
        //

        inline const_view_type view() const
        { return const_view_type(cptr(),M,N,Si,Sj,S,NonConj,M*N); }

        inline const_view_type transpose() const
        {
            return const_view_type(
                cptr(),N,M,Sj,Si,TMV_TransOf(S),NonConj,M*N); 
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),M,N,Si,Sj,S,isReal(T())?NonConj:Conj,M*N); 
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),N,M,Sj,Si,TMV_TransOf(S),
                isReal(T())?NonConj:Conj,M*N); 
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return const_uppertri_type(cptr(),N,Si,Sj,dt,S,NonConj); }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return const_lowertri_type(cptr(),M,Si,Sj,dt,S,NonConj); }

        inline const_vec_type constLinearView() const
        { return const_vec_type(cptr(),M*N,1,NonConj); }

        inline view_type view()
        { return view_type(ptr(),M,N,Si,Sj,S,NonConj,M*N); }

        inline view_type transpose() 
        {
            return view_type(ptr(),N,M,Sj,Si,
                           TMV_TransOf(S),NonConj,M*N); 
        }

        inline view_type conjugate()
        { 
            return view_type(
                ptr(),M,N,Si,Sj,S,isReal(T())?NonConj:Conj,M*N); 
        }

        inline view_type adjoint()
        { 
            return view_type(
                ptr(),N,M,Sj,Si,TMV_TransOf(S),
                isReal(T())?NonConj:Conj,M*N); 
        }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag) 
        { return uppertri_type(ptr(),N,Si,Sj,dt,S,NonConj); }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        { return lowertri_type(ptr(),M,Si,Sj,dt,S,NonConj); }

        inline vec_type linearView()
        { return vec_type(ptr(),M*N,1,NonConj); }

        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(const_vec_type ConstLinearView() const)
        { return constLinearView(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Transpose())
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint())
        { return adjoint(); }
        TMV_DEPRECATED(vec_type LinearView())
        { return linearView(); }


        //
        // I/O
        //

        inline void write(std::ostream& os) const
        { view().write(os); }

        inline void write(std::ostream& os, RT thresh) const
        { view().write(os,thresh); }

        inline size_t colsize() const { return M; }
        inline size_t rowsize() const { return N; }
        inline bool isSquare() const { return M == N; }
        TMV_DEPRECATED(bool IsSquare() const) { return isSquare(); }
        inline int stepi() const { return Si; }
        inline int stepj() const { return Sj; }
        inline bool isrm() const { return S == RowMajor; }
        inline bool iscm() const { return S == ColMajor; }
        inline bool isconj() const { return false; }
        inline StorageType stor() const { return isrm() ? RowMajor : ColMajor; }
        inline ConjType ct() const { return NonConj; }
        inline size_t ls() const { return M*N; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const
        { return S == RowMajor ? itsm[i*N+j] : itsm[j*M+i]; }

        inline T& ref(int i, int j)
        { return S == RowMajor ? itsm[i*N+j] : itsm[j*M+i]; }

    protected :

        StackArray<T,M*N> itsm;

    }; // SmallMatrix

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    template <class T2>
    inline void SmallMatrix<T,M,N,S,I>::Inverse(
        const MatrixView<T2>& minv) const
    { makeInverse(minv); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    template <class T2, StorageType S2, IndexStyle I2>
    inline void SmallMatrix<T,M,N,S,I>::Inverse(Matrix<T2,S2,I2>& minv) const
    { makeInverse(minv); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    template <StorageType S2, IndexStyle I2>
    inline void SmallMatrix<T,M,N,S,I>::InverseATA(Matrix<T,S2,I2>& ata) const
    { makeInverseATA(ata); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    template <class T2, StorageType S2, IndexStyle I2> 
    inline void SmallMatrix<T,M,N,S,I>::Inverse(
        SmallMatrix<T2,N,M,S2,I2>& minv) const
    { makeInverse(minv); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    template <class T2, StorageType S2, IndexStyle I2> 
    inline void SmallMatrix<T,M,N,S,I>::InverseATA(
        SmallMatrix<T2,N,N,S2,I2>& ata) const
    { makeInverseATA(ata); }

    //-------------------------------------------------------------------------

    //
    // Copy Matrices
    //

    template <int M, int N, StorageType S1, StorageType S2, class T1, class T2> 
    struct DoCopym {};
    template <int M, int N, StorageType S, class T1, class T2> 
    struct DoCopym<M,N,S,S,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(int i=0;i<M*N;++i) m2[i] = m1[i]; }
    };
    template <int M, int N, class T1, class T2> 
    struct DoCopym<M,N,ColMajor,RowMajor,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(int i=0;i<M;++i) for(int j=0;j<N;++j) m2[i*N+j] = m1[j*M+i]; }
    };
    template <int M, int N, class T1, class T2> 
    struct DoCopym<M,N,RowMajor,ColMajor,T1,T2>
    {
        inline DoCopym(const T1* m1, T2* m2)
        { for(int i=0;i<M;++i) for(int j=0;j<N;++j) m2[j*M+i] = m1[i*N+j]; }
    };

    template <int M, int N, StorageType S, class T>
    struct DoCopym<M,N,S,S,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };
    template <int M, int N, class T>
    struct DoCopym<M,N,ColMajor,RowMajor,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };
    template <int M, int N, class T>
    struct DoCopym<M,N,RowMajor,ColMajor,std::complex<T>,T>
    { inline DoCopym(const std::complex<T>*, T*) { TMVAssert(TMV_FALSE); } };

    template <int M, int N, StorageType S1, StorageType S2, class T1, class T2>
    inline void DoCopy(const T1* m1, T2* m2)
    { DoCopym<M,N,S1,S2,T1,T2>(m1,m2); }

    template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Copy(
        const SmallMatrix<T1,M,N,S1,I1>& m1, SmallMatrix<T2,M,N,S2,I2>& m2)
    { 
        TMVAssert(isComplex(T2()) || isReal(T1()));
        DoCopy<M,N,S1,S2>(m1.cptr(),m2.ptr()); 
    }

    //
    // Swap Matrices
    //

    template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
        SmallMatrix<T,M,N,S1,I1>& m1, SmallMatrix<T,M,N,S2,I2>& m2) 
    {
        if (S1==S2)
            for(int i=0;i<M*N;++i) TMV_SWAP(m1.ptr()[i],m2.ptr()[i]);
        else
            for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
                TMV_SWAP(m1.ref(i,j),m2.ref(i,j));
    }

    //
    // Functions of Matrices:
    //

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline T Det(const SmallMatrix<T,M,N,S,I>& m)
    { return m.det(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) LogDet(const SmallMatrix<T,M,N,S,I>& m)
    { return m.logDet(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline T Trace(const SmallMatrix<T,M,N,S,I>& m)
    { return m.trace(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) Norm(const SmallMatrix<T,M,N,S,I>& m)
    { return m.norm(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) NormSq(const SmallMatrix<T,M,N,S,I>& m)
    { return m.normSq(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) NormF(const SmallMatrix<T,M,N,S,I>& m)
    { return m.normF(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) Norm1(const SmallMatrix<T,M,N,S,I>& m)
    { return m.norm1(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) Norm2(const SmallMatrix<T,M,N,S,I>& m)
    { return m.norm2(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) NormInf(const SmallMatrix<T,M,N,S,I>& m)
    { return m.normInf(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) MaxAbsElement(const SmallMatrix<T,M,N,S,I>& m)
    { return m.maxAbsElement(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline TMV_RealType(T) MaxAbs2Element(const SmallMatrix<T,M,N,S,I>& m)
    { return m.maxAbs2Element(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Transpose(const SmallMatrix<T,M,N,S,I>& m)
    { return m.transpose(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Transpose(SmallMatrix<T,M,N,S,I>& m)
    { return m.transpose(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Conjugate(const SmallMatrix<T,M,N,S,I>& m)
    { return m.conjugate(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Conjugate(SmallMatrix<T,M,N,S,I>& m)
    { return m.conjugate(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Adjoint(const SmallMatrix<T,M,N,S,I>& m)
    { return m.adjoint(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Adjoint(SmallMatrix<T,M,N,S,I>& m)
    { return m.adjoint(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline QuotXm_1<T,T,M,N,S,I> Inverse(const SmallMatrix<T,M,N,S,I>& m)
    { return m.inverse(); }

    //
    // Matrix ==, != Matrix
    //

    template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline bool operator==(
        const SmallMatrix<T1,M,N,S1,I1>& m1, 
        const SmallMatrix<T2,M,N,S2,I2>& m2)
    { 
        if (S1==S2)
            for(int i=0;i<M*N;++i) {
                if (m1.cptr()[i] != m2.cptr()[i]) return false;
            }
        else
            for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
                if (m1(i,j) != m2(i,j)) return false;
            }
        return true;
    }

    template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline bool operator!=(
        const SmallMatrix<T1,M,N,S1,I1>& m1, 
        const SmallMatrix<T2,M,N,S2,I2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline std::ostream& operator<<(
        std::ostream& os, const SmallMatrix<T,M,N,S,I>& m)
    { m.write(os); return os; }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(
        std::istream& is, SmallMatrix<T,M,N,S,I>& m)
    { return is>>m.view(); }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    inline std::string TMV_Text(const SmallMatrix<T,M,N,S,I>& )
    { 
        std::ostringstream s;
        s << std::string("SmallMatrix<")<<TMV_Text(T())<<','<<M<<','<<N;
        s <<','<<TMV_Text(S)<<','<<TMV_Text(I)<<">"; 
        return s.str();
    }

} // namespace tmv

#endif
