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


#ifndef TMV_SmallMatrixArith_H
#define TMV_SmallMatrixArith_H

#include "tmv/TMV_SmallVectorArithFunc.h"
#include "tmv/TMV_SmallMatrixArithFunc.h"
#include "tmv/TMV_MatrixArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <typename T, ptrdiff_t M, ptrdiff_t N>
    class SmallMatrixComposite : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SmallMatrixComposite() {}
        inline SmallMatrixComposite(const SmallMatrixComposite<T,M,N>&) {}
        inline ~SmallMatrixComposite() {}

        inline ptrdiff_t colsize() const { return M; }
        inline ptrdiff_t rowsize() const { return N; }

        virtual void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m) const = 0;
        virtual void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m) const = 0;
        virtual void assignToM(MatrixView<real_type> m) const = 0;
        virtual void assignToM(MatrixView<complex_type> m) const = 0;
    };


    //
    // Scalar * Matrix
    //

    template <typename T, typename T1, ptrdiff_t M, ptrdiff_t N, int A>
    class ProdXm : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdXm(const T _x, const SmallMatrix<T1,M,N,A>& _m) :
            x(_x), m(_m) {}
        T getX() const { return x; }
        const SmallMatrix<T1,M,N,A>& getM() const { return m; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            m0=m;
            MultXV<M*N>(x,m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M && !SameStorage(m0,m)) {
                SmallMatrixCopy<M,N,S,ColMajor>(m.cptr(),m0.ptr());
                MultXV<M*N>(x,m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<M,N,S,RowMajor>(m.cptr(),m0.ptr());
                MultXV<M*N>(x,m0.ptr());
            } else {
                MultXM(x,m0=m);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M && !SameStorage(m0,m)) {
                SmallMatrixCopy<M,N,S,ColMajor>(m.cptr(),m0.ptr());
                MultXV<M*N>(x,m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<M,N,S,RowMajor>(m.cptr(),m0.ptr());
                MultXV<M*N>(x,m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                MultXM(x,m0=m);
            }
        }
    private:
        const T x;
        const SmallMatrix<T1,M,N,A>& m;
    };

    // m*=x
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator*=(SmallMatrix<T,M,N,A>& m, T x)
    { MultXV<M*N>(x,m.ptr()); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator*=(SmallMatrix<CT,M,N,A>& m, T x)
    { MultXV<M*N>(x,m.ptr()); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator*=(SmallMatrix<CT,M,N,A>& m, CCT x)
    { m *= CT(x); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator*=(SmallMatrix<CT,M,N,A>& m, VCT x)
    { m *= CT(x); return m; }

    // m/=x
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator/=(SmallMatrix<T,M,N,A>& m, T x)
    { m *= TMV_InverseOf(x); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator/=(SmallMatrix<CT,M,N,A>& m, T x)
    { m *= TMV_InverseOf(x); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator/=(SmallMatrix<CT,M,N,A>& m, CCT x)
    { m *= TMV_InverseOf(CT(x)); return m; }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator/=(SmallMatrix<CT,M,N,A>& m, VCT x)
    { m *= TMV_InverseOf(CT(x)); return m; }

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define X ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#include "tmv/TMV_AuxProdXM.h"
    // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef PRODXM
#undef GENMATRIX

    //
    // Matrix + Scalar
    //

    template <typename T, typename T1, ptrdiff_t N, int A>
    class SummX_1 : public SmallMatrixComposite<T,N,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        SummX_1(
            T TMV_DEBUGPARAM(_x1), const SmallMatrix<T1,N,N,A>& _m, T _x2) :
            m(_m), x2(_x2)
        { TMVAssert(_x1 == T(1)); }
        const SmallMatrix<T1,N,N,A>& getM() const
        { return m; }
        T getX2() const { return x2; }
        void assignTom(
            SmallMatrix<real_type,N,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<real_type,N,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<complex_type,N,N,RowMajor|CStyle>& m0) const
        { (m0=m) += x2; }
        void assignTom(
            SmallMatrix<complex_type,N,N,ColMajor|CStyle>& m0) const
        { (m0=m) += x2; }
        void assignTom(
            SmallMatrix<real_type,N,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<real_type,N,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<complex_type,N,N,RowMajor|FortranStyle>& m0) const
        { (m0=m) += x2; }
        void assignTom(
            SmallMatrix<complex_type,N,N,ColMajor|FortranStyle>& m0) const
        { (m0=m) += x2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
                real_type* m0p = m0.ptr();
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += TMV_REAL(x2);
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
                real_type* m0p = m0.ptr();
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += TMV_REAL(x2);
            } else {
                (m0=m) += TMV_REAL(x2);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
                complex_type* m0p = m0.ptr();
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += x2;
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
                complex_type* m0p = m0.ptr();
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += x2;
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                (m0=m) += x2;
            }
        }
    private:
        const SmallMatrix<T1,N,N,A>& m;
        const T x2;
    };

    template <typename T, typename T1, ptrdiff_t N, int A>
    class SummX : public SmallMatrixComposite<T,N,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        SummX(T _x1, const SmallMatrix<T1,N,N,A>& _m, T _x2) :
            x1(_x1), m(_m), x2(_x2) {}
        T getX1() const { return x1; }
        const SmallMatrix<T1,N,N,A>& getM() const
        { return m; }
        T getX2() const { return x2; }
        void assignTom(
            SmallMatrix<real_type,N,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=x1*m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<real_type,N,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=x1*m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<complex_type,N,N,RowMajor|CStyle>& m0) const
        { (m0=x1*m) += x2; }
        void assignTom(
            SmallMatrix<complex_type,N,N,ColMajor|CStyle>& m0) const
        { (m0=x1*m) += x2; }
        void assignTom(
            SmallMatrix<real_type,N,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=x1*m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<real_type,N,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); (m0=x1*m) += TMV_REAL(x2); }
        void assignTom(
            SmallMatrix<complex_type,N,N,RowMajor|FortranStyle>& m0) const
        { (m0=x1*m) += x2; }
        void assignTom(
            SmallMatrix<complex_type,N,N,ColMajor|FortranStyle>& m0) const
        { (m0=x1*m) += x2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
                real_type* m0p = m0.ptr();
                MultXV<N*N>(x1,m0p);
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += TMV_REAL(x2);
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
                real_type* m0p = m0.ptr();
                MultXV<N*N>(x1,m0p);
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += TMV_REAL(x2);
            } else {
                (m0=x1*m) += TMV_REAL(x2);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
                complex_type* m0p = m0.ptr();
                MultXV<N*N>(x1,m0p);
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += x2;
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m)) {
                SmallMatrixCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
                complex_type* m0p = m0.ptr();
                MultXV<N*N>(x1,m0p);
                for(ptrdiff_t i=0;i<N;i++) m0p[(N+1)*i] += x2;
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                (m0=x1*m) += x2;
            }
        }
    private:
        const T x1;
        const SmallMatrix<T1,N,N,A>& m;
        const T x2;
    };

    // m+=x
    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<T,N,N,A>& operator+=(SmallMatrix<T,N,N,A>& m, T x)
    { for(ptrdiff_t i=0;i<N;i++) m.ref(i,i) += x; return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator+=(SmallMatrix<CT,N,N,A>& m, T x)
    { for(ptrdiff_t i=0;i<N;i++) m.ref(i,i) += x; return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator+=(SmallMatrix<CT,N,N,A>& m, CCT x)
    { m += CT(x); return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator+=(SmallMatrix<CT,N,N,A>& m, VCT x)
    { m += CT(x); return m; }

    // m-=x
    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<T,N,N,A>& operator-=(SmallMatrix<T,N,N,A>& m, T x)
    { m += (-x); return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator-=(SmallMatrix<CT,N,N,A>& m, T x)
    { m += (-x); return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator-=(SmallMatrix<CT,N,N,A>& m, CCT x)
    { m += CT(-x); return m; }

    template <typename T, ptrdiff_t N, int A>
    inline SmallMatrix<CT,N,N,A>& operator-=(SmallMatrix<CT,N,N,A>& m, VCT x)
    { m += CT(-x); return m; }

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define SUMMX SummX
#define SUMMX_1 SummX_1
#define X1 ,N,N,A
#define X2 ,N,A
#define Y , ptrdiff_t N, int A
#include "tmv/TMV_AuxSumMX.h"
    // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

    //
    // Vector ^ Vector (OuterProduct)
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class OProdvv_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        OProdvv_1(
            const T TMV_DEBUGPARAM(_x), const SmallVector<T1,M,A1>& _v1,
            const SmallVector<T2,N,A2>& _v2) :
            v1(_v1), v2(_v2)
        { TMVAssert(_x == T(1)); }
        const SmallVector<T1,M,A1>& getV1() const { return v1; }
        const SmallVector<T2,N,A2>& getV2() const { return v2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M) {
                Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N) {
                Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
            } else {
                Rank1Update<false>(T(1),v1.view(),v2.view(),m0.view());
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            if (m0.stepi() == 1 && m0.stepj() == M && !m0.isconj()) {
                Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N && !m0.isconj()) {
                Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
            } else {
                Rank1Update<false>(T(1),v1.view(),v2.view(),m0.view());
            }
        }
    private:
        const SmallVector<T1,M,A1>& v1;
        const SmallVector<T2,N,A2>& v2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class OProdvv : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        OProdvv(
            const T _x, const SmallVector<T1,M,A1>& _v1,
            const SmallVector<T2,N,A2>& _v2) :
            x(_x), v1(_v1), v2(_v2) {}
        T getX() const { return x; }
        const SmallVector<T1,M,A1>& getV1() const { return v1; }
        const SmallVector<T2,N,A2>& getV2() const { return v2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M) {
                Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N) {
                Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
            } else {
                Rank1Update<false>(x,v1.view(),v2.view(),m0.view());
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            if (m0.stepi() == 1 && m0.stepj() == M && !m0.isconj()) {
                Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N && !m0.isconj()) {
                Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
            } else {
                Rank1Update<false>(x,v1.view(),v2.view(),m0.view());
            }
        }
    private:
        T x;
        const SmallVector<T1,M,A1>& v1;
        const SmallVector<T2,N,A2>& v2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, int A>
    class OProdVv : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        OProdVv(
            const T _x, const GenVector<T1>& _v1,
            const SmallVector<T2,N,A>& _v2) :
            x(_x), v1(_v1), v2(_v2) {}
        ptrdiff_t colsize() const { return v1.size(); }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const GenVector<T1>& getV1() const { return v1; }
        const SmallVector<T2,N,A>& getV2() const { return v2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<false>(x,v1,v2.view(),m0.view());
        }
        void assignToM(MatrixView<complex_type> m0) const
        { Rank1Update<false>(x,v1,v2.view(),m0.view()); }
    private:
        T x;
        const GenVector<T1>& v1;
        const SmallVector<T2,N,A>& v2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, int A>
    class OProdvV : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        OProdvV(
            const T _x, const SmallVector<T1,M,A>& _v1,
            const GenVector<T2>& _v2) :
            x(_x), v1(_v1), v2(_v2) {}
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return v2.size(); }
        T getX() const { return x; }
        const SmallVector<T1,M,A>& getV1() const { return v1; }
        const GenVector<T2>& getV2() const { return v2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            Rank1Update<false>(x,v1.view(),v2,m0.view());
        }
        void assignToM(MatrixView<complex_type> m0) const
        { Rank1Update<false>(x,v1.view(),v2,m0.view()); }
    private:
        T x;
        const SmallVector<T1,M,A>& v1;
        const GenVector<T2>& v2;
    };

    // m+=(v^v)
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<T,M,N,A>& operator+=(
        SmallMatrix<T,M,N,A>& m0, const OProdvv_1<T,T1,T2,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update_1<M,N,S>(
            opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<CT,M,N,A>& operator+=(
        SmallMatrix<CT,M,N,A>& m0, const OProdvv_1<T,T,T,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update_1<M,N,S>(
            opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m0, const OProdvv_1<T,T1,T2,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(T(1),opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m0, const OProdvv_1<T,T,T,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(CT(1),opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    // m-=(v^v)
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<T,M,N,A>& operator-=(
        SmallMatrix<T,M,N,A>& m0, const OProdvv_1<T,T1,T2,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update_m1<M,N,S>(
            opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<CT,M,N,A>& operator-=(
        SmallMatrix<CT,M,N,A>& m0, const OProdvv_1<T,T,T,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update_m1<M,N,S>(
            opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m0, const OProdvv_1<T,T1,T2,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(T(-1),opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m0, const OProdvv_1<T,T,T,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(CT(-1), opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    // m+=(x*v^v)
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<T,M,N,A>& operator+=(
        SmallMatrix<T,M,N,A>& m0, const OProdvv<T,T1,T2,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update<M,N,S>(
            opvv.getX(), opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<CT,M,N,A>& operator+=(
        SmallMatrix<CT,M,N,A>& m0, const OProdvv<T,T,T,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update<M,N,S>(
            opvv.getX(), opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator+=(
        SmallMatrix<T,M,N,A>& m0, const OProdVV<T,T1,T2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(opvv.getX(), opvv.getV1(), opvv.getV2(), m0.view());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator+=(
        SmallMatrix<CT,M,N,A>& m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(opvv.getX(), opvv.getV1(), opvv.getV2(), m0.view());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m0, const OProdvv<T,T1,T2,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(
            opvv.getX(), opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m0, const OProdvv<T,T,T,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(
            opvv.getX(), opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    // m-=(x*v^v)
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<T,M,N,A>& operator-=(
        SmallMatrix<T,M,N,A>& m0, const OProdvv<T,T1,T2,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update<M,N,S>(
            -opvv.getX(), opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline SmallMatrix<CT,M,N,A>& operator-=(
        SmallMatrix<CT,M,N,A>& m0, const OProdvv<T,T,T,M,N,A1,A2>& opvv)
    {
        const StorageType S = static_cast<StorageType>(A&AllStorageType);
        AddRank1Update<M,N,S>(
            -opvv.getX(), opvv.getV1().cptr(), opvv.getV2().cptr(), m0.ptr());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator-=(
        SmallMatrix<T,M,N,A>& m0, const OProdVV<T,T1,T2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), opvv.getV2(), m0.view());
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<CT,M,N,A>& operator-=(
        SmallMatrix<CT,M,N,A>& m0, const OProdVV<T,T,T>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(-opvv.getX(), opvv.getV1(), opvv.getV2(), m0.view());
        return m0;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m0, const OProdvv<T,T1,T2,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(
            -opvv.getX(), opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m0, const OProdvv<T,T,T,M,N,A1,A2>& opvv)
    {
        TMVAssert(m0.colsize() == opvv.colsize());
        TMVAssert(m0.rowsize() == opvv.rowsize());
        Rank1Update<true>(
            -opvv.getX(), opvv.getV1().view(), opvv.getV2().view(), m0);
        return m0;
    }

#define PRODMM OProdvv
#define PRODMM_1 OProdvv_1
#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallVector
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,M,A1
#define X2 ,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define OP operator^
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 OProdvv_1
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM OProdVv
#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallVector
#define PRODXM1 ProdXV
#define PRODXM2 ProdXv
#define X2 ,N,A
#define X3 ,N,A
#define Y , ptrdiff_t N, int A
#define OP operator^
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM OProdvV
#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXv
#define PRODXM2 ProdXV
#define X1 ,M,A
#define X3 ,M,A
#define Y , ptrdiff_t M, int A
#define OP operator^
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


    //
    // Matrix + Matrix
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm_1_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm_1_1(
            const T TMV_DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,A1>& _m1,
            const T TMV_DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
        const SmallMatrix<T1,M,N,A1>& getM1() const
        { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const
        { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
            } else {
                AddMM(T(1),m1.view(),T(1),m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(T(1),m1.view(),T(1),m2.view(),m0);
            }
        }
    private:
        const SmallMatrix<T1,M,N,A1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm_1_m1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm_1_m1(
            const T TMV_DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,A1>& _m1,
            const T TMV_DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
            } else {
                AddMM(T(1),m1.view(),T(-1),m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(T(1),m1.view(),T(-1),m2.view(),m0);
            }
        }
    private:
        const SmallMatrix<T1,M,N,A1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm_1_x : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm_1_x(
            const T TMV_DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,A1>& _m1,
            const T _x2, const SmallMatrix<T2,M,N,A2>& _m2) :
            m1(_m1), x2(_x2), m2(_m2)
        { TMVAssert(_x1 == T(1)); }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        T getX2() const { return x2; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x2*m2;
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x2*m2;
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            m0 = x2*m2;
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            m0 = x2*m2;
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x2*m2;
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x2*m2;
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            m0 = x2*m2;
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            m0 = x2*m2;
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                MultXV<M*N>(x2,m0.ptr());
                AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                MultXV<M*N>(x2,m0.ptr());
                AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
            } else {
                AddMM(T(1),m1.view(),x2,m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                MultXV<M*N>(x2,m0.ptr());
                AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                MultXV<M*N>(x2,m0.ptr());
                AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(T(1),m1.view(),x2,m2.view(),m0);
            }
        }
    private:
        const SmallMatrix<T1,M,N,A1>& m1;
        const T x2;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm_x_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm_x_1(
            const T _x1, const SmallMatrix<T1,M,N,A1>& _m1,
            const T TMV_DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,A2>& _m2) :
            x1(_x1), m1(_m1), m2(_m2)
        { TMVAssert(_x2 == T(1)); }
        T getX1() const { return x1; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
            } else {
                AddMM(x1,m1.view(),T(1),m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(x1,m1.view(),T(1),m2.view(),m0);
            }
        }
    private:
        const T x1;
        const SmallMatrix<T1,M,N,A1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm_x_m1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm_x_m1(
            const T _x1, const SmallMatrix<T1,M,N,A1>& _m1,
            const T TMV_DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,A2>& _m2) :
            x1(_x1), m1(_m1), m2(_m2)
        { TMVAssert(_x2 == T(-1)); }
        T getX1() const { return x1; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
            } else {
                AddMM(x1,m1.view(),T(1),m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(x1,m1.view(),T(1),m2.view(),m0);
            }
        }
    private:
        const T x1;
        const SmallMatrix<T1,M,N,A1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Summm : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Summm(
            const T _x1, const SmallMatrix<T1,M,N,A1>& _m1,
            const T _x2, const SmallMatrix<T2,M,N,A2>& _m2) :
            x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
        T getX1() const { return x1; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        T getX2() const { return x2; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            m0 = x1*m1;
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            m0 = x1*m1;
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                MultXV<M*N>(x1,m0.ptr());
                AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
    private:
        const T x1;
        const SmallMatrix<T1,M,N,A1>& m1;
        const T x2;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1>
    class SummM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        SummM(
            const T _x1, const SmallMatrix<T1,M,N,A1>& _m1,
            const T _x2, const GenMatrix<T2>& _m2) :
            x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return N; }
        T getX1() const { return x1; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        T getX2() const { return x2; }
        const GenMatrix<T2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            if (m2.stepi() == 1 && m2.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,ColMajor,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,ColMajor,ColMajor>(x2,m2.cptr(),m0.ptr());
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,ColMajor,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,ColMajor,RowMajor>(x2,m2.cptr(),m0.ptr());
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else if (m2.stepj() == 1 && m2.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,RowMajor,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,RowMajor,ColMajor>(x2,m2.cptr(),m0.ptr());
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,RowMajor,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,RowMajor,RowMajor>(x2,m2.cptr(),m0.ptr());
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            if (m2.stepi() == 1 && m2.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,ColMajor,ColMajor>(m2.cptr(),m0.ptr());
                    if (m2.isconj()) m0.conjugateSelf();
                    if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
                    if (x1 == T(1))
                        AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S1,ColMajor>(x1,m1.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,ColMajor,RowMajor>(m2.cptr(),m0.ptr());
                    if (m2.isconj()) m0.conjugateSelf();
                    if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
                    if (x1 == T(1))
                        AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S1,RowMajor>(x1,m1.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else if (m2.stepj() == 1 && m2.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,RowMajor,ColMajor>(m2.cptr(),m0.ptr());
                    if (m2.isconj()) m0.conjugateSelf();
                    if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
                    if (x1 == T(1))
                        AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S1,ColMajor>(x1,m1.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,RowMajor,RowMajor>(m2.cptr(),m0.ptr());
                    if (m2.isconj()) m0.conjugateSelf();
                    if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
                    if (x1 == T(1))
                        AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S1,RowMajor>(x1,m1.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
    private:
        const T x1;
        const SmallMatrix<T1,M,N,A1>& m1;
        const T x2;
        const GenMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
    class SumMm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        SumMm(
            const T _x1, const GenMatrix<T1>& _m1,
            const T _x2, const SmallMatrix<T2,M,N,A2>& _m2) :
            x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return N; }
        T getX1() const { return x1; }
        const GenMatrix<T1>& getM1() const { return m1; }
        T getX2() const { return x2; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m1.stepi() == 1 && m1.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,ColMajor,ColMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,ColMajor,RowMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else if (m2.stepj() == 1 && m2.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,RowMajor,ColMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,RowMajor,RowMajor>(m1.cptr(),m0.ptr());
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m1.stepi() == 1 && m1.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,ColMajor,ColMajor>(m1.cptr(),m0.ptr());
                    if (m1.isconj()) m0.conjugateSelf();
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,ColMajor,RowMajor>(m1.cptr(),m0.ptr());
                    if (m1.isconj()) m0.conjugateSelf();
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else if (m2.stepj() == 1 && m2.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                if (m0.stepi() == 1 && m0.stepj() == M) {
                    SmallMatrixCopy<M,N,RowMajor,ColMajor>(m1.cptr(),m0.ptr());
                    if (m1.isconj()) m0.conjugateSelf();
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else if (m0.stepj() == 1 && m0.stepi() == N) {
                    SmallMatrixCopy<M,N,RowMajor,RowMajor>(m1.cptr(),m0.ptr());
                    if (m1.isconj()) m0.conjugateSelf();
                    if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
                    if (x2 == T(1))
                        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
                    else
                        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
                    if (m0.isconj()) m0.conjugateSelf();
                } else {
                    AddMM(x1,m1.view(),x2,m2.view(),m0);
                }
            } else {
                AddMM(x1,m1.view(),x2,m2.view(),m0);
            }
        }
    private:
        const T x1;
        const GenMatrix<T1>& m1;
        const T x2;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    // m += m
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM_1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM_1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        TMVAssert(m1.colsize()==M);
        TMVAssert(m1.rowsize()==N);
        AddMM(T(1),m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        TMVAssert(m1.colsize()==M);
        TMVAssert(m1.rowsize()==N);
        AddMM(CT(1),m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==M);
        TMVAssert(m2.rowsize()==N);
        AddMM(T(1),m2,m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==M);
        TMVAssert(m2.rowsize()==N);
        AddMM(CT(1),m2,m1.view());
        return m1;
    }

    // m -= m
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM_m1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM_m1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        TMVAssert(m1.colsize()==M);
        TMVAssert(m1.rowsize()==N);
        AddMM(T(-1),m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m1, const SmallMatrix<T,M,N,A2>& m2)
    {
        TMVAssert(m1.colsize()==M);
        TMVAssert(m1.rowsize()==N);
        AddMM(CT(-1),m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==M);
        TMVAssert(m2.rowsize()==N);
        AddMM(T(-1),m2,m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.colsize()==M);
        TMVAssert(m2.rowsize()==N);
        AddMM(CT(-1),m2,m1.view());
        return m1;
    }

    // m += x*m
    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const ProdXm<T,T2,M,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM<M,N,S2,S1>(pxm.getX(),pxm.getM().cptr(),m.ptr()); return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const ProdXm<T,T,M,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM<M,N,S2,S1>(pxm.getX(),pxm.getM().cptr(),m.ptr()); return m;
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ProdXm<T,T2,M,N,A2>& pxm)
    {
        TMVAssert(m.colsize()==M);
        TMVAssert(m.rowsize()==N);
        AddMM(pxm.getX(),pxm.getM().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ProdXm<T,T,M,N,A2>& pxm)
    {
        TMVAssert(m.colsize()==M);
        TMVAssert(m.rowsize()==N);
        AddMM(pxm.getX(),pxm.getM().view(),m);
        return m;
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize()==M);
        TMVAssert(pxm.rowsize()==N);
        AddMM(pxm.getX(),pxm.getM(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const ProdXM<T,T>& pxm)
    {
        TMVAssert(pxm.colsize()==M);
        TMVAssert(pxm.rowsize()==N);
        AddMM(pxm.getX(),pxm.getM(),m.view());
        return m;
    }

    // m -= xm
    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const ProdXm<T,T2,M,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM<M,N,S2,S1>(-pxm.getX(),pxm.getM().cptr(),m.ptr()); return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const ProdXm<T,T,M,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        AddMM<M,N,S2,S1>(-pxm.getX(),pxm.getM().cptr(),m.ptr()); return m;
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ProdXm<T,T2,M,N,A2>& pxm)
    {
        TMVAssert(m.colsize()==M);
        TMVAssert(m.rowsize()==N);
        AddMM(-pxm.getX(),pxm.getM().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ProdXm<T,T,M,N,A2>& pxm)
    {
        TMVAssert(m.colsize()==M);
        TMVAssert(m.rowsize()==N);
        AddMM(-pxm.getX(),pxm.getM().view(),m);
        return m;
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize()==M);
        TMVAssert(pxm.rowsize()==N);
        AddMM(-pxm.getX(),pxm.getM(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const ProdXM<T,T>& pxm)
    {
        TMVAssert(pxm.colsize()==M);
        TMVAssert(pxm.rowsize()==N);
        AddMM(-pxm.getX(),pxm.getM(),m.view());
        return m;
    }


#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define SUMMM Summm
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X1 ,M,N,A1
#define X2 ,M,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#include "tmv/TMV_AuxSumMM.h"
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#include "tmv/TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define SUMMM SummM
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define X1 ,M,N,A1
#define X3 ,M,N,A1
#define Y , ptrdiff_t M, ptrdiff_t N, int A1
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define SUMMM SumMm
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define X2 ,M,N,A2
#define X3 ,M,N,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A2
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2


    //
    // Matrix * Matrix
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class Prodmm_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodmm_1(
            const T TMV_DEBUGPARAM(_x), const SmallMatrix<T1,M,K,A1>& _m1,
            const SmallMatrix<T2,K,N,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x == T(1)); }
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
            } else {
                MultMM<false>(T(1), m1.view(), m2.view(), m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                MultMM<false>(T(1), m1.view(), m2.view(), m0);
            }
        }
    private:
        const SmallMatrix<T1,M,K,A1>& m1;
        const SmallMatrix<T2,K,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class Prodmm : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodmm(
            const T _x, const SmallMatrix<T1,M,K,A1>& _m1,
            const SmallMatrix<T2,K,N,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2) {}
        T getX() const { return x; }
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,N,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        }
        void assignToM(MatrixView<real_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
            } else {
                MultMM<false>(x, m1.view(), m2.view(), m0);
            }
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
            const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
            if (m0.stepi() == 1 && m0.stepj() == M &&
                !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else if (m0.stepj() == 1 && m0.stepi() == N &&
                       !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
                MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
                if (m0.isconj()) m0.conjugateSelf();
            } else {
                MultMM<false>(x, m1.view(), m2.view(), m0);
            }
        }
    private:
        const T x;
        const SmallMatrix<T1,M,K,A1>& m1;
        const SmallMatrix<T2,K,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t K, ptrdiff_t N, int A2>
    class ProdMm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdMm(
            const T _x, const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,K,N,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == K); }
        ptrdiff_t colsize() const { return m1.colsize(); }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,N,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); MultMM<false>(x, m1, m2.view(), m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { MultMM<false>(x, m1, m2.view(), m0); }
    private:
        const T x;
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,K,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t K, int A1>
    class ProdmM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdmM(
            const T _x, const SmallMatrix<T1,M,K,A1>& _m1,
            const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m2.colsize() == K); }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return m2.rowsize(); }
        T getX() const { return x; }
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const GenMatrix<T2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); MultMM<false>(x, m1.view(), m2, m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { MultMM<false>(x, m1.view(), m2, m0); }
    private:
        const T x;
        const SmallMatrix<T1,M,K,A1>& m1;
        const GenMatrix<T2>& m2;
    };

    // m *= m
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        SmallMatrix<T,M,N,A1> m1_copy(m1);
        MultMM_1<M,N,N,S1,S2,S1>(m1_copy.cptr(),m2.cptr(),m1.ptr());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        SmallMatrix<CT,M,N,A1> m1_copy(m1);
        MultMM_1<M,N,N,S1,S2,S1>(m1_copy.cptr(),m2.cptr(),m1.ptr());
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        TMVAssert(m1.rowsize() == N);
        MultMM<false>(T(1),m1,m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        TMVAssert(m1.rowsize() == N);
        MultMM<false>(T(1),m1,m2.view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == N);
        MultMM<false>(T(1),m1.view(),m2,m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == N);
        MultMM<false>(T(1),m1.view(),m2,m1.view());
        return m1;
    }

    // m *= xm
    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const ProdXm<T,T2,N,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        SmallMatrix<T,M,N,A1> m1_copy(m1);
        MultMM<M,N,N,S1,S2,S1>(
            pxm.getX(),m1_copy.cptr(),pxm.getM().cptr(),m1.ptr());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const ProdXm<T,T,N,N,A2>& pxm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        SmallMatrix<CT,M,N,A1> m1_copy(m1);
        MultMM<M,N,N,S1,S2,S1>(
            pxm.getX(),m1_copy.cptr(),pxm.getM().cptr(),m1.ptr());
        return m1;
    }

    template <typename T, typename T2, ptrdiff_t N, int A2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const ProdXm<T,T2,N,N,A2>& pxm)
    {
        TMVAssert(m1.rowsize() == N);
        MultMM<false>(pxm.getX(),m1,pxm.getM().view(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const ProdXm<T,T,N,N,A2>& pxm)
    {
        TMVAssert(m1.rowsize() == N);
        MultMM<false>(pxm.getX(),m1,pxm.getM().view(),m1);
        return m1;
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const ProdXM<T,T2>& pxm)
    {
        TMVAssert(pxm.colsize() == N);
        TMVAssert(pxm.rowsize() == N);
        MultMM<false>(pxm.getX(),m1.view(),pxm.getM(),m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const ProdXM<T,T>& pxm)
    {
        TMVAssert(pxm.colsize() == N);
        TMVAssert(pxm.rowsize() == N);
        MultMM<false>(pxm.getX(),m1.view(),pxm.getM(),m1.view());
        return m1;
    }

    // m += mm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const Prodmm_1<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM_1<M,N,K,S2,S3,S1>(
            pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const Prodmm_1<T,T,T,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM_1<M,N,K,S2,S3,S1>(
            pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const Prodmm_1<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(T(1),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const Prodmm_1<T,T,T,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(T(1),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    // m -= mm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const Prodmm_1<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM_m1<M,N,K,S2,S3,S1>(
            pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const Prodmm_1<T,T,T,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM_m1<M,N,K,S2,S3,S1>(
            pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const Prodmm_1<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(T(-1),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const Prodmm_1<T,T,T,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(T(-1),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    // m += xmm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const Prodmm<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM<M,N,K,S2,S3,S1>(
            pmm.getX(),pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const Prodmm<T,T,T,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM<M,N,K,S2,S3,S1>(
            pmm.getX(),pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const Prodmm<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const Prodmm<T,T,T,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const ProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const ProdMM<T,T,T>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        MultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2(),m.view());
        return m;
    }

    // m -= xmm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const Prodmm<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM<M,N,K,S2,S3,S1>(
            -pmm.getX(),pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const Prodmm<T,T,T,M,N,K,A2,A3>& pmm)
    {
        const StorageType S1 = static_cast<StorageType>(A1&AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2&AllStorageType);
        const StorageType S3 = static_cast<StorageType>(A3&AllStorageType);
        AddMultMM<M,N,K,S2,S3,S1>(
            -pmm.getX(),pmm.getM1().cptr(),pmm.getM2().cptr(),m.ptr());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const Prodmm<T,T2,T3,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A2, int A3>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const Prodmm<T,T,T,M,N,K,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        MultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const ProdMM<T,T2,T3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const ProdMM<T,T,T>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        MultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2(),m.view());
        return m;
    }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM Prodmm
#define PRODMM_1 Prodmm_1
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define X1 ,M,K,A1
#define X2 ,K,N,A2
#define X3 ,M,N,K,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodmm_1
#define X3 ,M,N,K,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define PRODMM ProdmM
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define X1 ,M,K,A1
#define X3 ,M,K,A1
#define Y , ptrdiff_t M, ptrdiff_t K, int A1
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM ProdMm
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define X2 ,K,N,A2
#define X3 ,K,N,A2
#define Y , ptrdiff_t N, ptrdiff_t K, int A2
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2


    //
    // Element Product Matrix * Matrix
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class ElemProdmm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ElemProdmm(
            const T _x, const SmallMatrix<T1,M,N,A1>& _m1,
            const SmallMatrix<T2,M,N,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2)  {}
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }

        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(isReal(T()));
            ElemMultMM<false>(x, m1.view(), m2.view(), m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        { ElemMultMM<false>(x, m1.view(), m2.view(), m0); }
    private:
        const T x;
        const SmallMatrix<T1,M,N,A1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
    class ElemProdMm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ElemProdMm(
            const T _x, const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,M,N,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        {
            TMVAssert(m1.colsize() == M);
            TMVAssert(m1.rowsize() == N);
        }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,M,N,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); ElemMultMM<false>(x, m1, m2.view(), m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { ElemMultMM<false>(x, m1, m2.view(), m0); }
    private:
        const T x;
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,M,N,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1>
    class ElemProdmM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ElemProdmM(
            const T _x, const SmallMatrix<T1,M,N,A1>& _m1,
            const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        {
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
        }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const SmallMatrix<T1,M,N,A1>& getM1() const { return m1; }
        const GenMatrix<T2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); ElemMultMM<false>(x, m1.view(), m2, m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { ElemMultMM<false>(x, m1.view(), m2, m0); }
    private:
        const T x;
        const SmallMatrix<T1,M,N,A1>& m1;
        const GenMatrix<T2>& m2;
    };

    // m += xmm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m,
        const ElemProdmm<T,T2,T3,M,N,A2,A3>& pmm)
    {
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m,
        const ElemProdmm<T,T,T,M,N,A2,A3>& pmm)
    {
        ElemMultMM<true>(
            pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A3>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const ElemProdMm<T,T2,T3,M,N,A3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const ElemProdMm<T,T,T,M,N,A3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator+=(
        SmallMatrix<T,M,N,A1>& m, const ElemProdmM<T,T2,T3,M,N,A2>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator+=(
        SmallMatrix<CT,M,N,A1>& m, const ElemProdmM<T,T,T,M,N,A2>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ElemProdmm<T,T2,T3,M,N,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ElemProdmm<T,T,T,M,N,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A3>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ElemProdMm<T,T2,T3,M,N,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A3>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ElemProdMm<T,T,T,M,N,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator+=(
        MatrixView<T> m, const ElemProdmM<T,T2,T3,M,N,A2>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator+=(
        MatrixView<CT> m, const ElemProdmM<T,T,T,M,N,A2>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(pmm.getX(),pmm.getM1().view(),pmm.getM2(),m);
        return m;
    }

    // m -= xmm
    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m,
        const ElemProdmm<T,T2,T3,M,N,A2,A3>& pmm)
    {
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m,
        const ElemProdmm<T,T,T,M,N,A2,A3>& pmm)
    {
        ElemMultMM<true>(
            -pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A3>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const ElemProdMm<T,T2,T3,M,N,A3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A3>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const ElemProdMm<T,T,T,M,N,A3>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2().view(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator-=(
        SmallMatrix<T,M,N,A1>& m, const ElemProdmM<T,T2,T3,M,N,A2>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator-=(
        SmallMatrix<CT,M,N,A1>& m, const ElemProdmM<T,T,T,M,N,A2>& pmm)
    {
        TMVAssert(pmm.colsize() == M);
        TMVAssert(pmm.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2(),m.view());
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ElemProdmm<T,T2,T3,M,N,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ElemProdmm<T,T,T,M,N,A2,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A3>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ElemProdMm<T,T2,T3,M,N,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A3>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ElemProdMm<T,T,T,M,N,A3>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1(),pmm.getM2().view(),m);
        return m;
    }

    template <typename T, typename T2, typename T3, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<T> operator-=(
        MatrixView<T> m, const ElemProdmM<T,T2,T3,M,N,A2>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2(),m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator-=(
        MatrixView<CT> m, const ElemProdmM<T,T,T,M,N,A2>& pmm)
    {
        TMVAssert(m.colsize() == M);
        TMVAssert(m.rowsize() == N);
        ElemMultMM<true>(-pmm.getX(),pmm.getM1().view(),pmm.getM2(),m);
        return m;
    }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM ElemProdmm
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define OP ElemProd
#define X1 ,M,N,A1
#define X2 ,M,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define PRODMM ElemProdmM
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define OP ElemProd
#define X1 ,M,N,A1
#define X3 ,M,N,A1
#define Y , ptrdiff_t M, ptrdiff_t N, int A1
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM ElemProdMm
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define OP ElemProd
#define X2 ,M,N,A2
#define X3 ,M,N,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A2
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

    //
    // Matrix * Vector
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Prodmv_1 : public SmallVectorComposite<T,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodmv_1(
            const T TMV_DEBUGPARAM(_x), const SmallMatrix<T1,M,N,A1>& _m,
            const SmallVector<T2,N,A2>& _v) :
            m(_m), v(_v)
        { TMVAssert(_x == T(1)); }
        const SmallMatrix<T1,M,N,A1>& getM() const { return m; }
        const SmallVector<T2,N,A2>& getV() const { return v; }
        void assignTov(SmallVector<real_type,M,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,M,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<real_type,M,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,M,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
            else
                MultMV<false>(T(1),m.view(),v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            if (v0.step() == 1) {
                MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(T(1),m.view(),v.view(),v0);
        }
    private:
        const SmallMatrix<T1,M,N,A1>& m;
        const SmallVector<T2,N,A2>& v;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Prodmv : public SmallVectorComposite<T,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodmv(
            const T _x, const SmallMatrix<T1,M,N,A1>& _m,
            const SmallVector<T2,N,A2>& _v) :
            x(_x), m(_m), v(_v) {}
        T getX() const { return x; }
        const SmallMatrix<T1,M,N,A1>& getM() const { return m; }
        const SmallVector<T2,N,A2>& getV() const { return v; }
        void assignTov(SmallVector<real_type,M,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,M,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<real_type,M,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,M,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
            else
                MultMV<false>(x,m.view(),v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A1&AllStorageType);
            if (v0.step() == 1) {
                MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(x,m.view(),v.view(),v0);
        }
    private:
        const T x;
        const SmallMatrix<T1,M,N,A1>& m;
        const SmallVector<T2,N,A2>& v;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, int A>
    class ProdMv : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdMv(
            const T _x, const GenMatrix<T1>& _m,
            const SmallVector<T2,N,A>& _v) :
            x(_x), m(_m), v(_v)
        { TMVAssert(m.rowsize() == N); }
        ptrdiff_t size() const { return m.colsize(); }
        T getX() const { return x; }
        const GenMatrix<T1>& getM() const { return m; }
        const SmallVector<T2,N,A>& getV() const { return v; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(isReal(T()));
            MultMV<false>(x,m,v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        { MultMV<false>(x,m,v.view(),v0); }
    private:
        const T x;
        const GenMatrix<T1>& m;
        const SmallVector<T2,N,A>& v;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class ProdmV : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdmV(
            const T _x, const SmallMatrix<T1,M,N,A>& _m,
            const GenVector<T2>& _v) :
            x(_x), m(_m), v(_v)
        { TMVAssert(v.size() == N); }
        ptrdiff_t size() const { return M; }
        T getX() const { return x; }
        const SmallMatrix<T1,M,N,A>& getM() const { return m; }
        const GenVector<T2>& getV() const { return v; }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            TMVAssert(isReal(T()));
            if (v.step() == 1 && v0.step() == 1 && !SameStorage(v0,v))
                if (x == T(1))
                    MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
                else
                    MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
            else
                MultMV<false>(x,m.view(),v,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            if (v.step() == 1 && v0.step() == 1 && !SameStorage(v0,v) &&
                !v.isconj()) {
                if (x == T(1))
                    MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
                else
                    MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(x,m.view(),v,v0);
        }
    private:
        const T x;
        const SmallMatrix<T1,M,N,A>& m;
        const GenVector<T2>& v;
    };


    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Prodvm_1 : public SmallVectorComposite<T,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodvm_1(
            const T TMV_DEBUGPARAM(_x), const SmallVector<T1,M,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            v(_v), m(_m)
        { TMVAssert(_x == T(1)); }
        const SmallVector<T1,M,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
            else
                MultMV<false>(T(1),m.transpose(),v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            if (v0.step() == 1) {
                MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(T(1),m.transpose(),v.view(),v0);
        }
    private:
        const SmallVector<T1,M,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Prodvm : public SmallVectorComposite<T,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Prodvm(
            const T _x, const SmallVector<T1,M,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            x(_x), v(_v), m(_m) { }
        T getX() const { return x; }
        const SmallVector<T1,M,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        }
        void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
            else
                MultMV<false>(x,m.transpose(),v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A2&AllStorageType);
            if (v0.step() == 1) {
                MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(x,m.transpose(),v.view(),v0);
        }
    private:
        const T x;
        const SmallVector<T1,M,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class ProdVm : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdVm(
            const T _x, const GenVector<T1>& _v,
            const SmallMatrix<T2,M,N,A>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size() == M); }
        ptrdiff_t size() const { return N; }
        T getX() const { return x; }
        const GenVector<T1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            TMVAssert(isReal(T()));
            if (v0.step() == 1 && v.step() == 1 && !SameStorage(v0,v))
                if (x == T(1))
                    MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
                else
                    MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
            else
                MultMV<false>(x,m.transpose(),v,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            const StorageType S = static_cast<StorageType>(A&AllStorageType);
            if (v0.step() == 1 && v.step() == 1 && !SameStorage(v0,v) &&
                !v.isconj()) {
                if (x == T(1))
                    MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
                else
                    MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultMV<false>(x,m.transpose(),v,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const SmallMatrix<T2,M,N,A>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, int A>
    class ProdvM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        ProdvM(
            const T _x, const SmallVector<T1,M,A>& _v,
            const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(m.colsize() == M); }
        ptrdiff_t size() const { return m.rowsize(); }
        T getX() const { return x; }
        const SmallVector<T1,M,A>& getV() const { return v; }
        const GenMatrix<T2>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(isReal(T()));
            MultMV<false>(x,m.transpose(),v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        { MultMV<false>(x,m.transpose(),v.view(),v0); }
    private:
        const T x;
        const SmallVector<T1,M,A>& v;
        const GenMatrix<T2>& m;
    };

    // v *= m
    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        SmallVector<T,N> v_copy(v);
        MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        SmallVector<CT,N> v_copy(v);
        MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const SmallMatrix<T,N,N,A2>& m)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj()) {
            SmallVector<T,N> v_copy(v);
            MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
        } else {
            MultMV<false>(T(1),m.transpose(),v,v);
        }
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const SmallMatrix<T,N,N,A2>& m)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1) { // conj is ok
            SmallVector<CT,N> v_copy(v);
            MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
        } else {
            MultMV<false>(T(1),m.transpose(),v,v);
        }
        return v;
    }

    template <typename T, ptrdiff_t N, int A>
    inline SmallVector<T,N,A>& operator*=(
        SmallVector<T,N,A>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.colsize() == N);
        TMVAssert(m.rowsize() == N);
        MultMV<false>(T(1),m.transpose(),v.view(),v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A>
    inline SmallVector<CT,N,A>& operator*=(
        SmallVector<CT,N,A>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.colsize() == N);
        TMVAssert(m.rowsize() == N);
        MultMV<false>(T(1),m.transpose(),v.view(),v.view());
        return v;
    }

    // v *= xm
    template <typename T, typename T1, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const ProdXm<T,T1,N,N,A2>& pxm)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        SmallVector<T,N> v_copy(v);
        MultVM<N,N,S>(pxm.getX(),v_copy.cptr(),pxm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const ProdXm<T,T,N,N,A2>& pxm)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        SmallVector<CT,N> v_copy(v);
        MultVM<N,N,S>(pxm.getX(),v_copy.cptr(),pxm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, ptrdiff_t N, int A2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const ProdXm<T,T1,N,N,A2>& pxm)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj()) {
            SmallVector<T,N> v_copy(v);
            MultVM<N,N,S>(pxm.getX(),v_copy.cptr(),pxm.getM().cptr(),v.ptr());
        } else {
            MultMV<false>(pxm.getX(),pxm.getM().transpose(),v,v);
        }
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const ProdXm<T,T,N,N,A2>& pxm)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1) { // conj is ok
            SmallVector<CT,N> v_copy(v);
            MultVM<N,N,S>(pxm.getX(),v_copy.cptr(),pxm.getM().cptr(),v.ptr());
        } else {
            MultMV<false>(pxm.getX(),pxm.getM().transpose(),v,v);
        }
        return v;
    }

    template <typename T, typename T1, ptrdiff_t N, int A1>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const ProdXM<T,T1>& pxm)
    {
        TMVAssert(pxm.colsize() == N);
        TMVAssert(pxm.rowsize() == N);
        MultMV<false>(pxm.getX(),pxm.getM().transpose(),v.view(),v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const ProdXM<T,T>& pxm)
    {
        TMVAssert(pxm.colsize() == N);
        TMVAssert(pxm.rowsize() == N);
        MultMV<false>(pxm.getX(),pxm.getM().transpose(),v.view(),v.view());
        return v;
    }

    // v += mv
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,M,A1>& operator+=(
        SmallVector<T,M,A1>& v, const Prodmv_1<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV_1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,M,A1>& operator+=(
        SmallVector<CT,M,A1>& v, const Prodmv_1<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV_1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator+=(
        VectorView<T> v, const Prodmv_1<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1 && !v.isconj())
            AddMultMV_1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(T(1),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const Prodmv_1<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1)  // conj is ok
            AddMultMV_1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(T(1),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    // v -= mv
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,M,A1>& operator-=(
        SmallVector<T,M,A1>& v, const Prodmv_1<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV_m1<M,N,S>(pmv.getM().cptr(), pmv.getV().cptr(), v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,M,A1>& operator-=(
        SmallVector<CT,M,A1>& v, const Prodmv_1<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV_m1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator-=(
        VectorView<T> v, const Prodmv_1<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1 && !v.isconj())
            AddMultMV_m1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(T(-1),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const Prodmv_1<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1)  // conj is ok
            AddMultMV_m1<M,N,S>(pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(T(-1),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    // v += vm
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,N,A1>& operator+=(
        SmallVector<T,N,A1>& v, const Prodvm_1<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM_1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,N,A1>& operator+=(
        SmallVector<CT,N,A1>& v, const Prodvm_1<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM_1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator+=(
        VectorView<T> v, const Prodvm_1<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj())
            AddMultVM_1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(T(1),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const Prodvm_1<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1)  // conj is ok
            AddMultVM_1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(T(1),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    // v -= vm
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,N,A1>& operator-=(
        SmallVector<T,N,A1>& v, const Prodvm_1<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM_m1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,N,A1>& operator-=(
        SmallVector<CT,N,A1>& v, const Prodvm_1<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM_m1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator-=(
        VectorView<T> v, const Prodvm_1<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj())
            AddMultVM_m1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(T(-1),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const Prodvm_1<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1)  // conj is ok
            AddMultVM_m1<M,N,S>(pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(T(-1),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    // v += xmv
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,M,A1>& operator+=(
        SmallVector<T,M,A1>& v, const Prodmv<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV<M,N,S>(pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,M,A1>& operator+=(
        SmallVector<CT,M,A1>& v, const Prodmv<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV<M,N,S>(pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator+=(
        VectorView<T> v, const Prodmv<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1 && !v.isconj())
            AddMultMV<M,N,S>(
                pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(pmv.getX(),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const Prodmv<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1)  // conj is ok
            AddMultMV<M,N,S>(
                pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(pmv.getX(),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, int A1>
    inline SmallVector<T,M,A1>& operator+=(
        SmallVector<T,M,A1>& v, const ProdMV<T,T1,T2>& pmv)
    {
        TMVAssert(pmv.size() == M);
        MultMV<true>(pmv.getX(),pmv.getM().view(),pmv.getV().view(),v.view());
        return v;
    }

    template <typename T, ptrdiff_t M, int A1>
    inline SmallVector<CT,M,A1>& operator+=(
        SmallVector<CT,M,A1>& v, const ProdMV<T,T,T>& pmv)
    {
        TMVAssert(pmv.size() == M);
        MultMV<true>(pmv.getX(),pmv.getM().view(),pmv.getV().view(),v.view());
        return v;
    }

    // v -= xmv
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,M,A1>& operator-=(
        SmallVector<T,M,A1>& v, const Prodmv<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV<M,N,S>(-pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,M,A1>& operator-=(
        SmallVector<CT,M,A1>& v, const Prodmv<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        AddMultMV<M,N,S>(-pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator-=(
        VectorView<T> v, const Prodmv<T,T1,T2,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1 && !v.isconj())
            AddMultMV<M,N,S>(
                -pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(-pmv.getX(),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const Prodmv<T,T,T,M,N,A2,A3>& pmv)
    {
        const StorageType S = static_cast<StorageType>(A2&AllStorageType);
        TMVAssert(v.size() == M);
        if (v.step() == 1)  // conj is ok
            AddMultMV<M,N,S>(
                -pmv.getX(),pmv.getM().cptr(),pmv.getV().cptr(),v.ptr());
        else
            MultMV<true>(-pmv.getX(),pmv.getM().view(),pmv.getV().view(),v);
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, int A1>
    inline SmallVector<T,M,A1>& operator-=(
        SmallVector<T,M,A1>& v, const ProdMV<T,T1,T2>& pmv)
    {
        TMVAssert(pmv.size() == M);
        MultMV<true>(-pmv.getX(),pmv.getM().view(),pmv.getV().view(),v.view());
        return v;
    }

    template <typename T, ptrdiff_t M, int A1>
    inline SmallVector<CT,M,A1>& operator-=(
        SmallVector<CT,M,A1>& v, const ProdMV<T,T,T>& pmv)
    {
        TMVAssert(pmv.size() == M);
        MultMV<true>(-pmv.getX(),pmv.getM().view(),pmv.getV().view(),v.view());
        return v;
    }

    // v += xvm
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,N,A1>& operator+=(
        SmallVector<T,N,A1>& v, const Prodvm<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM<M,N,S>(pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,N,A1>& operator+=(
        SmallVector<CT,N,A1>& v, const Prodvm<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM<M,N,S>(pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator+=(
        VectorView<T> v, const Prodvm<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj())
            AddMultVM<M,N,S>(
                pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(pvm.getX(),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const Prodvm<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1)  // conj is ok
            AddMultVM<M,N,S>(
                pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(pvm.getX(),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    // v -= xvm
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<T,N,A1>& operator-=(
        SmallVector<T,N,A1>& v, const Prodvm<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM<M,N,S>(
            -pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2, int A3>
    inline SmallVector<CT,N,A1>& operator-=(
        SmallVector<CT,N,A1>& v, const Prodvm<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        AddMultVM<M,N,S>(
            -pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        return v;
    }

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<T> operator-=(
        VectorView<T> v, const Prodvm<T,T1,T2,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1 && !v.isconj())
            AddMultVM<M,N,S>(
                -pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(-pvm.getX(),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A2, int A3>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const Prodvm<T,T,T,M,N,A2,A3>& pvm)
    {
        const StorageType S = static_cast<StorageType>(A3&AllStorageType);
        TMVAssert(v.size() == N);
        if (v.step() == 1)  // conj is ok
            AddMultVM<M,N,S>(
                -pvm.getX(),pvm.getV().cptr(),pvm.getM().cptr(),v.ptr());
        else
            MultMV<true>(-pvm.getX(),pvm.getM().transpose(),pvm.getV().view(),v);
        return v;
    }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallVector
#define PRODMM Prodmv
#define PRODXM1 ProdXm
#define PRODXM2 ProdXv
#define PRODMM_1 Prodmv_1
#define X1 ,M,N,A1
#define X2 ,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodmv_1
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallMatrix
#define PRODMM Prodvm
#define PRODXM1 ProdXv
#define PRODXM2 ProdXm
#define PRODMM_1 Prodvm_1
#define X1 ,M,A1
#define X2 ,M,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodvm_1
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenVector
#define PRODMM ProdmV
#define PRODXM1 ProdXm
#define PRODXM2 ProdXV
#define X1 ,M,N,A1
#define X3 ,M,N,A1
#define Y , ptrdiff_t M, ptrdiff_t N, int A1
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,N,A1
#define Y , ptrdiff_t M, ptrdiff_t N, int A1
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallVector
#define PRODMM ProdMv
#define PRODXM1 ProdXM
#define PRODXM2 ProdXv
#define X2 ,N,A2
#define X3 ,N,A2
#define Y , ptrdiff_t N, int A2
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,N,A2
#define Y , ptrdiff_t N, int A2
#define GETM1 .getM()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenMatrix
#define PRODMM ProdvM
#define PRODXM1 ProdXv
#define PRODXM2 ProdXM
#define X1 ,M,A1
#define X3 ,M,A1
#define Y , ptrdiff_t M, int A1
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,A1
#define Y , ptrdiff_t M, int A1
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallMatrix
#define PRODMM ProdVm
#define PRODXM1 ProdXV
#define PRODXM2 ProdXm
#define X2 ,M,N,A2
#define X3 ,M,N,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,N,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2


    //
    // Scalar / Matrix
    //

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotXm_1 : public SmallMatrixComposite<T,N,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotXm_1(
            const T TMV_DEBUGPARAM(_x),
            const SmallMatrix<Tm,M,N,A>& _m) :
            m(_m)
        { TMVAssert(_x == T(1)); }
        const SmallMatrix<Tm,M,N,A>& getM() const { return m; }
        void assignTom(
            SmallMatrix<real_type,N,M,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<real_type,N,M,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<complex_type,N,M,RowMajor|CStyle>& m0) const
        { DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<complex_type,N,M,ColMajor|CStyle>& m0) const
        { DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<real_type,N,M,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<real_type,N,M,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<complex_type,N,M,RowMajor|FortranStyle>& m0) const
        { DoInverse(m,m0); }
        void assignTom(
            SmallMatrix<complex_type,N,M,ColMajor|FortranStyle>& m0) const
        { DoInverse(m,m0); }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); m.view().makeInverse(m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { m.view().makeInverse(m0); }
    private:
        const SmallMatrix<Tm,M,N,A>& m;
    };

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotXm : public SmallMatrixComposite<T,N,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotXm(const T _x, const SmallMatrix<Tm,M,N,A>& _m) :
            x(_x), m(_m) {}
        T getX() const { return x; }
        const SmallMatrix<Tm,M,N,A>& getM() const { return m; }
        void assignTom(
            SmallMatrix<real_type,N,M,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,N,M,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,N,M,RowMajor|CStyle>& m0) const
        { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,N,M,ColMajor|CStyle>& m0) const
        { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,N,M,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,N,M,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,N,M,RowMajor|FortranStyle>& m0) const
        { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,N,M,ColMajor|FortranStyle>& m0) const
        { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignToM(MatrixView<real_type> m0) const
        { TMVAssert(isReal(T())); m.view().makeInverse(m0); MultXM(x,m0); }
        void assignToM(MatrixView<complex_type> m0) const
        { m.view().makeInverse(m0); MultXM(x,m0); }
    private:
        const T x;
        const SmallMatrix<Tm,M,N,A>& m;
    };

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define QUOTXM QuotXm
#define X ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#include "tmv/TMV_AuxQuotXM.h"
#define X ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#define QUOTXM_1 QuotXm_1
#include "tmv/TMV_AuxQuotXMa.h"
#undef GENMATRIX
#undef PRODXM
#undef QUOTXM


    //
    // Vector / % Matrix
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Quotvm_1 : public SmallVectorComposite<T,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Quotvm_1(
            const T TMV_DEBUGPARAM(_x), const SmallVector<T1,M,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            v(_v), m(_m)
        { TMVAssert(_x == T(1)); }
        const SmallVector<T1,M,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); DoLDiv(m,v,v0); }
        void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { DoLDiv(m,v,v0); }
        void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); DoLDiv(m,v,v0); }
        void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { DoLDiv(m,v,v0); }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == N);
            TMVAssert(isReal(T()));
            m.view().LDiv(v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == N);
            m.view().LDiv(v.view(),v0);
        }
    private:
        const SmallVector<T1,M,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class Quotvm : public SmallVectorComposite<T,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Quotvm(
            const T _x, const SmallVector<T1,M,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            x(_x), v(_v), m(_m) {}
        T getX() const { return x; }
        const SmallVector<T1,M,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
        void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
        void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
        void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == N);
            TMVAssert(isReal(T()));
            m.view().LDiv(v.view(),v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == N);
            m.view().LDiv(v.view(),v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const SmallVector<T1,M,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotVm_1 : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotVm_1(
            const T TMV_DEBUGPARAM(_x), const GenVector<T1>& _v,
            const SmallMatrix<T2,M,N,A>& _m) :
            v(_v), m(_m)
        { TMVAssert(_x == T(1)); TMVAssert(v.size() == M); }
        ptrdiff_t size() const { return N; }
        const GenVector<T1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.view().LDiv(v,v0);
            MultXV(T(1),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.view().LDiv(v,v0);
            MultXV(T(1),v0);
        }
    private:
        const GenVector<T1>& v;
        const SmallMatrix<T2,M,N,A>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotVm : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotVm(
            const T _x, const GenVector<T1>& _v,
            const SmallMatrix<T2,M,N,A>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size() == M); }
        ptrdiff_t size() const { return N; }
        T getX() const { return x; }
        const GenVector<T1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.view().LDiv(v,v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.view().LDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const SmallMatrix<T2,M,N,A>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, int A>
    class QuotvM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotvM(
            const T _x, const SmallVector<T1,M,A>& _v,
            const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(m.colsize() == M); }
        ptrdiff_t size() const { return m.rowsize(); }
        T getX() const { return x; }
        const SmallVector<T1,M,A>& getV() const { return v; }
        const GenMatrix<T2>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.LDiv(v.view(),v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.LDiv(v.view(),v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const SmallVector<T1,M,A>& v;
        const GenMatrix<T2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class RQuotvm_1 : public SmallVectorComposite<T,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotvm_1(
            const T TMV_DEBUGPARAM(_x), const SmallVector<T1,N,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            v(_v), m(_m)
        { TMVAssert(_x == T(1)); }
        const SmallVector<T1,N,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,M,CStyle>& v0) const
        { TMVAssert(isReal(T())); DoRDiv(m,v,v0); }
        void assignTov(SmallVector<complex_type,M,CStyle>& v0) const
        { DoRDiv(m,v,v0); }
        void assignTov(SmallVector<real_type,M,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); DoRDiv(m,v,v0); }
        void assignTov(SmallVector<complex_type,M,FortranStyle>& v0) const
        { DoRDiv(m,v,v0); }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == M);
            TMVAssert(isReal(T()));
            m.view().RDiv(v.view(),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == M);
            m.view().RDiv(v.view(),v0);
        }
    private:
        const SmallVector<T1,N,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    class RQuotvm : public SmallVectorComposite<T,M>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotvm(
            const T _x, const SmallVector<T1,N,A1>& _v,
            const SmallMatrix<T2,M,N,A2>& _m) :
            x(_x), v(_v), m(_m) {}
        T getX() const { return x; }
        const SmallVector<T1,N,A1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A2>& getM() const { return m; }
        void assignTov(SmallVector<real_type,M,CStyle>& v0) const
        { TMVAssert(isReal(T())); DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
        void assignTov(SmallVector<complex_type,M,CStyle>& v0) const
        { DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
        void assignTov(SmallVector<real_type,M,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
        void assignTov(SmallVector<complex_type,M,FortranStyle>& v0) const
        { DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == M);
            TMVAssert(isReal(T()));
            m.view().RDiv(v.view(),v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == M);
            m.view().RDiv(v.view(),v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const SmallVector<T1,N,A1>& v;
        const SmallMatrix<T2,M,N,A2>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class RQuotVm_1 : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotVm_1(
            const T TMV_DEBUGPARAM(_x), const GenVector<T1>& _v,
            const SmallMatrix<T2,M,N,A>& _m) :
            v(_v), m(_m)
        { TMVAssert(_x==T(1)); TMVAssert(v.size() == N); }
        ptrdiff_t size() const { return M; }
        const GenVector<T1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.view().RDiv(v,v0);
            MultXV(T(1),v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.view().RDiv(v,v0);
            MultXV(T(1),v0);
        }
    private:
        const GenVector<T1>& v;
        const SmallMatrix<T2,M,N,A>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A>
    class RQuotVm : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotVm(
            const T _x, const GenVector<T1>& _v,
            const SmallMatrix<T2,M,N,A>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(v.size() == N); }
        ptrdiff_t size() const { return M; }
        T getX() const { return x; }
        const GenVector<T1>& getV() const { return v; }
        const SmallMatrix<T2,M,N,A>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.view().RDiv(v,v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.view().RDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const SmallMatrix<T2,M,N,A>& m;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, int A>
    class RQuotvM : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotvM(
            const T _x, const SmallVector<T1,N,A>& _v,
            const GenMatrix<T2>& _m) :
            x(_x), v(_v), m(_m)
        { TMVAssert(m.rowsize() == N); }
        ptrdiff_t size() const { return m.colsize(); }
        T getX() const { return x; }
        const SmallVector<T1,N,A>& getV() const { return v; }
        const GenMatrix<T2>& getM() const { return m; }
        void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            m.RDiv(v.view(),v0);
            MultXV(x,v0);
        }
        void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            m.RDiv(v.view(),v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const SmallVector<T1,N,A>& v;
        const GenMatrix<T2>& m;
    };

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator/=(
        SmallVector<T,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        DoLDivEq(m,v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator/=(
        SmallVector<CT,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        DoLDivEq(m,v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator%=(
        SmallVector<T,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        DoRDivEq(m,v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator%=(
        SmallVector<CT,N,A1>& v, const SmallMatrix<T,N,N,A2>& m)
    {
        DoRDivEq(m,v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<T,N,A1>& operator/=(
        SmallVector<T,N,A1>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.rowsize() == N);
        TMVAssert(m.colsize() == N);
        m.LDivEq(v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<CT,N,A1>& operator/=(
        SmallVector<CT,N,A1>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.rowsize() == N);
        TMVAssert(m.colsize() == N);
        m.LDivEq(v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<T,N,A1>& operator%=(
        SmallVector<T,N,A1>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.rowsize() == N);
        TMVAssert(m.colsize() == N);
        m.RDivEq(v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<CT,N,A1>& operator%=(
        SmallVector<CT,N,A1>& v, const GenMatrix<T>& m)
    {
        TMVAssert(m.rowsize() == N);
        TMVAssert(m.colsize() == N);
        m.RDivEq(v.view());
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<T> operator/=(
        VectorView<T> v, const SmallMatrix<T,N,N,A2>& m)
    {
        TMVAssert(v.size() == N);
        m.view().LDivEq(v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator/=(
        VectorView<CT> v, const SmallMatrix<T,N,N,A2>& m)
    {
        TMVAssert(v.size() == N);
        m.view().LDivEq(v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<T> operator%=(
        VectorView<T> v, const SmallMatrix<T,N,N,A2>& m)
    {
        TMVAssert(v.size() == N);
        m.view().RDivEq(v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator%=(
        VectorView<CT> v, const SmallMatrix<T,N,N,A2>& m)
    {
        TMVAssert(v.size() == N);
        m.view().RDivEq(v);
        return v;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const QuotXm_1<T,Tm,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const QuotXm_1<T,T,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),v);
        return v;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A1, int A2>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const QuotXm<T,Tm,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),v);
        MultXV<N>(qxm.getX(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1, int A2>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const QuotXm<T,T,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),v);
        MultXV<N>(qxm.getX(),v.ptr());
        return v;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const QuotXm_1<T,Tm,N,N,A2>& qxm)
    {
        TMVAssert(v.size() == N);
        qxm.getM().view().RDivEq(v);
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const QuotXm_1<T,T,N,N,A2>& qxm)
    {
        TMVAssert(v.size() == N);
        qxm.getM().view().RDivEq(v);
        return v;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const QuotXm<T,Tm,N,N,A2>& qxm)
    {
        TMVAssert(v.size() == N);
        qxm.getM().view().RDivEq(v);
        MultXV(qxm.getX(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const QuotXm<T,T,N,N,A2>& qxm)
    {
        TMVAssert(v.size() == N);
        qxm.getM().view().RDivEq(v);
        MultXV(qxm.getX(),v.ptr());
        return v;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A1>
    inline SmallVector<T,N,A1>& operator*=(
        SmallVector<T,N,A1>& v, const QuotXM<T,Tm>& qxm)
    {
        TMVAssert(qxm.getM().rowsize() == N);
        TMVAssert(qxm.getM().colsize() == N);
        qxm.getM().RDivEq(v.view());
        MultXV(qxm.getX(),v.ptr());
        return v;
    }

    template <typename T, ptrdiff_t N, int A1>
    inline SmallVector<CT,N,A1>& operator*=(
        SmallVector<CT,N,A1>& v, const QuotXM<T,T>& qxm)
    {
        TMVAssert(qxm.getM().rowsize() == N);
        TMVAssert(qxm.getM().colsize() == N);
        qxm.getM().RDivEq(v.view());
        MultXV(qxm.getX(),v.ptr());
        return v;
    }

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXv
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 Quotvm_1
#define QUOTMM Quotvm
#define RQUOTMM_1 RQuotvm_1
#define RQUOTMM RQuotvm
#define X1 ,M,A1
#define X1b ,N,A1
#define X2 ,M,N,A2
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM_1 Quotvm_1
#define PRODMM Quotvm
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM_1 RQuotvm_1
#define PRODMM RQuotvm
#define X3 ,M,N,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, int A1, int A2
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXv
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotvM
#define RQUOTMM RQuotvM
#define X1 ,M,A1
#define X3 ,M,A1
#define Y , ptrdiff_t M, int A1
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotvM
#define X3 ,M,A1
#define Y , ptrdiff_t M, int A1
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotvM
#define X3 ,M,A1
#define Y , ptrdiff_t M, int A1
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXV
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 QuotVm_1
#define QUOTMM QuotVm
#define RQUOTMM_1 RQuotVm_1
#define RQUOTMM RQuotVm
#define X2 ,M,N,A
#define X3 ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotVm
#define X3 ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotVm
#define X3 ,M,N,A
#define Y , ptrdiff_t M, ptrdiff_t N, int A
#define GETM1 .getV()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

    //
    // Matrix / % Matrix
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class Quotmm_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Quotmm_1(
            const T TMV_DEBUGPARAM(_x), const SmallMatrix<T1,K,N,A1>& _m1,
            const SmallMatrix<T2,K,M,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x == T(1)); }
        const SmallMatrix<T1,K,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,M,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { DoLDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { DoLDiv(m2,m1,m0); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().LDiv(m1.view(),m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            m2.view().LDiv(m1.view(),m0);
        }
    private:
        const SmallMatrix<T1,K,N,A1>& m1;
        const SmallMatrix<T2,K,M,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class Quotmm : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        Quotmm(
            const T _x, const SmallMatrix<T1,K,N,A1>& _m1,
            const SmallMatrix<T2,K,M,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2) {}
        T getX() const { return x; }
        const SmallMatrix<T1,K,N,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,M,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().LDiv(m1.view(),m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            m2.view().LDiv(m1.view(),m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const SmallMatrix<T1,K,N,A1>& m1;
        const SmallMatrix<T2,K,M,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t K, ptrdiff_t M, int A2>
    class QuotMm_1 : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotMm_1(
            const T TMV_DEBUGPARAM(_x), const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,K,M,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x==T(1)); TMVAssert(m1.colsize() == K); }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return m1.rowsize(); }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,M,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m2.view().LDiv(m1,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            m2.view().LDiv(m1,m0);
        }
    private:
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,K,M,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t K, ptrdiff_t M, int A2>
    class QuotMm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotMm(
            const T _x, const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,K,M,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.colsize() == K); }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return m1.rowsize(); }
        T getX() const { return x; }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,K,M,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m2.view().LDiv(m1,m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            m2.view().LDiv(m1,m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,K,M,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t K, ptrdiff_t N, int A1>
    class QuotmM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        QuotmM(
            const T _x, const SmallMatrix<T1,K,N,A1>& _m1,
            const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2) {}
        ptrdiff_t colsize() const { return m2.rowsize(); }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const SmallMatrix<T1,K,N,A1>& getM1() const { return m1; }
        const GenMatrix<T2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.LDiv(m1.view(),m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            m2.LDiv(m1.view(),m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const SmallMatrix<T1,K,N,A1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class RQuotmm_1 : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotmm_1(
            const T TMV_DEBUGPARAM(_x), const SmallMatrix<T1,M,K,A1>& _m1,
            const SmallMatrix<T2,N,K,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x == T(1)); }
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,N,K,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { DoRDiv(m2,m1,m0); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { DoRDiv(m2,m1,m0); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().RDiv(m1.view(),m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            m2.view().RDiv(m1.view(),m0);
        }
    private:
        const SmallMatrix<T1,M,K,A1>& m1;
        const SmallMatrix<T2,N,K,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2>
    class RQuotmm : public SmallMatrixComposite<T,M,N>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotmm(
            const T _x, const SmallMatrix<T1,M,K,A1>& _m1,
            const SmallMatrix<T2,N,K,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2) {}
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const SmallMatrix<T2,N,K,A2>& getM2() const { return m2; }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|CStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|CStyle>& m0) const
        { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|CStyle>& m0) const
        { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,ColMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<real_type,M,N,RowMajor|FortranStyle>& m0) const
        { TMVAssert(isReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,ColMajor|FortranStyle>& m0) const
        { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignTom(
            SmallMatrix<complex_type,M,N,RowMajor|FortranStyle>& m0) const
        { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().RDiv(m1.view(),m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == N);
            m2.view().RDiv(m1.view(),m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const SmallMatrix<T1,M,K,A1>& m1;
        const SmallMatrix<T2,N,K,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, ptrdiff_t K, int A2>
    class RQuotMm_1 : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotMm_1(
            const T TMV_DEBUGPARAM(_x), const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,N,K,A2>& _m2) :
            m1(_m1), m2(_m2)
        { TMVAssert(_x==T(1)); TMVAssert(m1.rowsize() == K); }
        ptrdiff_t colsize() const { return m1.colsize(); }
        ptrdiff_t rowsize() const { return N; }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,N,K,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().RDiv(m1,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            m2.view().RDiv(m1,m0);
        }
    private:
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,N,K,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, ptrdiff_t K, int A2>
    class RQuotMm : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotMm(
            const T _x, const GenMatrix<T1>& _m1,
            const SmallMatrix<T2,N,K,A2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m1.rowsize() == K); }
        ptrdiff_t colsize() const { return m1.colsize(); }
        ptrdiff_t rowsize() const { return N; }
        T getX() const { return x; }
        const GenMatrix<T1>& getM1() const { return m1; }
        const SmallMatrix<T2,N,K,A2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            TMVAssert(isReal(T()));
            m2.view().RDiv(m1,m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
            m2.view().RDiv(m1,m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const GenMatrix<T1>& m1;
        const SmallMatrix<T2,N,K,A2>& m2;
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t K, int A1>
    class RQuotmM : public MatrixComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        RQuotmM(
            const T _x, const SmallMatrix<T1,M,K,A1>& _m1,
            const GenMatrix<T2>& _m2) :
            x(_x), m1(_m1), m2(_m2)
        { TMVAssert(m2.rowsize() == K); }
        ptrdiff_t colsize() const { return M; }
        ptrdiff_t rowsize() const { return m2.colsize(); }
        T getX() const { return x; }
        const SmallMatrix<T1,M,K,A1>& getM1() const { return m1; }
        const GenMatrix<T2>& getM2() const { return m2; }
        void assignToM(MatrixView<real_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            m2.RDiv(m1.view(),m0);
            MultXM(x,m0);
        }
        void assignToM(MatrixView<complex_type> m0) const
        {
            TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
            m2.RDiv(m1.view(),m0);
            MultXM(x,m0);
        }
    private:
        const T x;
        const SmallMatrix<T1,M,K,A1>& m1;
        const GenMatrix<T2>& m2;
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator/=(
        SmallMatrix<T,M,N,A1>& m1, const SmallMatrix<T,M,M,A2>& m2)
    {
        DoLDivEq(m2,m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator/=(
        SmallMatrix<CT,M,N,A1>& m1, const SmallMatrix<T,M,M,A2>& m2)
    {
        DoLDivEq(m2,m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator%=(
        SmallMatrix<T,M,N,A1>& m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        DoRDivEq(m2,m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator%=(
        SmallMatrix<CT,M,N,A1>& m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        DoRDivEq(m2,m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator/=(
        SmallMatrix<T,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == M && m2.colsize() == M);
        m2.LDivEq(m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator/=(
        SmallMatrix<CT,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == M && m2.colsize() == M);
        m2.LDivEq(m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator%=(
        SmallMatrix<T,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == N && m2.colsize() == N);
        m2.RDivEq(m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator%=(
        SmallMatrix<CT,M,N,A1>& m1, const GenMatrix<T>& m2)
    {
        TMVAssert(m2.rowsize() == N && m2.colsize() == N);
        m2.RDivEq(m1.view());
        return m1;
    }

    template <typename T, ptrdiff_t M, int A2>
    inline MatrixView<T> operator/=(
        MatrixView<T> m1, const SmallMatrix<T,M,M,A2>& m2)
    {
        TMVAssert(m1.colsize() == M);
        m2.view().LDivEq(m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, int A2>
    inline MatrixView<CT> operator/=(
        MatrixView<CT> m1, const SmallMatrix<T,M,M,A2>& m2)
    {
        TMVAssert(m1.colsize() == M);
        m2.view().LDivEq(m1);
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<T> operator%=(
        MatrixView<T> m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        TMVAssert(m1.rowsize() == N);
        m2.view().RDivEq(m1);
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator%=(
        MatrixView<CT> m1, const SmallMatrix<T,N,N,A2>& m2)
    {
        TMVAssert(m1.rowsize() == N);
        m2.view().RDivEq(m1);
        return m1;
    }

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const QuotXm_1<T,Tm,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),m1);
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const QuotXm_1<T,T,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),m1);
        return m1;
    }

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const QuotXm<T,Tm,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),m1);
        MultXV<M*N>(qxm.getX(),m1.ptr());
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1, int A2>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const QuotXm<T,T,N,N,A2>& qxm)
    {
        DoRDivEq(qxm.getM(),m1);
        MultXV<M*N>(qxm.getX(),m1.ptr());
        return m1;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const QuotXm_1<T,Tm,N,N,A2>& qxm)
    {
        TMVAssert(m1.rowsize() == N);
        qxm.getM().view().RDivEq(m1);
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const QuotXm_1<T,T,N,N,A2>& qxm)
    {
        TMVAssert(m1.rowsize() == N);
        qxm.getM().view().RDivEq(m1);
        return m1;
    }

    template <typename T, typename Tm, ptrdiff_t N, int A2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m1, const QuotXm<T,Tm,N,N,A2>& qxm)
    {
        TMVAssert(m1.rowsize() == N);
        qxm.getM().view().RDivEq(m1);
        m1 *= qxm.getX();
        return m1;
    }

    template <typename T, ptrdiff_t N, int A2>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m1, const QuotXm<T,T,N,N,A2>& qxm)
    {
        TMVAssert(m1.rowsize() == N);
        qxm.getM().view().RDivEq(m1);
        m1 *= qxm.getX();
        return m1;
    }

    template <typename T, typename Tm, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<T,M,N,A1>& operator*=(
        SmallMatrix<T,M,N,A1>& m1, const QuotXM<T,Tm>& qxm)
    {
        TMVAssert(qxm.getM().rowsize() == N);
        TMVAssert(qxm.getM().colsize() == N);
        qxm.getM().RDivEq(m1.view());
        m1 *= qxm.getX();
        return m1;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A1>
    inline SmallMatrix<CT,M,N,A1>& operator*=(
        SmallMatrix<CT,M,N,A1>& m1, const QuotXM<T,T>& qxm)
    {
        TMVAssert(qxm.getM().rowsize() == N);
        TMVAssert(qxm.getM().colsize() == N);
        qxm.getM().RDivEq(m1.view());
        m1 *= qxm.getX();
        return m1;
    }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 Quotmm_1
#define QUOTMM Quotmm
#define RQUOTMM_1 RQuotmm_1
#define RQUOTMM RQuotmm
#define X1 ,K,N,A1
#define X1b ,M,K,A1
#define X2 ,K,M,A2
#define X2b ,N,K,A2
#define X3 ,M,N,K,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM_1 Quotmm_1
#define PRODMM Quotmm
#define X3 ,M,N,K,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM_1 RQuotmm_1
#define PRODMM RQuotmm
#define X3 ,M,N,K,A1,A2
#define Y , ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A1, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotmM
#define RQUOTMM RQuotmM
#define X1 ,K,N,A1
#define X1b ,M,K,A1
#define X3 ,K,N,A1
#define X3b ,M,K,A1
#define Y , ptrdiff_t K, ptrdiff_t N, int A1
#define Yb , ptrdiff_t M, ptrdiff_t K, int A1
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotmM
#define X3 ,K,N,A1
#define Y , ptrdiff_t K, ptrdiff_t N, int A1
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotmM
#define X3 ,M,K,A1
#define Y , ptrdiff_t M, ptrdiff_t K, int A1
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 QuotMm_1
#define QUOTMM QuotMm
#define RQUOTMM_1 RQuotMm_1
#define RQUOTMM RQuotMm
#define X2 ,K,M,A2
#define X2b ,N,K,A2
#define X3 ,K,M,A2
#define X3b ,N,K,A2
#define Y , ptrdiff_t K, ptrdiff_t M, int A2
#define Yb , ptrdiff_t N, ptrdiff_t K, int A2
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotMm
#define X3 ,K,M,A2
#define Y , ptrdiff_t K, ptrdiff_t M, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotMm
#define X3 ,N,K,A2
#define Y , ptrdiff_t N, ptrdiff_t K, int A2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

    //
    // P * m
    //

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class ProdPm : public MatrixComposite<T>
    {
    public:
        inline ProdPm(const Permutation& _p, const SmallMatrix<T,M,N,A>& _m) :
            p(_p), m(_m)
        { TMVAssert(m.colsize()==p.rowsize()); }
        inline ptrdiff_t colsize() const { return p.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline const Permutation& getP() const { return p; }
        inline const SmallMatrix<T,M,N,A>& getM() const { return m; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.applyOnLeft(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            p.applyOnLeft(m0=m);
        }
    private:
        const Permutation& p;
        const SmallMatrix<T,M,N,A>& m;
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class ProdmP : public MatrixComposite<T>
    {
    public:
        inline ProdmP(const SmallMatrix<T,M,N,A>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert(m.rowsize()==p.colsize()); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return p.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const SmallMatrix<T,M,N,A>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.applyOnRight(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            p.applyOnRight(m0=m);
        }
    private:
        const SmallMatrix<T,M,N,A>& m;
        const Permutation& p;
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline ProdmP<T,M,N,A> operator*(
        const SmallMatrix<T,M,N,A>& m, const Permutation& p)
    { return ProdmP<T,M,N,A>(m,p); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline ProdPm<T,M,N,A> operator*(
        const Permutation& p, const SmallMatrix<T,M,N,A>& m)
    { return ProdPm<T,M,N,A>(p,m); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator*=(
        SmallMatrix<T,M,N,A>& m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.applyOnRight(m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class QuotmP : public MatrixComposite<T>
    {
    public:
        inline QuotmP(const SmallMatrix<T,M,N,A>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert( m.colsize() == p.colsize() ); }
        inline ptrdiff_t colsize() const { return p.rowsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const SmallMatrix<T,M,N,A>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.inverse().applyOnLeft(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            p.inverse().applyOnLeft(m0=m);
        }
    protected:
        const SmallMatrix<T,M,N,A>& m;
        const Permutation& p;
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class RQuotmP : public MatrixComposite<T>
    {
    public:
        inline RQuotmP(const SmallMatrix<T,M,N,A>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert( m.rowsize() == p.rowsize() ); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return p.colsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const SmallMatrix<T,M,N,A>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.inverse().applyOnRight(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            p.inverse().applyOnRight(m0=m);
        }
    protected:
        const SmallMatrix<T,M,N,A>& m;
        const Permutation& p;
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline QuotmP<T,M,N,A> operator/(
        const SmallMatrix<T,M,N,A>& m, const Permutation& p)
    { return QuotmP<T,M,N,A>(m,p); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator/=(
        SmallMatrix<T,M,N,A>& m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.inverse().applyOnLeft(m);
        return m;
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline RQuotmP<T,M,N,A> operator%(
        const SmallMatrix<T,M,N,A>& m, const Permutation& p)
    { return RQuotmP<T,M,N,A>(m,p); }

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline SmallMatrix<T,M,N,A>& operator%=(
        SmallMatrix<T,M,N,A>& m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.inverse().applyOnRight(m);
        return m;
    }


} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
