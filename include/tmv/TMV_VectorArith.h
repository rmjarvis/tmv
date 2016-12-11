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


#ifndef TMV_VectorArith_H
#define TMV_VectorArith_H

#include "tmv/TMV_VectorArithFunc.h"
#include "tmv/TMV_VIt.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    // These are what we want to do no matter what type Tx is:
    template <typename T, int A, typename Tx>
    inline Vector<T,A>& operator+=(Vector<T,A>& v, const Tx& x)
    { v.view() += x; return v; }

    template <typename T, int A, typename Tx>
    inline Vector<T,A>& operator-=(Vector<T,A>& v, const Tx& x)
    { v.view() -= x; return v; }

    template <typename T, int A, typename Tx>
    inline Vector<T,A>& operator*=(Vector<T,A>& v, const Tx& x)
    { v.view() *= x; return v; }

    template <typename T, int A, typename Tx>
    inline Vector<T,A>& operator/=(Vector<T,A>& v, const Tx& x)
    { v.view() /= x; return v; }

    template <typename T, int A, typename Tx>
    inline Vector<T,A>& operator%=(Vector<T,A>& v, const Tx& x)
    { v.view() %= x; return v; }


    //
    // Vector * / Scalar
    //

    template <typename T, typename Tv>
    class ProdXV : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXV(T _x, const GenVector<Tv>& _v) : x(_x), v(_v) {}
        inline ptrdiff_t size() const { return v.size(); }
        inline T getX() const { return x; }
        inline const GenVector<Tv>& getV() const { return v; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(v0.size() == size());
            MultXV(x,v,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == size());
            MultXV(x,v,v0);
        }
    private:
        const T x;
        const GenVector<Tv>& v;
    };

    template <typename T>
    inline VectorView<T> operator*=(VectorView<T> v1, T x2)
    { MultXV(x2,v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator*=(VectorView<CT> v1, T x2)
    { MultXV(T(x2),v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator*=(VectorView<CT> v1, CCT x2)
    { MultXV(CT(x2),v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator*=(VectorView<CT> v1, VCT x2)
    { MultXV(CT(x2),v1); return v1; }

    template <typename T>
    inline VectorView<T> operator/=(VectorView<T> v1, T x2)
    { MultXV(TMV_InverseOf(x2),v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator/=(VectorView<CT> v1, T x2)
    { MultXV(TMV_InverseOf(x2),v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator/=(VectorView<CT> v1, CCT x2)
    { MultXV(TMV_InverseOf(CT(x2)),v1); return v1; }

    template <typename T>
    inline VectorView<CT> operator/=(VectorView<CT> v1, VCT x2)
    { MultXV(TMV_InverseOf(CT(x2)),v1); return v1; }

#define GENMATRIX GenVector
#define PRODXM ProdXV
#define GETM .getV()
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    //
    // Vector + Vector
    //

    template <typename T, typename T1, typename T2>
    class SumVV : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumVV(
            T _x1, const GenVector<T1>& _v1,
            T _x2, const GenVector<T2>& _v2) :
            x1(_x1),v1(_v1),x2(_x2), v2(_v2)
        { TMVAssert(v1.size() == v2.size()); }
        inline ptrdiff_t size() const { return v1.size(); }
        inline T getX1() const { return x1; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == v1.size());
            TMVAssert(isReal(T()));
            AddVV(x1,v1,x2,v2,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == v1.size());
            AddVV(x1,v1,x2,v2,v0);
        }
    private:
        const T x1;
        const GenVector<T1>& v1;
        const T x2;
        const GenVector<T2>& v2;
    };

    // v+=v
    template <typename T>
    inline VectorView<T> operator+=(VectorView<T> v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(1),v2,v1);
        return v1;
    }

    template <typename T>
    inline VectorView<CT> operator+=(VectorView<CT> v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(1),v2,v1);
        return v1;
    }

    // v-=v
    template <typename T>
    inline VectorView<T> operator-=(VectorView<T> v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(-1),v2,v1);
        return v1;
    }

    template <typename T>
    inline VectorView<CT> operator-=(VectorView<CT> v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(-1),v2,v1);
        return v1;
    }

    // v+=(x*v)
    template <typename T, typename T2>
    inline VectorView<T> operator+=(VectorView<T> v, const ProdXV<T,T2>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(pxv.getX(),pxv.getV(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator+=(VectorView<CT> v, const ProdXV<T,T>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(pxv.getX(),pxv.getV(),v);
        return v;
    }

    // v-=(x*v)
    template <typename T, typename T2>
    inline VectorView<T> operator-=(VectorView<T> v, const ProdXV<T,T2>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(-pxv.getX(),pxv.getV(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator-=(VectorView<CT> v, const ProdXV<T,T>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(-pxv.getX(),pxv.getV(),v);
        return v;
    }

#define SUMMM SumVV
#define GENMATRIX1 GenVector
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXV
#define PRODXM2 ProdXV
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxSumMM.h"
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxSumMMa.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

    //
    // Vector * Vector
    //

    template <typename T>
    inline T operator*(const GenVector<T>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1,v2);
    }

    template <typename T>
    inline CT operator*(const GenVector<CT>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1,v2);
    }

    template <typename T>
    inline CT operator*(const GenVector<T>& v1, const GenVector<CT>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return MultVV(v2,v1);
    }

    // v * (x*v)

    template <typename T, typename T2>
    inline T operator*(const GenVector<T>& v1, const ProdXV<T,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v1,v2.getV());
    }

    template <typename T>
    inline CT operator*(const GenVector<CT>& v1, const ProdXV<T,T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v1,v2.getV());
    }

    template <typename T, typename T2>
    inline CT operator*(const GenVector<T>& v1, const ProdXV<CT,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v2.getV(),v1);
    }

    // (x*v) * v

    template <typename T, typename T2>
    inline T operator*(const ProdXV<T,T2>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v1.getV(),v2);
    }

    template <typename T>
    inline CT operator*(const ProdXV<T,T>& v1, const GenVector<CT>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v2,v1.getV());
    }

    template <typename T, typename T2>
    inline CT operator*(const ProdXV<CT,T2>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v1.getV(),v2);
    }

    // (x*v) * (x*v)

    template <typename T, typename T1, typename T2>
    inline T operator*(const ProdXV<T,T1>& v1, const ProdXV<T,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV());
    }

    template <typename T, typename T1>
    inline CT operator*(const ProdXV<CT,T1>& v1, const ProdXV<T,T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV());
    }

    template <typename T, typename T2>
    inline CT operator*(const ProdXV<T,T>& v1, const ProdXV<CT,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV());
    }

    //
    // Element Product Vector * Vector
    //

    template <typename T, typename T1, typename T2>
    class ElemProdVV : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdVV(
            T _x, const GenVector<T1>& _v1, const GenVector<T2>& _v2) :
            x(_x), v1(_v1), v2(_v2)
        { TMVAssert(v1.size() == v2.size()); }
        inline ptrdiff_t size() const { return v1.size(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToV(VectorView<real_type> v0) const
        {
            TMVAssert(v0.size() == v1.size());
            TMVAssert(isReal(T()));
            ElemMultVV<false>(x,v1,v2,v0);
        }
        inline void assignToV(VectorView<complex_type> v0) const
        {
            TMVAssert(v0.size() == v1.size());
            ElemMultVV<false>(x,v1,v2,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v1;
        const GenVector<T2>& v2;
    };

    template <typename T, typename T2, typename T3>
    inline VectorView<T> operator+=(
        VectorView<T> v, const ElemProdVV<T,T2,T3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator+=(
        VectorView<CT> v, const ElemProdVV<T,T,T>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2(),v);
        return v;
    }

    template <typename T, typename T2, typename T3>
    inline VectorView<T> operator-=(
        VectorView<T> v, const ElemProdVV<T,T2,T3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator-=(
        VectorView<CT> v, const ElemProdVV<T,T,T>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2(),v);
        return v;
    }

#define PRODMM ElemProdVV
#define GENMATRIX1 GenVector
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXV
#define PRODXM2 ProdXV
#define OP ElemProd
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
