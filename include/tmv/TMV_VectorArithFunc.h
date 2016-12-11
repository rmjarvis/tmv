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


#ifndef TMV_VectorArithFunc_H
#define TMV_VectorArithFunc_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // v *= x
    template <typename T>
    void MultXV(const T x, VectorView<T> v2);
    // v2 = x * v1
    template <typename T, typename T1>
    void MultXV(const T x, const GenVector<T1>& v1, VectorView<T> v2);

    // v2 += x * v1
    template <typename T, typename T1>
    void AddVV(const T x, const GenVector<T1>& v1, VectorView<T> v2);
    // v3 = x1 * v1 + x2 * v2
    template <typename T, typename T1, typename T2>
    void AddVV(
        const T x1, const GenVector<T1>& v1,
        const T x2, const GenVector<T2>& v2, VectorView<T> v3);

    // v1 * v2 (dot product)
    // Note: the return type is the type of the first vector
    // This is important for mixing complex and real vectors
    template <typename T, typename T2>
    T MultVV(const GenVector<T>& v1, const GenVector<T2>& v2);

    template <bool add, typename T, typename Tx, typename Ty>
    void ElemMultVV(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, VectorView<T> z);

    template <typename T>
    class VectorComposite : public GenVector<T>
    {
    public:

        inline VectorComposite() {}
        inline VectorComposite(const VectorComposite<T>&) {}
        virtual inline ~VectorComposite() {}

        const T* cptr() const;
        inline ptrdiff_t step() const { return 1; }
        inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }

    private:
        mutable AlignedArray<T> _v;
    };

    // Specialize allowed complex combinations:
    template <typename T>
    inline void MultXV(const T x, VectorView<CT> v2)
    { MultXV(CT(x),v2); }
    template <typename T, typename T1>
    inline void MultXV(
        const T x, const GenVector<T1>& v1, VectorView<CT> v2)
    { MultXV(CT(x),v1,v2); }

    template <typename T, typename T1>
    inline void AddVV(
        const T x, const GenVector<T1>& v1, VectorView<CT> v2)
    { AddVV(CT(x),v1,v2); }
    template <typename T, typename T1, typename T2>
    inline void AddVV(
        const T x1, const GenVector<T1>& v1,
        const T x2, const GenVector<T2>& v2, VectorView<CT> v3)
    { AddVV(CT(x1),v1,CT(x2),v2,v3); }
    template <typename T>
    inline void AddVV(
        const CT x1, const GenVector<CT>& v1,
        const CT x2, const GenVector<T>& v2, VectorView<CT> v3)
    { AddVV(x2,v2,x1,v1,v3); }

    template <typename T>
    inline CT MultVV(const GenVector<T>& v1, const GenVector<CT>& v2)
    { return MultVV(v2,v1); }

    template <bool add, typename T, typename Tx, typename Ty>
    inline void ElemMultVV(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, VectorView<CT> z)
    { ElemMultVV<add>(CT(alpha),x,y,z); }

    // Specialize disallowed complex combinations:
    template <typename T>
    inline void MultXV(const CT , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void MultXV(
        const CT , const GenVector<Ta>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename T1>
    inline void AddVV(
        const CT , const GenVector<T1>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddVV(
        const CT , const GenVector<Ta>& ,
        const CT , const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultVV(
        const CT , const GenVector<Ta>& , const GenVector<Tb>& ,
        VectorView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
