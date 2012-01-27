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


#ifndef TMV_VectorArithFunc_H
#define TMV_VectorArithFunc_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // v *= x
    template <class T> 
    void MultXV(const T x, VectorView<T> v2);
    // v2 = x * v1
    template <class T, class T1> 
    void MultXV(const T x, const GenVector<T1>& v1, VectorView<T> v2);

    // v2 += x * v1
    template <class T, class T1> 
    void AddVV(const T x, const GenVector<T1>& v1, VectorView<T> v2);
    // v3 = x1 * v1 + x2 * v2
    template <class T, class T1, class T2> 
    void AddVV(
        const T x1, const GenVector<T1>& v1,
        const T x2, const GenVector<T2>& v2, VectorView<T> v3);

    // v1 * v2 (dot product)
    // Note: the return type is the type of the first vector
    // This is important for mixing complex and real vectors
    template <class T, class T2> 
    T MultVV(const GenVector<T>& v1, const GenVector<T2>& v2);

    template <bool add, class T, class Tx, class Ty> 
    void ElemMultVV(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, VectorView<T> z);

    template <class T> 
    class VectorComposite : public GenVector<T>
    {
    public:

        inline VectorComposite() {}
        inline VectorComposite(const VectorComposite<T>&) {}
        virtual inline ~VectorComposite() {}

        const T* cptr() const;
        inline int step() const { return 1; }
        inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }

    private:
        mutable AlignedArray<T> _v;
    };

    // Specialize allowed complex combinations:
    template <class T> 
    inline void MultXV(const T x, VectorView<CT> v2)
    { MultXV(CT(x),v2); }
    template <class T, class T1> 
    inline void MultXV(
        const T x, const GenVector<T1>& v1, VectorView<CT> v2)
    { MultXV(CT(x),v1,v2); }

    template <class T, class T1> 
    inline void AddVV(
        const T x, const GenVector<T1>& v1, VectorView<CT> v2)
    { AddVV(CT(x),v1,v2); }
    template <class T, class T1, class T2> 
    inline void AddVV(
        const T x1, const GenVector<T1>& v1,
        const T x2, const GenVector<T2>& v2, VectorView<CT> v3)
    { AddVV(CT(x1),v1,CT(x2),v2,v3); }
    template <class T> 
    inline void AddVV(
        const CT x1, const GenVector<CT>& v1,
        const CT x2, const GenVector<T>& v2, VectorView<CT> v3)
    { AddVV(x2,v2,x1,v1,v3); }

    template <class T> 
    inline CT MultVV(const GenVector<T>& v1, const GenVector<CT>& v2)
    { return MultVV(v2,v1); }

    template <bool add, class T, class Tx, class Ty> 
    inline void ElemMultVV(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, VectorView<CT> z)
    { ElemMultVV<add>(CT(alpha),x,y,z); }

    // Specialize disallowed complex combinations:
    template <class T> 
    inline void MultXV(const CT , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta> 
    inline void MultXV(
        const CT , const GenVector<Ta>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T, class T1> 
    inline void AddVV(
        const CT , const GenVector<T1>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddVV(
        const CT , const GenVector<Ta>& ,
        const CT , const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultVV(
        const CT , const GenVector<Ta>& , const GenVector<Tb>& ,
        VectorView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif 
