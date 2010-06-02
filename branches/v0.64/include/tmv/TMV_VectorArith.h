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


#ifndef TMV_VectorArith_H
#define TMV_VectorArith_H

#include "tmv/TMV_VectorArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    // These are what we want to do no matter what type Tx is:
    template <class T, IndexStyle I, class Tx> 
    inline Vector<T>& operator+=(Vector<T,I>& v, const Tx& x)
    { v.view() += x; return v; }

    template <class T, IndexStyle I, class Tx> 
    inline Vector<T>& operator-=(Vector<T,I>& v, const Tx& x)
    { v.view() -= x; return v; }

    template <class T, IndexStyle I, class Tx> 
    inline Vector<T>& operator*=(Vector<T,I>& v, const Tx& x)
    { v.view() *= x; return v; }

    template <class T, IndexStyle I, class Tx> 
    inline Vector<T>& operator/=(Vector<T,I>& v, const Tx& x)
    { v.view() /= x; return v; }

    template <class T, IndexStyle I, class Tx> 
    inline Vector<T>& operator%=(Vector<T,I>& v, const Tx& x)
    { v.view() %= x; return v; }


    //
    // Vector * / Scalar
    //

    template <class T, class Tv> 
    class ProdXV : public VectorComposite<T> 
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) xCT;
    public:
        inline ProdXV(T _x, const GenVector<Tv>& _v) : x(_x), v(_v) {}
        inline size_t size() const { return v.size(); }
        inline T getX() const { return x; }
        inline const GenVector<Tv>& getV() const { return v; }
        inline void assignToV(const VectorView<RT>& v0) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(v0.size() == size());
            MultXV(x,v,v0);
        }
        inline void assignToV(const VectorView<xCT>& v0) const
        {
            TMVAssert(v0.size() == size());
            MultXV(x,v,v0);
        }
    private:
        const T x;
        const GenVector<Tv>& v;
    };

    template <class T> 
    inline const VectorView<T>& operator*=(const VectorView<T>& v1, T x2)
    { MultXV(x2,v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator*=(const VectorView<CT>& v1, T x2)
    { MultXV(T(x2),v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator*=(const VectorView<CT>& v1, CCT x2)
    { MultXV(CT(x2),v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator*=(const VectorView<CT>& v1, VCT x2)
    { MultXV(CT(x2),v1); return v1; }

    template <class T> 
    inline const VectorView<T>& operator/=(const VectorView<T>& v1, T x2)
    { MultXV(T(1)/x2,v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator/=(const VectorView<CT>& v1, T x2)
    { MultXV(T(1)/x2,v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator/=(const VectorView<CT>& v1, CCT x2)
    { MultXV(T(1)/CT(x2),v1); return v1; }

    template <class T> 
    inline const VectorView<CT>& operator/=(const VectorView<CT>& v1, VCT x2)
    { MultXV(T(1)/CT(x2),v1); return v1; }

#define GENMATRIX GenVector
#define PRODXM ProdXV
#define GETM .getV()
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM
    //
    // Vector + Vector
    //

    template <class T, class T1, class T2> 
    class SumVV : public VectorComposite<T> 
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) xCT;
    public:
        inline SumVV(
            T _x1, const GenVector<T1>& _v1, 
            T _x2, const GenVector<T2>& _v2
        ) :
            x1(_x1),v1(_v1),x2(_x2), v2(_v2)
        { TMVAssert(v1.size() == v2.size()); }
        inline size_t size() const { return v1.size(); }
        inline T getX1() const { return x1; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToV(const VectorView<RT>& v0) const
        {
            TMVAssert(v0.size() == v1.size());
            TMVAssert(isReal(T()));
            AddVV(x1,v1,x2,v2,v0);
        }
        inline void assignToV(const VectorView<xCT>& v0) const
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
    template <class T> 
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v1, const GenVector<T>& v2) 
    { 
        TMVAssert(v1.size() == v2.size());
        AddVV(T(1),v2,v1); 
        return v1; 
    }

    template <class T> 
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(1),v2,v1); 
        return v1; 
    }

    // v-=v
    template <class T> 
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(-1),v2,v1); 
        return v1; 
    }

    template <class T> 
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v1.size() == v2.size());
        AddVV(T(-1),v2,v1); 
        return v1; 
    }

    // v+=(x*v)
    template <class T, class T2> 
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v, const ProdXV<T,T2>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(pxv.getX(),pxv.getV(),v);
        return v; 
    }

    template <class T> 
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v, const ProdXV<T,T>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(pxv.getX(),pxv.getV(),v);
        return v; 
    }

    // v-=(x*v)
    template <class T, class T2> 
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v, const ProdXV<T,T2>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(-pxv.getX(),pxv.getV(),v);
        return v; 
    }

    template <class T> 
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v, const ProdXV<T,T>& pxv)
    {
        TMVAssert(v.size() == pxv.size());
        AddVV(-pxv.getX(),pxv.getV(),v);
        return v; 
    }

#define GENMATRIX1 GenVector
#define GENMATRIX2 GenVector
#define SUMMM SumVV
#define PRODXM1 ProdXV
#define PRODXM2 ProdXV
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxSumMM.h"
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

    //
    // Vector * Vector
    //

    template <class T> 
    inline T operator*(const GenVector<T>& v1, const GenVector<T>& v2) 
    { 
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1,v2); 
    }

    template <class T> 
    inline CT operator*(const GenVector<CT>& v1, const GenVector<T>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return MultVV(v1,v2); 
    }

    template <class T> 
    inline CT operator*(const GenVector<T>& v1, const GenVector<CT>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return MultVV(v2,v1); 
    }

    // v * (x*v)

    template <class T, class T2> 
    inline T operator*(const GenVector<T>& v1, const ProdXV<T,T2>& v2) 
    { 
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v1,v2.getV()); 
    }

    template <class T> 
    inline CT operator*(const GenVector<CT>& v1, const ProdXV<T,T>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v1,v2.getV()); 
    }

    template <class T, class T2> 
    inline CT operator*(const GenVector<T>& v1, const ProdXV<CT,T2>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return v2.getX()*MultVV(v2.getV(),v1); 
    }

    // (x*v) * v

    template <class T, class T2> 
    inline T operator*(const ProdXV<T,T2>& v1, const GenVector<T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v1.getV(),v2); 
    }

    template <class T> 
    inline CT operator*(const ProdXV<T,T>& v1, const GenVector<CT>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v2,v1.getV()); 
    }

    template <class T, class T2> 
    inline CT operator*(const ProdXV<CT,T2>& v1, const GenVector<T>& v2)
    { 
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*MultVV(v1.getV(),v2); 
    }

    // (x*v) * (x*v)

    template <class T, class T1, class T2> 
    inline T operator*(const ProdXV<T,T1>& v1, const ProdXV<T,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV()); 
    }

    template <class T, class T1> 
    inline CT operator*(const ProdXV<CT,T1>& v1, const ProdXV<T,T>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV()); 
    }

    template <class T, class T2> 
    inline CT operator*(const ProdXV<T,T>& v1, const ProdXV<CT,T2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV()); 
    }

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif 
