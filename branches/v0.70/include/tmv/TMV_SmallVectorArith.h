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


#ifndef TMV_SmallVectorArith_H
#define TMV_SmallVectorArith_H

#include "tmv/TMV_SmallVectorArithFunc.h"
#include "tmv/TMV_VectorArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

    template <class T, int N> 
    class SmallVectorComposite : public VectorComposite<T>
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SmallVectorComposite() {}
        inline SmallVectorComposite(const SmallVectorComposite<T,N>&) {}
        virtual inline ~SmallVectorComposite() {}

        inline int size() const { return N; }
        virtual void assignTov(SmallVector<real_type,N,CStyle>& v) const = 0;
        virtual void assignTov(SmallVector<TMV_ComplexType(T),N,CStyle>& v) const = 0;
        virtual void assignTov(SmallVector<real_type,N,FortranStyle>& v) const = 0;
        virtual void assignTov(SmallVector<TMV_ComplexType(T),N,FortranStyle>& v) const = 0;
        virtual void assignToV(const VectorView<real_type>& v) const = 0;
        virtual void assignToV(const VectorView<TMV_ComplexType(T)>& v) const = 0;
    };


    //
    // Vector * / Scalar
    //

    template <class T, class T1, int N, IndexStyle I> 
    class ProdXv : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ProdXv(T _x, const SmallVector<T1,N,I>& _v) : x(_x), v(_v) {}
        inline T getX() const { return x; }
        inline const SmallVector<T1,N,I>& getV() const { return v; }
        inline void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); MultXV<N>(x,v.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { MultXV<N>(x,v.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); MultXV<N>(x,v.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { MultXV<N>(x,v.cptr(),v0.ptr()); }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v0.step() == 1) 
                MultXV<N>(x,v.cptr(),v0.ptr());
            else
                MultXV(x,v.view(),v0); 
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { 
            if (v0.step() == 1) {
                MultXV<N>(x,v.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                MultXV(x,v.view(),v0);
        }
    private:
        const T x;
        const SmallVector<T1,N,I>& v;
    };

    template <class T, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator*=(SmallVector<T,N,I>& v1, T x2)
    { MultXV<N>(x2,v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, T x2)
    { MultXV<N>(x2,v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, CCT x2)
    { MultXV<N>(CT(x2),v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, VCT x2)
    { MultXV<N>(CT(x2),v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator/=(SmallVector<T,N,I>& v1, T x2)
    { MultXV<N>(T(1)/x2,v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, T x2)
    { MultXV<N>(T(1)/x2,v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, CCT x2)
    { MultXV<N>(T(1)/CT(x2),v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, VCT x2)
    { MultXV<N>(T(1)/CT(x2),v1.ptr()); return v1; }

#define GENMATRIX SmallVector
#define PRODXM ProdXv
#define X ,N,I
#define Y ,int N, IndexStyle I
#define GETM .getV()
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    //
    // Vector + Vector
    //

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv_1_1 : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv_1_1(
            T TMV_DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
            T TMV_DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
            v1(_v1),v2(_v2)
        { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignToV(const VectorView<real_type>& v0) const
        {
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
            else
                AddVV(T(1),v1.view(),T(1),v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        {
            if (v0.step() == 1) {
                AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(T(1),v1.view(),T(1),v2.view(),v0);
        }
    private:
        const SmallVector<T1,N,I1>& v1;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv_1_m1 : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv_1_m1(
            T TMV_DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
            T TMV_DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
            v1(_v1),v2(_v2)
        { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
            else
                AddVV(T(1),v1.view(),T(-1),v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { 
            if (v0.step() == 1) {
                AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(T(1),v1.view(),T(-1),v2.view(),v0);
        }
    private:
        const SmallVector<T1,N,I1>& v1;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv_1_x : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv_1_x(
            T TMV_DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
            T _x2, const SmallVector<T2,N,I2>& _v2) :
            v1(_v1),x2(_x2),v2(_v2)
        { TMVAssert(_x1 == T(1)); }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
            else
                AddVV(T(1),v1.view(),x2,v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        {
            if (v0.step() == 1) {
                AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(T(1),v1.view(),x2,v2.view(),v0);
        }
    private:
        const SmallVector<T1,N,I1>& v1;
        const T x2;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv_x_1 : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv_x_1(
            T _x1, const SmallVector<T1,N,I1>& _v1, 
            T TMV_DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
            x1(_x1),v1(_v1),v2(_v2)
        { TMVAssert(_x2 == T(1)); }
        inline T getX1() const { return x1; }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(SmallVector<real_type,N,CStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,CStyle>& v0) const
        { AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<real_type,N,FortranStyle>& v0) const
        { TMVAssert(isReal(T())); AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignTov(SmallVector<complex_type,N,FortranStyle>& v0) const
        { AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
            else
                AddVV(x1,v1.view(),T(1),v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { 
            if (v0.step() == 1) {
                AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(x1,v1.view(),T(1),v2.view(),v0);
        }
    private:
        const T x1;
        const SmallVector<T1,N,I1>& v1;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv_x_m1 : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv_x_m1(
            T _x1, const SmallVector<T1,N,I1>& _v1, 
            T TMV_DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
            x1(_x1),v1(_v1),v2(_v2)
        { TMVAssert(_x2 == T(-1)); }
        inline T getX1() const { return x1; }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(
            SmallVector<real_type,N,CStyle>& v0) const
        { 
            TMVAssert(isReal(T()));
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<complex_type,N,CStyle>& v0) const
        { 
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<real_type,N,FortranStyle>& v0) const
        { 
            TMVAssert(isReal(T()));
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<complex_type,N,FortranStyle>& v0) const
        {
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        }
        inline void assignToV(const VectorView<real_type>& v0) const
        {
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
            else
                AddVV(x1,v1.view(),T(-1),v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        {
            if (v0.step() == 1) {
                AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(x1,v1.view(),T(-1),v2.view(),v0);
        }
    private:
        const T x1;
        const SmallVector<T1,N,I1>& v1;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class Sumvv : public SmallVectorComposite<T,N> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline Sumvv(
            T _x1, const SmallVector<T1,N,I1>& _v1, 
            T _x2, const SmallVector<T2,N,I2>& _v2) :
            x1(_x1),v1(_v1),x2(_x2), v2(_v2) {}
        inline T getX1() const { return x1; }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignTov(
            SmallVector<real_type,N,CStyle>& v0) const
        { 
            TMVAssert(isReal(T()));
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<complex_type,N,CStyle>& v0) const
        {
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<real_type,N,FortranStyle>& v0) const
        {
            TMVAssert(isReal(T()));
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        }
        inline void assignTov(
            SmallVector<complex_type,N,FortranStyle>& v0) const
        {
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        }
        inline void assignToV(const VectorView<real_type>& v0) const
        {
            TMVAssert(isReal(T()));
            if (v0.step() == 1)
                AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
            else
                AddVV(x1,v1.view(),x2,v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        {
            if (v0.step() == 1) {
                AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else
                AddVV(x1,v1.view(),x2,v2.view(),v0);
        }
    private:
        const T x1;
        const SmallVector<T1,N,I1>& v1;
        const T x2;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I> 
    class SumVv : public VectorComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumVv(
            T _x1, const GenVector<T1>& _v1, 
            T _x2, const SmallVector<T2,N,I>& _v2) :
            x1(_x1),v1(_v1),x2(_x2), v2(_v2) { TMVAssert(v1.size() == N); }
        inline int size() const { return N; }
        inline T getX1() const { return x1; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const SmallVector<T2,N,I>& getV2() const { return v2; }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v1.step() == 1 && v0.step() == 1 && !SameStorage(v0,v1)) {
                if (x1 == T(1)) 
                    if (x2 == T(1))
                        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
                else 
                    if (x2 == T(1))
                        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
            }
            else AddVV(x1,v1,x2,v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { 
            if (v1.step() == 1 && v0.step() == 1 && !v1.isconj() &&
                !SameStorage(v0,v1)) {
                if (x1 == T(1)) 
                    if (x2 == T(1))
                        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
                else 
                    if (x2 == T(1))
                        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else AddVV(x1,v1,x2,v2.view(),v0);
        }
    private:
        const T x1;
        const GenVector<T1>& v1;
        const T x2;
        const SmallVector<T2,N,I>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I> 
    class SumvV : public VectorComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline SumvV(
            T _x1, const SmallVector<T1,N,I>& _v1, 
            T _x2, const GenVector<T2>& _v2) :
            x1(_x1),v1(_v1),x2(_x2), v2(_v2) { TMVAssert(v2.size() == N); }
        inline int size() const { return N; }
        inline T getX1() const { return x1; }
        inline const SmallVector<T1,N,I>& getV1() const { return v1; }
        inline T getX2() const { return x2; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            if (v2.step() == 1 && v0.step() == 1 && !SameStorage(v0,v2)) {
                if (x1 == T(1)) 
                    if (x2 == T(1))
                        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
                else 
                    if (x2 == T(1))
                        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
            }
            else AddVV(x1,v1.view(),x2,v2,v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { 
            if (v2.step() == 1 && v0.step() == 1 && !v2.isconj() &&
                !SameStorage(v0,v2) ) {
                if (x1 == T(1)) 
                    if (x2 == T(1))
                        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
                else 
                    if (x2 == T(1))
                        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else if (x2 == T(-1))
                        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
                    else
                        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
                if (v0.isconj()) v0.conjugateSelf();
            }
            else AddVV(x1,v1.view(),x2,v2,v0);
        }
    private:
        const T x1;
        const SmallVector<T1,N,I>& v1;
        const T x2;
        const GenVector<T2>& v2;
    };

    // v+=v
    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<T,N,I1>& operator+=(
        SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
    { AddVV_1<N>(v2.cptr(),v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<CT,N,I1>& operator+=(
        SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
    { AddVV_1<N>(v2.cptr(),v1.ptr()); return v1; }

    // v-=v
    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<T,N,I1>& operator-=(
        SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2)
    { AddVV_m1<N>(v2.cptr(),v1.ptr()); return v1; }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<CT,N,I1>& operator-=(
        SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
    { AddVV_m1<N>(v2.cptr(),v1.ptr()); return v1; }

    // v+=(x*v)
    template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<T,N,I1>& operator+=(
        SmallVector<T,N,I1>& v, const ProdXv<T,T2,N,I2>& v2)
    { AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr()); return v; }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<CT,N,I1>& operator+=(
        SmallVector<CT,N,I1>& v, const ProdXv<T,T,N,I2>& v2)
    { AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr()); return v; }

    // v-=(x*v)
    template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<T,N,I1>& operator-=(
        SmallVector<T,N,I1>& v, const ProdXv<T,T2,N,I2>& v2)
    { AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr()); return v; }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline SmallVector<CT,N,I1>& operator-=(
        SmallVector<CT,N,I1>& v, const ProdXv<T,T,N,I2>& v2)
    { AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr()); return v; }

    // Mix with Vector
    template <class T, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.step() == 1 && !v2.isconj())
            AddVV_1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(1),v2,v1.view()); 
        return v1; 
    }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.step() == 1)
            AddVV_1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(1),v2,v1.view()); 
        return v1; 
    }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.step() == 1 && !v2.isconj())
            AddVV_m1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(-1),v2,v1.view()); 
        return v1; 
    }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v1, const GenVector<T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.step() == 1)
            AddVV_m1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(-1),v2,v1.view()); 
        return v1; 
    }

    template <class T, class Tv> 
    class ProdXV;

    template <class T, class T2, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v, const ProdXV<T,T2>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.getV().step() == 1 && !v2.getV().isconj())
            AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(v2.getX(),v2.getV(),v.view()); 
        return v; 
    }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v, const ProdXV<T,T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.getV().step() == 1 && !v2.getV().isconj())
            AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(v2.getX(),v2.getV(),v.view()); 
        return v; 
    }

    template <class T, class T2, int N, IndexStyle I> 
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v, const ProdXV<T,T2>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.getV().step() == 1 && !v2.getV().isconj())
            AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(-v2.getX(),v2.getV(),v.view()); 
        return v; 
    }

    template <class T, int N, IndexStyle I> 
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v, const ProdXV<T,T>& v2) 
    {
        TMVAssert(v2.size() == N);
        if (v2.getV().step() == 1 && !v2.getV().isconj())
            AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(-v2.getX(),v2.getV(),v.view()); 
        return v; 
    }

    template <class T, int N, IndexStyle I> 
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v1, const SmallVector<T,N,I>& v2) 
    { 
        if (v1.step() == 1 && !v1.isconj())
            AddVV_1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(1),v2.view(),v1);
        return v1; 
    }

    template <class T, int N, IndexStyle I> 
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v1, const SmallVector<T,N,I>& v2) 
    {
        if (v1.step() == 1)
            AddVV_1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(1),v2.view(),v1);
        return v1;
    }

    // v-=v
    template <class T, int N, IndexStyle I> 
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v1, const SmallVector<T,N,I>& v2)
    { 
        if (v1.step() == 1 && !v1.isconj())
            AddVV_m1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(-1),v2.view(),v1);
        return v1;
    }

    template <class T, int N, IndexStyle I> 
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v1, const SmallVector<T,N,I>& v2) 
    { 
        if (v1.step() == 1)
            AddVV_m1<N>(v2.cptr(),v1.ptr());
        else
            AddVV(T(-1),v2.view(),v1);
        return v1;
    }

    // v+=(x*v)
    template <class T, class T2, int N, IndexStyle I> 
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v, const ProdXv<T,T2,N,I>& v2)
    {
        if (v.step() == 1 && !v.isconj())
            AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(v2.getX(),v2.getV().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I> 
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v, const ProdXv<T,T,N,I>& v2)
    { 
        if (v.step() == 1)
            AddVV<N>(v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(v2.getX(),v2.getV().view(),v);
        return v;
    }

    // v-=(x*v)
    template <class T, class T2, int N, IndexStyle I> 
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v, const ProdXv<T,T2,N,I>& v2)
    { 
        if (v.step() == 1 && !v.isconj())
            AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(-v2.getX(),v2.getV().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I> 
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v, const ProdXv<T,T,N,I>& v2)
    {
        if (v.step() == 1)
            AddVV<N>(-v2.getX(),v2.getV().cptr(),v.ptr());
        else
            AddVV(-v2.getX(),v2.getV().view(),v);
        return v;
    }

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallVector
#define SUMMM Sumvv
#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,N,I1
#define X2 ,N,I2
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxSumMM.h"
    // These get undef'ed in TMV_AuxSumMM.h
#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenVector
#define SUMMM SumvV
#define PRODXM1 ProdXv
#define PRODXM2 ProdXV
#define X1 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallVector
#define SUMMM SumVv
#define PRODXM1 ProdXV
#define PRODXM2 ProdXv
#define X2 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

    //
    // Element Product Vector * Vector
    //

    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    class ElemProdvv : public VectorComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdvv(
            T _x, const SmallVector<T1,N,I1>& _v1, 
            const SmallVector<T2,N,I2>& _v2) :
            x(_x),v1(_v1),v2(_v2) {}
        inline int size() const { return N; }
        inline T getX() const { return x; }
        inline const SmallVector<T1,N,I1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I2>& getV2() const { return v2; }
        inline void assignToV(const VectorView<real_type>& v0) const
        {
            TMVAssert(isReal(T()));
            ElemMultVV<false>(x,v1.view(),v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { ElemMultVV<false>(x,v1.view(),v2.view(),v0); }
    private:
        const T x;
        const SmallVector<T1,N,I1>& v1;
        const SmallVector<T2,N,I2>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I> 
    class ElemProdVv : public VectorComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdVv(
            T _x, const GenVector<T1>& _v1, const SmallVector<T2,N,I>& _v2) :
            x(_x),v1(_v1),v2(_v2) { TMVAssert(v1.size() == N); }
        inline int size() const { return N; }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV1() const { return v1; }
        inline const SmallVector<T2,N,I>& getV2() const { return v2; }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            ElemMultVV<false>(x,v1,v2.view(),v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { ElemMultVV<false>(x,v1,v2.view(),v0); }
    private:
        const T x;
        const GenVector<T1>& v1;
        const SmallVector<T2,N,I>& v2;
    };

    template <class T, class T1, class T2, int N, IndexStyle I> 
    class ElemProdvV : public VectorComposite<T> 
    {
    public:
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        inline ElemProdvV(
            T _x, const SmallVector<T1,N,I>& _v1, const GenVector<T2>& _v2) :
            x(_x),v1(_v1),v2(_v2) { TMVAssert(v2.size() == N); }
        inline int size() const { return N; }
        inline T getX() const { return x; }
        inline const SmallVector<T1,N,I>& getV1() const { return v1; }
        inline const GenVector<T2>& getV2() const { return v2; }
        inline void assignToV(const VectorView<real_type>& v0) const
        { 
            TMVAssert(isReal(T()));
            ElemMultVV<false>(x,v1.view(),v2,v0);
        }
        inline void assignToV(const VectorView<complex_type>& v0) const
        { ElemMultVV<false>(x,v1.view(),v2,v0); }
    private:
        const T x;
        const SmallVector<T1,N,I>& v1;
        const GenVector<T2>& v2;
    };

    template <class T, int N, IndexStyle I, class T2, class T3>
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v, const ElemProdVV<T,T2,T3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I>
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v, const ElemProdVV<T,T,T>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, class T3>
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v, const ElemProdVV<T,T2,T3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I>
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v, const ElemProdVV<T,T,T>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, IndexStyle I2, class T3>
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v, const ElemProdvV<T,T2,T3,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I2>
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v, const ElemProdvV<T,T,T,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, IndexStyle I2, class T3>
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v, const ElemProdvV<T,T2,T3,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I2>
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v, const ElemProdvV<T,T,T,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, class T3, IndexStyle I3>
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v, const ElemProdVv<T,T2,T3,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I3>
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v, const ElemProdVv<T,T,T,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, class T3, IndexStyle I3>
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v, const ElemProdVv<T,T2,T3,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I3>
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v, const ElemProdVv<T,T,T,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, IndexStyle I2, class T3, IndexStyle I3>
    inline SmallVector<T,N,I>& operator+=(
        SmallVector<T,N,I>& v, const ElemProdvv<T,T2,T3,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I2, IndexStyle I3>
    inline SmallVector<CT,N,I>& operator+=(
        SmallVector<CT,N,I>& v, const ElemProdvv<T,T,T,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, class T2, IndexStyle I2, class T3, IndexStyle I3>
    inline SmallVector<T,N,I>& operator-=(
        SmallVector<T,N,I>& v, const ElemProdvv<T,T2,T3,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, IndexStyle I, IndexStyle I2, IndexStyle I3>
    inline SmallVector<CT,N,I>& operator-=(
        SmallVector<CT,N,I>& v, const ElemProdvv<T,T,T,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v.view());
        return v;
    }

    template <class T, int N, class T2, IndexStyle I2, class T3>
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v, const ElemProdvV<T,T2,T3,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2(),v);
        return v;
    }

    template <class T, int N, IndexStyle I2>
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v, const ElemProdvV<T,T,T,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2(),v);
        return v;
    }

    template <class T, int N, class T2, IndexStyle I2, class T3>
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v, const ElemProdvV<T,T2,T3,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2(),v);
        return v;
    }

    template <class T, int N, IndexStyle I2>
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v, const ElemProdvV<T,T,T,N,I2>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2(),v);
        return v;
    }

    template <class T, int N, class T2, class T3, IndexStyle I3>
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v, const ElemProdVv<T,T2,T3,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I3>
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v, const ElemProdVv<T,T,T,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, class T2, class T3, IndexStyle I3>
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v, const ElemProdVv<T,T2,T3,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I3>
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v, const ElemProdVv<T,T,T,N,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, class T2, IndexStyle I2, class T3, IndexStyle I3>
    inline const VectorView<T>& operator+=(
        const VectorView<T>& v, const ElemProdvv<T,T2,T3,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I2, IndexStyle I3>
    inline const VectorView<CT>& operator+=(
        const VectorView<CT>& v, const ElemProdvv<T,T,T,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, class T2, IndexStyle I2, class T3, IndexStyle I3>
    inline const VectorView<T>& operator-=(
        const VectorView<T>& v, const ElemProdvv<T,T2,T3,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v);
        return v;
    }

    template <class T, int N, IndexStyle I2, IndexStyle I3>
    inline const VectorView<CT>& operator-=(
        const VectorView<CT>& v, const ElemProdvv<T,T,T,N,I2,I3>& pvv)
    {
        TMVAssert(v.size() == pvv.size());
        ElemMultVV<true>(-pvv.getX(),pvv.getV1().view(),pvv.getV2().view(),v);
        return v;
    }


#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallVector
#define PRODMM ElemProdvv
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,N,I1
#define X2 ,N,I2
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM1 .getV()
#define GETM2 .getV()
#define OP ElemProd
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM1 .getV1()
#define GETM2 .getV2()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenVector
#define PRODMM ElemProdvV
#define PRODXM1 ProdXv
#define PRODXM2 ProdXV
#define OP ElemProd
#define X1 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallVector
#define PRODMM ElemProdVv
#define PRODXM1 ProdXV
#define PRODXM2 ProdXv
#define OP ElemProd
#define X2 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM1 .getV()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

    //
    // Vector * Vector
    //

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline T operator*(
        const SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
    { return MultVV<N>(v1.cptr(),v2.cptr()); }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2)
    { return MultVV<N>(v1.cptr(),v2.cptr()); }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const SmallVector<T,N,I1>& v1, const SmallVector<CT,N,I2>& v2)
    { return MultVV<N>(v2.cptr(),v1.cptr()); }

    // v * (x*v)
    template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline T operator*(
        const SmallVector<T,N,I1>& v1, const ProdXv<T,T2,N,I2>& v2) 
    { return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr()); }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const SmallVector<CT,N,I1>& v1, const ProdXv<T,T,N,I2>& v2)
    { return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr()); }

    template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const SmallVector<T,N,I1>& v1, const ProdXv<CT,T2,N,I2>& v2)
    { return v2.getX()*MultVV<N>(v2.getV().cptr(),v1.cptr()); }

    // (x*v) * v
    template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
    inline T operator*(
        const ProdXv<T,T1,N,I1>& v1, const SmallVector<T,N,I2>& v2)
    { return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr()); }

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const ProdXv<T,T,N,I1>& v1, const SmallVector<CT,N,I2>& v2)
    { return v1.getX()*MultVV<N>(v2.cptr(),v1.getV().cptr()); }

    template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const ProdXv<CT,T1,N,I1>& v1, const SmallVector<T,N,I2>& v2)
    { return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr()); }

    // (x*v) * (x*v)
    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline T operator*(
        const ProdXv<T,T1,N,I1>& v1, const ProdXv<T,T2,N,I2>& v2)
    { return v1.getX()*v2.getX()*MultVV<N>(v1.getV().cptr(),v2.getV().cptr()); }

    template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const ProdXv<CT,T1,N,I1>& v1, const ProdXv<T,T,N,I2>& v2)
    { return v1.getX()*v2.getX()*MultVV<N>(v1.getV().cptr(),v2.getV().cptr()); }

    template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline CT operator*(
        const ProdXv<T,T,N,I1>& v1, const ProdXv<CT,T2,N,I2>& v2)
    { return v1.getX()*v2.getX()*MultVV<N>(v1.getV().cptr(),v2.getV().cptr()); }


    // Mix with Vector:

    // v * v
    template <class T, int N, IndexStyle I> 
    inline T operator*(const GenVector<T>& v1, const SmallVector<T,N,I>& v2) 
    { 
        if (v1.step() == 1 && !v1.isconj())
            return MultVV<N>(v1.cptr(),v2.cptr());
        else 
            return MultVV(v1,v2.view());
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const GenVector<CT>& v1, const SmallVector<T,N,I>& v2)
    { 
        if (v1.step() == 1 && !v1.isconj())
            return MultVV<N>(v1.cptr(),v2.cptr());
        else 
            return MultVV(v1,v2.view()); 
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const GenVector<T>& v1, const SmallVector<CT,N,I>& v2)
    {
        if (v1.step() == 1 && !v1.isconj())
            return MultVV<N>(v2.cptr(),v1.cptr());
        else 
            return MultVV(v2.view(),v1);
    }

    // v * (x*v)
    template <class T, class T2, int N, IndexStyle I> 
    inline T operator*(const GenVector<T>& v1, const ProdXv<T,T2,N,I>& v2) 
    { 
        if (v1.step() == 1 && !v1.isconj())
            return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr());
        else 
            return v2.getX()*MultVV(v1,v2.getV().view());
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const GenVector<CT>& v1, const ProdXv<T,T,N,I>& v2)
    {
        if (v1.step() == 1 && !v1.isconj())
            return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr());
        else 
            return v2.getX()*MultVV(v1,v2.getV().view());
    }

    template <class T, class T2, int N, IndexStyle I> 
    inline CT operator*(const GenVector<T>& v1, const ProdXv<CT,T2,N,I>& v2)
    {
        if (v1.step() == 1 && !v1.isconj())
            return v2.getX()*MultVV<N>(v2.getV().cptr(),v1.cptr());
        else 
            return v2.getX()*MultVV(v2.getV().view(),v1);
    }

    // (x*v) * v
    template <class T, class T1, int N, IndexStyle I> 
    inline T operator*(const ProdXV<T,T1>& v1, const SmallVector<T,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr());
        else 
            return v1.getX()*MultVV(v1.getV(),v2.view());
    }

    template <class T, class T1, int N, IndexStyle I> 
    inline CT operator*(const ProdXV<CT,T1>& v1, const SmallVector<T,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr());
        else 
            return v1.getX()*MultVV(v1.getV(),v2.view());
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const ProdXV<T,T>& v1, const SmallVector<CT,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*MultVV<N>(v2.cptr(),v1.getV().cptr());
        else 
            return v1.getX()*MultVV(v2.view(),v1.getV());
    }

    // (x*v) * (x*v)
    template <class T, class T1, class T2, int N, IndexStyle I> 
    inline T operator*(const ProdXV<T,T1>& v1, const ProdXv<T,T2,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v1.getV().cptr(),v2.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV().view());
    }

    template <class T, class T1, int N, IndexStyle I> 
    inline CT operator*(const ProdXV<CT,T1>& v1, const ProdXv<T,T,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v1.getV().cptr(),v2.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v1.getV(),v2.getV().view());
    }

    template <class T, class T2, int N, IndexStyle I> 
    inline CT operator*(const ProdXV<T,T>& v1, const ProdXv<CT,T2,N,I>& v2)
    {
        if (v1.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v2.getV().cptr(),v1.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v2.getV().view(),v1.getV());
    }

    // v * v 
    template <class T, int N, IndexStyle I> 
    inline T operator*(const SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
    {
        if (v2.step() == 1 && !v2.step())
            return MultVV<N>(v1.cptr(),v2.cptr());
        else 
            return MultVV(v1.view(),v2);
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const SmallVector<CT,N,I>& v1, const GenVector<T>& v2)
    {
        if (v2.step() == 1 && !v2.step())
            return MultVV<N>(v1.cptr(),v2.cptr());
        else 
            return MultVV(v1.view(),v2);
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const SmallVector<T,N,I>& v1, const GenVector<CT>& v2)
    { 
        if (v2.step() == 1 && !v2.step())
            return MultVV<N>(v2.cptr(),v1.cptr());
        else 
            return MultVV(v2,v1.view());
    }

    // v * (x*v)
    template <class T, class T2, int N, IndexStyle I> 
    inline T operator*(const SmallVector<T,N,I>& v1, const ProdXV<T,T2>& v2) 
    { 
        if (v2.getV().step() == 1)
            return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr());
        else 
            return v2.getX()*MultVV(v1.view(),v2.getV());
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const SmallVector<CT,N,I>& v1, const ProdXV<T,T>& v2)
    {
        if (v2.getV().step() == 1)
            return v2.getX()*MultVV<N>(v1.cptr(),v2.getV().cptr());
        else 
            return v2.getX()*MultVV(v1.view(),v2.getV());
    }

    template <class T, class T2, int N, IndexStyle I> 
    inline CT operator*(const SmallVector<T,N,I>& v1, const ProdXV<CT,T2>& v2)
    { 
        if (v2.getV().step() == 1)
            return v2.getX()*MultVV<N>(v2.getV().cptr(),v1.cptr());
        else 
            return v2.getX()*MultVV(v2.getV(),v1.view());
    }

    // (x*v) * v
    template <class T, class T1, int N, IndexStyle I> 
    inline T operator*(const ProdXv<T,T1,N,I>& v1, const GenVector<T>& v2)
    {
        if (v2.step() == 1 && !v2.step())
            return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr());
        else 
            return v1.getX()*MultVV(v1.getV().view(),v2);
    }

    template <class T, class T1, int N, IndexStyle I> 
    inline CT operator*(const ProdXv<CT,T1,N,I>& v1, const GenVector<T>& v2)
    {
        if (v2.step() == 1 && !v2.step())
            return v1.getX()*MultVV<N>(v1.getV().cptr(),v2.cptr());
        else 
            return v1.getX()*MultVV(v1.getV().view(),v2);
    }

    template <class T, int N, IndexStyle I> 
    inline CT operator*(const ProdXv<T,T,N,I>& v1, const GenVector<CT>& v2)
    {
        if (v2.step() == 1 && !v2.step())
            return v1.getX()*MultVV<N>(v2.cptr(),v1.getV().cptr());
        else 
            return v1.getX()*MultVV(v2,v1.getV().view());
    }

    // (x*v) * (x*v)

    template <class T, class T1, class T2, int N, IndexStyle I> 
    inline T operator*(const ProdXv<T,T1,N,I>& v1, const ProdXV<T,T2>& v2)
    {
        if (v2.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v1.getV().cptr(),v2.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v1.getV().view(),v2.getV());
    }

    template <class T, class T1, int N, IndexStyle I> 
    inline CT operator*(const ProdXv<CT,T1,N,I>& v1, const ProdXV<T,T>& v2)
    {
        if (v2.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v1.getV().cptr(),v2.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v1.getV().view(),v2.getV());
    }

    template <class T, class T2, int N, IndexStyle I> 
    inline CT operator*(const ProdXv<T,T,N,I>& v1, const ProdXV<CT,T2>& v2)
    {
        if (v2.getV().step() == 1)
            return v1.getX()*v2.getX()*
                MultVV<N>(v2.getV().cptr(),v1.getV().cptr());
        else 
            return v1.getX()*v2.getX()*MultVV(v2.getV(),v1.getV().view());
    }


} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif 
