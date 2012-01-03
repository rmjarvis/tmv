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


#ifndef TMV_SmallVectorArithFunc_H
#define TMV_SmallVectorArithFunc_H

#define CT std::complex<T>

namespace tmv {

    template <class T1, class T2>
    struct OKTypeHelper { enum { ok=true }; };
    template <class T>
    struct OKTypeHelper<CT,T> { enum { ok=false }; };

    template <class T1, class T2, class T3>
    struct OKTypeHelper2 { enum { ok=true }; };
    template <class T>
    struct OKTypeHelper2<T,CT,T> { enum { ok=false }; };
    template <class T>
    struct OKTypeHelper2<CT,T,T> { enum { ok=false }; };
    template <class T>
    struct OKTypeHelper2<CT,CT,T> { enum { ok=false }; };

#define OKTypes(T1,T2) (OKTypeHelper<T1,T2>::ok)
#define OKTypes2(T1,T2,T3) (OKTypeHelper2<T1,T2,T3>::ok)


    // v *= x
    template <int N, class T, class Tv, bool ok> 
    struct MultXVHelper 
    { MultXVHelper(const T, Tv*) {} };
    template <int N, class T, class Tv>
    struct MultXVHelper<N,T,Tv,true>
    {
        MultXVHelper(const T x, Tv* v)
        { for(int i=0;i<N;++i) v[i] *= x; }
    };
    template <int N, class T, class Tv> 
    inline void MultXV(const T x, Tv* v)
    { MultXVHelper<N,T,Tv,OKTypes(T,Tv)>(x,v); }

    // v2 = x * v1
    template <int N, class T, class T1, class T2, bool ok>
    struct MultXVHelper2
    { MultXVHelper2(const T, const T1*, T2*) {} };
    template <int N, class T, class T1, class T2>
    struct MultXVHelper2<N,T,T1,T2,true>
    {
        MultXVHelper2(const T x, const T1* v1, T2* v2)
        { for(int i=0;i<N;++i) v2[i] = x * v1[i]; }
    };
    template <int N, class T, class T1, class T2> 
    inline void MultXV( const T x, const T1* v1, T2* v2)
    { MultXVHelper2<N,T,T1,T2,OKTypes2(T,T1,T2)>(x,v1,v2); }

    // v2 += x * v1
    template <int N, class T1, class T2, bool ok>
    struct AddVV_1Helper
    { AddVV_1Helper(const T1*, T2*) {} };
    template <int N, class T1, class T2>
    struct AddVV_1Helper<N,T1,T2,true>
    {
        inline AddVV_1Helper(const T1* v1, T2* v2)
        { for(int i=0;i<N;++i) v2[i] += v1[i]; }
    };
    template <int N, class T1, class T2> 
    inline void AddVV_1(const T1* v1, T2* v2)
    { AddVV_1Helper<N,T1,T2,OKTypes(T1,T2)>(v1,v2); }

    template <int N, class T1, class T2, bool ok>
    struct AddVV_m1Helper
    { AddVV_m1Helper(const T1*, T2*) {} };
    template <int N, class T1, class T2>
    struct AddVV_m1Helper<N,T1,T2,true>
    {
        inline AddVV_m1Helper(const T1* v1, T2* v2)
        { for(int i=0;i<N;++i) v2[i] -= v1[i]; }
    };
    template <int N, class T1, class T2> 
    inline void AddVV_m1(const T1* v1, T2* v2)
    { AddVV_m1Helper<N,T1,T2,OKTypes(T1,T2)>(v1,v2); }

    template <int N, class T, class T1, class T2, bool ok>
    struct AddVVHelper
    { AddVVHelper(const T, const T1*, T2*) {} };
    template <int N, class T, class T1, class T2>
    struct AddVVHelper<N,T,T1,T2,true>
    {
        inline AddVVHelper(const T x, const T1* v1, T2* v2) 
        { for(int i=0;i<N;++i) v2[i] += x*v1[i]; }
    };
    template <int N, class T, class T1, class T2> 
    inline void AddVV(const T x, const T1* v1, T2* v2)
    { AddVVHelper<N,T,T1,T2,OKTypes2(T,T1,T2)>(x,v1,v2); }

    // v3 = x1 * v1 + x2 * v2
    template <int N, class T1, class T2, class T3> 
    inline void AddVV_1_1(const T1* v1, const T2* v2, T3* v3)
    { DoCopy<N>(v1,v3); AddVV_1<N>(v2,v3); }
    template <int N, class T1, class T2, class T3> 
    inline void AddVV_1_m1(const T1* v1, const T2* v2, T3* v3)
    { DoCopy<N>(v1,v3); AddVV_m1<N>(v2,v3); }
    template <int N, class T, class T1, class T2, class T3>
    inline void AddVV_1_x(const T1* v1, const T x2, const T2* v2, T3* v3)
    { MultXV<N>(x2,v2,v3); AddVV_1<N>(v1,v3); }
    template <int N, class T, class T1, class T2, class T3>
    inline void AddVV_x_1(const T x1, const T1* v1, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV_1<N>(v2,v3); }
    template <int N, class T, class T1, class T2, class T3>
    inline void AddVV_x_m1(const T x1, const T1* v1, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV_m1<N>(v2,v3); }
    template <int N, class T, class T1, class T2, class T3> 
    inline void AddVV(const T x1, const T1* v1, const T x2, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV<N>(x2,v2,v3); }

    // v1 * v2 (dot product)
    template <class T1, class T2> 
    struct ProdType { typedef T1 Tprod; };
    template <class T> 
    struct ProdType<T,std::complex<T> > { typedef std::complex<T> Tprod; };
#define ProductType(T1,T2) typename ProdType<T1,T2>::Tprod

    template <int N, class T, class T1, class T2>
    struct MultVVHelper
    {
        inline MultVVHelper(T& res, const T1* v1, const T2* v2) 
        { for(int i=0;i<N;++i) res += v1[i]*v2[i]; }
    };
    template <int N, class T1, class T2> 
    inline ProductType(T1,T2) MultVV(const T1* v1, const T2* v2)
    { 
        ProductType(T1,T2) res(0);
        MultVVHelper<N,ProductType(T1,T2),T1,T2>(res,v1,v2);
        return res;
    }
#undef ProductType

    template <class T, class T1, class T2, class T3, int N, IndexStyle I1, IndexStyle I2, IndexStyle I3>
    inline void ElemMultVV(
        const T alpha, const SmallVector<T1,N,I1>& v1, 
        const SmallVector<T2,N,I2>& v2, SmallVector<T3,N,I3>& v3)
    {
        if (alpha == T(1))
            for(int i=0;i<N;++i) v3.ref(i) = v1.cref(i) * v2.cref(i);
        else
            for(int i=0;i<N;++i) v3.ref(i) = alpha * v1.cref(i) * v2.cref(i);
    }
    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2, IndexStyle I3> 
    inline void ElemMultVV(
        const CT , const SmallVector<T1,N,I1>& ,
        const SmallVector<T2,N,I2>& , SmallVector<T,N,I3>& )
    { TMVAssert(TMV_FALSE); }


    template <class T, class T1, class T2, class T3, int N, IndexStyle I1, IndexStyle I2, IndexStyle I3>
    inline void AddElemMultVV(
        const T alpha, const SmallVector<T1,N,I1>& v1, 
        const SmallVector<T2,N,I2>& v2, SmallVector<T3,N,I3>& v3)
    {
        if (alpha == T(1))
            for(int i=0;i<N;++i) v3.ref(i) += v1.cref(i) * v2.cref(i);
        else if (alpha == T(-1))
            for(int i=0;i<N;++i) v3.ref(i) -= v1.cref(i) * v2.cref(i);
        else
            for(int i=0;i<N;++i) v3.ref(i) += alpha * v1.cref(i) * v2.cref(i);
    }
    template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2, IndexStyle I3> 
    inline void AddElemMultVV(
        const CT , const SmallVector<T1,N,I1>& ,
        const SmallVector<T2,N,I2>& , SmallVector<T,N,I3>& )
    { TMVAssert(TMV_FALSE); }

#undef OKTypes
#undef OKTypes2

} // namespace tmv

#undef CT

#endif 
