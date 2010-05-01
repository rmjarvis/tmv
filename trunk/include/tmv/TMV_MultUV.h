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


#ifndef TMV_MultUV_H
#define TMV_MultUV_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Vector.h"
#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_Prefetch.h"

#ifdef PRINTALGO_UV
#include <iostream>
#endif

namespace tmv {

    // Defined below:
    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void NoAliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void AliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3);

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, int ix, class T, class M2>
    inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class V1, int ix, class T, class M2>
    inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class V1, int ix, class T, class M2>
    inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // Defined in TMV_MultUV.cpp
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstUpperTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
    template <class T1, DiagType D, bool C1, class T2, bool C2, class T3>
    void InstAddMultMV(
        const T3 x,
        const ConstLowerTriMatrixView<T1,D,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

    //
    // Matrix * Vector
    //

    // Q1 is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_Q1 200 
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

    // Q2 is the minimum size to copy a vector if its step == UNKNOWN.
#define TMV_Q2 4

    // Q3 is the crossover memory size to start using prefetch commands.
    // This is undoubtedly a function of the L1 (and L2?) cache size,
    // but 2KBytes is probably not too bad for most machines.
    // (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_Q3 2048


    // Note: all algorithms here are designed to work if v2 and v3 are 
    // the same storage.  So it works for things like v *= U without 
    // any need for a temporary copy.
    template <int algo, int s, bool add, 
              int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper;

    // algo 0: s = 0, so nothing to do
    template <bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<0,0,add,ix,T,M1,V2,V3>
    {
        static void call(const Scaling<ix,T>& , const M1& , const V2& , V3& ) 
        {} 
    };

    // algo 1: s == 1, so simplifies to a scalar product
    template <bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<1,1,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            Maybe<add>::add( 
                v3.ref(0) , 
                ZProd<false,false>::prod(
                    x, Maybe<!M1::_unit>::prod(m1.cref(0,0) , v2.cref(0)))); 
        }
    };

    // algo 11: The basic column major loop for UpperTri
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<11,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int N = (s == UNKNOWN ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_diag_type M1d;
            PT2 Xj;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M1d::const_nonconj_type::const_iterator IT1d;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;
            IT1 A0j = m1.get_col(0,0,1).nonConj().begin();
            IT2 X = v2.nonConj().begin();
            IT3 Y0 = v3.begin();
            IT3 Y = Y0;
            const int Astepj = m1.stepj();

            const bool dopref = N * sizeof(T1) >= TMV_Q3;

            Prefetch_Read(m1.cptr());
            Prefetch_Read(v2.cptr());
            Prefetch_MultiWrite(v3.ptr());

            for(int j=0;N--;++j) {
                // loop from j = 0 .. N-1
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X++);
                    // y.subVector(0,j) += x(j) * A.col(j,0,j);
                    MultXV_Helper<-4,UNKNOWN,true,0,PT2,M1c,V3>::call2(
                        j,Scaling<0,PT2>(Xj),A0j,Y0);
                    A0j.shiftP(Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                    Maybe<add>::add(*Y++,Maybe<!unit>::prod(m1.cref(j,j),Xj));
                } else {
                    A0j.shiftP(Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                    ++X; 
                    Maybe<!add>::set(*Y++,T3(0));
                }
            }
        }
    };

    // algo 12: The basic row major loop for UpperTri
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<12,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        { 
            int N = (s == UNKNOWN ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename M1::const_row_sub_type M1r;
            PT Yi;

            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            // Actually Aii is the address of A(i,i+1)
            IT1 Aii = m1.get_row(0,1,N).nonConj().begin();
            IT2 X = v2.nonConj().begin();
            IT3 Y = v3.begin();
            const int Adiagstep = m1.diagstep();
            const int Astepi = m1.stepi();

            const bool dopref = N * sizeof(T1) >= TMV_Q3;

            Prefetch_Read(m1.cptr());
            Prefetch_MultiRead(v2.cptr());
            Prefetch_Write(v3.ptr());

            // [ A00 A01 A02 ] [  0 ]   [ A01 B2 ]
            // [  0  A11 A12 ] [ B2 ] = [ A11 B2 ]
            // [  0   0  A22 ] [  0 ]   [    0   ]
            while (N > 0 && v2(N-1) == T2(0)) Maybe<!add>::set(v3(--N),T3(0));
            int i1=0;
            while (N > 0 && *X == T2(0)) { ++i1; --N; ++X; ++Aii; }
            for(int i=0;i<i1;++i) {
                Yi = MultVV_Helper<-4,UNKNOWN,M1r,V2>::call2(N,Aii-1,X);
                Aii.shiftP(Astepi);
                if (dopref) Prefetch_Read(Aii.get()-1);
                Maybe<add>::add(*Y++, x * Yi);
            }

            for(int i=i1;N--;++i) {
                // loop from i = 0 .. N-1
                // Yi = A.row(i,i,N) * x.subVector(i,N)
                //    = A(i,i) * x(i)
                //      + A.row(i,i+1,N) * x.subVector(i+1,N)
                Yi = Maybe<!unit>::template zprod<false,c2>(m1.cref(i,i),*X++);
                Yi += MultVV_Helper<-4,UNKNOWN,M1r,V2>::call2(N,Aii,X);
                Aii.shiftP(Adiagstep);
                if (dopref) Prefetch_Read(Aii.get()-1);
                Maybe<add>::add(*Y++, x * Yi);
            }
        }
    };

    // algo 15: UpperTri: unroll by rows
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<15,s,add,ix,T,M1,V2,V3>
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                Unroller<I,N/2>::unroll(x,m1,v2,v3);
                Unroller<I+N/2,N-N/2>::unroll(x,m1,v2,v3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                typedef typename M1::const_row_sub_type M1r;
                typedef typename V2::const_subvector_type V2s;
                const bool unit = M1::_unit;
                Maybe<add>::add(
                    v3.ref(I) , x * ( 
                        Maybe<!unit>::prod(m1.cref(I,I) , v2.cref(I)) +
                        MultVV_Helper<-4,s-I-1,M1r,V2s>::call(
                            m1.get_row(I,I+1,s),v2.cSubVector(I+1,s))));
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix,T>& , const M1& , const V2& , V3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 15: N,s,x = "<<s<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,v2,v3); 
        }
    };

    // algo 21: The basic column major loop for LowerTri
    // This is designed to work if v2 and v3 are the same storage.
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<21,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = (s == UNKNOWN ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_diag_type M1d;
            PT2 Xj;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M1d::const_nonconj_type::const_iterator IT1d;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            // Actually Ajj is the address of A(j,j+1)
            IT1 Ajj = m1.get_col(N-1,N,N).nonConj().begin();
            IT2 X = v2.nonConj().begin() + N-1;
            IT3 Y = v3.begin() + N;
            const int Adiagstep = m1.diagstep();

            const bool dopref = N * Adiagstep * sizeof(T1) >= TMV_Q3;

            Prefetch_Read(v2.cptr());
            Prefetch_MultiWrite(v3.ptr());
            if (dopref) Prefetch_Read(Ajj.get());
            else Prefetch_Read(m1.cptr());

            for(int j=N,len=0;j--;++len) {
                // loop from j = N-1 .. 0
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X--);
                    // y.subVector(j+1,N) += x(j) * A.col(j,j+1,N);
                    MultXV_Helper<-4,UNKNOWN,true,0,PT2,M1c,V3>::call2(
                        len,Scaling<0,PT2>(Xj),Ajj,Y--);
                    Ajj.shiftP(-Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                    // y(j) (+)= x(j) * A(j,j);
                    Maybe<add>::add(*Y , Maybe<!unit>::prod(m1.cref(j,j),Xj));
                } else {
                    Ajj.shiftP(-Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                    --X;
                    Maybe<!add>::set(*--Y,T3(0));
                }
            }
        }
    };

    // algo 22: The basic row major loop for LowerTri
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<22,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        { 
            int N = (s == UNKNOWN ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename M1::const_row_sub_type M1r;
            PT Yi;

            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            IT1 Ai0 = m1.get_row(N-1,0,N-1).nonConj().begin();
            IT2 X0 = v2.nonConj().begin();
            IT2 X = X0 + N-1;
            IT3 Y = v3.begin() + N-1;
            const int Astepi = m1.stepi();

            const bool dopref = N * m1.diagstep() * sizeof(T1) >= TMV_Q3;

            Prefetch_MultiRead(v2.cptr());
            Prefetch_Write(v3.ptr());
            if (dopref) Prefetch_Read(Ai0.get());
            else Prefetch_Read(m1.cptr());

            // [ A00  0   0  ] [  0 ]   [    0   ]
            // [ A10 A11  0  ] [ B2 ] = [ A11 B2 ]
            // [ A20 A21 A22 ] [  0 ]   [ A21 B2 ]
            int i=N-1;
            int i1=0;
            while (N > 0 && *X0 == T2(0)) { 
                Maybe<!add>::set(v3(i1++),T3(0));
                --N; ++Ai0; ++X0;
            }
            int i2=i;
            while (N > 0 && *X == T2(0)) { --i2; --N; --X; }
            for(;i>i2;--i) {
                Yi = MultVV_Helper<-4,UNKNOWN,M1r,V2>::call2(N,Ai0,X0);
                Ai0.shiftP(-Astepi);
                if (dopref) Prefetch_Read(Ai0.get());
                Maybe<add>::add(*Y--, x * Yi);
            }

            for(;N--;--i) {
                // loop from i = N-1 .. 0
                // Yi = A.row(i,0,i+1) * x.subVector(0,i+1)
                //    = A(i,i) * x(i)
                //      + A.row(i,0,i) * x.subVector(0,i)
                Yi = MultVV_Helper<-4,UNKNOWN,M1r,V2>::call2(N,Ai0,X0);
                Ai0.shiftP(-Astepi);
                if (dopref) Prefetch_Read(Ai0.get());
                Yi += Maybe<!unit>::template zprod<false,c2>(m1.cref(i,i),*X--);
                Maybe<add>::add(*Y--, x * Yi);
            } 
        }
    };

    // algo 25: LowerTri: unroll by rows
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<25,s,add,ix,T,M1,V2,V3>
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                Unroller<I+N/2,N-N/2>::unroll(x,m1,v2,v3);
                Unroller<I,N/2>::unroll(x,m1,v2,v3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                typedef typename M1::const_row_sub_type M1r;
                typedef typename V2::const_subvector_type V2s;
                const bool unit = M1::_unit;
                Maybe<add>::add(
                    v3.ref(I) , x * ( 
                        Maybe<!unit>::prod(m1.cref(I,I) , v2.cref(I)) +
                        MultVV_Helper<-4,I,M1r,V2s>::call(
                            m1.get_row(I,0,I),v2.cSubVector(0,I))));
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix,T>& , const M1& , const V2& , V3& ) {}
        };

        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 25: N,s,x = "<<s<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,v2,v3); 
        }
    };

    // algo 43: colmajor, v3.step == UNKNOWN, so maybe copy v3
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<43,s,add,ix,T,M1,V2,V3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = s == UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 43: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
#ifdef TMV_OPT_SCALE
            if (N > TMV_Q2) {
#endif
                MultUV_Helper<85,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef TMV_OPT_SCALE
            } else 
                MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
        }
    };

    // algo 53: rowmajor, v2.step == UNKNOWN, so maybe copy v2
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<53,s,add,ix,T,M1,V2,V3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = s == UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 53: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
#ifdef TMV_OPT_SCALE
            if (N > TMV_Q2) {
#endif
                MultUV_Helper<81,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef TMV_OPT_SCALE
            } else 
                MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
        }
    };

    // algo 81: copy v2
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<81,s,add,ix,T,M1,V2,V3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = s == UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 81: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename V2::value_type T2;
            typedef typename VCopyHelper<T2,s,false>::type V2c;
            V2c v2c(N);
            typedef typename V2c::view_type V2cv;
            typedef typename V2c::const_view_type V2ccv;
            V2cv v2cv = v2c.view();
            V2ccv v2ccv = v2c.view();
            CopyV_Helper<-2,s,V2,V2cv>::call(v2,v2cv);
            MultUV_Helper<-2,s,add,ix,T,M1,V2ccv,V3>::call(x,m1,v2ccv,v3);
        }
    };

    // algo 85: v3c = x*m1*v2, v3 (+)= v3c
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<85,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = s == UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 85: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename V3::value_type T3;
            typedef typename VCopyHelper<T3,s,false>::type V3c;
            V3c v3c(N);
            typedef typename V3c::view_type V3cv;
            typedef typename V3c::const_view_type V3ccv;
            V3cv v3cv = v3c.view();
            V3ccv v3ccv = v3c.view();
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultUV_Helper<-2,s,false,ix,T,M1,V2,V3cv>::call(x,m1,v2,v3cv);
            MultXV_Helper<-2,s,add,1,RT,V3ccv,V3>::call(one,v3ccv,v3);
        }
    };

    // algo -4: No branches or copies
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                M1::_upper ? (
                    unroll ? 15 :
                    M1::_colmajor ? 11 :
                    M1::_rowmajor ? 12 :
                    V2::_step == 1 ? 12 : V3::_step == 1 ? 11 : 12 ) :
                ( // lowertri
                    unroll ? 25 :
                    M1::_colmajor ? 21 :
                    M1::_rowmajor ? 22 :
                    V2::_step == 1 ? 22 : V3::_step == 1 ? 21 : 22 );
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-3,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product
            //
            // UpperTri:
            // 11 = column major, simple for loop
            // 12 = row major, simple for loop
            // 15 = fully unroll by rows
            //
            // LowerTri:
            // 21 = column major, simple for loop
            // 22 = row major, simple for loop
            // 25 = fully unroll by rows
            //
            // Copy a vector to new storage:
            // 81 = copy v2
            // 85 = temp v3 = x*m1*v2

#if TMV_OPT == 0
            const int algo = 
                M1::_upper ?
                ( M1::_colmajor ? 11 : 12 ) :
                ( M1::_colmajor ? 21 : 22 );
#else
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                ( s == 0 ) ? 0 : // trivial - nothing to do
                ( s == 1 ) ? 1 : // trivial - s = 1
                M1::_upper ? (
                    unroll ? 15 :
                    M1::_colmajor ? (
                        V3::_step == UNKNOWN ? (
                            s == UNKNOWN ? 43 :
                            s > TMV_Q2 ? 85 : 11 ) :
                        11 ) :
                    M1::_rowmajor ? (
                        V2::_step == UNKNOWN ? (
                            s == UNKNOWN ? 53 :
                            s > TMV_Q2 ? 81 : 12 ) :
                        12 ) :
                    V2::_step == 1 ? 12 : V3::_step == 1 ? 11 : 12 ) :
                ( // lowertri
                    unroll ? 25 :
                    M1::_colmajor ? (
                        V3::_step == UNKNOWN ? (
                            s == UNKNOWN ? 43 :
                            s > TMV_Q2 ? 85 : 21 ) :
                        21 ) :
                    M1::_rowmajor ? (
                        V2::_step == UNKNOWN ? (
                            s == UNKNOWN ? 53 :
                            s > TMV_Q2 ? 81 : 22 ) :
                        22 ) :
                    V2::_step == 1 ? 22 : V3::_step == 1 ? 21 : 22 );
#endif
#ifdef PRINTALGO_UV
            std::cout<<"InlineMultUV: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, 
              int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<97,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultUV_Helper<-2,s,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 98: call InstMultMV
    template <int s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<98,s,false,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typename V3::value_type xx(x);
            InstMultMV(xx,m1.xdView(),v2.xView(),v3.xView());
        }
    };
    template <int s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<98,s,true,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typename V3::value_type xx(x);
            InstAddMultMV(xx,m1.xdView(),v2.xView(),v3.xView());
        }
    };

    // algo -2: Check for inst
    template <int s, bool add, 
              int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-2,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst = 
                M1::unknownsizes &&
                V2::unknownsizes &&
                V3::unknownsizes &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                V3::_conj ? 97 :
                inst ? 98 : 
                -3;
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 99: Check for aliases
    template <int s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<99,s,true,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 ( (M1::_upper && v2.step()>=v3.step()) ||
                   (M1::_lower && v2.step()<=v3.step()) ||
                   !SameStorage(v2.vec(),v3.vec())) ) {
                // No aliasing (or no clobbering)
                MultUV_Helper<-2,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultUV_Helper<85,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else { 
                // SameStorage(v2,v3)
                // Use temporary for v2
                MultUV_Helper<81,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // if !add, then don't need the temporary if v2,v3 are aliased.
    template <int s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<99,s,false,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 ( (M1::_upper && v2.step()>=v3.step()) ||
                   (M1::_lower && v2.step()<=v3.step()) ||
                   !SameStorage(v2.vec(),v3.vec())) ) {
                // No aliasing (or no clobbering)
                MultUV_Helper<-2,s,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultUV_Helper<85,s,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else { 
                AliasCopy(v2,v3);
                NoAliasMultMV<false>(x,m1,v3,v3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, 
              int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-1,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const bool checkalias =
                M1::_size == UNKNOWN &&
                V2::_size == UNKNOWN && 
                V3::_size == UNKNOWN;
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                checkalias ? 99 : 
                -2;
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    template <int algo, bool add, int ix, class T, class M1, class V2, class V3>
    inline void DoMultUV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V2::_size,V3::_size>::same));
        TMVAssert(m1.size() == v3.size());
        TMVAssert(m1.size() == v2.size());
        TMVAssert(v2.size() == v3.size());
        const int s = Sizes<Sizes<M1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        M1v m1v = m1.cView();
        V2v v2v = v2.cView();
        V3v v3v = v3.cView();
        MultUV_Helper<algo,s,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3)
    { DoMultUV<-1,add>(x,m1,v2,v3); }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void NoAliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3)
    { DoMultUV<-2,add>(x,m1,v2,v3); }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3)
    { DoMultUV<-3,add>(x,m1,v2,v3); }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void AliasMultMV(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseVector_Calc<V2>& v2, 
        BaseVector_Mutable<V3>& v3)
    { DoMultUV<99,add>(x,m1,v2,v3); }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,V3::_size>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(m2.rowsize() == v3.size());
        MultMV<add>(x,m2.transpose(),v1.vec(),v3.vec());
    }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,V3::_size>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(m2.rowsize() == v3.size());
        NoAliasMultMV<add>(x,m2.transpose(),v1.vec(),v3.vec());
    }

    template <bool add, int ix, class T, class V1, class M2, class V3>
    inline void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,V3::_size>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(m2.rowsize() == v3.size());
        AliasMultMV<add>(x,m2.transpose(),v1.vec(),v3.vec());
    }

    template <class V1, int ix, class T, class M2>
    inline void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { MultVM<false>(x,v1.vec(),m2.mat(),v1.vec()); }

    template <class V1, int ix, class T, class M2>
    inline void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { NoAliasMultVM<false>(x,v1.vec(),m2.mat(),v1.vec()); }

    template <class V1, int ix, class T, class M2>
    inline void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { AliasMultVM<false>(x,v1.vec(),m2.mat(),v1.vec()); }

} // namespace tmv

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3

#endif 
