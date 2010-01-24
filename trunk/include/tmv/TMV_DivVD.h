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


#ifndef TMV_DivVD_H
#define TMV_DivVD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_BaseVector.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class V1, class M2, class V3>
    inline void LDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void InlineLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <int ix, class T, class V1, class V2, class V3>
    inline void ElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void NoAliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void AliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);


    // Defined in TMV_DivVD.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstLDivVM(
        const T3 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstElemDivVV(
        const T3 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

    //
    // Element-wise division:
    // v3(i) = x * v1(i) / v2(i)
    //

    template <int algo, int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper;

    // algo 0: size == 0, nothing to do
    template <int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<0,0,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(const Scaling<ix,T>&, const V1&, const V2&, V3&) {}
        static void call2(int, const Scaling<ix,T>&, IT1, IT2, IT3) {}
    };

    // algo 11: simple for loop
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<1,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;
            if (n) do {
                *C++ = ZProd<false,false>::prod(
                    x, ZProd<c1,c2>::quot(*A++,*B++));
            } while (--n);
        }
    };  

    // TODO: Write SSE algo's.

    // algo -4: No branches or copies
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-4,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        enum { algo = (
                1 ) };
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::vconj);
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            TMVStaticAssert(!V3::vconj);
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call2(n,x,A,B,C);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-3,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { ElemDivVV_Helper<-4,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); }
    };

    // algo 97: Conjugate
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<97,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { 
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemDivVV_Helper<-2,size,ix,T,V1c,V2c,V3c>::call(x,v1c,v2c,v3c);
        }
    };

    // algo 98: Call inst
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<98,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { InstElemDivVV(v1.xView(),v2.xView(),v3.xView()); }
    };

    // algo -2: Check for inst
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-2,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::vconj ? 97 : 
                inst ? 98 : 
                -4;
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    // algo 99: Check for aliases
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<99,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber1 =
                !SameStorage(v1,v3) ||
                ExactSameStorage(v1,v3) || 
                v1.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v1.step());
            const bool noclobber2 =
                !SameStorage(v2,v3) ||
                ExactSameStorage(v2,v3) || 
                v2.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v2.step());
            if (noclobber1) {
                if (noclobber2) {
                    // No aliasing (or no clobering)
                    ElemDivVV_Helper<-2,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
                } else { 
                    // Need a temporary for v2
                    NoAliasElemDivVV(x,v1,v2.copy(),v3);
                }
            } else {
                if (noclobber2) {
                    // Need a temporary for v1
                    NoAliasElemDivVV(v1.copy(),v2,v3);
                } else {
                    // Need a temporary for v3
                    typename V3::copy_type v3c(v3.size());
                    NoAliasElemDivVV(x,v1,v2,v3c);
                    NoAliasCopy(v3c,v3);
                }
            }
        }
    };

    // algo -1: Check for aliases?
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-1,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber =
                VStepHelper<V1,V3>::noclobber &&
                VStepHelper<V2,V3>::noclobber;
            const bool checkalias =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
                !noclobber;
            const int algo =
                checkalias ? 99 : 
                -2;
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    template <int ix, class T, class V1, class V2, class V3>
    inline void ElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-1,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void NoAliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-2,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-4,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void AliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<99,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    // TODO: Check for singular DiagMatrix
    template <int ix, class T, class V1, class M2, class V3>
    inline void LDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-1,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-2,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void InlineLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-4,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasLDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<99,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    //
    // v1 /= m2
    //

    template <class V1, class M2>
    inline void LDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { LDivVM(v1,m2,v1); }
    template <class V1, class M2>
    inline void NoAliasLDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasLDivVM(v1,m2,v1); }
    template <class V1, class M2>
    inline void AliasLDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { AliasLDivVM(v1,m2,v1); }

    //
    // v3 = v1 % m2
    //

    template <int ix, class T, class V1, class M2, class V3>
    inline void RDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDivVM(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasRDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasLDivVM(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasRDivVM(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { AliasLDivVM(x,v1,m2,v3); }

    //
    // v1 %= m2
    //

    template <class V1, class M2>
    inline void RDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { LDivVM(v1,m2,v1); }
    template <class V1, class M2>
    inline void NoAliasRDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasLDivVM(v1,m2,v1); }
    template <class V1, class M2>
    inline void AliasRDivEqVM(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { AliasLDivVM(v1,m2,v1); }

} // namespace tmv

#endif 
