

#ifndef TMV_MultBV_H
#define TMV_MultBV_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_BaseVector.h"
#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_ElemMultVV.h"
#include "TMV_MultMV_Funcs.h"

#ifdef PRINTALGO_BV
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#endif

#ifdef XDEBUG_BV
#include <iostream>
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_CopyM.h"
#endif

//
// BandMatrix * Vector
//

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_MV_ZeroIX (ix==0)
//#define TMV_MV_ZeroIX (ix!=1)

namespace tmv {

    // Defined in TMV_MultBV.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper;

    // algo 0: cs or rs = 0, so nothing to do
    // Correction: if rs = 0, cs != 0 and !add, then we need to do v3.setZero().
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<0,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const V2& , V3& v3) 
        {
#ifdef PRINTALGO_BV
            const int M = cs == TMV_UNKNOWN ? int(v3.size()) : cs;
            std::cout<<"BV algo 0: M,N,cs,rs,x = "<<M<<','<<0<<
                ','<<cs<<','<<rs<<','<<T(0)<<std::endl;
#endif
            Maybe<(!add && rs == 0 && cs != 0)>::zero(v3); 
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<11,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_BV
            std::cout<<"BV algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;

            typedef typename V2::const_subvector_type V2s;
            typedef typename V2s::const_nonconj_type::const_iterator IT2;
            const bool c2 = V2::_conj;

            typedef typename V2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            PT2 Xj;

            typedef typename V3::subvector_type V3s;
            typedef typename V3s::iterator IT3;

            const int Astepj = m1.stepj();
            const int Adiagstep = m1.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m1.nhi();
            const int j2 = TMV_MIN(N,M-m1.nlo());
            const int j3 = TMV_MIN(N,M+m1.nhi());
            int len = m1.nlo()+1;

            IT1 A0j = m1.get_col(0,0,len).begin().nonConj();
            IT2 X = v2.begin().nonConj();
            IT3 Y0 = v3.begin();

            Maybe<!add>::zero(v3);
            int j=0;
            for(;j<j1;++j) {
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X++);
                    MultXV_Helper<-4,xx,true,0,PT2,M1c,V3s>::call2(
                        len,Scaling<0,PT2>(Xj),A0j,Y0);
                } else {
                    ++X;
                }
                A0j.shiftP(Astepj);
                if (len < M) ++len;
            }
            if (j1 < j2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(;j<j2;++j) {
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X++);
                    MultXV_Helper<-4,lh,true,0,PT2,M1c,V3s>::call2(
                        len,Scaling<0,PT2>(Xj),A0j,Y0);
                } else {
                    ++X;
                }
                A0j.shiftP(Adiagstep);
                ++Y0;
            }
            if (j1 >= j2) ++len;
            for(;j<j3;++j) {
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X++);
                    MultXV_Helper<-4,xx,true,0,PT2,M1c,V3s>::call2(
                        --len,Scaling<0,PT2>(Xj),A0j,Y0);
                } else {
                    ++X; --len;
                }
                A0j.shiftP(Adiagstep);
                ++Y0;
            }
        }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<21,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_BV
            std::cout<<"BV algo 21: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;

            typedef typename V2::const_subvector_type V2s;
            typedef typename V2s::const_nonconj_type::const_iterator IT2;

            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT;
            PT Yi;

            typedef typename V3::subvector_type V3s;
            typedef typename V3s::iterator IT3;

            const int Astepi = m1.stepi();
            const int Adiagstep = m1.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m1.nlo();
            const int i2 = TMV_MIN(M,N-m1.nhi());
            const int i3 = TMV_MIN(M,N+m1.nlo());
            int len = m1.nhi()+1;

            IT1 Ai0 = m1.get_row(0,0,len).begin().nonConj();
            IT2 X0 = v2.begin().nonConj();
            IT3 Y = v3.begin();

            int i=0;
            for(;i<i1;++i) {
                Yi = MultVV_Helper<-4,xx,M1r,V2s>::call2(len,Ai0,X0);
                Maybe<add>::add(*Y++, ZProd<false,false>::prod(x,Yi));
                Ai0.shiftP(Astepi);
                if (len < N) ++len;
            }
            if (i1 < i2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(;i<i2;++i) {
                Yi = MultVV_Helper<-4,lh,M1r,V2s>::call2(len,Ai0,X0);
                Maybe<add>::add(*Y++, ZProd<false,false>::prod(x,Yi));
                Ai0.shiftP(Adiagstep);
                ++X0;
            }
            if (i1 >= i2) ++len;
            for(;i<i3;++i) {
                Yi = MultVV_Helper<-4,xx,M1r,V2s>::call2(--len,Ai0,X0);
                Maybe<add>::add(*Y++, ZProd<false,false>::prod(x,Yi));
                Ai0.shiftP(Adiagstep);
                ++X0;
            }
            if (i3 < M && !add) for(int i=i3;i<M;++i) *Y++ = T3(0);
        }
    };

    // algo 31: Loop over diagonals
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<31,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_BV
            std::cout<<"BV algo 31: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_diag_sub_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;

            typedef typename V2::const_subvector_type V2s;
            typedef typename V2s::const_nonconj_type::const_iterator IT2;

            typedef typename V3::subvector_type V3s;
            typedef typename V3s::iterator IT3;

            const int Astepj = m1.stepj();
            const int Astepi = m1.stepi();
            int len = TMV_MIN(M-m1.nlo(),N);

            IT1 A0 = m1.get_diag(-m1.nlo()).begin().nonConj();
            IT2 X0 = v2.begin().nonConj();
            IT3 Y0 = v3.begin(); Y0 += m1.nlo();

            Maybe<!add>::zero(v3);
            for(int k=m1.nlo();k;--k) {
                ElemMultVV_Helper<-4,xx,true,ix,T,M1c,V2s,V3s>::call2(
                    len,x,A0,X0,Y0);
                A0.shiftP(-Astepi);
                --Y0;
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const int ds = IntTraits2<cs,rs>::min;
            ElemMultVV_Helper<-4,ds,true,ix,T,M1c,V2s,V3s>::call2(
                len,x,A0,X0,Y0);
            for(int k=1;k<=m1.nhi();++k) {
                A0.shiftP(Astepj);
                ++X0;
                if (k+len > N) --len;
                ElemMultVV_Helper<-4,xx,true,ix,T,M1c,V2s,V3s>::call2(
                    len,x,A0,X0,Y0);
            }
        }
    };

    // algo 51: ix == 0, !add, so might want to use algo 53
    // to do the scaling at the end.
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<51,cs,rs,false,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_BV
            std::cout<<"BV algo 51: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M1::_colmajor || M1::_diagmajor);
            if (M >= N)
                MultBV_Helper<-4,cs,rs,false,ix,T,M1,V2,V3>::call(
                    x,m1,v2,v3);
            else 
                MultBV_Helper<53,cs,rs,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 53: column major, !add, apply x at the end
    template <int cs, int rs, class T, class M1, class V2, class V3>
    struct MultBV_Helper<53,cs,rs,false,0,T,M1,V2,V3>
    {
        static void call(
            const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_BV
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"BV algo 53: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultBV_Helper<-4,cs,rs,false,1,RT,M1,V2,V3>::call(one,m1,v2,v3);
            ScaleV_Helper<-3,cs,0,T,V3>::call(x,v3);
        }
    };

    // algo 81: copy v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<81,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_BV
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"BV algo 81: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typename V3::noalias_type v3na = v3.noAlias();
            MultMV<add>(x,m1,v2.copy(),v3na);
        }
    };

    // algo 82: copy x*v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<82,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_BV
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"BV algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename V3::noalias_type v3na = v3.noAlias();
            MultMV<add>(one,m1,(x*v2).calc(),v3na);
        }
    };

    // algo 83: Copy v2, figure out where to put x
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<83,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            if (N >= M) {
                MultBV_Helper<81,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                MultBV_Helper<82,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 81
    template <int cs, int rs, bool add, class T, class M1, class V2, class V3>
    struct MultBV_Helper<83,cs,rs,add,1,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<1,T>& x, const M1& m1, const V2& v2, V3& v3)
        { MultBV_Helper<81,cs,rs,add,1,T,M1,V2,V3>::call(x,m1,v2,v3); }
    };

    // algo 84: v3c = m1*v2, v3 (+)= x*v3c
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<84,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_BV
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"BV algo 84: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typename V3::noalias_type v3na = v3.noAlias();
            MultXV<add>(x,(m1*v2).calc(),v3na);
        }
    };

    // algo 85: v3c = x*m1*v2, v3 (+)= v3c
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<85,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_BV
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"BV algo 85: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename V3::noalias_type v3na = v3.noAlias();
            MultXV<add>(one,(x*m1*v2).calc(),v3na);
        }
    };

    // algo 86: Use temporary for v3, figure out where to put x
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<86,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            if (N >= M) {
                MultBV_Helper<84,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                MultBV_Helper<85,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 84
    template <int cs, int rs, bool add, class T, class M1, class V2, class V3>
    struct MultBV_Helper<86,cs,rs,add,1,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<1,T>& x, const M1& m1, const V2& v2, V3& v3)
        { MultBV_Helper<84,cs,rs,add,1,T,M1,V2,V3>::call(x,m1,v2,v3); }
    };

    // algo 90: call inst
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<90,cs,rs,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<90,cs,rs,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 91: call inst alias
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<91,cs,rs,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<91,cs,rs,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<97,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultBV_Helper<-2,cs,rs,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<197,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultBV_Helper<99,cs,rs,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<98,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 !SameStorage(v2,v3) ) {
                // No aliasing
                MultBV_Helper<-2,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultBV_Helper<86,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                // SameStorage(v2,v3)
                // Use temporary for v2
                MultBV_Helper<83,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<99,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                V3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultBV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<-4,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ? (
                    M1::_colmajor ? 11 : M1::_rowmajor ? 21 : 31 ) :
                M1::_colmajor ? (
                    ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 11 :
                    ( TMV_MV_ZeroIX && !add && rs > cs ) ? 53 : 
                    11 ) :
                M1::_rowmajor ? 21 :
                M1::_diagmajor ? (
                    ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 31 :
                    ( TMV_MV_ZeroIX && !add && rs > cs ) ? 53 : 
                    31 ) :
                V2::_step == 1 ? 21 : V3::_step == 1 ? 11 : 31;
#ifdef PRINTALGO_BV
            const int M = cs == TMV_UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"BV algo -4\n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
            std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"v2 = "<<v2<<std::endl;
            std::cout<<"v3 = "<<v3<<std::endl;
#endif
            MultBV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef PRINTALGO_BV
            std::cout<<"v3 => "<<v3<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<-3,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = cs or rs == 0, so nothing to do
            //
            // Regular algorithm:
            // 11 = column major
            // 21 = row major
            // 31 = diag major
            //
            // Meta algorithms
            // 51 = ix==0 && !add, so might want algo 53
            // 53 = !add, apply x at the end
            //
            // Copy a vector to new storage:
            // 81 = copy v2
            // 82 = copy x*v2
            // 83 = copy v2, figure out where to put x
            // 84 = temp v3 = m1*v2
            // 85 = temp v3 = x*m1*v2
            // 86 = temp v3, figure out where to put x
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
                TMV_OPT == 0 ? (
                    M1::_colmajor ? 11 : M1::_rowmajor ? 21 : 31 ) :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ? (
                    M1::_colmajor ? 11 : M1::_rowmajor ? 21 : 31 ) :
                M1::_colmajor ? (
                    ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? (
                        V3::_step == TMV_UNKNOWN ? 86 : 
                        ( TMV_MV_ZeroIX && !add ) ? 51 : 
                        11 ) :
                    V3::_step == TMV_UNKNOWN ? ( rs > cs ? 84 : 85 ) :
                    ( TMV_MV_ZeroIX && !add && rs > cs ) ? 53 : 
                    11 ) :
                M1::_rowmajor ? (
                    ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? (
                        V2::_step == TMV_UNKNOWN ? 83 : 21 ) :
                    V2::_step == TMV_UNKNOWN ? ( rs > cs ? 81 : 82 ) :
                    21 ) :
                M1::_diagmajor ? (
                    ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? (
                        V2::_step == TMV_UNKNOWN ? 83 : 
                        V3::_step == TMV_UNKNOWN ? 86 : 
                        ( TMV_MV_ZeroIX && !add ) ? 51 : 
                        31 ) :
                    V2::_step == TMV_UNKNOWN ? ( rs > cs ? 81 : 82 ) : 
                    V3::_step == TMV_UNKNOWN ? ( rs > cs ? 84 : 85 ) :
                    ( TMV_MV_ZeroIX && !add && rs > cs ) ? 53 : 
                    31 ) :
                V2::_step == 1 ? 21 : V3::_step == 1 ? 11 : 31;
#ifdef PRINTALGO_BV
            std::cout<<"InlineMultMV: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
            std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
#ifdef XDEBUG_BV
            typedef typename V3::real_type RT;
            typedef typename V3::value_type T3;
            Matrix<T3> m1c = m1;
            Vector<T3> v2c = v2;
            Vector<T3> v3i = v3;
            Vector<T3> v3c = v3;
            MultMV<add>(x,m1c,v2c,v3c);
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"v2 = "<<v2<<std::endl;
            if (add) std::cout<<"v3 = "<<v3<<std::endl;
            std::cout<<"v3c = "<<v3c<<std::endl;
#endif
            MultBV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef XDEBUG_BV
            std::cout<<"v3 => "<<v3<<std::endl;
            if (Norm(v3-v3c) > 1.e-3*(Norm(m1c)*Norm(v2c)+(add?Norm(v3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"v2 = "<<v2c<<std::endl;
                std::cout<<"v3 = "<<v3i<<std::endl;
                std::cout<<"v3 => "<<v3<<std::endl;
                std::cout<<"Correct v3 = "<<v3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<-2,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                V3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultBV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultBV_Helper<-1,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                V3::_checkalias ? 99 : 
                -2;
            MultBV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultBV_Helper<-1,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultBV_Helper<-3,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineAliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultBV_Helper<98,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

} // namespace tmv

#undef TMV_MV_SCALE
#undef TMV_MV_ZeroIX

#endif 
