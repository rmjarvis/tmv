

#ifndef TMV_MultBD_H
#define TMV_MultBD_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_ElemMultVV.h"
#include "TMV_MultXV.h"
#include "TMV_SmallBandMatrix.h"

//#define PRINTALGO_BD

#ifdef PRINTALGO_BD
#include <iostream>
#endif

namespace tmv {

    // Defined in TMV_MultBD.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, BandMatrixView<T3> m3);


    //
    // BandMatrix * DiagMatrix
    //

    // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
    // It doesn't really seem to matter much either way.
#define TMV_BD_ZeroIX (ix == 0)

    template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper;

    // algo 0: cs or rs = 0, so nothing to do
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<0,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const M2& , M3& ) 
        {} 
    };

    // algo 11: Loop over columns
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<11,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = (cs == Unknown ? m3.colsize() : cs);
            const int N = (rs == Unknown ? m3.rowsize() : rs);
#ifdef PRINTALGO_BD
            std::cout<<"BD algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3c::iterator IT3;
            const int rowstep1 = m1.stepj();
            const int rowstep3 = m3.stepj();
            const int diagstep1 = m1.diagstep();
            const int diagstep3 = m3.diagstep();

            const int xx = Unknown;
            const int c2 = M2::_conj;
            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m1.nhi();
            const int j2 = TMV_MIN(N,M-m1.nlo());
            const int j3 = TMV_MIN(N,M+m1.nhi());

            int len = m1.nlo()+1;
            IT1 it1 = m1.get_col(0,0,len).begin().nonConj();
            IT2 it2 = m2.diag().begin().nonConj();
            IT3 it3 = m3.get_col(0,0,len).begin();

            PT2 dj;
            int j=0;
            for(;j<j1;++j) {
                dj = ZProd<false,c2>::prod(x,*it2++);
                MultXV_Helper<-4,xx,add,0,PT2,M1c,M3c>::call2(len,dj,it1,it3);
                it1.shiftP(rowstep1);
                it3.shiftP(rowstep3);
                if (len < M) ++len;
            } 
            if (j1 < j2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(;j<j2;++j) {
                dj = ZProd<false,c2>::prod(x,*it2++);
                MultXV_Helper<-4,lh,add,0,PT2,M1c,M3c>::call2(len,dj,it1,it3);
                it1.shiftP(diagstep1);
                it3.shiftP(diagstep3);
            }
            if (j1 >= j2) ++len;
            for(;j<j3;++j) {
                dj = ZProd<false,c2>::prod(x,*it2++);
                MultXV_Helper<-4,xx,add,0,PT2,M1c,M3c>::call2(--len,dj,it1,it3);
                it1.shiftP(diagstep1);
                it3.shiftP(diagstep3);
            }
        }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<21,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = (cs == Unknown ? m3.colsize() : cs);
            const int N = (rs == Unknown ? m3.rowsize() : rs);
#ifdef PRINTALGO_BD
            std::cout<<"BD algo 21: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3r::iterator IT3;
            const int colstep1 = m1.stepi();
            const int colstep3 = m3.stepi();
            const int diagstep1 = m1.diagstep();
            const int diagstep3 = m3.diagstep();
 
            const int xx = Unknown;
            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m1.nlo();
            const int i2 = TMV_MIN(M,N-m1.nhi());
            const int i3 = TMV_MIN(M,N+m1.nlo());

            int len = m1.nhi()+1;
            IT1 it1 = m1.get_row(0,0,len).begin().nonConj();
            IT2 it2 = m2.diag().begin().nonConj();
            IT3 it3 = m3.get_row(0,0,len).begin();

            int i=0;
            for (;i<i1;++i) {
                ElemMultVV_Helper<-4,xx,add,ix,T,M1r,M2d,M3r>::call2(
                    len,x,it1,it2,it3);
                it1.shiftP(colstep1);
                it3.shiftP(colstep3);
                if (len < N) ++len;
            }
            if (i1 < i2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for (;i<i2;++i) {
                ElemMultVV_Helper<-4,lh,add,ix,T,M1r,M2d,M3r>::call2(
                    len,x,it1,it2++,it3);
                it1.shiftP(diagstep1);
                it3.shiftP(diagstep3);
            }
            if (i1 >= i2) ++len;
            for (;i<i3;++i) {
                ElemMultVV_Helper<-4,xx,add,ix,T,M1r,M2d,M3r>::call2(
                    --len,x,it1,it2++,it3);
                it1.shiftP(diagstep1);
                it3.shiftP(diagstep3);
            }
        }
    };

    // algo 31: Loop over diagonals
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<31,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = (cs == Unknown ? m3.colsize() : cs);
            const int N = (rs == Unknown ? m3.rowsize() : rs);
#ifdef PRINTALGO_BD
            std::cout<<"BD algo 31: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_diag_sub_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M3::diag_sub_type M3r;
            typedef typename M3r::iterator IT3;
            const int colstep1 = m1.stepi();
            const int colstep3 = m3.stepi();
            const int rowstep1 = m1.stepj();
            const int rowstep3 = m3.stepj();
 
            const int xx = Unknown;
            const int ds = IntTraits2<cs,rs>::min;
            int len = TMV_MIN(M-m1.nlo(),N);
            IT1 it1 = m1.get_diag(-m1.nlo()).begin().nonConj();
            IT2 it2 = m2.diag().begin().nonConj();
            IT3 it3 = m3.get_diag(-m1.nlo()).begin();

            for (int k=m1.nlo();k;--k) {
                ElemMultVV_Helper<-4,xx,add,ix,T,M1r,M2d,M3r>::call2(
                    len,x,it1,it2,it3);
                it1.shiftP(-colstep1);
                it3.shiftP(-colstep3);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            ElemMultVV_Helper<-4,ds,add,ix,T,M1r,M2d,M3r>::call2(
                len,x,it1,it2,it3);
            for (int k=1;k<=m1.nhi();++k) {
                it1.shiftP(rowstep1);
                it3.shiftP(rowstep3);
                if (k+len > N) --len;
                ElemMultVV_Helper<-4,xx,add,ix,T,M1r,M2d,M3r>::call2(
                    len,x,it1,++it2,it3);
            }
        }
    };

    template <int ix, class T, class M> class ProdXM;

    // algo 82: copy x*m2
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<82,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_BD
            const int N = (rs == Unknown ? m3.rowsize() : rs);
            const int M = (cs == Unknown ? m3.colsize() : cs);
            std::cout<<"BD algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::real_type RT;
            Scaling<1,RT> one;
            typedef typename Traits2<T,typename M2::value_type>::type PT;
            typedef typename MCopyHelper<PT,Diag,rs,rs>::type M2c;
            M2c m2c = ProdXM<ix,T,M2>(x,m2);
            MultBD_Helper<-2,cs,rs,add,1,RT,M1,M2c,M3>::call(one,m1,m2c,m3);
        }
    };

    // algo 90: call inst
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<90,cs,rs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<90,cs,rs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 91: call inst alias
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<91,cs,rs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<91,cs,rs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<97,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultBD_Helper<-2,cs,rs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<197,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultBD_Helper<99,cs,rs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<98,cs,rs,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultBD_Helper<-2,cs,rs,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultBD_Helper<82,cs,rs,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // SameStorage(m1,m3)
                // Need a temporary
                typedef typename M1::copy_type M1c;
                M1c m1c = m1;
                MultBD_Helper<-2,cs,rs,true,ix,T,M1c,M2,M3>::call(x,m1c,m2,m3);
            }
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<98,cs,rs,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultBD_Helper<-2,cs,rs,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultBD_Helper<82,cs,rs,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // SameStorage(m1,m3)
                // Let Copy handle the aliasing
                Copy(m1,m3);
                typedef typename M3::const_view_type M3c;
                M3c m3c = m3.constView();
                MultBD_Helper<-2,cs,rs,false,ix,T,M3c,M2,M3>::call(x,m3c,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<99,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
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
                M3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultBD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<-4,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
            const bool bothdm = M1::_diagmajor && M3::_diagmajor;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                TMV_OPT == 0 ? 11 :
                bothcm ? 11 :
                bothrm ? 21 :
                bothdm ? 31 :
                11;
#ifdef PRINTALGO_BD
            std::cout<<"XB: algo -4\n";
            std::cout<<"add = "<<add<<"  x = "<<ix<<" "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            MultBD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#ifdef PRINTALGO_BD
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<-3,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms to choose from:
            //  0 = cs or rs == 0, so nothing to do
            //  1 = cs == 1: reduces to trivial MultDV function
            //  2 = rs == 1: reduces to trivial MultXV function
            //
            // 11 = column major loop: MultXV on each column
            // 21 = row major loop: MultDV on each row
            // 31 = diag major loop: MultDV on each diagonal
            //
            // 82 = copy x*m2

            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
            const bool bothdm = M1::_diagmajor && M3::_diagmajor;
            const bool docopy = 
                TMV_OPT == 0 ? false :
                TMV_BD_ZeroIX || M2::_diagstep != 1;

            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                TMV_OPT == 0 ? 11 :
                bothcm ? 11 :
                bothrm ? (docopy ? 82 : 21) :
                bothdm ? (docopy ? 82 : 31) :
                11;
#ifdef PRINTALGO_BD
            std::cout<<"InlineMultBD: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
#endif
#ifdef XDEBUG_BD
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            MultMM<add>(x,m1c,m2c,m3c);
#endif
            MultBD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            if (m3.nlo() > m1.nlo()) 
                Maybe<!add>::zero2(m3.diagRange(-m3.nlo(),-m1.nlo()));
            if (m3.nhi() > m1.nhi()) 
                Maybe<!add>::zero2(m3.diagRange(m1.nhi()+1,m3.nhi()+1));
#ifdef PRINTALGO_BD
            if (m3.nlo() > m1.nlo() || m3.nhi() > m1.nhi()) 
                std::cout<<"After zeros: m2 => "<<m2<<std::endl;
#endif
#ifdef XDEBUG_BD
            if (Norm(m3-m3c) > 1.e-3*(Norm(m1c)*Norm(m2c)+(add?Norm(m3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"m2 = "<<m2c<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 => "<<m3<<std::endl;
                std::cout<<"Correct m3 = "<<m3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<-2,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
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
                M3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultBD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultBD_Helper<-1,cs,rs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                M3::_checkalias ? 99 : 
                -2;
            MultBD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rowsize() == m2.size());
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs1 = Sizes<M3::_rowsize,M1::_rowsize>::size;
        const int rs = Sizes<rs1,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBD_Helper<-1,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rowsize() == m2.size());
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs1 = Sizes<M3::_rowsize,M1::_rowsize>::size;
        const int rs = Sizes<rs1,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBD_Helper<-3,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rowsize() == m2.size());
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs1 = Sizes<M3::_rowsize,M1::_rowsize>::size;
        const int rs = Sizes<rs1,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultBD_Helper<98,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    // 
    // B = D * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    TMV_INLINE void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Band_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        MultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    // 
    // B *= D
    //
    
    template <class M1, int ix, class T, class M2>
    TMV_INLINE void MultEqMM(
        BaseMatrix_Band_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { MultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    //
    // M = B * D or D * B
    //
    
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        typename BMVO<M3>::b b3 = BandMatrixViewOf(m3,m1.nlo(),m1.nhi());
        MultMM<add>(x,m1,m2,b3);
        if (m1.nlo() < m1.colsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cRowRange(m1.nlo()+1,m1.colsize()),
                    m1.colsize()-m1.nlo()-2,0));
        if (m1.nhi() < m1.rowsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cColRange(m1.nhi()+1,m1.rowsize()),
                    0,m1.rowsize()-m1.nhi()-2));
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_rowsize>::same));
        TMVAssert(m3.colsize() == m2.colsize());
        TMVAssert(m3.rowsize() == m2.rowsize());
        typename BMVO<M3>::b b3 = BandMatrixViewOf(m3,m2.nlo(),m2.nhi());
        MultMM<add>(x,m1,m2,b3);
        if (m2.nlo() < m2.colsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cRowRange(m2.nlo()+1,m2.colsize()),
                    m2.colsize()-m2.nlo()-2,0));
        if (m2.nhi() < m2.rowsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m3.cColRange(m2.nhi()+1,m2.rowsize()),
                    0,m2.rowsize()-m2.nhi()-2));
    }

    //
    // U = B * D or D * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        const int lo = Maybe<M3::_upper>::select(m1.nlo(),m1.nhi());
        const int hi = Maybe<M3::_upper>::select(m1.nhi(),m1.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b3 = BandMatrixViewOf(m3,hi);
        MultMM<add>(x,m1,m2,b3);
        if (hi < m3.size()-1)
            Maybe<!add>::zero2(m3.offDiag(hi+1));
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m3.colsize() == m2.colsize());
        TMVAssert(m3.rowsize() == m2.rowsize());
        const int lo = Maybe<M3::_upper>::select(m2.nlo(),m2.nhi());
        const int hi = Maybe<M3::_upper>::select(m2.nhi(),m2.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b3 = BandMatrixViewOf(m3,hi);
        MultMM<add>(x,m1,m2,b3);
        if (hi < m3.size()-1)
            Maybe<!add>::zero2(m3.offDiag(hi+1));
    }


    //
    // D = B * D or D * B
    //

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3)
    {
        TMVAssert(m1.nlo() == 0);
        TMVAssert(m1.nhi() == 0);
        MultMM<add>(x,DiagMatrixViewOf(m1.diag()),m2,m3);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Band<M2>& m2, BaseMatrix_Diag_Mutable<M3>& m3)
    {
        TMVAssert(m2.nlo() == 0);
        TMVAssert(m2.nhi() == 0);
        MultMM<add>(x,m1,DiagMatrixViewOf(m2.diag()),m3);
    }

} // namespace tmv

#undef TMV_BD_ZeroIX

#endif 
