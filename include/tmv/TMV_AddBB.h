

#ifndef TMV_AddBB_H
#define TMV_AddBB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_Scaling.h"
#include "TMV_AddVV.h"
#include "TMV_MultXB.h"
#include "TMV_MultXM_Funcs.h"

#ifdef PRINTALGO_AddBB
#include <iostream>
#include "TMV_BandMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_AddBB.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMM(
        const T3 x1, const ConstBandMatrixView<T1,C1>& m1,
        const T3 x2, const ConstBandMatrixView<T2,C2>& m2,
        BandMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMM(
        const T3 x1, const ConstBandMatrixView<T1,C1>& m1,
        const T3 x2, const ConstBandMatrixView<T2,C2>& m2,
        BandMatrixView<T3> m3);


    //
    // BandMatrix + BandMatrix
    //

    template <int algo, int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper;

    // algo 0: size == 0, nothing to do
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<0,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& , const M1& , 
            const Scaling<ix2,T2>& , const M2& , M3& )
        {
#ifdef PRINTALGO_AddBB
            std::cout<<"AddBB algo 0\n";
#endif
        }
    };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddBB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"AddBB algo 1: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::const_linearview_type M2l;
            typedef typename M3::linearview_type M3l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            M3l m3l = m3.linearView();
            InlineAddVV(x1,m1l,x2,m2l,m3l);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<11,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_AddBB
            std::cout<<"AddBB algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            TMVAssert(m1.nlo() == m2.nlo());
            TMVAssert(m1.nhi() == m2.nhi());
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;
            typedef typename M3c::iterator IT3;
            TMVStaticAssert(!M3c::_conj);
            const int rowstep1 = m1.stepj();
            const int rowstep2 = m2.stepj();
            const int rowstep3 = m3.stepj();
            const int diagstep1 = m1.diagstep();
            const int diagstep2 = m2.diagstep();
            const int diagstep3 = m3.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m1.nhi();
            const int j2 = TMV_MIN(N,M-m1.nlo());
            const int j3 = TMV_MIN(N,M+m1.nhi());
            int len = m1.nlo()+1;
            IT1 it1 = m1.get_col(0,0,len).begin().nonConj();
            IT2 it2 = m2.get_col(0,0,len).begin().nonConj();
            IT3 it3 = m3.get_col(0,0,len).begin();
            for(int j=0;j<j1;++j) {
                AddVV_Helper<-4,xx,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    len,x1,it1,x2,it2,it3);
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                it3.shiftP(rowstep3);
                if (len < M) ++len;
            }
            if (j1 < j2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(int j=j1;j<j2;++j) {
                AddVV_Helper<-4,lh,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    len,x1,it1,x2,it2,it3);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
                it3.shiftP(diagstep3);
            }
            for(int j=j2;j<j3;++j) {
                AddVV_Helper<-4,xx,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    --len,x1,it1,x2,it2,it3);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
                it3.shiftP(diagstep3);
            }
        }
    };

    // algo 12: Loop over rows
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<12,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_AddBB
            std::cout<<"AddBB algo 12: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            TMVAssert(m1.nlo() == m2.nlo());
            TMVAssert(m1.nhi() == m2.nhi());
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;
            typedef typename M3r::iterator IT3;
            TMVStaticAssert(!M3r::_conj);
            const int colstep1 = m1.stepi();
            const int colstep2 = m2.stepi();
            const int colstep3 = m3.stepi();
            const int diagstep1 = m1.diagstep();
            const int diagstep2 = m2.diagstep();
            const int diagstep3 = m3.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m1.nlo();
            const int i2 = TMV_MIN(M,N-m1.nhi());
            const int i3 = TMV_MIN(M,N+m1.nlo());
            int len = m1.nhi()+1;
            IT1 it1 = m1.get_row(0,0,len).begin().nonConj();
            IT2 it2 = m2.get_row(0,0,len).begin().nonConj();
            IT3 it3 = m3.get_row(0,0,len).begin();
            for(int i=0;i<i1;++i) {
                AddVV_Helper<-4,xx,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    len,x1,it1,x2,it2,it3);
                it1.shiftP(colstep1);
                it2.shiftP(colstep2);
                it3.shiftP(colstep3);
                if (len < N) ++len;
            }
            if (i1 < i2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(int i=i1;i<i2;++i) {
                AddVV_Helper<-4,lh,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    len,x1,it1,x2,it2,it3);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
                it3.shiftP(diagstep3);
            }
            for(int i=i2;i<i3;++i) {
                AddVV_Helper<-4,xx,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    --len,x1,it1,x2,it2,it3);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
                it3.shiftP(diagstep3);
            }
        }
    };

    // algo 13: Loop over diagonals
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<13,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
#ifdef PRINTALGO_AddBB
            std::cout<<"AddBB algo 13: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            TMVAssert(m1.nlo() == m2.nlo());
            TMVAssert(m1.nhi() == m2.nhi());
            typedef typename M1::const_diag_sub_type M1d;
            typedef typename M2::const_diag_sub_type M2d;
            typedef typename M3::diag_sub_type M3d;
            typedef typename M1d::const_nonconj_type::const_iterator IT1;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M3d::iterator IT3;
            TMVStaticAssert(!M3d::_conj);
            const int colstep1 = m1.stepi();
            const int colstep2 = m2.stepi();
            const int colstep3 = m3.stepi();
            const int rowstep1 = m1.stepj();
            const int rowstep2 = m2.stepj();
            const int rowstep3 = m3.stepj();
            IT1 it1 = m1.get_diag(-m1.nlo()).begin().nonConj();
            IT2 it2 = m2.get_diag(-m1.nlo()).begin().nonConj();
            IT3 it3 = m3.get_diag(-m1.nlo()).begin();
            int len = TMV_MIN(M-m1.nlo(),N);
            for(int k=m1.nlo();k;--k) {
                AddVV_Helper<-4,xx,ix1,T1,M1d,ix2,T2,M2d,M3d>::call2(
                    len,x1,it1,x2,it2,it3);
                it1.shiftP(-colstep1);
                it2.shiftP(-colstep2);
                it3.shiftP(-colstep3);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const int ds = IntTraits2<cs,rs>::min;
            AddVV_Helper<-4,ds,ix1,T1,M1d,ix2,T2,M2d,M3d>::call2(
                len,x1,it1,x2,it2,it3);
            for(int k=0;k<m1.nhi();++k) {
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                it3.shiftP(rowstep3);
                if (k+len >= N) --len;
                AddVV_Helper<-4,xx,ix1,T1,M1d,ix2,T2,M2d,M3d>::call2(
                    len,x1,it1,x2,it2,it3);
            }
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<90,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<91,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAliasAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<97,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddBB_Helper<-2,cs,rs,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<197,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddBB_Helper<99,cs,rs,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<98,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddBB
            std::cout<<"AddBB Check alias\n";
            std::cout<<"s1,2 = "<<s1<<"  "<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing 
#ifdef PRINTALGO_AddBB
                std::cout<<"No alias\n";
#endif
                AddBB_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            } else if (!s2) {
#ifdef PRINTALGO_AddBB
                std::cout<<"Do m1 first\n";
#endif
                // Alias with m1 only, do m1 first
                AliasMultXM<false>(x1,m1,m3);
                NoAliasMultXM<true>(x2,m2,m3);
            } else if (!s1) {
#ifdef PRINTALGO_AddBB
                std::cout<<"Do m2 first\n";
#endif
                // Alias with m2 only, do m2 first
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1,m3);
            } else {
#ifdef PRINTALGO_AddBB
                std::cout<<"Need temporary\n";
#endif
                // Need a temporary
                typename M1::copy_type m1c = m1;
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1c,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<99,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                conj ? 197 :
                inst ? 91 :
                98;
            AddBB_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<-4,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M3::value_type T3;
            const bool allrm = M1::_rowmajor && M2::_rowmajor && M3::_rowmajor;
            const bool allcm = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? 13 :
                allcm ? 11 :
                allrm ? 12 :
                13;
#ifdef PRINTALGO_AddBB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"AddBB algo -4: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddBB_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<-3,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            const int lo = TMV_MIN(m1.nlo(),m2.nlo());
            const int hi = TMV_MIN(m1.nhi(),m2.nhi());
            typedef typename M1::const_diagrange_type M1x;
            typedef typename M2::const_diagrange_type M2x;
            typedef typename M3::diagrange_type M3x;
            M1x m1x = m1.cDiagRange(-lo,hi+1);
            M2x m2x = m2.cDiagRange(-lo,hi+1);
            AddBB_Helper<-4,cs,rs,ix1,T1,M1x,ix2,T2,M2x,M3>::call(
                x1,m1,x2,m2,m3);
            if (m1.nlo() > m2.nlo()) {
                M1x m1a = m1.cDiagRange(-m1.nlo(),-lo);
                M3x m3a = m3.cDiagRange(-m1.nlo(),-lo);
                NoAliasMultXM<false>(x1,m1a,m3a);
            }
            if (m2.nlo() > m1.nlo()) {
                M2x m2b = m2.cDiagRange(-m2.nlo(),-lo);
                M3x m3b = m3.cDiagRange(-m2.nlo(),-lo);
                NoAliasMultXM<false>(x2,m2b,m3b);
            }
            if (m1.nhi() > m2.nhi()) {
                M1x m1c = m1.cDiagRange(hi+1,m1.nhi()+1);
                M3x m3c = m3.cDiagRange(hi+1,m1.nhi()+1);
                NoAliasMultXM<false>(x1,m1c,m3c);
            }
            if (m2.nhi() > m1.nhi()) {
                M2x m2d = m2.cDiagRange(hi+1,m2.nhi()+1);
                M3x m3d = m3.cDiagRange(hi+1,m2.nhi()+1);
                NoAliasMultXM<false>(x2,m2d,m3d);
            }
            const int maxlo = TMV_MAX(m1.nlo(),m2.nlo());
            const int maxhi = TMV_MAX(m1.nhi(),m2.nhi());
            if (m3.nlo() > maxlo)
                m3.diagRange(-m3.nlo(),-maxlo).setZero();
            if (m3.nhi() > maxhi)
                m3.diagRange(maxhi+1,m3.nhi()+1).setZero();
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                conj ? 97 :
                inst ? 90 :
                -3;
#ifdef PRINTALGO_AddBB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"AddBB algo -2: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            AddBB_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
#ifdef PRINTALGO_AddBB
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddBB_Helper<-1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = 
                M3::_checkalias ? 99 : 
                -2;
#ifdef PRINTALGO_AddBB
            const int M = cs==TMV_UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"AddBB algo -1: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddBB_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() <= m3.nlo());
        TMVAssert(m2.nlo() <= m3.nlo());
        TMVAssert(m1.nhi() <= m3.nhi());
        TMVAssert(m2.nhi() <= m3.nhi());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddBB_Helper<-1,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() <= m3.nlo());
        TMVAssert(m2.nlo() <= m3.nlo());
        TMVAssert(m1.nhi() <= m3.nhi());
        TMVAssert(m2.nhi() <= m3.nhi());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddBB_Helper<-2,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void InlineAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() <= m3.nlo());
        TMVAssert(m2.nlo() <= m3.nlo());
        TMVAssert(m1.nhi() <= m3.nhi());
        TMVAssert(m2.nhi() <= m3.nhi());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddBB_Helper<-3,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void InlineAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() <= m3.nlo());
        TMVAssert(m2.nlo() <= m3.nlo());
        TMVAssert(m1.nhi() <= m3.nhi());
        TMVAssert(m2.nhi() <= m3.nhi());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddBB_Helper<98,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Band<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Band<M2>& m2, 
        BaseMatrix_Band_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m1.nlo() <= m3.nlo());
        TMVAssert(m2.nlo() <= m3.nlo());
        TMVAssert(m1.nhi() <= m3.nhi());
        TMVAssert(m2.nhi() <= m3.nhi());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddBB_Helper<99,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

} // namespace tmv

#endif 
