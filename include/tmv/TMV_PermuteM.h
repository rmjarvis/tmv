
#ifndef TMV_PermuteM_H
#define TMV_PermuteM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_SwapV.h"
#include "TMV_Prefetch.h"

//#define PRINTALGO_PERM

#ifdef PRINTALGO_PERM
#include <iostream>
#include "TMV_MatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_Matrix.cpp
    template <class T>
    void InstPermuteRows(
        MatrixView<T> m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2);
    template <class T>
    void InstReversePermuteRows(
        MatrixView<T> m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2);

    //
    // PermuteRows
    //
    
    // Defined in TMV_Matrix.cpp
    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper;

    // algo 11: Simple loop over columns
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<11,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"PermuteRows algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            for (ptrdiff_t j=0;j<N;++j) 
                m.get_col(j).permute(p,i1,i2);
        }
    };

    // algo 12: Simple loop over rows
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<12,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"PermuteRows algo 12: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            p += i1;
            for(ptrdiff_t i=i1;i<i2;++i,++p) {
                TMVAssert(*p < m.colsize());
                m.cSwapRows(i,*p);
            }
        }
    };

    // algo 13: Loop over columns in blocks, then over rows with a block
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<13,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"PermuteRows algo 13: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::row_type M1r;
            typedef typename M1r::iterator IT;

            ptrdiff_t N_32 = (N>>5); // N_32 = N/32
            const ptrdiff_t Nx = N - (N_32<<5); // Nx = N % 32
            const ptrdiff_t rsx = rs == Unknown ? Unknown : (rs % 32);
            IT it1 = m.get_row(0).begin();
            const ptrdiff_t stepi = m.stepi();
            it1.shiftP(i1*stepi);
            p += i1;
            if (N_32) do {
                const ptrdiff_t* pi = p;
                IT it1i = it1;
                for(ptrdiff_t i=i1;i<i2;++i,++pi) {
                    TMVAssert(*pi < m.colsize());
                    if (*pi != i) {
                        IT it2 = it1i;
                        it2.shiftP((*pi-i)*stepi);
                        SwapV_Helper<-3,32,M1r,M1r>::call2(32,it1i,it2);
                    }
                    it1i.shiftP(stepi);
                }
                it1 += 32;
            } while (--N_32);
            if (Nx) {
                for(ptrdiff_t i=i1;i<i2;++i,++p) {
                    TMVAssert(*p < m.colsize());
                    if (*p != i) {
                        IT it2 = it1;
                        it2.shiftP((*p-i)*stepi);
                        SwapV_Helper<-3,rsx,M1r,M1r>::call2(Nx,it1,it2);
                    }
                    it1.shiftP(stepi);
                }
            }
        }
    };

    // algo 14: Loop over rows with iterators
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<14,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"PermuteRows algo 14: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::row_type M1r;
            typedef typename M1r::iterator IT;
            IT it1 = m.get_row(0).begin();
            const ptrdiff_t stepi = m.stepi();
            it1.shiftP(i1*stepi);
            p += i1;
            for(ptrdiff_t i=i1;i<i2;++i,++p) {
                TMVAssert(*p < m.colsize());
                if (*p != i) {
                    IT it2 = it1;
                    it2.shiftP((*p-i)*stepi);
                    SwapV_Helper<-3,rs,M1r,M1r>::call2(N,it1,it2);
                }
                it1.shiftP(stepi);
            }
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<90,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"PermuteRows algo 90: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            InstPermuteRows(m.xView(),p,i1,i2); 
        } 
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<97,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"PermuteRows algo 97: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            PermuteRows_Helper<-2,cs,rs,Mc>::call(mc,p,i1,i2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<-3,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const int algo = 
                TMV_OPT == 0 ? 12 :
                (M1::_colmajor && (rs != Unknown && rs <= 32)) ? 11 :
                M1::_rowmajor ? 14 : M1::_colmajor ? 13 :
                12;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"InlinePermuteRows: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"i1,i2 = "<<i1<<','<<i2<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            PermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
#ifdef PRINTALGO_PERM
            //std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<-2,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"PermuteRows: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"i1,i2 = "<<i1<<','<<i2<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            PermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
#ifdef PRINTALGO_PERM
            //std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct PermuteRows_Helper<-1,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        { PermuteRows_Helper<-2,cs,rs,M1>::call(m,p,i1,i2); }
    };

    template <class M>
    inline void PermuteRows(
        BaseMatrix_Rec_Mutable<M>& m,
        const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        PermuteRows_Helper<-2,cs,rs,Mv>::call(mv,p,i1,i2); 
    }

    template <class M>
    inline void InlinePermuteRows(
        BaseMatrix_Rec_Mutable<M>& m, 
        const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        PermuteRows_Helper<-3,cs,rs,Mv>::call(mv,p,i1,i2); 
    }


    //
    // ReversePermuteRows
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper;

    // algo 11: Simple loop over columns
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<11,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"ReversePermuteRows algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            for (ptrdiff_t j=0;j<N;++j) 
                m.get_col(j).reversePermute(p,i1,i2);
        }
    };

    // algo 12: Simple loop over rows
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<12,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"ReversePermuteRows algo 12: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            p += i2-1;
            for(ptrdiff_t i=i2-1;i>=i1;--i,--p) {
                TMVAssert(*p < m.colsize());
                m.cSwapRows(i,*p);
            }
        }
    };

    // algo 13: Loop over columns in blocks, then over rows with a block
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<13,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"ReversePermuteRows algo 13: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::row_type M1r;
            typedef typename M1r::iterator IT;

            ptrdiff_t N_32 = (N>>5); // N_32 = N/32
            const ptrdiff_t Nx = N - (N_32<<5); // Nx = N % 32
            const ptrdiff_t rsx = rs == Unknown ? Unknown : (rs % 32);
            IT it1 = m.get_row(0).begin();
            const ptrdiff_t stepi = m.stepi();
            it1.shiftP((i2-1)*stepi);
            p += i2-1;

            if (N_32) do {
                const ptrdiff_t* pi = p;
                IT it1i = it1;
                for(ptrdiff_t i=i2-1;i>=i1;--i,--pi) {
                    TMVAssert(*pi < m.colsize());
                    if (*pi != i) {
                        IT it2 = it1i;
                        it2.shiftP((*pi-i)*stepi);
                        SwapV_Helper<-3,32,M1r,M1r>::call2(32,it1i,it2);
                    }
                    it1i.shiftP(-stepi);
                }
                it1 += 32;
            } while (--N_32);
            if (Nx) {
                for(ptrdiff_t i=i2-1;i>=i1;--i,--p) {
                    TMVAssert(*p < m.colsize());
                    if (*p != i) {
                        IT it2 = it1;
                        it2.shiftP((*p-i)*stepi);
                        SwapV_Helper<-3,rsx,M1r,M1r>::call2(Nx,it1,it2);
                    }
                    it1.shiftP(-stepi);
                }
            }
        }
    };

    // algo 14: Loop over rows with iterators
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<14,cs,rs,M1>
    {
        static void call(M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            std::cout<<"ReversePermuteRows algo 14: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::row_type M1r;
            typedef typename M1r::iterator IT;
            IT it1 = m.get_row(0).begin();
            const ptrdiff_t stepi = m.stepi();
            Prefetch_Write(it1.get()+i1*stepi);
            it1.shiftP((i2-1)*stepi);
            p += i2-1;
            for(ptrdiff_t i=i2-1;i>=i1;--i,--p) {
                TMVAssert(*p < m.colsize());
                if (*p != i) {
                    IT it2 = it1;
                    it2.shiftP((*p-i)*stepi);
                    SwapV_Helper<-3,rs,M1r,M1r>::call2(N,it1,it2);
                }
                it1.shiftP(-stepi);
            }
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<90,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        { 
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"ReversePermuteRows algo 90: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            InstReversePermuteRows(m.xView(),p,i1,i2); 
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<97,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"ReversePermuteRows algo 97: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReversePermuteRows_Helper<-2,cs,rs,Mc>::call(mc,p,i1,i2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<-3,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            const int algo = 
                TMV_OPT == 0 ? 12 :
                (M1::_rowmajor && (rs != Unknown && rs <= 32)) ? 14 :
                M1::_colmajor ? 11 : M1::_rowmajor ? 13 :
                12;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"InlineReversePermuteRows: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"i1,i2 = "<<i1<<','<<i2<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            ReversePermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
#ifdef PRINTALGO_PERM
            //std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<-2,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
#ifdef PRINTALGO_PERM
            const ptrdiff_t M = cs == Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m.rowsize() : rs;
            std::cout<<"ReversePermuteRows: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"i1,i2 = "<<i1<<','<<i2<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            ReversePermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
#ifdef PRINTALGO_PERM
            //std::cout<<"m = "<<m<<std::endl;
#endif
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct ReversePermuteRows_Helper<-1,cs,rs,M1>
    {
        static TMV_INLINE void call(
            M1& m, const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
        { ReversePermuteRows_Helper<-2,cs,rs,M1>::call(m,p,i1,i2); }
    };

    template <class M>
    inline void ReversePermuteRows(
        BaseMatrix_Rec_Mutable<M>& m,
        const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReversePermuteRows_Helper<-2,cs,rs,Mv>::call(mv,p,i1,i2); 
    }

    template <class M>
    inline void InlineReversePermuteRows(
        BaseMatrix_Rec_Mutable<M>& m, 
        const ptrdiff_t* p, const ptrdiff_t i1, const ptrdiff_t i2)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReversePermuteRows_Helper<-3,cs,rs,Mv>::call(mv,p,i1,i2); 
    }

} // namespace tmv

#endif
