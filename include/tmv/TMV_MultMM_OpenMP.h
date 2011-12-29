
#ifndef TMV_MultMM_OPENMP_H
#define TMV_MultMM_OPENMP_H
#ifdef _OPENMP

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_MultMM.h"

namespace tmv {

    // Defined in TMV_MultMM_OpenMP.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_OpenMP(
        const T3 x, 
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_OpenMP(
        const T3 x, 
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    //
    // Algo 69: OpenMPMultMM
    //

    template <int algo, int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper;

    // algo 69: Split the output matrix, m3, into strips.  
    // Use either row strips or columns strips according to whether there
    // are more rows or columns.  Also, round each strip up to a multiple 
    // of 16 rows or columns to maximize the efficiency of blocking in 
    // each section.
    // Then call the normal block algorithm for each section.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<69,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
#ifdef PRINTALGO_MM
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 69: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

#ifdef TMV_MM_USE_RECURSIVE_BLOCK
            const int Mb = cs == TMV_UNKNOWN ? TMV_UNKNOWN : (cs >> 6);
            const int Nb = rs == TMV_UNKNOWN ? TMV_UNKNOWN : (rs >> 6);
            const int Kb = xs == TMV_UNKNOWN ? TMV_UNKNOWN : (xs >> 6);
            const int Kb2 = IntTraits2<Kb,Kb>::prod;
            const int MbNbKb2 = IntTraits2<IntTraits2<Mb,Nb>::prod,Kb2>::prod;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;

            // If we are in this function, then we pretty much know that 
            // we want to do one of the large matrix algorithms (63,64,68)
            // for the sub-problems.  So we call algo 72 to determine which
            // one to use.  
            // The algo1 selection here mimics that selection
            // when the sizes are known.
            // It's not perfect, since it uses the unthreaded M,N,K, rather
            // than either M/nthreads or N/nthreads, but it should usually
            // select a pretty good algorithm, and it keeps the 
            // compiler from instantiating the three possible algorithms
            // due to the if statements in algo 72.
            const int algo1 = 
                inst ? -2 : 
                (cs == TMV_UNKNOWN || rs == TMV_UNKNOWN || xs == TMV_UNKNOWN) ? 72 :
#ifdef TMV_MM_USE_WINOGRAD
                (cs >= TMV_MM_MIN_WINOGRAD && rs >= TMV_MM_MIN_WINOGRAD && 
                 xs >= TMV_MM_MIN_WINOGRAD) ? 68 :
#endif
#ifdef TMV_MM_USE_RECURSIVE_BLOCK
                (MbNbKb2 >= TMV_MM_MIN_RECURSIVE) ? 63 :
#endif
                64;

#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int mythread = omp_get_thread_num();
                if (num_threads == 1) {
                    MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                        x,m1,m2,m3);
                } else if (M > N) {
                    int Mx = M / num_threads;
                    Mx = ((((Mx-1)>>4)+1)<<4); // round up to mult of 16
                    int i1 = mythread * Mx;
                    int i2 = (mythread+1) * Mx;
                    if (i2 > M || mythread == num_threads-1) i2 = M;
                    if (i1 < M) {
                        // Need to make sure, since we rounded up Mx!
                        typedef typename M1::const_rowrange_type M1r;
                        typedef typename M3::rowrange_type M3r;
                        const int csx = TMV_UNKNOWN; 
                        M1r m1r = m1.cRowRange(i1,i2);
                        M3r m3r = m3.cRowRange(i1,i2);

                        MultMM_Helper<
                            algo1,csx,rs,xs,add,ix,T,M1r,M2,M3r>::call(
                                x,m1r,m2,m3r);
                    }
                } else {
                    int Nx = N / num_threads;
                    Nx = ((((Nx-1)>>4)+1)<<4); 
                    int j1 = mythread * Nx;
                    int j2 = (mythread+1) * Nx;
                    if (j2 > N || mythread == num_threads-1) j2 = N;
                    if (j1 < N)  {
                        typedef typename M2::const_colrange_type M2c;
                        typedef typename M3::colrange_type M3c;
                        const int rsx = TMV_UNKNOWN; 
                        M2c m2c = m2.cColRange(j1,j2);
                        M3c m3c = m3.cColRange(j1,j2);
                        MultMM_Helper<
                            algo1,cs,rsx,xs,add,ix,T,M1,M2c,M3c>::call(
                                x,m1,m2c,m3c);
                    }
                }
            }
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<90,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM_OpenMP(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<90,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM_OpenMP(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -3: Only one algorithm here, so do it.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<-3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 69;
            MultMM_OpenMP_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                inst ? 90 :
                -3;
            MultMM_OpenMP_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_OpenMP_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T> x, const M1& m1, const M2& m2, M3& m3)
        {
            MultMM_OpenMP_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_OpenMP(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultMM_OpenMP_Helper<-2,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM_OpenMP(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultMM_OpenMP_Helper<-3,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(
            x,m1v,m2v,m3v);
    }



} // namespace tmv

#endif
#endif 
