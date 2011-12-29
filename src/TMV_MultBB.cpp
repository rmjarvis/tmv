
#include "tmv/TMV_MultBB.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultXB.h"
#include "tmv/TMV_ScaleB.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"

namespace tmv {

    template <bool add, int ix, class Tx, class T1, int C1, class T2, int C2, class T3>
    static void DoMultMM(
        const Scaling<ix,Tx>& x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3)
    {
        TMVAssert(m3.iscm() || m3.isdm());
        if (m3.iscm()) {
            BandMatrixView<T3,ColMajor> m3cm = m3.cmView();
            TMVAssert(m2.iscm() || m2.isrm() || m2.isdm());
            if (m2.iscm()) { // xcc
                TMVAssert(m1.iscm() || m1.isrm() || m1.isdm());
                if (m1.iscm())
                    InlineMultMM<add>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m1.isrm())
                    InlineMultMM<add>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m1.isdm())
                    InlineMultMM<add>(x,m1.dmView(),m2.cmView(),m3cm);
            } else if (m2.isrm()) { // crc
                TMVAssert(m1.iscm());
                InlineMultMM<add>(x,m1.cmView(),m2.rmView(),m3cm);
            } else {
                TMVAssert(m1.iscm());
                InlineMultMM<add>(x,m1.cmView(),m2.dmView(),m3cm);
            }
        } else if (m3.isdm()) {
            BandMatrixView<T3,DiagMajor> m3dm = m3.dmView();
            TMVAssert(m2.isdm());
            TMVAssert(m1.isdm());
            InlineMultMM<add>(x,m1.dmView(),m2.dmView(),m3dm);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3)
    {
#if TMV_OPT <= 2 
        m3.setZero();
        InstAddMultMM(T3(1),m1,m2,m3);
        if (x != T3(1)) InstScale(x,m3);
#else
        typedef typename Traits<T3>::real_type RT;
        const Scaling<1,RT> one;
        if (m3.colsize() > 0 && m3.rowsize() > 0) {
            if (m1.rowsize() > 0) {
                if (m3.iscm()) {
                    if (m2.iscm()) {
                        if (m1.iscm() || m1.isrm() || m1.isdm()) {
                            if (x == T3(1)) {
                                DoMultMM<false>(one,m1,m2,m3);
                            } else {
                                DoMultMM<false>(Scaling<0,T3>(x),m1,m2,m3);
                            }
                        } else {
                            BandMatrix<T3,RowMajor|NoDivider> m1c = x*m1;
                            DoMultMM<false>(one,m1c.constView().xView(),m2,m3);
                        }
                    } else if (m2.isrm()) {
                        if (m1.iscm()) {
                            if (x == T3(1)) {
                                DoMultMM<false>(one,m1,m2,m3);
                            } else {
                                DoMultMM<false>(Scaling<0,T3>(x),m1,m2,m3);
                            }
                        } else {
                            BandMatrix<T3,ColMajor|NoDivider> m1c = x*m1;
                            DoMultMM<false>(one,m1c.constView().xView(),m2,m3);
                        }
                    } else if (m2.isdm()) {
                        if (m1.iscm()) {
                            if (x == T3(1)) {
                                DoMultMM<false>(one,m1,m2,m3);
                            } else {
                                DoMultMM<false>(Scaling<0,T3>(x),m1,m2,m3);
                            }
                        } else {
                            BandMatrix<T3,ColMajor|NoDivider> m1c = x*m1;
                            DoMultMM<false>(one,m1c.constView().xView(),m2,m3);
                        }
                    } else {
                        BandMatrix<T3,ColMajor|NoDivider> m2c = x*m2;
                        InstMultMM(T3(1),m1,m2c.constView().xView(),m3);
                    }
                } else if (m3.isdm()) {
                    if (m1.isdm() && m2.isdm()) {
                        if (x == T3(1)) {
                            DoMultMM<false>(one,m1,m2,m3);
                        } else {
                            DoMultMM<false>(Scaling<0,T3>(x),m1,m2,m3);
                        }
                    } else if (m1.isdm()) {
                        BandMatrix<T3,DiagMajor|NoDivider> m2c = x*m2;
                        DoMultMM<false>(one,m1,m2c.constView().xView(),m3);
                    } else if (m2.isdm()) {
                        BandMatrix<T3,DiagMajor|NoDivider> m1c = x*m1;
                        DoMultMM<false>(one,m1c.constView().xView(),m2,m3);
                    } else {
                        const int lo = TMV_MIN(m3.colsize()-1,m1.nlo()+m2.nlo());
                        const int hi = TMV_MIN(m3.rowsize()-1,m1.nhi()+m2.nhi());
                        BandMatrix<T3,ColMajor|NoDivider> m3c(
                            m3.colsize(),m3.rowsize(),lo,hi);
                        InstMultMM(T3(1),m1,m2,m3c.xView());
                        InstMultXM(x,m3c.constView().xView(),m3);
                    }
                } else if (m3.isrm()) {
                    InstMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
                } else  {
                    const int lo = TMV_MIN(m3.colsize()-1,m1.nlo()+m2.nlo());
                    const int hi = TMV_MIN(m3.rowsize()-1,m1.nhi()+m2.nhi());
                    BandMatrix<T3,ColMajor|NoDivider> m3c(
                        m3.colsize(),m3.rowsize(),lo,hi);
                    InstMultMM(T3(1),m1,m2,m3c.xView());
                    InstMultXM(x,m3c.constView().xView(),m3);
                }
            } else {
                m3.setZero();
            }
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3)
    {
        typedef typename Traits<T3>::real_type RT;
        const Scaling<1,RT> one;
        if (m3.colsize() > 0 && m3.rowsize() > 0 && m1.rowsize() > 0) {
            if (m3.iscm()) {
                if (m2.iscm()) {
                    if (m1.iscm() || m1.isrm() || m1.isdm()) {
                        if (x == T3(1)) {
                            DoMultMM<true>(one,m1,m2,m3);
                        } else {
#if TMV_OPT <= 2
                            // To avoid instantiating x!=1 paths, we need to 
                            // copy something.  Arbitrarily pick m1.
                            BandMatrix<T3,RowMajor|NoDivider> m1c = x*m1;
                            DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
#else
                            DoMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
#endif
                        }
                    } else {
                        BandMatrix<T3,RowMajor|NoDivider> m1c = x*m1;
                        DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                    }
                } else if (m2.isrm()) {
                    if (m1.iscm()) {
                        if (x == T3(1)) {
                            DoMultMM<true>(one,m1,m2,m3);
                        } else {
#if TMV_OPT <= 2
                            // Copy m2, since CCC is more efficient than CRC
                            BandMatrix<T3,ColMajor|NoDivider> m2c = x*m2;
                            DoMultMM<true>(one,m1,m2c.constView().xView(),m3);
#else
                            DoMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
#endif
                        }
                    } else {
                        BandMatrix<T3,ColMajor|NoDivider> m1c = x*m1;
                        DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                    }
                } else if (m2.isdm()) {
                    if (m1.iscm()) {
                        if (x == T3(1)) {
                            DoMultMM<true>(one,m1,m2,m3);
                        } else {
#if TMV_OPT <= 2
                            // Copy m2, since CCC is more efficient than CDC
                            BandMatrix<T3,ColMajor|NoDivider> m2c = x*m2;
                            DoMultMM<true>(one,m1,m2c.constView().xView(),m3);
#else
                            DoMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
#endif
                        }
                    } else {
                        BandMatrix<T3,ColMajor|NoDivider> m1c = x*m1;
                        DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                    }
                } else {
                    BandMatrix<T3,ColMajor|NoDivider> m2c = x*m2;
                    InstAddMultMM(T3(1),m1,m2c.constView().xView(),m3);
                }
            } else if (m3.isdm()) {
                if (m1.isdm() && m2.isdm()) {
                    if (x == T3(1)) {
                        DoMultMM<true>(one,m1,m2,m3);
                    } else {
#if TMV_OPT <= 2
                        // Arbitrarily pick m1 to copy.
                        BandMatrix<T3,DiagMajor|NoDivider> m1c = x*m1;
                        DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
#else
                        DoMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
#endif
                    }
                } else if (m1.isdm()) {
                    BandMatrix<T3,DiagMajor|NoDivider> m2c = x*m2;
                    DoMultMM<true>(one,m1,m2c.constView().xView(),m3);
                } else if (m2.isdm()) {
                    BandMatrix<T3,DiagMajor|NoDivider> m1c = x*m1;
                    DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                } else {
                    const int lo = TMV_MIN(m3.colsize()-1,m1.nlo()+m2.nlo());
                    const int hi = TMV_MIN(m3.rowsize()-1,m1.nhi()+m2.nhi());
                    BandMatrix<T3,ColMajor|NoDivider> m3c(
                        m3.colsize(),m3.rowsize(),lo,hi);
                    InstMultMM(T3(1),m1,m2,m3c.xView());
                    InstAddMultXM(x,m3c.constView().xView(),m3);
                }
            } else if (m3.isrm()) {
                InstAddMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
            } else  {
                const int lo = TMV_MIN(m3.colsize()-1,m1.nlo()+m2.nlo());
                const int hi = TMV_MIN(m3.rowsize()-1,m1.nhi()+m2.nhi());
                BandMatrix<T3,ColMajor|NoDivider> m3c(
                    m3.colsize(),m3.rowsize(),lo,hi);
                InstMultMM(T3(1),m1,m2,m3c.xView());
                InstAddMultXM(x,m3c.constView().xView(),m3);
            }
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstBandMatrixView<T2,C2>& m2, BandMatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultBB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


