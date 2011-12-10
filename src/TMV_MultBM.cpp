
//#undef NDEBUG
//#include "TMV.h"

#include "TMV_Blas.h"
#include "tmv/TMV_MultBM.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultXB.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ScaleB.h"

#ifdef BLAS
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"
#endif

namespace tmv {

    template <bool add, int ix, class Tx, class T1, int C1, class T2, int C2, class T3>
    static void DoMultMM(
        const Scaling<ix,Tx>& x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        TMVAssert(m3.iscm() || m3.isrm());
        TMVAssert(m2.iscm() || m2.isrm());
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            if (m2.iscm()) { // xcc
                if (m1.iscm())
                    InlineMultMM<add>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m1.isrm())
                    InlineMultMM<add>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m1.isdm())
                    InlineMultMM<add>(x,m1.dmView(),m2.cmView(),m3cm);
                else {
                    BandMatrix<T1,ColMajor|NoDivider> m1c = m1;
                    InlineMultMM<add>(x,m1c.cmView(),m2.cmView(),m3cm);
                }
            } else { // crc
                if (m1.iscm()) 
                    InlineMultMM<add>(x,m1.cmView(),m2.rmView(),m3cm);
                else {
                    //BandMatrix<T1,ColMajor|NoDivider> m1cm = m1;
                    //InlineMultMM<add>(x,m1cm.cmView(),m2.rmView(),m3cm);
                    InlineMultMM<add>(x,m1,m2.rmView(),m3cm);
                }
            }
        } else {
            MatrixView<T3,RowMajor> m3rm = m3.rmView();
            if (m2.iscm()) { // rcr
                if (m1.isrm()) 
                    InlineMultMM<add>(x,m1.rmView(),m2.cmView(),m3rm);
                else {
                    //BandMatrix<T1,RowMajor|NoDivider> m1rm = m1;
                    //InlineMultMM<add>(x,m1rm.rmView(),m2.cmView(),m3rm);
                    InlineMultMM<add>(x,m1,m2.cmView(),m3rm);
                }
            } else { // rrr, crr
                if (m1.iscm())
                    InlineMultMM<add>(x,m1.cmView(),m2.rmView(),m3rm);
                else if (m1.isrm())
                    InlineMultMM<add>(x,m1.rmView(),m2.rmView(),m3rm);
                else {
                    //BandMatrix<T1,RowMajor|NoDivider> m1rm = m1;
                    //InlineMultMM<add>(x,m1rm.rmView(),m2.rmView(),m3rm);
                    InlineMultMM<add>(x,m1,m2.rmView(),m3rm);
                }
            }
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
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
                if (m3.iscm() || m3.isrm()) {
                    if (m2.iscm() || m2.isrm()) {
                        DoMultMM<false>(one,m1,m2,m3);
                        if (x != T3(1)) InstScale(x,m3);
                    } else {
                        Matrix<T3,ColMajor|NoDivider> m2c = x*m2;
                        DoMultMM<false>(one,m1,m2c.constView().xView(),m3);
                    }
                } else  {
                    Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
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
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        typedef typename Traits<T3>::real_type RT;
        const Scaling<1,RT> one;
        if (m3.colsize() > 0 && m3.rowsize() > 0 && m1.rowsize() > 0) {
            if (m3.iscm() || m3.isrm()) {
                if (m2.iscm() || m2.isrm()) {
                    if (x == T3(1)) {
                        DoMultMM<true>(one,m1,m2,m3);
                    } else {
#if TMV_OPT <= 2
                        // To avoid instantiating x!=1 paths, we need to 
                        // copy something.  m1 is probably the smallest, 
                        // since it is banded, so use that.
                        if (m3.iscm()) {
                            BandMatrix<T3,ColMajor|NoDivider> m1c = x*m1;
                            DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                        } else {
                            BandMatrix<T3,RowMajor|NoDivider> m1c = x*m1;
                            DoMultMM<true>(one,m1c.constView().xView(),m2,m3);
                        }
#else
                        DoMultMM<true>(Scaling<0,T3>(x),m1,m2,m3);
#endif
                    }
                } else {
                    Matrix<T3,ColMajor|NoDivider> m2c = x*m2;
                    DoMultMM<true>(one,m1,m2c.constView().xView(),m3);
                }
            } else  {
                Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
                InstMultMM(T3(1),m1,m2,m3c.xView());
                InstAddMultXM(x,m3c.constView().xView(),m3);
            }
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstBandMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultBM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


