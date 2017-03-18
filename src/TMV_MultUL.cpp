
//#define PRINTALGO_UL
//#define XDEBUG_UL

#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_ConjugateV.h"


// The most common reason to do this function is to basically undo an
// LU decomposition.  So something like:
// m2 = m.unitLowerTri() * m.upperTri();
// This operation can be done in place.  So we can copy the L and U
// triangles into the final storage and then do the operation there.
//
// Thus, we would like to convert all possible calls of a MultMM(L,U,m)
// or MultMM(U,L,m) type into this in-place version.
// The only real problem with doing that is that we could have 
// neither L or U be UnitDiag, and they could have different data on
// their diagonals.  In this case, it is impossible to do the in-place
// algorithm, because the diagonals clash.
//
// So we do a slightly inefficient thing for this.  We copy either 
// U or L to ColMajor storage to match m so that the call to InlineMultMM
// will be of the same signature as the in-place calls.

namespace tmv {

    template <class M1, class M2, class T>
    static inline void DoMultUL(
        const M1& m1, const M2& m2, MatrixView<T,ColMajor> m3)
    {
        Scaling<1,typename Traits<T>::real_type> one;

        if (!m1.isunit() && !m2.isunit()) {
            // Then we might have a clash of the diagonals...
            if (m1.iscm() && m2.iscm() && !m1.isconj() && !m2.isconj()) {
                // Then it's ok, since all cm anyway.
                InlineMultMM<false>(
                    one,m1.nonConj().cmView(),m2.nonConj().cmView(),m3);
            } else if (m1.iscm() && !SameStorage(m1,m3) && !m1.isconj()) {
                // Then can copy m2
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrixView<T,ColMajor>,
                        LowerTriMatrixView<T,ColMajor> >::type M2x;
                M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
                InstCopy(m2,m2x.xView());
                InlineMultMM<false>(one,m1.nonConj().cmView(),m2x.cmView(),m3);
            } else if (m2.iscm() && !SameStorage(m2,m3) && !m2.isconj()) {
                // Then can copy m1
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrixView<T,ColMajor>,
                        LowerTriMatrixView<T,ColMajor> >::type M1x;
                M1x m1x = Maybe<M1::_upper>::uppertri(m3,m1.dt());
                InstCopy(m1,m1x.xView());
                InlineMultMM<false>(one,m1x.cmView(),m2.nonConj().cmView(),m3);
            } else if (m1.iscm() && !m1.isconj()) {
                // Need temporary storage for m2.
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrix<T,NonUnitDiag|ColMajor>,
                        LowerTriMatrix<T,NonUnitDiag|ColMajor> >::type M2c;
                M2c m2c(m2.size());
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrixView<T,ColMajor>,
                        LowerTriMatrixView<T,ColMajor> >::type M2cv;
                M2cv m2cv = m2c.viewAsUnknownDiag(m2.dt());
                InstCopy(m2,m2cv.xView());
                InlineMultMM<false>(one,m1.nonConj().cmView(),m2cv,m3);
            } else {
                // Need temporary storage for m1, and can copy m2
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrix<T,NonUnitDiag|ColMajor>,
                        LowerTriMatrix<T,NonUnitDiag|ColMajor> >::type M1c;
                M1c m1c(m1.size());
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrixView<T,ColMajor>,
                        LowerTriMatrixView<T,ColMajor> >::type M1cv;
                M1cv m1cv = m1c.viewAsUnknownDiag(m1.dt());
                InstCopy(m1,m1cv.xView());

                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrixView<T,ColMajor>,
                        LowerTriMatrixView<T,ColMajor> >::type M2x;
                M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
                if (!SameStorage(m2,m3)) {
                    InstCopy(m2,m2x.xView());
                } else if (m2.isconj()) {
                    m2x.conjugateSelf();
                }

                InlineMultMM<false>(one,m1cv,m2x.cmView(),m3);
            }
        } else {
            // No clash of the diagonals.  Copy both to m3, and do
            // the calculation in place.
            typedef typename TypeSelect<M1::_upper,
                    UpperTriMatrixView<T,ColMajor>,
                    LowerTriMatrixView<T,ColMajor> >::type M1x;
            M1x m1x = Maybe<M1::_upper>::uppertri(m3,m1.dt());
            InstCopy(m1,m1x.xView());

            typedef typename TypeSelect<M2::_upper,
                    UpperTriMatrixView<T,ColMajor>,
                    LowerTriMatrixView<T,ColMajor> >::type M2x;
            M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
            InstCopy(m2,m2x.xView());

            InlineMultMM<false>(one,m1x.cmView(),m2x.cmView(),m3);
        }
    }

    template <class T, class M1, class M2>
    static inline void DoInstMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
#ifdef XDEBUG_LU
        std::cout<<"Start DoInstMultMM:\n";
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
        std::cout<<"ptrs = "<<m1.cptr()<<"  "<<m2.cptr()<<"  "<<m3.cptr()<<std::endl;
        Matrix<T> m3c = x * Matrix<T>(m1) * Matrix<T>(m2);
#endif
        if (m3.iscm()) {
            DoMultUL(m1,m2,m3.cmView());
            InstScale(x,m3);
        } else if (m3.isrm()) {
            DoMultUL(m2.transpose(),m1.transpose(),m3.transpose().cmView());
            InstScale(x,m3);
        } else {
            Matrix<T,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
            DoMultUL(m1,m2,m3c.cmView());
            InstMultXM(x,m3c.constView().xView(),m3.xView());
        }
#ifdef XDEBUG_LU
        std::cout<<"In DoInstMultMM: m3 => "<<m3<<std::endl;
        std::cout<<"m3c = "<<m3c<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(m3-m3c)<<std::endl;
        if (Norm(m3-m3c) > 1.e-3 * Norm(m1)*Norm(m2)) abort();
#endif
    }

    template <class T, class M1, class M2>
    static inline void DoInstAddMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        Matrix<T,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
        DoMultUL(m1,m2,m3c.cmView());
        InstAddMultXM(x,m3c.constView().xView(),m3.xView());
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstLowerTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstUpperTriMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultUL.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


