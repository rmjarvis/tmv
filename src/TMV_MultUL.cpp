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

//#define PRINTALGO_UL

#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SimpleMatrix.h"

// The most common reason to do this function is to basically undo an
// LU decomposition.  So something like:
// m2 = m.unitLowerTri() * m.UpperTri();
// This operation can be done in place.  So we can copy the L and U
// triangles into the final storage and then do the operation there.
//
// Thus, we would like to convert all possible calls of a MultMM(L,U,m)
// of MultMM(U,L,m) type into this in-place version.
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
    static inline void DoMultUL(const M1& m1, const M2& m2, MatrixView<T,1> m3)
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
                        UpperTriMatrixView<T,UnknownDiag,1>,
                        LowerTriMatrixView<T,UnknownDiag,1> >::type M2x;
                M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
                InstCopy(m2,m2x.xdView());
                InlineMultMM<false>(one,m1.nonConj().cmView(),m2x.cmView(),m3);
            } else if (m2.iscm() && !SameStorage(m2,m3) && !m2.isconj()) {
                // Then can copy m1
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrixView<T,UnknownDiag,1>,
                        LowerTriMatrixView<T,UnknownDiag,1> >::type M1x;
                M1x m1x = Maybe<M1::_upper>::uppertri(m3,m1.dt());
                InstCopy(m1,m1x.xdView());
                InlineMultMM<false>(one,m1x.cmView(),m2.nonConj().cmView(),m3);
            } else if (m1.iscm() && !m1.isconj()) {
                // Need temporary storage for m2.
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrix<T,NonUnitDiag,ColMajor>,
                        LowerTriMatrix<T,NonUnitDiag,ColMajor> >::type M2c;
                M2c m2c(m2.size());
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrixView<T,UnknownDiag,1>,
                        LowerTriMatrixView<T,UnknownDiag,1> >::type M2cv;
                M2cv m2cv = m2c.viewAsUnknownDiag(m2.dt());
                InstCopy(m2,m2cv.xdView());
                InlineMultMM<false>(one,m1.nonConj().cmView(),m2cv,m3);
            } else {
                // Need temporary storage for m1, and can copy m2
                typedef typename TypeSelect<M2::_upper,
                        UpperTriMatrixView<T,UnknownDiag,1>,
                        LowerTriMatrixView<T,UnknownDiag,1> >::type M2x;
                M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
                InstCopy(m2,m2x.xdView());

                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrix<T,NonUnitDiag,ColMajor>,
                        LowerTriMatrix<T,NonUnitDiag,ColMajor> >::type M1c;
                M1c m1c(m1.size());
                typedef typename TypeSelect<M1::_upper,
                        UpperTriMatrixView<T,UnknownDiag,1>,
                        LowerTriMatrixView<T,UnknownDiag,1> >::type M1cv;
                M1cv m1cv = m1c.viewAsUnknownDiag(m1.dt());
                InstCopy(m1,m1cv.xdView());
                InlineMultMM<false>(one,m1cv,m2x.cmView(),m3);
            }
        } else {
            // No clash of the diagonals.  Copy both to m3, and do
            // the calculation in place.
            typedef typename TypeSelect<M1::_upper,
                    UpperTriMatrixView<T,UnknownDiag,1>,
                    LowerTriMatrixView<T,UnknownDiag,1> >::type M1x;
            M1x m1x = Maybe<M1::_upper>::uppertri(m3,m1.dt());
            InstCopy(m1,m1x.xdView());

            typedef typename TypeSelect<M2::_upper,
                    UpperTriMatrixView<T,UnknownDiag,1>,
                    LowerTriMatrixView<T,UnknownDiag,1> >::type M2x;
            M2x m2x = Maybe<M2::_upper>::uppertri(m3,m2.dt());
            InstCopy(m2,m2x.xdView());

            InlineMultMM<false>(one,m1x.cmView(),m2x.cmView(),m3);
        }
    }

    template <class T, class M1, class M2>
    static inline void DoInstMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        if (m3.iscm()) {
            DoMultUL(m1,m2,m3.cmView());
            InstScale(x,m3);
        } else if (m3.isrm()) {
            DoMultUL(m2.transpose(),m1.transpose(),m3.transpose().cmView());
            InstScale(x,m3);
        } else {
            SimpleMatrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
            DoMultUL(m1,m2,m3c.cmView());
            InstMultXM(x,m3c.constView().xView(),m3.xView());
        }
    }

    template <class T, class M1, class M2>
    static inline void DoInstAddMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        SimpleMatrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
        DoMultUL(m1,m2,m3c.cmView());
        InstAddMultXM(x,m3c.constView().xView(),m3.xView());
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { DoInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { DoInstAddMultMM(x,m1,m2,m3); }


#define InstFile "TMV_MultUL.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


