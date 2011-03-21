//////////////////////////////////////////////////////////////////////////////
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

//#define PRINTALGO_UU

#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_ScaleU.h"

namespace tmv {

    template <class M2, class M3>
    static inline void DoMultEqUU(const M2& m2, M3& m3)
    {
        typedef typename M2::value_type T2;
        typedef typename M3::real_type RT;
        Scaling<1,RT> one;
        if (m2.iscm()) 
            InlineMultMM<false>(one,m3.constView(),m2.cmView(),m3);
        else if (m2.isrm())
            InlineMultMM<false>(one,m3.constView(),m2.rmView(),m3);
        else {
            DoMultEqUU(m2.copy().xdView(),m3);
        }
    }

    template <class M1, class M2, class M3>
    static inline void DoMultUU(const M1& m1, const M2& m2, M3& m3)
    {
        typedef typename M2::value_type T2;
        typedef typename M3::real_type RT;

        // Can save some compiled code by copying m1 to m3.
        // This means the InlineMultMM call will always have m1 and m3
        // with the same majority.  And m1 will never be conjugated.
        //
        // Given the alias checking that has been done already,
        // there are some cases that were allowed that won't work
        // with the copy.
        //
        // Specifically, if diagonal of m2 has the same storage as m3,
        // then we need to make sure the copy won't change it.
        // (Remember that the only way m2 can have the same storage as m3
        // is if it is in the opposite triangle, so the diagonal is all 
        // we have to worry about.)
        if ( SameStorage(m2,m3) && !m2.isunit() &&
             (M1::_conj || m1.isunit() || !SameStorage(m1,m3)) ) {
            // Then the copy will clobber m2.diag().  
            // Need to use a temporary for m2.
            DoMultUU(m1,m2.copy().xdView(),m3);
        } else {
            // OK to copy m1 to m3.  However, the copy might be trivial.
            // Check if m1 is same storage as m3.
            if ( SameStorage(m1,m3) ) {
                if (ExactSameStorage(m1,m3)) {
                    // Just check diagonal and conj.
                    Maybe<M1::_conj>::conjself(m3);
                    if (m1.isunit() && !m3.isunit()) m3.diag().setAllTo(RT(1));
                } else { // Opposite storage
                    if (m1.size() > 1)
                        InstCopy(m1.offDiag().xdView(),m3.offDiag().xdView());
                    if (!m3.isunit()) {
                        if (m1.isunit()) m3.diag().setAllTo(RT(1));
                        else Maybe<M1::_conj>::conjself2(m3.diag());
                    }
                }
            } else {
                InstCopy(m1,m3.xdView());
            }
            DoMultEqUU(m2,m3);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,UnknownDiag> m3)
    {
        if (m3.iscm()) {
            UpperTriMatrixView<T3,UnknownDiag,1> m3cm = m3.cmView();
            DoMultUU(m1,m2,m3cm);
            if (x != T3(1)) InstScale(x,m3.viewAsNonUnitDiag());
        } else if (m3.isrm()) {
            UpperTriMatrixView<T3,UnknownDiag,UNKNOWN,1> m3rm = m3.rmView();
            DoMultUU(m1,m2,m3rm);
            if (x != T3(1)) InstScale(x,m3.viewAsNonUnitDiag());
        } else {
            SimpleMatrix<T3,ColMajor> m3c(m3.size(),m3.size(),T3(0));
            UpperTriMatrixView<T3,UnknownDiag,1> m3ct = m3c.upperTri(m3.dt());
            DoMultUU(m1,m2,m3ct);
            if (m3.isunit())
                InstCopy(m3ct.constView().xdView(),m3);
            else
                InstMultXM(x,m3ct.constView().xdView(),m3.viewAsNonUnitDiag());
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,NonUnitDiag> m3)
    {
        typedef typename Traits2<T1,T2>::type T12;
        SimpleMatrix<T12,ColMajor> m1c(m3.size(),m3.size());
        UpperTriMatrixView<T12,UnknownDiag,1> m1ct = m1c.upperTri(
            (m1.isunit() && m2.isunit()) ? UnitDiag : NonUnitDiag);
        InstCopy(m1,m1ct.xdView());
        DoMultUU(m1,m2,m1ct);
        InstAddMultXM(x,m1ct.constView().xdView(),m3);
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,UnknownDiag> m3)
    {
        if (m3.iscm()) {
            LowerTriMatrixView<T3,UnknownDiag,1> m3cm = m3.cmView();
            DoMultUU(m1,m2,m3cm);
            if (x != T3(1)) InstScale(x,m3.viewAsNonUnitDiag());
        } else if (m3.isrm()) {
            LowerTriMatrixView<T3,UnknownDiag,UNKNOWN,1> m3rm = m3.rmView();
            DoMultUU(m1,m2,m3rm);
            if (x != T3(1)) InstScale(x,m3.viewAsNonUnitDiag());
        } else {
            SimpleMatrix<T3,ColMajor> m3c(m3.size(),m3.size(),T3(0));
            LowerTriMatrixView<T3,UnknownDiag,1> m3ct = m3c.lowerTri(m3.dt());
            DoMultUU(m1,m2,m3ct);
            if (m3.isunit())
                InstCopy(m3ct.constView().xdView(),m3);
            else
                InstMultXM(x,m3ct.constView().xdView(),m3.viewAsNonUnitDiag());
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,NonUnitDiag> m3)
    {
        typedef typename Traits2<T1,T2>::type T12;
        SimpleMatrix<T12,ColMajor> m1c(m3.size(),m3.size());
        LowerTriMatrixView<T12,UnknownDiag,1> m1ct = m1c.lowerTri(
            (m1.isunit() && m2.isunit()) ? UnitDiag : NonUnitDiag);
        InstCopy(m1,m1ct.xdView());
        DoMultUU(m1,m2,m1ct);
        InstAddMultXM(x,m1ct.constView().xdView(),m3);
    }


#define InstFile "TMV_MultUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


