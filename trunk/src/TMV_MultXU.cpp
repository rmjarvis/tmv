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

#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ScaleU.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ConjugateV.h"

namespace tmv {

    template <bool add, class T, class M1, class M2>
    static void DoMultXM(const T x, const M1& m1, M2& m2)
    {
        if (x == T(1))
            InlineMultXM<add>(Scaling<1,T>(x),m1,m2); 
        else if (x == T(-1))
            InlineMultXM<add>(Scaling<-1,T>(x),m1,m2); 
        else if (x == T(0))
            Maybe<!add>::zero(m2);
        else
            InlineMultXM<add>(Scaling<0,T>(x),m1,m2); 
    }

    template <bool add, class T, class M1, class M2>
    static void DoMultXM(const std::complex<T> x, const M1& m1, M2& m2)
    { 
        if (imag(x) == T(0)) {
            if (real(x) == T(1))
                InlineMultXM<add>(Scaling<1,T>(real(x)),m1,m2); 
            else if (real(x) == T(-1))
                InlineMultXM<add>(Scaling<-1,T>(real(x)),m1,m2); 
            else if (real(x) == T(0))
                Maybe<!add>::zero(m2);
            else
                InlineMultXM<add>(Scaling<0,T>(real(x)),m1,m2); 
        } else 
            InlineMultXM<add>(Scaling<0,std::complex<T> >(x),m1,m2); 
    }

    template <class T1, int C1, class T2>
    void InstMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        UpperTriMatrixView<T2,NonUnitDiag> m2)
    {
        if (m1.iscm() && m2.iscm()) {
            UpperTriMatrixView<T2,NonUnitDiag|ColMajor> m2cm = m2.cmView();
            DoMultXM<false>(x,m1.cmView(),m2cm);
        } else if (m1.isrm() && m2.isrm()) {
            UpperTriMatrixView<T2,NonUnitDiag|RowMajor> m2rm = m2.rmView();
            DoMultXM<false>(x,m1.rmView(),m2rm);
        } else {
            InstCopy(m1,m2.xView());
            InstScale(x,m2);
        }
    }

    template <class T1, int C1, class T2>
    void InstAddMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        UpperTriMatrixView<T2,NonUnitDiag> m2)
    {
        if (m2.iscm()) {
            UpperTriMatrixView<T2,NonUnitDiag|ColMajor> m2cm = m2;
            if (m1.iscm())
                DoMultXM<true>(x,m1.cmView(),m2cm);
            else 
                DoMultXM<true>(x,m1,m2cm);
        } else if (m2.isrm()) {
            UpperTriMatrixView<T2,NonUnitDiag|RowMajor> m2rm = m2;
            if (m1.isrm())
                DoMultXM<true>(x,m1.rmView(),m2rm);
            else 
                DoMultXM<true>(x,m1,m2rm);
        } else {
            DoMultXM<true>(x,m1,m2);
        }
    }



#define InstFile "TMV_MultXU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


