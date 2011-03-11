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


#ifndef TMV_LUInverse_H
#define TMV_LUInverse_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultUL.h"
#include "TMV_InvertU.h"
#include "TMV_Permutation.h"

#ifdef PRINTALGO_LU
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined below:
    template <class M1> 
    void LU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P);
    template <class M1>
    void InlineLU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P);

    template <class M1, class M2> 
    void LU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2> 
    void InlineLU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2);

    // Defined in TMV_LUInverse.cpp
    template <class T1>
    void InstLU_Inverse(MatrixView<T1> m1, const Permutation& P);

    template <class T1, class T2, bool C1> 
    void InstLU_InverseATA(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P, 
        const bool trans, MatrixView<T2> m2);


    template <int algo, int cs, int rs, class M2>
    struct LU_Inverse_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<0,cs,rs,M1>
    { static void call(M1& , const Permutation& ) {} };

    // algo 11: Normal case
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<11,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUInverse algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // m1 = (PLU)^-1
            //    = U^-1 L^-1 Pt
            typename M1::uppertri_type U = m1.upperTri();
            typename M1::unit_lowertri_type L = m1.unitLowerTri();
            U.invertSelf();
            L.invertSelf();
            const Scaling<1,typename M1::real_type> one;
            NoAliasMultMM<false>(one,U,L,m1);
            P.inverse().applyOnRight(m1);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<-3,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                11;
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUInverse\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LU_Inverse_Helper<algo,cs,rs,M1>::call(m1,P);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<97,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        { 
            typedef typename M1::conjugate_type M1c;
            M1c m1c = m1.conjugate();
            LU_Inverse_Helper<-1,cs,rs,M1c>::call(m1c,P);
        }
    };

    // algo 98: call InstLU_Inverse
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<98,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        { InstLU_Inverse(m1.xView(),P); }
    };

    // algo -1: Check for inst
    template <int cs, int rs, class M1>
    struct LU_Inverse_Helper<-1,cs,rs,M1>
    {
        static void call(M1& m1, const Permutation& P)
        {
            typedef typename M1::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
                Traits<T2>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M1::_conj ? 97 :
                inst ? 98 :
                -3;
            LU_Inverse_Helper<algo,cs,rs,M1>::call(m1,P);
        }
    };

    template <class M1> 
    static void InlineLU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P)
    {
        const int cs = M1::_colsize;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        M1v m1v = m1.cView();
        LU_Inverse_Helper<-3,cs,rs,M1v>::call(m1v,P);
    }

    template <class M1> 
    static void LU_Inverse(
        BaseMatrix_Rec_Mutable<M1>& m1, const Permutation& P)
    {
        const int cs = M1::_colsize;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        M1v m1v = m1.cView();
        LU_Inverse_Helper<-1,cs,rs,M1v>::call(m1v,P);
    }

    template <int algo, int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<0,cs,rs,M1,M2>
    { static void call(const M1& , const Permutation& , const bool, M2& ) {} };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<11,cs,rs,M1,M2>
    {
        static void call(const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUInverseATA algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"m1 (=LU) = "<<m1<<std::endl;
            std::cout<<"trans = "<<trans<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // (At A)^-1 = A^-1 (A^-1)t
            // = (U^-1 L^-1 Pt) (P L^-1t U^-1t)
            // = U^-1 L^-1 L^-1t U^-1t
            //
            // if PLU is really AT, then
            // A^-1 = P L^-1T U^-1T
            // (At A)^-1 = P L^-1T U^-1T U^-1* L^-1* Pt

            typename M1::const_unit_lowertri_type L = m1.unitLowerTri();
            typename M1::const_uppertri_type U = m1.upperTri();

            if (trans) {
                typename M2::uppertri_type uinv = m2.upperTri();
#ifdef PRINTALGO_LU
                std::cout<<"u = "<<U<<std::endl;
#endif
                uinv = U.inverse();
#ifdef PRINTALGO_LU
                std::cout<<"uinv = "<<uinv<<std::endl;
#endif
                //m2 = uinv.transpose() * uinv.conjugate();
                // This has to be done this way becuase only the second
                // matrix can be in the opposite storage for m=L*U,
                // and if these are known sized matrices, the alias check 
                // won't get called.
                m2.transpose() = uinv.adjoint() * uinv;
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                m2 /= L.transpose();
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                m2 %= L.conjugate();
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                P.inverse().applyOnRight(m2);
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                P.applyOnLeft(m2);
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
            } else {
                typename M2::unit_lowertri_type linv = m2.unitLowerTri();
#ifdef PRINTALGO_LU
                std::cout<<"l = "<<L<<std::endl;
#endif
                linv = L.inverse();
#ifdef PRINTALGO_LU
                std::cout<<"linv = "<<linv<<std::endl;
#endif
                m2 = linv * linv.adjoint();
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                m2 /= U;
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
                m2 %= U.adjoint();
#ifdef PRINTALGO_LU
                std::cout<<"m2 => "<<m2<<std::endl;
#endif
            }
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<-3,cs,rs,M1,M2>
    {
        static void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                11;
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUInverseATA\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LU_InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,P,trans,m2);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<97,cs,rs,M1,M2>
    {
        static void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LU_InverseATA_Helper<-1,cs,rs,M1c,M2c>::call(m1c,P,trans,m2c);
        }
    };

    // algo 98: call InstLU_Inverse
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<98,cs,rs,M1,M2>
    {
        static void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        { 
            TMVAssert(m1.iscm());
            InstLU_InverseATA(m1.xView().cmView(),P,trans,m2.xView()); 
        }
    };

    // algo -1: Check for inst
    template <int cs, int rs, class M1, class M2>
    struct LU_InverseATA_Helper<-1,cs,rs,M1,M2>
    {
        static void call(
            const M1& m1, const Permutation& P, const bool trans, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                M2::_conj ? 97 :
                inst ? 98 :
                -3;
            LU_InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,P,trans,m2);
        }
    };

    template <class M1, class M2> 
    static void InlineLU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        LU_InverseATA_Helper<-3,cs,rs,M1v,M2v>::call(m1v,P,trans,m2v);
    }

    template <class M1, class M2> 
    static void LU_InverseATA(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        const bool trans, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        LU_InverseATA_Helper<-1,cs,rs,M1v,M2v>::call(m1v,P,trans,m2v);
    }

} // namespace tmv

#endif

