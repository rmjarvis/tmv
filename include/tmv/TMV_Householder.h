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


//
// This file contains the code to implement Householder reflections.
//
// A Householder reflection is a unitary matrix H which transforms
// a given vector x into y e0 = (y,0,0,0...0) where |y| = |x|
//
// H is of the form: 
//
// H = I - beta u ut = I - beta (x - y e0) (xt - y* e0t)
//
// H x = (1 - beta (|x|^2 - y* x0)) x + beta (y |x|^2 - |y|^2 x0) e0
//
// Since H x is supposed to equal y e0, the coefficient of x must be 0, so
//
// beta = 1 / (|x|^2 - y* x0)
//
// With this definition for beta, the coefficient of e0 is:
//
// (y |x|^2 - |y|^2 x0) / (|x|^2 - y* x0)
// = |y|^2 (y - x0) / y* (y - x0) = y
//
// Check that H is unitary:
//
// HtH = (I - beta* (x - y e0) (xt - y* e0t))(I - beta (x - y e0) (xt - y* e0t))
//     = I - (beta*+beta) (x xt - y e0 xt - y* x e0t + |y|^2 e0 e0t)
//         + |beta|^2 (x - y e0) (|x|^2 - y* x0 - y x0* + |y|^2) (xt - y* e0t)
//     = I - [beta + beta* - |beta|^2 (1/beta + 1/beta*)] (x-ye0)(xt-y* e0t)
//     = I - (beta + beta* - beta* - beta)(x-y e0)(xt-y* e0t)
//     = I
//
// Note that we still have a choice on the value for y.  The only constraint
// is that |y| = |x| (the vector norm).  So y = |x| exp(i theta)
//
// There are two options which make some sense: 
// choosing a real y ( +-|x| ) or choosing y so that beta is real
// ( 1/(|x|^2 +- |x||x0|) ).
//
// The better choice is to make beta real, since then the matrix is
// Hermitian as well as unitary (ie. H = Ht, H*H = I), which also
// leads to the determinant being -1, which is actually required for
// it to be a "reflection" matrix.  It also has the advantage of (slightly)
// quicker multiplies, since the multiplication of H involves a multiplication
// by beta.  This step takes half the time if beta is real.  Of course, it's
// not the bulk of the calculation for large matrices, hence the "(slightly)".
//
// Unfortunately, LaPack made the wrong choice here.  They chose to make
// y real.  The _only_ advantage (that I can see) for this choice is 
// in calculating the SVD of a matrix.  Then the bidiagonalization
// using Householder matrices goes directly to a real bidiagonal, so
// it doesn't require a second step using Givens matrices to make it real.
// This is a pretty minor advantage which doesn't offset the speed loss
// from the complex multiplies for the betas, and makes keeping track of 
// the determinant much harder.
//
// Prior to TMV v0.70, I chose the LaPack convention to keep compatibility
// with their routines.  However, with TMV v0.70, I decided to switch to 
// using the true reflection convention, so now I have real beta values.
//
// We still have two possible real values for beta:
// ( 1/(|x|^2 +- |x||x0|) ).
//
// We want to minimize the effects of rounding errors, which means maximize
// 1/beta = |x|^2 +- |x| |x0|
// Clearly the + option is larger, so that is the choice we use.
//
// The corresponding value of y is then: y = -|x| (x0/|x0|) 
// (If x0 == 0, use y = |x|.)
//
//
// The determinant of this Householder non-reflection is:
//
// det(H) = det( I - beta (x+|x|e0) (x+|x|e0)t )
//
// A simple, but not (to me) obvious identity is:
//
// det(I-alpha u ut) = 1 - alpha |u|^2
//
// Aside: The proof is pretty cool and comes from the following identity:
// (the so-called Matrix Determinant Lemma)
//   ( I   0 ) ( I+uvT  u ) (  I   0 ) = ( I    u   )
//   ( vT  1 ) (   0    1 ) ( -vT  1 )   ( 0  1+vTu )
// The determinant of the left side is det(I+uvT).
// The determinant of the right side is 1+vTu.
// Our special case comes from v = -alpha u.conj().
//
//
// Thus, det(H) = 1 - beta |x - y e0|^2 
// |x - y e0|^2 = |x0-y|^2 + |x_1..N|^2
//              = |x0|^2 + |x_1..N|^2 - x0* y - x0 y* + |y|^2
//              = |x|^2 - x0* y - x0 y* + |x|^2
//              = 2|x|^2 - 2 Re( x0* (-|x| x0/|x0|) )
//              = 2|x|^2 + 2 Re( |x0| |x| )
//              = 2|x|^2 + 2|x0||x|
// det(H) = 1 - 2 beta (|x|^2 + |x0||x|)
//        = 1 - 2 (|x|^2 + |x0||x|) / (|x|^2 + |x0||x|)
//        = -1
//
// There are two further complications to deal with in practice.
//
// First, we don't introduce new storage to store the u part of the 
// Householder matrix definition.  Instead, we store it in the same location
// as the original vector that we reflect into y e0.  However, this means
// that we don't have space to store the first element of u, since we
// need to keep the resulting y value there.  So the solution is to 
// recognize that there is another degeneracy that we can take advantage of.
// Specifically, if u is scaled by a value w, and beta is scaled by 1/w^2,
// then the Householder matrix is unchanged.
//
// This means that we can scale u so that the first element is some given
// value.  This is valid, since u0 can never be zero:  
// If x0 != 0:
//     u0 = x0 - y = x0 + |x| x0 / |x0| = x0 * (1 + |x|/|x0|) != 0
// If x0 == 0:
//     u0 = -y = -|x| != 0  (assuming |x| != 0, see below)
//
// So we can scale u by 1/u0, and get a version that has u0 == 1.
// In addition we also have to scale beta by u0^2.  
// Then all further usage of the Householder matrix can proceed
// with the implicit value that u(0) == 1, rather than the value that
// is actually stored in that location.
//
// Second, the above prescription runs into a problem in the case of an
// exactly zero input vector, since x0 == 0 and |x| == 0, so u0 = 0.
// Furthermore, if x0 != 0, but the rest of x is zero, then there is 
// no point in doing the reflection, since the vector is already in the
// right form.  So we can save some calculation by using H = I instead
// of the normal calculation.  (In this case, H is not a reflection, so 
// det = 1, not -1.)  We mark this case by setting beta = 0.
//

#ifndef TMV_Householder_H
#define TMV_Householder_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_Rank1VVM.h"

// TODO: Put functions with arithmetic into .cpp file with Inst, etc.
// so don't need these:
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_ProdVV.h"
#include "tmv/TMV_OProdVV.h"
#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_SumMM.h"

//#define XDEBUG_HOUSE

#ifdef XDEBUG_HOUSE
#include <iostream>
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_VectorIO.h"
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class V>
    static inline void HouseholderReflect(
        typename V::reference x0, V& u, typename V::real_type& beta)
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
#ifdef XDEBUG_HOUSE
        Vector<T> xx(u.size()+1);
        xx(0) = x0;
        xx.subVector(1,xx.size()) = u;
#endif

        // Finds the Householder matrix H which rotates v into y e0.
        // The vector v of the Householder matrix is stored in v,
        // except for the first element.
        // Beta is the return value.
        
        // Since this routine involves squares of elements, we risk overflow
        // and underflow problems if done naively.  The simplest (although
        // probably not the most efficient) solution is to scale all the 
        // intermediate values by the maximum abs value in the vector.

        RT scale = u.maxAbs2Element();
        RT absx0 = TMV_ABS(x0);
        if (absx0 > scale) scale = absx0;

        if (TMV_Underflow(scale)) {
            // Then the situation is hopeless, and we should just zero out
            // the whole vector.
            u.setZero();
            x0 = T(0);
            beta = RT(0);
            return;
        }

        // Determine normx = |x|
        RT invscale = RT(1)/scale;
        RT normsqx = u.normSq(invscale);

        // if all of x other than first element are 0, H is identity
        if (normsqx == RT(0)) {
            // Set them all to x explicitly in case underflow led to the 0.
            u.setZero();
            beta = RT(0);
            return;
        }

        // Finish calculating normx in usual case
        absx0 *= invscale;
        x0 *= invscale;
        RT normsqx0 = absx0*absx0;
        normsqx += normsqx0;
        RT normx = TMV_SQRT(normsqx);

        // y = -|x| (x0/|x0|) 
        T y = absx0 == RT(0) ? normx : -normx*x0/absx0;

        // beta = 1 / (|x|^2 + |x||x0|)
        // H = I - beta v vt
        // with u = x - y e0 in first column
        // Renormalize beta,u so that u(0) = 1
        T u0 = x0-y;
        RT normu0 = TMV_NORM(u0);
        beta = normu0 / (normsqx + normx * absx0);
        T invu0 = RT(1)/u0;

        // Sometimes this combination can underflow, so check.
        T scaled_invu0 = invu0 * invscale;
        if ((TMV_REAL(invu0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invu0))) ||
            (TMV_IMAG(invu0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invu0)))) {
            u *= invscale;
            u *= invu0;
        } else {
            u *= scaled_invu0;
        }
                        
        x0 = y*scale;

#ifdef XDEBUG_HOUSE
        Vector<T> uu(xx.size());
        uu(0) = T(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T> H = T(1)-beta*(uu^uu.conjugate());
        // Check the following:
        // H * xx = (y,0,0.0...)
        // |y| = Norm(xx)
        Vector<T> Hxx = H * xx;
        if (TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(x0)) > 0.0001*TMV_ABS(x0) ||
            TMV_ABS(Norm(xx)-TMV_ABS(x0)) > 0.0001*TMV_ABS(x0) ||
            Norm(Hxx.subVector(1,Hxx.size())) > 0.0001*TMV_ABS(x0)) {
            cerr<<"Householder Reflect:\n";
            cerr<<"Input: x = "<<xx<<endl;
            cerr<<"x0 = "<<x0<<endl;
            cerr<<"Norm(x) = "<<Norm(xx)<<endl;
            cerr<<"Output: uu= "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"xx = "<<xx<<endl;
            cerr<<"Hxx = "<<Hxx<<endl;
            cerr<<"abs(abs(hxx(0))-abs(x0)) = "<<
                TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(x0))<<endl;
            cerr<<"abs(Norm(xx)-abs(x0)) = "<<
                TMV_ABS(Norm(xx))-TMV_ABS(x0)<<endl;
            cerr<<"Norm(Hxx(1,N)) = "<<Norm(Hxx.subVector(1,Hxx.size()))<<endl;
            abort();
        }
#endif
    }

    template <class V>
    static inline void HouseholderReflect(V& x, typename V::real_type& beta)
    {
        typename V::reference x0 = x.ref(0);
        typename V::subvector_type u = x.subVector(1,x.size());
        HouseholderReflect(x0,u,beta);
    }

    template <class V>
    static inline bool HouseholderUnReflect(
        typename V::reference y, V& u, typename V::real_type& beta)
    {
        // This is similar, except that the roles of y and x0 are swapped.
        // That is, the final rotated vector y e0 is presumed known, as is 
        // the bulk of the input vector.  
        // The unreflected value of the first element in the vector is returned 
        // as y.
        // The rest of the vector is transformed into the Householder vector.
        // This is used for downdating QR decompositions.
        // The return value is true if successful, false if |y| < |x|.

        typedef typename V::value_type T;
        typedef typename V::real_type RT;

        RT scale = u.maxAbs2Element();
        RT absy = TMV_ABS(y);
        if (absy > scale) scale = absy;

        if (TMV_Underflow(scale)) {
            // Then the situation is hopeless, and we should just zero out
            // the whole vector.
            u.setZero();
            y = T(0);
            beta = RT(0);
            return true; // But still declared successful.
        }

        // Determine |y| = |x| now, but need to find x0 from |x1| and |y|.
        RT invscale = RT(1)/scale;
        RT normsqx1 = u.normSq(invscale);

        // if all of x other than first element are 0, H is identity
        if (normsqx1 == RT(0)) {
            // Set them all to x explicitly in case underflow led to the 0.
            u.setZero();
            beta = RT(0);
            return true;
        }

        // Finish calculating x0 in usual case
        absy *= invscale;
        y *= invscale;
        RT normsqx = absy*absy;
        RT normsqx0 = normsqx - normsqx1;
        if (normsqx0 < RT(0)) return false;

        // y = -|x| (x0/|x0|) 
        // x0 = -|x0| (y/|y|)
        RT absx0 = TMV_SQRT(normsqx0);
        T x0 = -absx0*y/absy;

        // beta = 1 / (|x|^2 + |x||x0|)
        // H = I - beta v vt
        // with u = x - y e0 in first column
        // Renormalize beta,u so that u(0) = 1
        T u0 = x0-y;
        RT normu0 = TMV_NORM(u0);
        beta = normu0 / (normsqx + absy * absx0);
        T invu0 = RT(1)/u0;

        // Sometimes this combination can underflow, so check.
        T scaled_invu0 = invu0 * invscale;
        if ((TMV_REAL(invu0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invu0))) ||
            (TMV_IMAG(invu0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invu0)))) {
            u *= invscale;
            u *= invu0;
        } else {
            u *= scaled_invu0;
        }
                        
        y = x0*scale;
    }

    template <class V, class M0, class Mx, class Vt>
    static inline void HouseholderMultEq(
        const V& u, const typename V::real_type beta,
        BaseVector_Mutable<M0>& m0, BaseMatrix_Rec_Mutable<Mx>& mx,
        BaseVector_Mutable<Vt>& temp)
    {
        // This routine takes 
        // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 ut ) ( m0 ) ]
        // ( mx )      ( mx )   [ ( mx )        ( u )          ( mx ) ]
        // 
        // ( m0 ) -= beta (   m0 + ut mx   )
        // ( mx )         ( u (m0 + ut mx) )

        TMVAssert(u.size() == mx.colsize());
        TMVAssert(m0.size() == mx.rowsize());
        TMVAssert(temp.size() == mx.rowsize());
#ifdef XDEBUG_HOUSE
        //std::cout<<"Start Householder::multEq"<<std::endl;
        //std::cout<<"u = "<<u<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m0 = "<<m0.vec()<<std::endl;
        //std::cout<<"mx = "<<mx.mat()<<std::endl;
        typedef typename V::value_type T;
        typedef typename Mx::value_type T2;
        Matrix<T2> mm(mx.colsize()+1,m0.size());
        mm.row(0) = m0.vec();
        mm.rowRange(1,mm.colsize()) = mx.mat();
#endif
        typedef typename V::real_type RT;

        if (m0.size() > 0 && beta != RT(0)) {
            //std::cout<<"Not trivial.\n";
            temp = u.conjugate() * mx.mat();
            //std::cout<<"temp = "<<temp<<std::endl;
            temp += m0.vec();
            //std::cout<<"temp => "<<temp<<std::endl;
            temp *= -beta;
            //std::cout<<"temp => "<<temp<<std::endl;

            m0.vec() += temp;
            //std::cout<<"m0 => "<<m0<<std::endl;
            mx.mat() += u^temp;
            //std::cout<<"m0 => "<<mx<<std::endl;
        }
#ifdef XDEBUG_HOUSE
        Vector<T> uu(u.size()+1);
        uu(0) = T(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T> H = T(1) - beta*(uu^uu.conjugate());
        Matrix<T2> Hm = H * mm;
        Matrix<T2> Hm2(mx.colsize()+1,m0.size());
        Hm2.row(0) = m0.vec();
        Hm2.rowRange(1,Hm2.colsize()) = mx.mat();
        if (Norm(Hm-Hm2) > 0.001*Norm(Hm)) {
            cerr<<"Householder::multEq\n";
            cerr<<"Input: m = "<<mm<<endl;
            cerr<<"uu = "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Output: m = "<<Hm2<<endl;
            abort();
        }
        //std::cout<<"mx => "<<mx.mat()<<std::endl;
#endif
    }

    template <class V, class M, class V2>
    static inline void HouseholderMultEq(
        V& v, const typename V::real_type beta,
        BaseMatrix_Rec_Mutable<M>& m,
        BaseVector_Mutable<V2>& temp)
    {
        // This routine takes 
        // m <- H m = m - beta v vt m
        // Where v is the full (1 u) vector including a spot for the 1.
        // However, the actual value of that location is not used.

        TMVAssert(v.size() == m.colsize());
        TMVAssert(temp.size() == m.rowsize());
#ifdef XDEBUG_HOUSE
        //std::cout<<"Start Householder::multEq"<<std::endl;
        //std::cout<<"v = "<<v<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m = "<<m.mat()<<std::endl;
        typedef typename V::value_type T;
        typedef typename Mx::value_type T2;
        Matrix<T2> mm = m;
#endif
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        typedef typename M::transpose_type Mt;
        typedef typename V::const_conjugate_type Vc;
        const int cs = Sizes<M::_colsize,V::_size>::size;
        const int rs = Sizes<M::_rowsize,V2::_size>::size;

        if (m.rowsize() > 0 && beta != RT(0)) {
            T v0 = v.cref(0);
            v.ref(0) = 1.;

            //temp = -beta * v.conjugate() * m.mat();
            const Scaling<0,RT> x(-beta);
            Mt mt = m.mat().transpose();
            Vc vc = v.conjugate();
            MultMV_Helper<-4,rs,cs,false,0,RT,Mt,Vc,V2>::call(x,mt,vc,temp.vec());

            //m.mat() += v^temp;
            const Scaling<1,RT> one;
            Rank1VVM_Helper<-4,cs,rs,true,1,RT,V,V2,M>::call(one,v.vec(),temp.vec(),m.mat());

            v.ref(0) = v0;
        }
#ifdef XDEBUG_HOUSE
        Matrix<T> H = T(1) - beta*(v^v.conjugate());
        Matrix<T2> Hm = H * mm;
        if (Norm(Hm-m) > 0.001*Norm(Hm)) {
            cerr<<"Householder::multEq\n";
            cerr<<"Input: m = "<<mm<<endl;
            cerr<<"v = "<<v<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Output: m = "<<m<<endl;
            abort();
        }
        //std::cout<<"mx => "<<mx.mat()<<std::endl;
#endif
    }

    template <class V, class Vx>
    static inline void HouseholderMultEq(
        const V& u, const typename V::real_type beta,
        typename Vx::reference v0, BaseVector_Mutable<Vx>& vx)
    {
        // This routine takes 
        // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 ut ) ( m0 ) ]
        // ( mx )      ( mx )   [ ( mx )        ( u )          ( mx ) ]
        // 
        // ( m0 ) -= beta (   m0 + ut mx   )
        // ( mx )         ( u (m0 + ut mx) )

        TMVAssert(u.size() == vx.size());
        typedef typename V::value_type T;
#ifdef XDEBUG_HOUSE
        //std::cout<<"Start Householder::multEq"<<std::endl;
        //std::cout<<"u = "<<u<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"v0 = "<<v0<<std::endl;
        //std::cout<<"vx = "<<vx.vec()<<std::endl;
        typedef typename Vx::value_type T2;
        Vector<T2> vv(vx.size()+1);
        vv.ref(0) = v0;
        vv.subVector(1,vv.size()) = vx.vec();
#endif
        typedef typename V::real_type RT;

        if (beta != RT(0)) {
            //std::cout<<"Not trivial.\n";
            T temp = -beta * (v0 + u.conjugate() * vx.vec());
            //std::cout<<"temp = "<<temp<<std::endl;

            v0 += temp;
            //std::cout<<"m0 => "<<m0<<std::endl;
            vx.vec() += u*temp;
            //std::cout<<"m0 => "<<mx<<std::endl;
        }
#ifdef XDEBUG_HOUSE
        Vector<T> uu(u.size()+1);
        uu(0) = T(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T> H = T(1) - beta*(uu^uu.conjugate());
        Matrix<T2> Hv = H * vv;
        Matrix<T2> Hv2(vx.size()+1);
        Hv2(0) = v0;
        Hv2.subVector(1,Hv2.size()) = vx.vec();
        if (Norm(Hv-Hv2) > 0.001*Norm(Hv)) {
            cerr<<"Householder::multEq\n";
            cerr<<"Input: v = "<<vv<<endl;
            cerr<<"uu = "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hv = "<<Hv<<endl;
            cerr<<"Output: v = "<<Hv2<<endl;
            abort();
        }
        //std::cout<<"vx => "<<vx.vec()<<std::endl;
#endif
    }

#if 0
    template <class V, class M1, class M2>
    static inline void HouseholderMult(
        const V& u, const typename V::real_type beta,
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) 
    {
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(u.size()+1 == m1.colsize());
        TMVAssert(u.size()+1 == m2.colsize());
#ifdef XDEBUG_HOUSE
        typedef typename V::value_type T;
        typedef typename M2::value_type T2;
        //std::cout<<"Start Householder::mult"<<std::endl;
        //std::cout<<"u = "<<u<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m1 = "<<m1.mat()<<std::endl;
        //std::cout<<"m2 = "<<m2.mat()<<std::endl;
#endif
        typedef typename V::real_type RT;

        if (m1.rowsize() > 0) {
            if (beta != RT(0)) {
                const int N = m1.colsize();
                typename M1::const_row_type m0 = m1.mat().row(0);
                typename M1::const_rowrange_type mx = m1.mat().rowRange(1,N);

                typename M2::row_type m2_0 = m2.mat().row(0);
                typename M2::rowrange_type m2_x = m2.mat().rowRange(1,N);

                m2_0 = u.conjugate() * mx;
                m2_0 += m0;
                m2_0 *= -beta;

                m2_x = mx;
                m2_x += u^m2_0;
                m2_0 += m0;
            } else {
                m2.mat() = m1.mat();
            }
        }
#ifdef XDEBUG_HOUSE
        Vector<T> uu(u.size()+1);
        uu(0) = T(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T> H = T(1) - beta*(uu^uu.conjugate());
        Matrix<T2> Hm = H * m1.mat();
        if (Norm(Hm-m2.mat()) > 0.001*Norm(Hm)) {
            cerr<<"Householder::mult\n";
            cerr<<"Input: m = "<<m1.mat()<<endl;
            cerr<<"uu = "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Output: m = "<<m2.mat()<<endl;
            abort();
        }
        //std::cout<<"mx => "<<mx<<std::endl;
#endif
    }

    template <class V, class V1, class V2>
    static inline void HouseholderMult(
        const V& u, const typename V::real_type beta,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V::real_type RT;
        TMVAssert(u.size()+1 == v1.size());
        TMVAssert(u.size()+1 == v2.size());
        if (beta != RT(0)) {
            const int N = v1.size();
            typename V2::value_type temp = u.conjugate() * v1.subVector(1,N);
            temp += v1.cref(0);
            temp *= beta;
            v2.ref(0) = v1.cref(0) - temp;
            v2.subVector(1,N) = v1.subVector(1,N) - temp*u;
        }
    }
#endif

    template <class V>
    static inline void HouseholderUnpack(
        V& u, const typename V::real_type beta, typename V::reference u0)
    {
        typedef typename V::value_type T;
        // This routine takes u <- H*e0

        if (beta == T(0)) {
            // then all but u(0) is already 0
            u0 = T(1);
        } else {
            // u <- (I-beta (1 u) (1 u)t) e0 = e0 - beta (1 u)
            u0 = T(1)-beta;
            u *= -beta;
        }
    }

    template <class M1, class M2, class RT>
    void BlockHouseholderAugment(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, RT beta)
    {
        // All but the last columns of the input matrices, Y,Z are such that
        // I - Y'Z'Y't is a Block Householder matrix (the product of several
        // individual Householder matrices).
        // The last column of Y has the vector for the next Householder matrix.
        //
        // Z is output such that:
        // I - Y Z Yt = (I - Y' Z' Y't) (I - beta v vt)
        //
        // Y Z Yt = Y' Z' Y't + beta v vt - beta Y' Z' Y't v vt
        //
        // Blocking Y as [ Y' v ] and Z as [ Z'  z   ]
        //                                 [ 0  beta ]
        // Y Z Yt = [ Y' v ] [ Z' Y't + z vt ]
        //                   [    beta vt    ]
        //        = Y' Z' Y't + beta v vt + Y' z vt
        //
        // Comparing these equations, we find:
        // z = -beta Z' Y't v
        //
        // Remember that the first element of v is not stored, but rather
        // is assumed to be 1.

        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        int M = Y.colsize();
        int N = Y.rowsize()-1; // # of cols already computed
#ifdef XDEBUG_HOUSE
        typedef typename M1::value_type T;
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z);
        Matrix<T> H0 = T(1) - 
            Y0.colRange(0,N)*Z.subTriMatrix(0,N)*Y0.colRange(0,N).adjoint();
        Matrix<T> H1(Y.colsize(),Y.colsize());
        H1.setToIdentity();
        H1.subMatrix(N,M,N,M) -= beta * (Y0.col(N,N,M)^Y0.col(N,N,M).conjugate());
        Matrix<T> H2 = H0*H1;
#endif

        if (beta == RT(0)) {
            Z.col(N,0,N+1).setZero();
        } else if (N == 0) {
            Z.ref(0,0) = beta;
        } else {
            typename M1::const_col_sub_type v = Y.col(N,N+1,M);
            typename M2::col_sub_type z = Z.col(N,0,N);
            z = Y.subMatrix(N+1,M,0,N).adjoint()*v;
            z += Y.row(N,0,N).conjugate();
            z = -beta * Z.subTriMatrix(0,N) * z;
            Z.ref(N,N) = beta;
        }
#ifdef XDEBUG_HOUSE
        Matrix<T> YY(Y);
        YY.upperTri().setToIdentity();
        Matrix<T> HH = T(1) - YY * Z * YY.adjoint();
        if (Norm(HH-H2) > 0.001*Norm(H2)) {
            cerr<<"BlockHouseholderAugment\n";
            cerr<<"Input: Y = "<<Y0<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"Z = "<<Z0<<endl;
            cerr<<"H0 = "<<H0<<endl;
            cerr<<"H1 = "<<H1<<endl;
            cerr<<"H0*H1 = "<<H2<<endl;
            cerr<<"Output: Y = "<<Y<<endl;
            cerr<<"Z = "<<Z<<endl;
            cerr<<"H = "<<HH<<endl;
            abort();
        }
#endif
    }

    template <class M1, class M2, class V>
    void BlockHouseholderMakeZ(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, 
        const BaseVector_Calc<V>& beta)
    {
        // This routine calculates the Z component of the BlockHouseholder
        // formulation for Q.  Y contains the v's for the Householder matrices,
        // and beta contains the beta's.  
        //
        // The output BlockHouseholder matrix I-YZYt is the product of 
        // the adjoints, H0t H1t ... HNt, since this is the product which 
        // we actually use for most calculations.  
        // If you want the product of the H's, input beta.conjugate instead.
        //
        // Note that the Y matrix is really just the unit lower trapezoidal 
        // component of the input Y.
        TMVAssert(Y.rowsize() == Z.rowsize());
        TMVAssert(Y.rowsize() == beta.size());
        TMVAssert(Y.colsize() >= Y.rowsize());

        typedef typename M1::value_type T;
        typedef typename M1::real_type RT;
#ifdef XDEBUG_HOUSE
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z);
        Vector<T> beta0 = beta;
        Matrix<T> Htot(Y.colsize(),Y.colsize());
        Htot.setToIdentity();
        for(int i=0;i<int(Y.rowsize());i++) {
            Matrix<T> H(Y.colsize(),Y.colsize());
            H.setToIdentity();
            H.subMatrix(i,Y.colsize(),i,Y.colsize()) -= beta(i) * 
                (Y0.col(i,i,Y0.colsize()) ^ 
                 Y0.col(i,i,Y0.colsize()).conjugate());
            Htot *= H.adjoint();
        }
#endif

        const int M = Y.colsize();
        const int N = Y.rowsize();

        if (N > 2) {
            int j1 = (N+1)/2;
            typename M1::const_colrange_type Y1 = Y.colRange(0,j1);
            typename M2::subtrimatrix_type Z1 = Z.subTriMatrix(0,j1);
            typename V::const_subvector_type beta1 = beta.subVector(0,j1);
            BlockHouseholderMakeZ(Y1,Z1,beta1);

            typename M1::const_submatrix_type Y2 = Y.subMatrix(j1,M,j1,N);
            typename M2::subtrimatrix_type Z2 = Z.subTriMatrix(j1,N);
            typename V::const_subvector_type beta2 = beta.subVector(j1,N);
            BlockHouseholderMakeZ(Y2,Z2,beta2);

            // (I-Y1 Z1 Y1t)(I-Y2 Z2 Y2t) = 
            // I - Y1 Z1 Y1t - Y2 Z2 Y2t + Y1 Z1 Y1t Y2 Z2 Y2t
            // Y = [ Y1 Y2 ]
            // Z = [ Z1 Z3 ]
            //     [ 0  Z2 ]
            // Y Z Yt = [ Y1 Y2 ] [ Z1 Y1t + Z3 Y2t ]
            //                    [      Z2 Y2t     ]
            //        = Y1 Z1 Y1t + Y2 Z2 Y2t + Y1 Z3 Y2t
            // So, Z3 = -Z1 Y1t Y2 Z2
            // Remember that the Y's are Unit Lower Trapezoidal, so do the 
            // rectangle and triangle parts separately.
            //
            typename M2::submatrix_type Z3 = Z.subMatrix(0,j1,j1,N);
            Z3 = Y1.rowRange(j1,N).adjoint() *
                Y.subMatrix(j1,N,j1,N).unitLowerTri();
            Z3 += Y1.rowRange(N,M).adjoint() * Y.subMatrix(N,M,j1,N);
            Z3 = -Z1*Z3;
            Z3 *= Z2;
        } else if (N==2) {
            // Y = ( Y00  0  )   Z = ( Z00  Z01 )
            //     ( Yx0 Yx1 )       (  0   Z11 )
            //
            // I - Y Z Yt = I - ( Z00 Y00   Z01 Y00           ) ( Y00* Yx0t )
            //                  ( Z00 Yx0   Z01 Yx0 + Z11 Yx1 ) (  0   Yx1t )
            // = ( 1 - Z00 Y00 Y00*   -Z00 Y00 Yx0t - Z01 Y00 Yx1t      )
            //   ( -Z00 Y00* Yx0      I - Z00 Yx0 Yx0t - Z01 Yx0 Yx1t 
            //                               - Z11 Yx1 Yx1t             )
            //
            //
            // (I - b0* v0 v0t)*( 1            0     )
            //                  ( 0   I - b1* v1 v1t )
            // [ Let v0 = ( v00 v0x ) ]
            //  = ( 1 - b0* v00 v00*   -b0* v00 v0xt    ) ( 1        0        )
            //    ( -b0* v00* v0x      I - b0* v0x v0xt ) ( 0  I - b1* v1 v1t )
            //  = ( 1 - b0* v00 v00*   -b0* v00 v0xt + b0* b1* v00 v0xt v1 v1t )
            //    ( -b0* v00* v0x      I - b0* v0x v0xt - b1* v1 v1t 
            //                             + b0* b1* v0x v0xt v1 v1t           )
            //
            // Matching the two results, we find:
            //
            // Y00 = v00
            // Yx0 = v0x
            // Yx1 = v1
            // Z00 = b0*
            // Z11 = b1*
            // Z01 = -b0* b1* v0xt v1
            const RT b0 = beta.cref(0);
            const RT b1 = beta.cref(1);
            Z.ref(0,0) = b0;
            Z.ref(1,1) = b1;
            T temp = Y.col(0,2,M).conjugate()*Y.col(1,2,M);
            temp += TMV_CONJ(Y.cref(1,0));
            Z.ref(0,1) = -b0*b1*temp;
        } else { // N == 1
            // I - Y Z Yt = I - beta* v vt
            // Therefore:
            // Y.col(0) = v
            // Z(0,0) = beta*
            Z.ref(0,0) = beta.cref(0);
        }
#ifdef XDEBUG_HOUSE
        Matrix<T> Y2(Y);
        Y2.upperTri().setToIdentity();
        Matrix<T> Hnet = T(1) - Y2*Z*Y2.adjoint();
        if (Norm(Htot-Hnet) > 0.001*Norm(Htot)) {
            cerr<<"BlockHouseholderMakeZ\n";
            cerr<<"Input: Y = "<<Y0<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<Htot<<endl;
            cerr<<"Output: Y = "<<Y2<<endl;
            cerr<<"Z = "<<Z<<endl;
            cerr<<"H = "<<Hnet<<endl;
            cerr<<"Norm(H-Htot) = "<<Norm(Htot-Hnet)<<endl;
            abort();
        }
#endif
    }

    template <class M1, class M2, class M3, class M4>
    void BlockHouseholderLMult(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        // The input Y,Z are such that (I - YZYt) is a Block Householder matrix.
        // The upper square portion of Y is taken to be unit lower triangular.
        // ie. the diagonal and upper triangular portion are not referenced.
        // The routine then finds m <- (I - YZYt) m

        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());

#ifdef XDEBUG_HOUSE
        typedef typename M1::value_type T;
        typedef typename M3::value_type T3;
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z);
        Matrix<T3> m0(m);
        Matrix<T> H = T(1) - Y0*Z0*Y0.adjoint();
        Matrix<T3> Hm = H*m0;
#endif

        int M = Y.colsize();
        int N = Y.rowsize();
        const int cs = Sizes<M1::_colsize,M3::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_size>::size;
        const int xs = M3::_rowsize;
        const int cs1 = IntTraits2<cs,rs>::diff;
        const int Si1 = M1::_stepi;
        const int Sj1 = M1::_stepj;
        const int C1 = M1::_conj;
        const int Si3 = M3::_stepi;
        const int Sj3 = M3::_stepj;
        const int C3 = M3::_conj;
        typedef typename M1::value_type T1;
        typedef typename M3::value_type T3;
        typedef typename MViewHelper<T1,Rec,rs,rs,Si1,Sj1,C1>::ctype M1a;
        typedef typename MViewHelper<T1,Rec,cs1,rs,Si1,Sj1,C1>::ctype M1b;
        typedef typename MViewHelper<T3,Rec,rs,xs,Si3,Sj3,C3>::type M3a;
        typedef typename MViewHelper<T3,Rec,cs1,xs,Si3,Sj3,C3>::type M3b;

        M1a Ya = Y.rowRange(0,N);
        M1b Yb = Y.rowRange(N,M);
        M3a ma = m.rowRange(0,N);
        M3b mb = m.rowRange(N,M);
 
        // temp = ZYtm
        temp = Ya.unitLowerTri().adjoint() * ma;
        temp += Yb.adjoint() * mb;
        temp = -Z * temp;
        ma += Ya.unitLowerTri() * temp;
        mb += Yb * temp;

#ifdef XDEBUG_HOUSE
        if (Norm(Hm-m) > 0.001*Norm(m0)*Norm(H)) {
            cerr<<"BlockHouseholderLMult\n";
            cerr<<"Input: Y = "<<Y0<<endl;
            cerr<<"Z = "<<Z0<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Ya = "<<Ya<<endl;
            cerr<<"Yb = "<<Yb<<endl;
            cerr<<"temp = "<<temp<<endl;
            cerr<<"Output: m = "<<m<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            abort();
        }
#endif
    }

    template <class M1, class M2, class M3, class M4>
    void BlockHouseholderLDiv(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        // The routine finds m <- (I - YZYt)^-1 m
        // = (I - YZtYt) m

        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());
#ifdef XDEBUG_HOUSE
        typedef typename M1::value_type T;
        typedef typename M3::value_type T3;
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z.adjoint());
        Matrix<T3> m0(m);
        Matrix<T> Hinv = T(1) - Y0*Z0*Y0.adjoint();
        Matrix<T3> Hm = Hinv*m0;
#endif

        int M = Y.colsize();
        int N = Y.rowsize();
        const int cs = Sizes<M1::_colsize,M3::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_size>::size;
        const int xs = M3::_rowsize;
        const int cs1 = IntTraits2<cs,rs>::diff;
        const int Si1 = M1::_stepi;
        const int Sj1 = M1::_stepj;
        const int C1 = M1::_conj;
        const int Si3 = M3::_stepi;
        const int Sj3 = M3::_stepj;
        const int C3 = M3::_conj;
        typedef typename M1::value_type T1;
        typedef typename M3::value_type T3;
        typedef typename MViewHelper<T1,Rec,rs,rs,Si1,Sj1,C1>::ctype M1a;
        typedef typename MViewHelper<T1,Rec,cs1,rs,Si1,Sj1,C1>::ctype M1b;
        typedef typename MViewHelper<T3,Rec,rs,xs,Si3,Sj3,C3>::type M3a;
        typedef typename MViewHelper<T3,Rec,cs1,xs,Si3,Sj3,C3>::type M3b;

        M1a Ya = Y.rowRange(0,N);
        M1b Yb = Y.rowRange(N,M);
        M3a ma = m.rowRange(0,N);
        M3b mb = m.rowRange(N,M);

        // temp = ZtYtm
        temp = Ya.unitLowerTri().adjoint() * ma;
        temp += Yb.adjoint() * mb;
        temp = -Z.adjoint() * temp;
        ma += Ya.unitLowerTri() * temp;
        mb += Yb * temp;

#ifdef XDEBUG_HOUSE
        if (Norm(Hm-m) > 0.001*Norm(m0)*Norm(Hinv)) {
            cerr<<"BlockHouseholderLDiv\n";
            cerr<<"Input: Y = "<<Y0<<endl;
            cerr<<"Z = "<<Z0<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"Hinv = "<<Hinv<<endl;
            cerr<<"Output: m = "<<m<<endl;
            abort();
        }
#endif
    }

    template <class M1, class M2, class M3>
    void BlockHouseholderUnpack(
        BaseMatrix_Rec_Mutable<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Tri_Mutable<M3>& temp)
    {
        // The BlockHouseholder matrix H is defined by Y,Z.
        // Y is unpacked in place.
        TMVAssert(Y.colsize() > 0);
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(Y.rowsize() == Z.size());
        TMVStaticAssert(M3::_upper);
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(!M3::_unknowndiag);
        TMVAssert(temp.size() == Y.rowsize());

        int M = Y.colsize();
        int N = Y.rowsize();
        const int cs = M1::_colsize;
        const int rs = Sizes<M1::_rowsize,M2::_size>::size;
        const int cs1 = IntTraits2<cs,rs>::diff;
        const int Si1 = M1::_stepi;
        const int Sj1 = M1::_stepj;
        const int C1 = M1::_conj;
        typedef typename M1::value_type T;
        typedef typename MViewHelper<T,Rec,rs,rs,Si1,Sj1,C1>::type M1a;
        typedef typename MViewHelper<T,Rec,cs1,rs,Si1,Sj1,C1>::type M1b;

        M1a Ya = Y.rowRange(0,N);
        M1b Yb = Y.rowRange(N,M);
 
        // Make Y equal to:
        // Ht [ I ] = (I - YZYt) [ I ]
        //    [ 0 ]              [ 0 ]
        temp = -Z * Ya.unitLowerTri().adjoint();
        Yb *= temp;
        Ya = Ya.unitLowerTri() * temp;
        Ya.diag().addToAll(T(1));
    }


} // namespace tmv

#endif
