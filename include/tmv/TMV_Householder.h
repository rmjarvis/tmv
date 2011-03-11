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

#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

    template <class V> class Householder;

    template <class V>
    struct Traits<Householder<V> >
    {
        typedef typename V::value_type value_type;
        typedef typename V::real_type real_type;
        typedef typename V::complex_type complex_type;
        enum { isreal = V::isreal };
        enum { iscomplex = V::iscomplex };

        typedef Householder<V> type;
        typedef Matrix<value_type> copy_type;
        typedef copy_type calc_type;
        typedef type eval_type;
        typedef Householder<V> inverse_type;

        enum { _colsize = IntTraits<V::_size>::Sp1 };
        enum { _rowsize = _colsize };
        enum { _shape = SquareRec };
        enum { _fort = V::_fort };
        enum { _calc = false };
    };

    // The template parameter is the vector type used to store u_1..N.
    // Note: V itself might be const, in which case many of the 
    // non-const methods won't work.
    // Also, the normal usage would have V be a VectorView, so the 
    // V copies are trivial.  However, you can have the Householder matrix
    // copy all the elements into its own storage by having V be 
    // a deep copying class like Vector.
    template <class V>
    class Householder : public BaseMatrix<Householder<V> >
    {
    public :
        typedef typename V::value_type T;
        typedef typename V::real_type RT;

        typedef Householder<V> type;
        typedef const Householder<typename V::const_conjugate_type> 
            const_conjugate_type;
        typedef const_conjugate_type const_transpose_type;
        typedef const type& const_adjoint_type;
        typedef const type& const_inverse_type;

        //
        // Constructors
        //
        
        // The first constructor performs the reflection on the given
        // vector and ends up with the Householder vector in the location
        // of x.
        //
        // Sometimes (e.g. for QR_Update), the first element of the 
        // reflected vector is separated from the rest of the vector.
        // So we pass in the first element and the rest as separate parameters.
        //
        // x is converted to u on output.
        Householder(T& x0, const V& x) : u(x)
        { doReflection(x0); }

        // In case the first element is some other reference type:
        template <class RefT>
        Householder(RefT x0, const V& x) : u(x)
        { T temp=x0; doReflection(temp); x0 = temp; }

        // This is a convenience constructor that takes x0 and x as
        // a single vector:
        template <class Vfull>
        Householder(Vfull& xfull) : 
            u(xfull.subVector(1,xfull.size()))
        { T temp=xfull.cref(0); doReflection(temp); xfull.ref(0) = temp; }

        // The next constructor assumes that the reflection is already
        // done, so u and beta are passed in, and can then be used.
        Householder(const V& _u, const RT& _beta) : u(_u), beta(_beta) {}

        // The next constructor reverses the roles of x0 and y.  
        // The extra bool parameter would normally be set to true, 
        // meaning unreflect = true, but in fact it is ignored, as the
        // different functionality is actually set by the signature. 
        // In this case, the x0 parameter is actually y on input and 
        // x0 on output (rather than the usual x0 on input and y on output).
        // Again, x is the unreflected vector on input, and is
        // converted to u on output.  
        // N.B. This version is used for QR downdating.
        Householder(T& x0, const V& x, bool) : u(x)
        { doUnReflection(x0); }
        template <class RefT>
        Householder(RefT x0, const V& x, bool) : u(x)
        { T temp=x0; doReflection(temp); x0 = temp; }
        template <class Vfull>
        Householder(Vfull& xfull, bool) : u(xfull.subVector(1,xfull.size()))
        { T temp=xfull.cref(0); doUnReflection(temp); xfull.ref(0) = temp; }

        //
        // Special Accessors: 
        //
        
        // beta
        RT getBeta() const { return beta; }

        // u except for the first element (which is 1)
        const V& getU() const { return u; }

        // determinant
        RT det() const 
        { return beta == RT(0) ? RT(1) : RT(-1); }


        //
        // Perform Householder multiplication:
        // 
        
        // This first one allows the first row and the rest in different
        // places:
        template <class M0, class Mx> 
        void multEq(M0& m0, Mx& mx) const;

        // Here m0 and mx are a single matrix:
        template <class M>
        void multEq(M& m) const
        {
            typename M::row_type m0 = m.row(0);
            typename M::rowrange_type mx = m.rowRange(1,m.colsize());
            multEq(m0,mx);
        }

        // If result is new storage, then there is no need for a temporary
        // vector during the calculation.  The temporary isn't very
        // much, but still worth having the separate function.
        template <class M1, class M2>
        void mult(const M1& m1, M2& m2) const;

        // Repeat for vectors:
        template <class Vx> 
        void multEqV(typename Vx::reference v0, Vx& vx) const;
        template <class V2>
        void multEqV(V2& v) const
        {
            typename Vx::reference v0 = v.ref(0);
            typename Vx::subvector_type vx = v.rowRange(1,v.size());
            multEqV(v0,vx);
        }
        template <class V1, class V2>
        void multV(const V1& v1, V2& v2) const;

 
        
        // 
        // Unpack the Householder vector in place
        //

        // Convert the vector used for the Householder matrix along
        // with the v0 location into H*e0.
        // Note: after doing this function, the Householder matrix
        // shouldn't be used further, since it isn't correct.
        void unPack(T& v0);

        // In case v0 is some other reference type
        template <class RefT>
        void unPack(RefT v0)
        { T temp = v0; unPack(temp); v0 = temp; }


        //
        // Now all the stuff to make it work as a BaseMatrix:
        // (Including overloads of some functions that are trivial to 
        //  calculate.)
        //

        size_t colsize() const { return u.size()+1; }
        size_t rowsize() const { return u.size()+1; }
        T cref(size_t i, size_t j) const 
        {
            // H = [ I - beta ( 1 ) ( 1 ut ) ]
            //     [          ( u )          ]
            //   = [ 1-beta      -beta ut   ]
            //     [ -beta u    I-beta u ut ]
            return 
                i==j ? 
                ( i==0 ? T(1)-beta : T(1)-beta*TMV_NORM(u.cref(i-1)) ) : 
                ( i==0 ? -beta*TMV_CONJ(u.cref(j-1)) : 
                  j==0 ? -beta*u.cref(j-1) :
                  -beta*u.cref(i-1)*TMV_CONJ(u.cref(j-1)) );
        }

        RT trace() const
        { return RT(colsize()) - beta*(RT(1)+u.NormSq()); }

        template <class M2>
        void assignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            const int N = u.size()+1;
            m2.mat().ref(0,0) = RT(1)-beta;
            m2.mat().col(0,1,N) = -beta*u;
            m2.mat().row(0,1,N) = -beta*u.conjugate();
            m2.mat().subMatrix(1,1,N,N) = -beta*u^u.conjugate();
            m2.mat().subMatrix(1,1,N,N) += RT(1);
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            // TODO use NoAlias versions:
            assignTo(m2);
        }

        const_conjugate_type conjugate() const
        { return const_conjugate_type(u,beta); }
        const_transpose_type transpose() const
        { return const_transpose_type(u,beta); }
        const_adjoint_type conjugate() const
        { return *this; }
        const_inverse_type inverse() const
        { return *this; }

    private :
        void doReflection(T& x0);
        void doUnReflection(T& x0);

        V u; 
        RT beta;
    };

    template <class V>
    static void Householder<V>::doReflection(typename Householder<V>::T& x0)
    {
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

        RT scale = x.maxAbs2Element();
        RT absx0 = TMV_Abs(x0);
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
        if (normsqx == RT(0) && TMV_IMAG(x0) == RT(0)) {
            // Set them all to x explicitly in case underflow led to the 0.
            u.setZero();
            x0 = T(0);
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
        RT absx0 = TMV_ABS(x0);
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
        if ((TMV_REAL(invv0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invv0))) ||
            (TMV_IMAG(invv0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invv0)))) {
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
    static void Householder<V>::doUnReflection(typename Householder<V>::T& y)
    {
        // This is similar, except that the roles of y and x0 are swapped.
        // That is, the final rotated vector y e0 is presumed known, as is 
        // the bulk of the input vector.  
        // The unreflected value of the first element in the vector is returned 
        // as y.
        // The rest of the vector is transformed into the Householder vector.
        // This is used for downdating QR decompositions.
        // The return value is true if successful, false if |y| < |x|.

        TMVAssert(u.size() > 0);
#ifdef XDEBUG_HOUSE
        Vector<T> xx(u.size()+1);
        xx.subVector(1,xx.size()) = u;
        T yy = y;
#endif
        // abs(y) = |x|, so |x0|^2 is |y|^2 - |x_1..N|^2
        // normsqx1 = |x_1..N|^2
        RT normsqx1 = NormSq(u);

        // if all of x other than first element are 0, H is identity
        if (normsqx1 == RT(0)) {
            // Set them all to x explicitly in case underflow led to the 0.
            u.setZero();
            beta = RT(0);
            return;
        }

        // Finish calculating x0 in usual case
        RT normsqx = TMV_NORM(y);
        RT normx = TMV_SQRT(normsqx);

        if (normsqx * TMV_Epsilon<T>() == RT(0)) {
            // Then we need to rescale, since underflow will cause 
            // rounding errors
            const RT eps = TMV_Epsilon<T>();
            // Epsilon is a pure power of 2, so no rounding errors from 
            // rescaling.
            u /= eps;
            y /= eps;
            doUnReflection(y);
            y *= eps;
        } else if (RT(1)/normsqx == RT(0)) {
            // Then we have overflow, so we need to rescale:
            const RT eps = TMV_Epsilon<T>();
            RT scale = eps; normx *= eps;
            while (normx > RT(1)) { scale *= eps; normx *= eps; }
            u *= scale;
            y *= scale;
            doUnReflection(y);
            y /= scale;
        } else {
            // The usual case with no rescaling required:
            RT normsqx0 = normsqx - normsqx1;
            if (normsqx0 < RT(0)) {
                // TODO: Do something more with this?
                throw NonPosDef();
            }

            // y = -|x| (x0/|x0|) 
            // x0 = -|x0| (y/|x|)
            RT absx0 = TMV_SQRT(normsqx0);
            T x0 = absx0 == RT(0) ? T(0) : -absx0*y/normx;

            // beta = 1 / (|x|^2 + |x||x0|)
            // H = I - beta v vt
            // with u = x - y e0 in first column
            // Renormalize beta,u so that u(0) = 1
            T u0 = x0-y;
            RT normu0 = TMV_NORM(u0);
            beta = normu0 / (normsqx + normx * absx0);
            u /= u0;
            y = x0;
        }
#ifdef XDEBUG_HOUSE
        xx(0) = y;
        Vector<T> uu(xx.size());
        uu(0) = T(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T> H = T(1)-beta*(uu^uu.conjugate());
        // Check the following:
        // H * xx = (yy,0,0.0...)
        // |xx| = |yy|
        Vector<T> Hxx = H * xx;
        if (TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(yy)) > 0.0001*TMV_ABS(yy) ||
            TMV_ABS(Norm(xx)-TMV_ABS(yy)) > 0.0001*TMV_ABS(yy) ||
            Norm(Hxx.subVector(1,Hxx.size())) > 0.0001*TMV_ABS(yy)) {
            cerr<<"Householder UnReflect:\n";
            cerr<<"Input: x = "<<xx<<endl;
            cerr<<"       y = "<<yy<<endl;
            cerr<<"x0 => "<<y<<endl;
            cerr<<"Norm(x) = "<<Norm(xx)<<endl;
            cerr<<"Output: uu= "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"xx = "<<xx<<endl;
            cerr<<"Hxx = "<<Hxx<<endl;
            cerr<<"abs(abs(hxx(0))-abs(x0)) = "<<
                TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(yy))<<endl;
            cerr<<"abs(Norm(xx)-abs(x0)) = "<<
                TMV_ABS(Norm(xx))-TMV_ABS(yy)<<endl;
            cerr<<"Norm(Hxx(1,N)) = "<<Norm(Hxx.subVector(1,Hxx.size()))<<endl;
            abort();
        }
#endif
    }

    template <class V> template <class M0, class Mx>
    static void Householder<V>::multEq(M0& m0, Mx& mx)
    {
        // This routine takes 
        // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 ut ) ( m0 ) ]
        // ( mx )      ( mx )   [ ( mx )        ( u )          ( mx ) ]
        // 
        // ( m0 ) -= beta (   m0 + ut mx   )
        // ( mx )         ( u (m0 + ut mx) )

        TMVAssert(u.size() == mx.colsize());
        TMVAssert(m0.size() == mx.rowsize());
#ifdef XDEBUG
        //std::cout<<"Start Householder::multEq"<<std::endl;
        //std::cout<<"u = "<<u<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m0 = "<<m0<<std::endl;
        //std::cout<<"mx = "<<mx<<std::endl;
        Matrix<T2> mm(mx.colsize()+1,m0.size());
        mm.row(0) = m0;
        mm.rowRange(1,mm.colsize()) = mx;
#endif

        if (m0.size() > 0 && beta != RT(0)) {
            typedef typename M0::copy_type Vtemp;
            Vtemp temp = u.conjugate() * mx;
            temp += m0;
            temp *= beta;

            m0 -= temp;
            mx -= u^temp;
        }
#ifdef XDEBUG
        Vector<T1> uu(u.size()+1);
        uu(0) = T1(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T1> H = T1(1) - beta*(uu^uu.conjugate());
        Matrix<T2> Hm = H * mm;
        Matrix<T2> Hm2(mx.colsize()+1,m0.size());
        Hm2.row(0) = m0;
        Hm2.rowRange(1,Hm2.colsize()) = mx;
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
        //std::cout<<"mx => "<<mx<<std::endl;
#endif
    }

    template <class V> template <class M1, class M2>
    static void Householder<V>::mult(const M1& m1, M2& m2)
    {
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(u.size()+1 == m1.colsize());
        TMVAssert(u.size()+1 == m2.colsize());
#ifdef XDEBUG
        //std::cout<<"Start Householder::mult"<<std::endl;
        //std::cout<<"u = "<<u<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        //std::cout<<"m1 = "<<m1<<std::endl;
        //std::cout<<"m2 = "<<m2<<std::endl;
#endif

        if (m1.rowsize() > 0) {
            if (beta != RT(0)) {
                const int N = m1.colsize();
                typename M1::const_row_type m0 = m1.row(0);
                typename M1::const_rowrange_type mx = m1.rowRange(1,N);

                typename M2::row_type m2_0 = u.conjugate() * mx;
                m2_0 += m0;
                m2_0 *= -beta;

                typename M2::rowrange_type m2_x = m2.rowRange(1,N);
                m2_x = mx;
                m2_x += u^m2_0;
                m2_0 += m0;
            } else {
                m2 = m1;
            }
        }
#ifdef XDEBUG
        Vector<T1> uu(u.size()+1);
        uu(0) = T1(1);
        uu.subVector(1,uu.size()) = u;
        Matrix<T1> H = T1(1) - beta*(uu^uu.conjugate());
        Matrix<T2> Hm = H * m1;
        if (Norm(Hm-m2) > 0.001*Norm(Hm)) {
            cerr<<"Householder::mult\n";
            cerr<<"Input: m = "<<m1<<endl;
            cerr<<"uu = "<<uu<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Output: m = "<<m2<<endl;
            abort();
        }
        //std::cout<<"mx => "<<mx<<std::endl;
#endif
    }

    template <class V> template <class Vx>
    static void Householder<V>::multEqV(typename Vx::reference v0, Vx& vx)
    {
        // ( v0 ) -= beta (   v0 + ut vx   )
        // ( vx )         ( u (v0 + ut vx) )
        TMVAssert(u.size() == vx.size());
        if (beta != RT(0)) {
            typename Vx::value_type temp = u.conjugate() * vx;
            temp += v0;
            temp *= beta;
            v0 -= temp;
            vx -= temp*u;
        }
    }

    template <class V> template <class V1, class V2>
    static void Householder<V>::multV(const V1& v1, V2& v2)
    {
        TMVAssert(u.size()+1 == v1.size());
        TMVAssert(u.size()+1 == v2.size());
        if (beta != RT(0)) {
            const int N = m1.colsize();
            typename V2::value_type temp = u.conjugate() * v1.subVector(1,N);
            temp += v1.cref(0);
            temp *= -beta;
            v2.ref(0) = v2.cref(0) - temp;
            v2.subVector(1,N) = v1.subVector(1,N) - temp*u;
        }
    }

    template <class V> 
    static void Householder<V>::unpack(T& u0)
    {
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


    //
    // H * v
    // TODO: break out add option at compile time rather than runtime if
    //

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void MultMV(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else {
            m1.multV(v2,v3);
            Scale(x,v3);
        }
    }
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void NoAliasMultMV(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m1,v2,v3); }
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void AliasMultMV(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m1,v2,v3); }

    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasMultMV<add>(x,m2.transpose(),v1,v3); }
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { AliasMultMV<add>(x,m2.transpose(),v1,v3); }
 
    template <int ix, class T, class V1, class V2>
    static void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { m2.transpose().multEqV(v1); Scale(x,v1); }
    template <int ix, class T, class V1, class V2>
    static void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { MultEqVM(v1,x,m2); }
    template <int ix, class T, class V1, class V2>
    static void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { MultEqVM(v1,x,m2); }

    //
    // v / H:
    // v3 = H^-1 * v1 = H*v1
    //

    template <int ix, class T, class V1, class V2, class V3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { m2.multV(v1,v3); Scale(x,v3); }
    template <int ix, class T, class V1, class V2, class V3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class V2, class V3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }

    template <class V1, class V2>
    static void LDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { m2.multEqV(v1); }
    template <class V1, class V2>
    static void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { LDivEq(v1,m2); }
    template <class V1, class V2>
    static void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { LDivEq(v1,m2); }

    //
    // v % H:
    // v3 = v1 H^-1 = v1 H
    // v3 = HT v1
    //

    template <int ix, class T, class V1, class V3>
    static void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { m2.transpose().multV(v1,v3); Scale(x,v3); }
    template <int ix, class T, class V1, class V3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class V3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const Householder<V2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }

    template <class V1, class V2>
    static void RDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { m2.transpose().multEqV(v1); }
    template <class V1, class V2>
    static void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { RDivEq(v1,m2); }
    template <class V1, class V2>
    static void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const Householder<V2>& m2)
    { RDivEq(v1,m2); }



    //
    // H * m
    // TODO: break out add option at compile time rather than runtime if
    //

    template <bool add, int ix, class T, class V1, class M2, class M3>
    static void MultMM(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else {
            m1.mult(m2,m3);
            Scale(x,m3);
        }
    }
    template <bool add, int ix, class T, class V1, class M2, class M3>
    static void NoAliasMultMM(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }
    template <bool add, int ix, class T, class V1, class M2, class M3>
    static void AliasMultMM(
        const Scaling<ix,T>& x, const Householder<V1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class V2, class M3>
    static void MultVM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(x,m2.transpose(),m1,m3); }
    template <bool add, int ix, class T, class M1, class V2, class M3>
    static void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasMultMM<add>(x,m2.transpose(),m1,m3); }
    template <bool add, int ix, class T, class M1, class V2, class M3>
    static void AliasMultVM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasMultMM<add>(x,m2.transpose(),m1,m3); }
 
    template <int ix, class T, class M1, class V2>
    static void MultEqMM(
        BaseMatrix_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { m2.transpose().multEq(m1); Scale(x,m1); }
    template <int ix, class T, class M1, class V2>
    static void NoAliasMultEqMM(
        BaseMatrix_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { MultEqMM(m1,x,m2); }
    template <int ix, class T, class M1, class V2>
    static void AliasMultEqMM(
        BaseMatrix_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const Householder<V2>& m2)
    { MultEqMM(m1,x,m2); }

    //
    // m / H:
    //

    template <int ix, class T, class M1, class V2, class M3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { m2.mult(m1,m3); Scale(x,m3); }
    template <int ix, class T, class M1, class V2, class M3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }
    template <int ix, class T, class M1, class V2, class M3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }

    template <class M1, class V2>
    static void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { m2.multEq(m1); }
    template <class M1, class V2>
    static void NoAliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { LDivEq(m1,m2); }
    template <class M1, class V2>
    static void AliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { LDivEq(m1,m2); }

    //
    // m % H:
    //

    template <int ix, class T, class M1, class V2, class M3>
    static void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { m2.transpose().mult(m1,m3); Scale(x,m3); }
    template <int ix, class T, class M1, class V2, class M3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { RDiv(x,m1,m2,m3); }
    template <int ix, class T, class M1, class V2, class M3>
    static void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const Householder<V2>& m2, BaseMatrix_Mutable<M3>& m3)
    { RDiv(x,m1,m2,m3); }

    template <class M1, class V2>
    static void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { m2.transpose().multEq(m1); }
    template <class M1, class V2>
    static void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { RDivEq(m1,m2); }
    template <class M1, class V2>
    static void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const Householder<V2>& m2)
    { RDivEq(m1,m2); }


#if 0
    // I'll deal with these next...
    template <class T> 
    void BlockHouseholderAugment(
        const GenMatrix<T>& Y, const UpperTriMatrixView<T>& Z, T beta)
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
        TMVAssert(!Z.isconj());
        int M = Y.colsize();
        int N = Y.rowsize()-1; // # of cols already computed
#ifdef XDEBUG
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

        if (beta == T(0)) {
            Z.col(N,0,N+1).setZero();
        } else if (N == 0) {
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr() >= Z.first);
            TMVAssert(Z.ptr() < Z.last);
#endif
            *Z.ptr() = beta;
        } else {
            ConstVectorView<T> v = Y.col(N,N+1,M);
            VectorView<T> z = Z.col(N,0,N);
            z = Y.subMatrix(N+1,M,0,N).adjoint()*v;
            z += Y.row(N,0,N).conjugate();
            z = -beta * Z.subTriMatrix(0,N) * z;
            // Z(N,N) = beta
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) >= Z.first);
            TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) < Z.last);
#endif
            *(Z.ptr() + N*(Z.stepi()+Z.stepj())) = beta;
        }
#ifdef XDEBUG
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

    template <class T> 
    void BlockHouseholderMakeZ(
        const GenMatrix<T>& Y, const UpperTriMatrixView<T>& Z, 
        const GenVector<T>& beta)
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
        TMVAssert(!Z.isconj());
        TMVAssert(!Y.isconj());
        TMVAssert(!beta.isconj());
#ifdef XDEBUG
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

        if (N==1) {
            // I - Y Z Yt = I - beta* v vt
            // Therefore:
            // Y.col(0) = v
            // Z(0,0) = beta*
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr() >= Z.first);
            TMVAssert(Z.ptr() < Z.last);
#endif
            *Z.ptr() = TMV_CONJ(*beta.cptr());
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
            T* Z00 = Z.ptr();
            T* Z01 = Z00 + Z.stepj();
            T* Z11 = Z01 + Z.stepi();
            const T cb0 = TMV_CONJ(*beta.cptr());
            const T cb1 = TMV_CONJ(*(beta.cptr() + beta.step()));
            const T cY10 = TMV_CONJ(*(Y.cptr() + Y.stepi()));
            // Z(0,0) = TMV_CONJ(beta(0))
#ifdef TMVFLDEBUG
            TMVAssert(Z00 >= Z.first);
            TMVAssert(Z00 < Z.last);
            TMVAssert(Z01 >= Z.first);
            TMVAssert(Z01 < Z.last);
            TMVAssert(Z11 >= Z.first);
            TMVAssert(Z11 < Z.last);
#endif
            *Z00 = cb0;
            // Z(1,1) = TMV_CONJ(beta(1))
            *Z11 = cb1;
            T temp = Y.col(0,2,M).conjugate()*Y.col(1,2,M);
            // temp += TMV_CONJ(Y(1,0))
            temp += cY10;
            // Z(0,1) = -Z(0,0)*Z(1,1)*temp;
            *Z01 = -cb0*cb1*temp;
        } else {
            int j1 = (N+1)/2;
            ConstMatrixView<T> Y1 = Y.colRange(0,j1);
            UpperTriMatrixView<T> Z1 = Z.subTriMatrix(0,j1);
            ConstVectorView<T> beta1 = beta.subVector(0,j1);
            BlockHouseholderMakeZ(Y1,Z1,beta1);

            ConstMatrixView<T> Y2 = Y.subMatrix(j1,M,j1,N);
            UpperTriMatrixView<T> Z2 = Z.subTriMatrix(j1,N);
            ConstVectorView<T> beta2 = beta.subVector(j1,N);
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
            MatrixView<T> Z3 = Z.subMatrix(0,j1,j1,N);
            Z3 = Y1.rowRange(j1,N).adjoint() *
                Y.subMatrix(j1,N,j1,N).lowerTri(UnitDiag);
            Z3 += Y1.rowRange(N,M).adjoint() * Y.subMatrix(N,M,j1,N);
            Z3 = -Z1*Z3;
            Z3 *= Z2;
        }
#ifdef XDEBUG
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

    template <class T, class T2> 
    void BlockHouseholderLMult(
        const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
        const MatrixView<T2>& m)
    {
        // The input Y,Z are such that (I - YZYt) is a Block Householder matrix.
        // The upper square portion of Y is taken to be unit lower triangular.
        // ie. the diagonal and upper triangular portion are not referenced.
        // The routine then finds m <- (I - YZYt) m

        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
#ifdef XDEBUG
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z);
        Matrix<T2> m0(m);
        Matrix<T> H = T(1) - Y0*Z0*Y0.adjoint();
        Matrix<T2> Hm = H*m0;
#endif

        int M = Y.colsize();
        int N = Y.rowsize();

        if (m.iscm()) {
            Matrix<T2,ColMajor> ZYtm = 
                Y.rowRange(0,N).lowerTri(UnitDiag).adjoint() * m.rowRange(0,N);
            ZYtm += Y.rowRange(N,M).adjoint() * m.rowRange(N,M);
            ZYtm = Z * ZYtm;
            m.rowRange(0,N) -= Y.rowRange(0,N).lowerTri(UnitDiag) * ZYtm;
            m.rowRange(N,M) -= Y.rowRange(N,M) * ZYtm;
        } else {
            Matrix<T2,RowMajor> ZYtm = 
                Y.rowRange(0,N).lowerTri(UnitDiag).adjoint() * m.rowRange(0,N);
            ZYtm += Y.rowRange(N,M).adjoint() * m.rowRange(N,M);
            ZYtm = Z * ZYtm;
            m.rowRange(0,N) -= Y.rowRange(0,N).lowerTri(UnitDiag) * ZYtm;
            m.rowRange(N,M) -= Y.rowRange(N,M) * ZYtm;
        }
#ifdef XDEBUG
        if (Norm(Hm-m) > 0.001*Norm(m0)*Norm(H)) {
            cerr<<"BlockHouseholderLMult\n";
            cerr<<"Input: Y = "<<Y0<<endl;
            cerr<<"Z = "<<Z0<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Output: m = "<<m<<endl;
            abort();
        }
#endif
    }

    template <class T, class T2> 
    void BlockHouseholderLDiv(
        const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
        const MatrixView<T2>& m)
    {
        // The routine finds m <- (I - YZYt)^-1 m
        // = (I - YZtYt) m

        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
#ifdef XDEBUG
        Matrix<T> Y0(Y);
        Y0.upperTri().setToIdentity();
        Matrix<T> Z0(Z.adjoint());
        Matrix<T2> m0(m);
        Matrix<T> Hinv = T(1) - Y0*Z0*Y0.adjoint();
        Matrix<T2> Hm = Hinv*m0;
#endif

        int M = Y.colsize();
        int N = Y.rowsize();

        if (m.isrm()) {
            Matrix<T2,RowMajor> ZtYtm = 
                Y.rowRange(0,N).lowerTri(UnitDiag).adjoint() * m.rowRange(0,N);
            ZtYtm += Y.rowRange(N,M).adjoint() * m.rowRange(N,M);
            ZtYtm = Z.adjoint() * ZtYtm;
            m.rowRange(0,N) -= Y.rowRange(0,N).lowerTri(UnitDiag) * ZtYtm;
            m.rowRange(N,M) -= Y.rowRange(N,M) * ZtYtm;
        } else {
            Matrix<T2,ColMajor> ZtYtm = 
                Y.rowRange(0,N).lowerTri(UnitDiag).adjoint() * m.rowRange(0,N);
            ZtYtm += Y.rowRange(N,M).adjoint() * m.rowRange(N,M);
            ZtYtm = Z.adjoint() * ZtYtm;
            m.rowRange(0,N) -= Y.rowRange(0,N).lowerTri(UnitDiag) * ZtYtm;
            m.rowRange(N,M) -= Y.rowRange(N,M) * ZtYtm;
        }
#ifdef XDEBUG
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

    template <class T> 
    void BlockHouseholderUnpack(
        const MatrixView<T>& Y, const GenUpperTriMatrix<T>& Z,
        const MatrixView<T>& m)
    {
        // This routine multiplies the rest of the matrix by the 
        // BlockHouseholder matrix Ht, defined by Y,Z.
        // Then Y is unpacked in place.
        TMVAssert(Y.colsize() > 0);
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(Y.colsize() == m.colsize());
        TMVAssert(Y.rowsize() == Z.size());

        int M = Y.colsize();
        int N = Y.rowsize();

        // Multiply the rest of m by Ht
        BlockHouseholderLMult(Y,Z,m);
        // Make the first N columns equal to 
        // Ht [ I ] = (I - YZYt) [ I ]
        //    [ 0 ]              [ 0 ]
        UpperTriMatrix<T,NonUnitDiag,RowMajor> temp = 
            -Z * Y.rowRange(0,N).lowerTri(UnitDiag).adjoint();
        Y.rowRange(N,M) *= temp;
        Y.rowRange(0,N) = Y.rowRange(0,N).lowerTri(UnitDiag) * temp;
        Y.rowRange(0,N).diag().addToAll(T(1));
    }
#endif

#if 0
    // Expand this out in the calling routine.
    template <class T> T HouseholderReflect(const MatrixView<T>& m, T& det)
    {
        // Multiplies m by a Householder matrix H which rotates
        // the first column into y e0.
        // The vector v of the  Householder matrix is stored
        // in the first column of m, except for the first element.
        // Beta is the return value.
        // The rest of the matrix is multiplied by H.
        // (For the first column, this means that y is the new first element.)
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.rowsize() > 0);
#ifdef XDEBUG
        Matrix<T> m0(m);
#endif

        const VectorView<T> v = m.col(0,1,m.colsize());
        T beta;
        if (m.isconj()) {
            T m00 = TMV_CONJ(*m.cptr());
            beta = HouseholderReflect(m00,v,det);
#ifdef TMVFLDEBUG
            TMVAssert(m.ptr() >= m.first);
            TMVAssert(m.ptr() < m.last);
#endif
            *m.ptr() = TMV_CONJ(m00);
        } else {
            beta = HouseholderReflect(*m.ptr(),v,det);
        }
        if (beta != T(0)) HouseholderLMult(v,beta,m.colRange(1,m.rowsize()));
#ifdef XDEBUG
        Vector<T> vv(m.colsize());
        vv(0) = T(1);
        vv.subVector(1,vv.size()) = m.col(0,1,vv.size());
        Matrix<T> H = T(1)-beta*(vv^vv.conjugate());
        Matrix<T> Hm = H * m0;
        Matrix<T> Hm2 = m;
        Hm2.col(0,1,vv.size()).setZero();
        if (Norm(Hm-Hm2) > 0.001*Norm(m0)) {
            cerr<<"Householder Reflect\n";
            cerr<<"Input: m0 = "<<m0<<endl;
            cerr<<"Output: vv = "<<vv<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"m = "<<m<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"vv^vvt = "<<(vv^vv.conjugate())<<endl;
            cerr<<"beta*vv^vvt = "<<beta*(vv^vv.conjugate())<<endl;
            cerr<<"1-beta*vv^vvt = "<<(T(1)-beta*(vv^vv.conjugate()))<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Hm2 = "<<Hm2<<endl;
            abort();
        }
#endif
        return beta;
    }
#endif

#if 0
    // Expand this out in the calling routine.
    template <class T> 
    void HouseholderUnpack(const MatrixView<T>& m, T beta)
    {
        // The input matrix is taken to have a Householder vector
        // stored in the first column (not including the first element.   
        // This routine multiplies the rest of the matrix by the Householder   
        // matrix Ht.
        // The First column is then set to Ht times e0.

        TMVAssert(m.colsize() > 0);
        TMVAssert(m.rowsize() > 0);

        // Multiply the rest of m by Ht = I - beta* v vt
        HouseholderLMult(
            m.col(0,1,m.colsize()),TMV_CONJ(beta),
            m.subMatrix(0,m.colsize(),1,m.rowsize()));
        // Finally, make the first column equal to Ht times e0
        if (m.isconj()) {
            T m00 = TMV_CONJ(*m.cptr());
            HouseholderUnpack(m00,m.col(0,1,m.colsize()),beta);
#ifdef TMVFLDEBUG
            TMVAssert(m.ptr() >= m.first);
            TMVAssert(m.ptr() < m.last);
#endif
            *m.ptr() = TMV_CONJ(m00);
        } else {
            HouseholderUnpack(*m.ptr(),m.col(0,1,m.colsize()),beta);
        }
    }
#endif

} // namespace tmv

#endif
