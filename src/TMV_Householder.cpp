///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG

#include "tmv/TMV_Householder.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_VIt.h"

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    T HouseholderReflect(T& x0, VectorView<T> x, T& det)
    {
#ifdef XDEBUG
        Vector<T> xx(x.size()+1);
        xx(0) = x0;
        xx.subVector(1,xx.size()) = x;
        std::cout<<"Start reflect: x0 = "<<x0<<", x = "<<x<<std::endl;
#endif

        // Since this routine involves squares of elements, we risk overflow
        // and underflow problems if done naively.  The simplest (although
        // probably not the most efficient) solution is to scale all the 
        // intermediate values by the maximum abs value in the vector.
        //std::cout<<"Start HouseholderReflect"<<std::endl;
        //std::cout<<"x = "<<x<<std::endl;
        RT scale = x.maxAbs2Element();
        //std::cout<<"scale = "<<scale<<std::endl;
        RT absx0 = TMV_ABS(x0);
        //std::cout<<"x0,absx0 = "<<x0<<"  "<<absx0<<std::endl;
        if (absx0 > scale) scale = absx0;
        //std::cout<<"scale => "<<scale<<std::endl;
        if (TMV_Underflow(scale)) {
            // Then the situation is hopeless, and we should just zero out
            // the whole vector.
            //std::cout<<"scale is essentially 0\n";
            x.setZero();
            x0 = T(0);
            return T(0);
        }

        // Finds the Householder matrix H which rotates v into y e0.
        // The vector v of the Householder matrix is stored in v,
        // except for the first element.
        // Beta is the return value.

        // Determine normx = |x|
        RT invscale = TMV_InverseOf(scale);
        RT normsqx = x.normSq(invscale);
        //std::cout<<"invscale = "<<invscale<<std::endl;
        //std::cout<<"new normsqx = "<<normsqx<<std::endl;

        // if all of x other than first element are 0, H is identity
        if (normsqx == RT(0) && TMV_IMAG(x0) == RT(0)) {
            //std::cout<<"no reflection necessary\n";
            // Set them all to x explicitly in case underflow let to the 0.
            x.setZero();
            return T(0); 
            // Determinant in this case is 1 (H = I), so leave det alone.
        }

        // Finish calculating normx in usual case
        absx0 *= invscale;
        x0 *= invscale;
        //std::cout<<"absx0 => "<<absx0<<", x0 => "<<x0<<std::endl;
        RT normsqx0 = absx0*absx0;
        normsqx += normsqx0;
        RT normx = TMV_SQRT(normsqx);
        //std::cout<<"normx = "<<normx<<std::endl;

        // y = +- |x|
        RT y =  TMV_REAL(x0) > 0 ? -normx : normx;
        //std::cout<<"y = "<<y<<std::endl;

        // beta = 1 / (|x|^2 + |x| x0)
        // H = I - beta v vt
        // with v = x - y e0 in first column
        // Renormalize beta,v so that v(0) = 1
        T v0 = x0-y;
        RT normv0 = TMV_NORM(v0);
        T beta = TMV_Divide(normv0 , (normsqx - y * x0));
        T invv0 = TMV_InverseOf(v0);
        //std::cout<<"v0 = "<<v0<<", normv0 = "<<normv0;
        //std::cout<<", beta = "<<beta<<", invv0 = "<<invv0<<std::endl;

        // Sometimes this combination can underflow, so check.
        T scaled_invv0 = invv0 * invscale;
        //std::cout<<"scaled_invv0 = "<<scaled_invv0<<std::endl;
        if ((TMV_REAL(invv0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invv0))) ||
            (TMV_IMAG(invv0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invv0)))) {
            //std::cout<<"Two steps:\n";
            x *= invscale;
            //std::cout<<"x -> "<<x<<std::endl;
            x *= invv0;
        } else {
            x *= scaled_invv0;
        }
        //std::cout<<"x => "<<x<<std::endl;

        x0 = y*scale;
        //std::cout<<"x0 = "<<x0<<std::endl;

        // Determinant of H = -beta^2/|beta|^2
        // But we are actually keeping track of the determinant of 
        // Q which is now multiplied by Ht.
        // The determinant of Ht = -conj(beta)^2/|beta|^2
        if (det != T(0)) {
            if (TMV_IMAG(beta) == RT(0))
                det = -det;
            else 
                det *= -TMV_CONJ(beta*beta)/TMV_NORM(beta);
        }
        //std::cout<<"det = "<<det<<std::endl;
#ifdef XDEBUG
        Vector<T> vv(xx.size());
        vv(0) = T(1);
        vv.subVector(1,vv.size()) = x;
        Matrix<T> H = T(1)-beta*(vv^vv.conjugate());
        // Check the following:
        // H * xx = (Norm(xx),0,0.0...)
        // x0 = Norm(xx)
        Vector<T> Hxx = H * xx;
        if (!(TMV_ABS(TMV_ABS(Hxx(0)/scale)-TMV_ABS(x0/scale)) <=0.0001*TMV_ABS(x0/scale)) ||
            !(TMV_ABS(Norm(xx/scale)-TMV_ABS(x0/scale)) <=0.0001*TMV_ABS(x0/scale)) ||
            !(Norm(Hxx.subVector(1,Hxx.size())/scale) <=0.0001*TMV_ABS(x0/scale))) {
            cerr<<"Householder Reflect:\n";
            cerr<<"Input: x = "<<xx<<endl;
            cerr<<"Norm(x) = "<<Norm(xx)<<endl;
            cerr<<"Output: v = "<<vv<<endl;
            cerr<<"x0 = "<<x0<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"xx = "<<xx<<endl;
            cerr<<"Hxx = "<<Hxx<<endl;
            cerr<<"abs(hxx(0)) = "<<TMV_ABS(Hxx(0))<<std::endl;
            cerr<<"abs(x0) = "<<TMV_ABS(x0)<<std::endl;
            cerr<<"abs(hxx(0))-abs(x0) = "<<TMV_ABS(Hxx(0))-TMV_ABS(x0)<<endl;
            cerr<<"abs(abs(hxx(0))-abs(x0)) = "<<
                TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(x0))<<endl;
            cerr<<"abs(Norm(xx)-abs(x0)) = "<<
                TMV_ABS(Norm(xx))-TMV_ABS(x0)<<endl;
            cerr<<"Norm(Hxx(1,N)) = "<<Norm(Hxx.subVector(1,Hxx.size()))<<endl;
            abort();
        }
#endif
        return beta;
    }

    template <class T> T HouseholderReflect(MatrixView<T> m, T& det)
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

        VectorView<T> v = m.col(0,1,m.colsize());
        T beta;
        if (m.isconj()) {
            T m00 = TMV_CONJ(*m.cptr());
            beta = HouseholderReflect(m00,v,det);
#ifdef TMVFLDEBUG
            TMVAssert(m.ptr() >= m._first);
            TMVAssert(m.ptr() < m._last);
#endif
            *m.ptr() = TMV_CONJ(m00);
        } else {
            beta = HouseholderReflect(*m.ptr(),v,det);
        }
        //std::cout<<"Before LMult: m = "<<m<<std::endl;
        //std::cout<<"v = "<<v<<std::endl;
        //std::cout<<"beta = "<<beta<<std::endl;
        if (beta != T(0)) HouseholderLMult(v,beta,m.colRange(1,m.rowsize()));
        //std::cout<<"m => "<<m<<std::endl;
#ifdef XDEBUG
        if (TMV_Underflow(m0.col(0).maxAbsElement())) return beta;
        Vector<T> vv(m.colsize());
        vv(0) = T(1);
        vv.subVector(1,vv.size()) = m.col(0,1,vv.size());
        Matrix<T> H = T(1)-beta*(vv^vv.conjugate());
        Matrix<T> Hm = H * m0;
        Matrix<T> Hm2 = m;
        Hm2.col(0,1,vv.size()).setZero();
        if (!(Norm(Hm-Hm2) <=0.001*Norm(m0))) {
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

#if 0
    template <class T> 
    T HouseholderReflect(ConjRef<T> x0, VectorView<T> x, T& det)
    {
        T& x0r = x0.getRef();
        x0r = TMV_CONJ(x0r);
        return HouseholderReflect(x0r,x,det);
        // x0r ends up real, so we don't need to care about the conjugation
        // of the output value.
    }

    template <class T> 
    T HouseholderReflect(VarConjRef<T> x0, VectorView<T> x, T& det)
    {
        T& x0r = x0.getRef();
        if (x0.isconj()) x0r = TMV_CONJ(x0r);
        return HouseholderReflect(x0r,x,det);
    }
#endif

    template <class T> 
    T HouseholderReflect(VectorView<T> x, T& det)
    {
        // Same as above, but takes (x0,x) to be contiguous
        return HouseholderReflect(x(0),x.subVector(1,x.size()),det);
    }

    template <class T> 
    bool HouseholderUnReflect(T& y, VectorView<T> x, T& beta)
    {
        // This is similar, except that the roles of y and x0 are swapped.
        // That is, the final rotated vector y e0 is presumed known, as is 
        // the bulk of the input vector.  
        // The unreflected value of the first element in the vector is returned 
        // as y.
        // The rest of the vector is transformed into the Householder vector.
        // This is used for downdating QR decompositions.
        // The return value is true if successful, false if |y| < |x|.

        TMVAssert(x.size() > 0);
        TMVAssert(TMV_IMAG(y) == RT(0));

        RT normsqx1 = NormSq(x);

        // if all of x other than first element are 0, H is identity
        if (normsqx1 == RT(0)) {
            beta = 0; 
            return true; 
            // Determinant in this case is 1 (H = I), so leave det alone.
        }

        RT normsqx = TMV_SQR(TMV_REAL(y));
        RT normsqx0 = normsqx - normsqx1;
        if (normsqx0 < RT(0)) return false;

        // Same consideration on the +-: Maximize 1/beta.
        RT x0 =  TMV_REAL(y) > 0 ? -TMV_SQRT(normsqx0) : TMV_SQRT(normsqx0);

        // beta = 1 / (|x|^2 + |x| x0)
        // H = I - beta v vt
        // with v = x - y e0 in first column
        // Renormalize beta,v so that v(0) = 1
        T v0 = x0-y;
        RT normv0 = TMV_NORM(v0);
        beta = normv0 / (normsqx - y * x0);

        x /= v0;
        y = x0;

        return true;
    }

#if 0
    template <class T> 
    bool HouseholderUnReflect(ConjRef<T> x0, VectorView<T> x, T& beta)
    {
        T& x0r = x0.getRef();
        TMVAssert(TMV_IMAG(x0r) == RT(0));
        return HouseholderUnReflect(x0r,x,beta);
    }

    template <class T> 
    bool HouseholderUnReflect(VarConjRef<T> x0, VectorView<T> x, T& beta)
    {
        TMVAssert(TMV_IMAG(x0.getRef()) == RT(0));
        return HouseholderUnReflect(x0.getRef(),x,beta);
    }
#endif

    template <class T1, class T2> 
    void HouseholderLMult(
        const GenVector<T1>& v, T1 beta, VectorView<T2> m0,
        MatrixView<T2> mx)
    {
        // The input vector, v, is taken to be the vector for a  
        // Householder matrix, H.  This routine takes 
        // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 vt ) ( m0 ) ]
        // ( mx )      ( mx )   [ ( mx )        ( v )          ( mx ) ]
        // 
        // ( m0 ) -= beta (   m0 + vt mx   )
        // ( mx )         ( v (m0 + vt mx) )

        TMVAssert(v.size() == mx.colsize());
        TMVAssert(m0.size() == mx.rowsize());
#ifdef XDEBUG
        Matrix<T2> mm(mx.colsize()+1,m0.size());
        mm.row(0) = m0;
        mm.rowRange(1,mm.colsize()) = mx;
        std::cout<<"Start HouseholderLMult"<<std::endl;
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"beta = "<<beta<<std::endl;
        std::cout<<"m0 = "<<m0<<std::endl;
        std::cout<<"mx = "<<mx<<std::endl;
#endif

        if (m0.size() == 0) return;
        else if (beta != T1(0)) {
            Vector<T2> temp(mx.rowsize());
            temp = v.conjugate() * mx;
            temp += m0;
            temp *= beta;

            m0 -= temp;
            mx -= v^temp;
        }
#ifdef XDEBUG
        std::cout<<"m0 => "<<m0<<std::endl;
        std::cout<<"mx => "<<mx<<std::endl;
        Vector<T1> vv(v.size()+1);
        vv(0) = T1(1);
        vv.subVector(1,vv.size()) = v;
        Matrix<T1> H = T1(1) - beta*(vv^vv.conjugate());
        Matrix<T2> Hm = H * mm;
        Matrix<T2> Hm2(mx.colsize()+1,m0.size());
        Hm2.row(0) = m0;
        Hm2.rowRange(1,Hm2.colsize()) = mx;
        if (!(Norm(Hm-Hm2) <=0.001*Norm(Hm))) {
            cerr<<"HouseholderLMult\n";
            cerr<<"Input: m = "<<mm<<endl;
            cerr<<"vv = "<<vv<<endl;
            cerr<<"H = "<<H<<endl;
            cerr<<"Hm = "<<Hm<<endl;
            cerr<<"Output: m = "<<Hm2<<endl;
            abort();
        }
#endif
    }

    template <class T1, class T2> 
    void HouseholderLMult(
        const GenVector<T1>& v, T1 beta, MatrixView<T2> m)
    {
        // The same as above, except m0 and mx are a single contiguous
        // matrix, m.
        TMVAssert(m.colsize() > 0);
        HouseholderLMult(v,beta,m.row(0),m.rowRange(1,m.colsize()));
    }

    template <class T> 
    void HouseholderUnpack(T& v0, VectorView<T> v, T beta)
    {
        // The input matrix is taken to have a Householder vector
        // stored in the first column (not including the first element.   
        // This routine multiplies the rest of the matrix by the Householder   
        // matrix Ht.
        // The First column is then set to Ht times e0.

        if (beta == T(0)) {
            // then all but v(0) is already 0
            v0 = T(1);
        } else {
            // v <- (I-beta* (1 v)T (1 v*)) e0 = e0 - beta* v 
            v0 = T(1)-TMV_CONJ(beta);
            v *= -TMV_CONJ(beta);
        }
    }

#if 0
    template <class T> 
    void HouseholderUnpack(ConjRef<T> v0, VectorView<T> v, T beta)
    {
        T vv = v0;
        HouseholderUnpack(vv,v,beta);
        v0 = vv;
    }

    template <class T> 
    void HouseholderUnpack(VarConjRef<T> v0, VectorView<T> v, T beta)
    {
        if (v0.isconj()) {
            T vv = v0;
            HouseholderUnpack(vv,v,beta);
            v0 = vv;
        } else {
            HouseholderUnpack(v0.getRef(),v,beta);
        }
    }
#endif


    template <class T> 
    void HouseholderUnpack(MatrixView<T> m, T beta)
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
            TMVAssert(m.ptr() >= m._first);
            TMVAssert(m.ptr() < m._last);
#endif
            *m.ptr() = TMV_CONJ(m00);
        } else {
            HouseholderUnpack(*m.ptr(),m.col(0,1,m.colsize()),beta);
        }
    }

    template <class T> 
    void BlockHouseholderAugment(
        const GenMatrix<T>& Y, UpperTriMatrixView<T> Z, T beta)
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
        ptrdiff_t M = Y.colsize();
        ptrdiff_t N = Y.rowsize()-1; // # of cols already computed
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
            TMVAssert(Z.ptr() >= Z._first);
            TMVAssert(Z.ptr() < Z._last);
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
            TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) >= Z._first);
            TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) < Z._last);
#endif
            *(Z.ptr() + N*(Z.stepi()+Z.stepj())) = beta;
        }
#ifdef XDEBUG
        Matrix<T> YY(Y);
        YY.upperTri().setToIdentity();
        Matrix<T> HH = T(1) - YY * Z * YY.adjoint();
        if (!(Norm(HH-H2) <=0.001*Norm(H2))) {
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
        const GenMatrix<T>& Y, UpperTriMatrixView<T> Z, 
        const GenVector<T>& beta)
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
    {
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
        for(ptrdiff_t i=0;i<Y.rowsize();i++) {
            Matrix<T> H(Y.colsize(),Y.colsize());
            H.setToIdentity();
            H.subMatrix(i,Y.colsize(),i,Y.colsize()) -= beta(i) * 
                (Y0.col(i,i,Y0.colsize()) ^ 
                 Y0.col(i,i,Y0.colsize()).conjugate());
            Htot *= H.adjoint();
        }
#endif

        const ptrdiff_t M = Y.colsize();
        const ptrdiff_t N = Y.rowsize();

        if (N==1) {
            // I - Y Z Yt = I - beta* v vt
            // Therefore:
            // Y.col(0) = v
            // Z(0,0) = beta*
#ifdef TMVFLDEBUG
            TMVAssert(Z.ptr() >= Z._first);
            TMVAssert(Z.ptr() < Z._last);
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
            TMVAssert(Z00 >= Z._first);
            TMVAssert(Z00 < Z._last);
            TMVAssert(Z01 >= Z._first);
            TMVAssert(Z01 < Z._last);
            TMVAssert(Z11 >= Z._first);
            TMVAssert(Z11 < Z._last);
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
            ptrdiff_t j1 = (N+1)/2;
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
        if (!(Norm(Htot-Hnet) <=0.001*Norm(Htot))) {
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
        MatrixView<T2> m)
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

        ptrdiff_t M = Y.colsize();
        ptrdiff_t N = Y.rowsize();

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
        if (!(Norm(Hm-m) <=0.001*Norm(m0)*Norm(H))) {
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
        MatrixView<T2> m)
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

        ptrdiff_t M = Y.colsize();
        ptrdiff_t N = Y.rowsize();

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
        if (!(Norm(Hm-m) <=0.001*Norm(m0)*Norm(Hinv))) {
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
        MatrixView<T> Y, const GenUpperTriMatrix<T>& Z,
        MatrixView<T> m)
    {
        // This routine multiplies the rest of the matrix by the 
        // BlockHouseholder matrix Ht, defined by Y,Z.
        // Then Y is unpacked in place.
        TMVAssert(Y.colsize() > 0);
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(Y.colsize() == m.colsize());
        TMVAssert(Y.rowsize() == Z.size());

        ptrdiff_t M = Y.colsize();
        ptrdiff_t N = Y.rowsize();

        // Multiply the rest of m by Ht
        BlockHouseholderLMult(Y,Z,m);
        // Make the first N columns equal to 
        // Ht [ I ] = (I - YZYt) [ I ]
        //    [ 0 ]              [ 0 ]
        UpperTriMatrix<T,NonUnitDiag|RowMajor> temp = 
            -Z * Y.rowRange(0,N).lowerTri(UnitDiag).adjoint();
        Y.rowRange(N,M) *= temp;
        Y.rowRange(0,N) = Y.rowRange(0,N).lowerTri(UnitDiag) * temp;
        Y.rowRange(0,N).diag().addToAll(T(1));
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_Householder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


