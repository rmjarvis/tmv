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


#include "TMV_SymLDLDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <bool herm, class T> 
    typename Sym2x2_Helper<herm,T>::d_type SymInvert_2x2(T& a, T& b, T& c)
    {
        typedef typename Traits<T>::real_type RT;
        typedef typename Sym2x2_Helper<herm,T>::d_type d_type;
        // Rescale to help avoid underflow/overflow:
        RT max = TMV_MAX(TMV_ABS2(a),TMV_MAX(TMV_ABS2(b),TMV_ABS2(c)));
        a /= max;
        b /= max;
        c /= max;

        // Invert matrix [ a  c* ] -->  1/(ab-|c|^2) [ b -c* ]
        //               [ c  b  ]                   [ -c a  ]
        d_type d = Sym2x2_Helper<herm,T>::calculateD(a,b,c); 
        TMV_SWAP(a,b);
        d *= max;
        a /= d; 
        b /= d;
        c /= -d;
        return d*max;
    }

    template <class T, class T1ab, class T1c> 
    static void LMultEq_2x2(T1ab a, T1ab b, T1c c, T1c cc, T& x, T& y)
    {
        // Solve [ x ] <- [ a  cc ] [ x ]
        //       [ y ]    [ c  b  ] [ y ]
        T tempx = x;
        x = (a*x+cc*y);
        y = (b*y+c*tempx);
    }

    template <class T, class T1ab, class T1c> 
    static void LMultEq_2x2(
        T1ab a, T1ab b, T1c c, T1c cc, MatrixView<T> m)
    {
        // Solve m <- [ a  cc ] m 
        //            [ c  b  ]
        TMVAssert(m.colsize() == 2);
        //cout<<"LMultEq_2x2: m = "<<TMV_Text(m)<<"  "<<m<<endl;
        TMVAssert(m.ct() == NonConj);
        if (m.isrm()) {
            T* m0 = m.ptr();
            T* m1 = m0 + m.stepi();
            for (ptrdiff_t k=m.rowsize();k>0;--k,++m0,++m1) {
#ifdef TMVFLDEBUG
                TMVAssert(m0 >= m._first);
                TMVAssert(m0 < m._last);
                TMVAssert(m1 >= m._first);
                TMVAssert(m1 < m._last);
#endif
                LMultEq_2x2(a,b,c,cc,*m0,*m1);
            }
        } else {
            const ptrdiff_t sj = m.stepj();
            T* m0 = m.ptr();
            T* m1 = m0 + m.stepi();
            for (ptrdiff_t k=m.rowsize();k>0;--k,m0+=sj,m1+=sj) {
#ifdef TMVFLDEBUG
                TMVAssert(m0 >= m._first);
                TMVAssert(m0 < m._last);
                TMVAssert(m1 >= m._first);
                TMVAssert(m1 < m._last);
#endif
                LMultEq_2x2(a,b,c,cc,*m0,*m1);
            }
        }
    }

    template <bool herm, class T, class T1> 
    static void SymLMultEq_2x2(T1 a, T1 b, T1 c, MatrixView<T> m)
    {
        //cout<<"SymLMultEq_2x2: m = "<<TMV_Text(m)<<"  "<<m<<endl;
        TMVAssert(m.colsize() == 2);
        if (herm)
            if (m.isconj())
                LMultEq_2x2(TMV_REAL(a),TMV_REAL(b),TMV_CONJ(c),c,
                            m.conjugate());
            else
                LMultEq_2x2(TMV_REAL(a),TMV_REAL(b),c,TMV_CONJ(c),m);
        else
            if (m.isconj())
                LMultEq_2x2(TMV_CONJ(a),TMV_CONJ(b),TMV_CONJ(c),TMV_CONJ(c),
                            m.conjugate());
            else
                LMultEq_2x2(a,b,c,c,m);
    }

    template <bool herm, class T, class T1> 
    static void SymRMultEq_2x2(T1 a, T1 b, T1 c, MatrixView<T> m)
    {
        TMVAssert(m.rowsize() == 2);
        if (herm)
            if (m.isconj())
                LMultEq_2x2(TMV_REAL(a),TMV_REAL(b),c,TMV_CONJ(c),m.adjoint());
            else
                LMultEq_2x2(TMV_REAL(a),TMV_REAL(b),TMV_CONJ(c),c,
                            m.transpose());
        else
            if (m.isconj())
                LMultEq_2x2(TMV_CONJ(a),TMV_CONJ(b),TMV_CONJ(c),TMV_CONJ(c),
                            m.adjoint());
            else
                LMultEq_2x2(a,b,c,c,m.transpose());
    }

    template <bool herm, class T, class T1> 
    void PseudoDiag_LDivEq(
        const GenVector<T1>& D, const GenVector<T1>& xD,
        MatrixView<T> m)
    {
        TMVAssert(D.size() == m.colsize());
        TMVAssert(xD.size()+1 == m.colsize());
        TMVAssert(D.ct() == NonConj);
        TMVAssert(xD.ct() == NonConj);

#ifdef XDEBUG
        //cout<<"Start PseudoDiag_LDivEq\n";
        //cout<<"xD = "<<xD<<endl;
        Matrix<T> m0(m);
        Matrix<T1> DD(D.size(),D.size(),T1(0));
        DD.diag() = D;
        DD.diag(-1) = xD;
        DD.diag(1) = herm ? xD.conjugate() : xD.view();
        //cout<<"xD = "<<xD<<endl;
#endif

        const T1* Di = D.cptr();
        const T1* xDi = xD.cptr();

        const ptrdiff_t N = D.size();
        const ptrdiff_t sd = D.step();
        const ptrdiff_t sx = xD.step();

        for(ptrdiff_t i=0;i<N;) {
            if (i==N-1 || *xDi == T1(0)) {
                if (herm) m.row(i) /= TMV_REAL(*Di);
                else m.row(i) /= *Di;
                Di+=sd; xDi+=sx; ++i;
            } else {
                T1 x = *Di;
                T1 y = *(Di+=sd);
                T1 z = *xDi;
                SymInvert_2x2<herm>(x,y,z);
                SymLMultEq_2x2<herm>(x,y,z,m.rowRange(i,i+2));
                Di+=sd,xDi+=2*sx,i+=2;
            }
        }
#ifdef XDEBUG
        //cout<<"xD = "<<xD<<endl;
        Matrix<T> m2 = DD * m;
        if (Norm(m2-m0) > 0.001*Norm(DD)*Norm(m0)) {
            cerr<<"PseudoDiag_LDivEq: m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"xD = "<<xD<<endl;
            cerr<<"DD = "<<DD<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"m2 = DD*m "<<m2<<endl;
            abort();
        }
        //cout<<"Done PseudoDiag_LDivEq: xD = "<<xD<<endl;
#endif
    }

    template <bool herm, class T, class T1> 
    void PseudoDiag_LMultEq(
        const GenVector<T1>& D, const GenVector<T1>& xD,
        MatrixView<T> m)
    {
        TMVAssert(D.size() == m.colsize());
        TMVAssert(xD.size()+1 == m.colsize());
        TMVAssert(D.ct() == NonConj);
        TMVAssert(xD.ct() == NonConj);

#ifdef XDEBUG
        //cout<<"Start PseudoDiag_LMultEq: xD = "<<xD<<endl;
        Matrix<T> m0(m);
        Matrix<T1> DD(D.size(),D.size(),T1(0));
        DD.diag() = D;
        DD.diag(-1) = xD;
        DD.diag(1) = herm ? xD.conjugate() : xD.view();
        Matrix<T> m2 = DD * m;
        //cout<<"In PseudoLMultEq:\n";
        //cout<<"D = "<<D<<endl;
        //cout<<"xD = "<<xD<<endl;
        //cout<<"DD = "<<DD<<endl;
        //cout<<"m = "<<m<<endl;
        //cout<<"m0 = "<<m0<<endl;
        //cout<<"DD*m = "<<m2<<endl;
#endif

        const T1* Di = D.cptr();
        const T1* xDi = xD.cptr();

        const ptrdiff_t N = D.size();
        const ptrdiff_t sd = D.step();
        const ptrdiff_t sx = xD.step();

        for(ptrdiff_t i=0;i<N;) {
            if (i==N-1 || *xDi == T1(0)) {
                if (herm) 
                    m.row(i) *= TMV_REAL(*Di);
                else
                    m.row(i) *= *Di;
                Di+=sd; xDi+=sx; ++i;
            } else {
                T1 x = *Di;
                T1 y = *(Di+=sd);
                T1 z = *xDi;
                SymLMultEq_2x2<herm>(x,y,z,m.rowRange(i,i+2));
                Di+=sd,xDi+=2*sx,i+=2;
            }
        }
#ifdef XDEBUG
        //cout<<"xD = "<<xD<<endl;
        //cout<<"Done: m = "<<m<<endl;
        if (Norm(m2-m) > 0.00001*Norm(DD)*Norm(m0)) {
            cerr<<"PseudoDiag_LMultEq: m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"xD = "<<xD<<endl;
            cerr<<"DD = "<<DD<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"-> m = "<<m<<endl;
            abort();
        }
        //cout<<"Done PseudoDiag_LMultEq: xD = "<<xD<<endl;
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymLDLPseudo.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


