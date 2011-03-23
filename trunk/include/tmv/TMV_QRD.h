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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QR Decomposition.
//
// The basic idea of an QR decomposition is that any 
// matrix A can be decomposed into a unitary matrix
// times an upper triangle matrix.
//
// A = Q R
//
// We do this using Householder transformations, which can
// be stored in a lower triangle matrix, thus other than the 
// diagonal, which they both need, we can store them in the 
// place of the original matrix. It is more convenient to 
// keep the diagonal of Q in place and take the diagonal of R
// separately.
//
// If R is not singular, the solution to A x = b is found from
// Q R x = b
// R x = Qt b
// which can be solved by back-substitutaion.
//
// If m > n, this does not actually give a solution to A x = b,
// since Q Qt != I  (Q is only column orthogonal if m > n.)
// But it does give the value of x which minimizes the 2-norm
// of (A x - b)
//
// The 2-norm of v is the square root of vt v, so
// |A x - b|^2 =
// (A x - b)t (A x - b) = 
// (xt At - bt) (A x - b) =
// xt At A x - bt A x - xt At b + bt b =
// xt Rt Qt Q R x - bt Q R x - xt Rt Qt b + bt b =
// |R x - Qt b|^2 + |b|^2 - |Qt b|^2
// Clearly the x which minimizes this is the solution of R x = Qt b.
//
// If R is singular, then you need QRP Decomposition (see TMV_QRPDiv.h).
//


#ifndef TMV_QRD_H
#define TMV_QRD_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseVector.h"
#include "TMV_Divider.h"
//#include "TMV_Matrix.h"
//#include "TMV_NormU.h"
//#include "TMV_MultUL.h"
//#include "TMV_Permutation.h"

//#include "TMV_PackedQ.h"
//#include "TMV_QRDiv.h"
//#include "TMV_QRDecompose.h"
//#include "TMV_QRInverse.h"

namespace tmv {

    template <bool small, class M>
    struct QRD_Impl;

    template <class M>
    class QRD : 
        public Divider<typename M::value_type>
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;
        typedef typename M::complex_type CT;

        // This next bit finds the storage type to use for the lu matrix
        // regardless of what kind of matrix M is.  e.g. this should
        // work even if M is a TriMatrix or a BandMatrix, etc.
        enum { cs = IntTraits2<M::_colsize,M::_rowsize>::max };
        enum { rs = IntTraits2<M::_colsize,M::_rowsize>::min };
        typedef typename MCopyHelper<T,Rec,cs,rs,false,false>::type qr_type;

        typedef typename qr_type::const_view_type getqr_type;
        typedef typename PackedQ<qr_type,beta_type> getq_type;
        typedef typename qr_type::const_uppertri_type getr_type;
        typedef const beta_type& getbeta_type;

        //
        // Constructors
        //

        // Sets up the internal storage and does the decomposition.
        QRD(const M& A, bool _inplace=false);

        // The copy constructor has transfer of ownership semantics.
        // This way an QRD object can be returned by value, and the 
        // copy is cheap.  I don't think there is any reason to use
        // a more sophisticated technique like shared_ptr or something
        // similar, since there are no non-const methods.
        QRD(const QRD<M>& rhs);

        // Clean up the internal storage
        ~QRD();


        //
        // Division: (not in place)
        // 
        
        template <class M1, class M2>
        void solve(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V1, class V2>
        void solve(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const;

        template <class M1, class M2>
        void solveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V1, class V2>
        void solveTranspose(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const;

        // These are the virtual functions from the Divider base class.
        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<RT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<CT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solve(m1,m2); } 
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<RT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<CT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT> v2) const 
        { solve(v1,v2); } 
        void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solve(v1,v2); } 

        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<RT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<CT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solveTranspose(m1,m2); }
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<RT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<CT>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT> v2) const 
        { solveTranspose(v1,v2); }
        void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& v1, 
            VectorView<CT,UNKNOWN,true> v2) const 
        { solveTranspose(v1,v2); }

        //
        // Perform the division in place
        // (This requires a temporary for QR division.)
        //
        
        template <class M2>
        void solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V2>
        void solveInPlace(BaseVector_Mutable<V2>& v2) const;

        template <class M2>
        void solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const;
        template <class V2>
        void solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const;

        // These are the virtual functions from the Divider base class.
        void doSolveInPlace(MatrixView<RT> m2) const 
        { solveInPlace(m2);  }
        void doSolveInPlace(MatrixView<CT> m2) const 
        { solveInPlace(m2);  }
        void doSolveInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solveInPlace(m2);  }
        void doSolveInPlace(VectorView<RT> v2) const 
        { solveInPlace(v2);  }
        void doSolveInPlace(VectorView<CT> v2) const 
        { solveInPlace(v2);  }
        void doSolveInPlace(VectorView<CT,UNKNOWN,true> v2) const 
        { solveInPlace(v2);  }

        void doSolveTransposeInPlace(MatrixView<RT> m2) const 
        { solveTransposeInPlace(m2); } 
        void doSolveTransposeInPlace(MatrixView<CT> m2) const 
        { solveTransposeInPlace(m2); } 
        void doSolveTransposeInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
        { solveTransposeInPlace(m2); } 
        void doSolveTransposeInPlace(VectorView<RT> v2) const 
        { solveTransposeInPlace(v2); } 
        void doSolveTransposeInPlace(VectorView<CT> v2) const 
        { solveTransposeInPlace(v2); } 
        void doSolveTransposeInPlace(
            VectorView<CT,UNKNOWN,true> v2) const 
        { solveTransposeInPlace(v2); } 


        //
        // Determinant
        //
        
        T det() const;
        RT logDet(T* sign) const;
        bool isSingular() const;

        
        //
        // Inverse
        //
        
        template <class M2>
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const;

        void doMakeInverse(MatrixView<RT> minv) const 
        { makeInverse(minv); } 
        void doMakeInverse(MatrixView<CT> minv) const 
        { makeInverse(minv); } 
        void doMakeInverse(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> minv) const 
        { makeInverse(minv); } 


        //
        // InverseATA
        //
        
        template <class M2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const;

        void doMakeInverseATA(MatrixView<RT> ata) const
        { makeInverseATA(ata); }
        void doMakeInverseATA(MatrixView<CT> ata) const
        { makeInverseATA(ata); }
        void doMakeInverseATA(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> ata) const
        { makeInverseATA(ata); }


        // 
        // Condition (kappa_inf)
        //

        RT condition(RT normInf) const;


        //
        // Access Decomposition
        //

        bool isTrans() const;
        getq_type getQ() const;
        getr_type getR() const;
        getqr_type getQR() const;
        getbeta_type getBeta() const;

        bool preferInPlace() const { return false; }

    private :

        enum { small = (
                M::_colsize != UNKNOWN && M::_rowsize != UNKNOWN &&
                M::_colsize <= 8 && M::_rowsize <= 8 ) };

        mutable std::auto_ptr<QRD_Impl<small,M> > pimpl;

        size_t colsize() const;
        size_t rowsize() const;

        // op= not allowed.
        QRD<M>& operator=(const QRD<M>&);
    };
    
    template <bool isvalid, bool istrans>
    struct QRHelper;

    template <>
    struct QRHelper<true,false>
    {
        template <class M0, class V0, class M2>
        static void makeInverse(const M0& QRx, const V0& beta, M2& m2)
        { QR_Inverse(QRx,beta,m2); }
        template <class M0, class V0, class M2>
        static void solveInPlace(const M0& QRx, const V0& beta, M2& m2)
        { QR_SolveInPlace(QRx,beta,m2); }
        template <class M0, class V0, class M1, class M2>
        static void solve(const M0& QRx, const V0& beta, const M1& m1, M2& m2)
        { QR_Solve(QRx,beta,m1,m2); }
    };
    template <>
    struct QRHelper<true,true>
    {
        template <class M0, class V0, class M2>
        static void makeInverse(const M0& QRx, const V0& beta, M2& m2)
        {
            typename M2::transpose_type m2t = m2.transpose();
            QR_Inverse(QRx,beta,m2t);
        }
        template <class M0, class V0, class M2>
        static void solveInPlace(const M0& QRx, const V0& beta, M2& m2)
        { QR_SolveTransposeInPlace(QRx,beta,m2); }
        template <class M0, class V0, class M1, class M2>
        static void solve(const M0& QRx, const V0& beta, const M1& m1, M2& m2)
        { QR_SolveTranspose(QRx,beta,m1,m2); }
    };
    template <bool istrans>
    struct QRHelper<false,istrans>
    {
        template <class M0, class V0, class M2>
        static void makeInverse(const M0& , const V0& , M2& ) {}
        template <class M0, class V0, class M2>
        static void solveInPlace(const M0& , const V0& , M2& ) {}
        template <class M0, class V0, class M1, class M2>
        static void solve(const M0& , const V0& , const M1& , M2& ) {}
    };

    template <class M>
    struct QRD_Impl<true,M>
    {
        enum { istrans = M::_colmajor < M::_rowmajor };
        typedef typename QRD<M>::qr_type::view_type qrx_type;
        enum { cs = typename QRD<M>::qr_type::_colsize };
        enum { rs = typename QRD<M>::qr_type::_colsize };

        QRD_Impl(const M& A, bool ) : QRx( SmallQRx.view() )
        {
            TMVStaticAssert(M::_colsize != UNKNOWN);
            TMVStaticAssert(M::_rowsize != UNKNOWN);
            TMVStaticAssert(M::_colsize == istrans ? rs : cs);
            TMVStaticAssert(M::_rowsize == istrans ? cs : rs);
            QRx = Maybe<istrans>::transpose(A);
        }

        template <class M1, class M2>
        void solve(const M1& m1, M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::solve(QRx,beta,m1,m2); 
        }

        template <class M1, class M2>
        void solveTranspose(const M1& m1, M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,!istrans>::solve(QRx,beta,m1,m2); 
        }

        template <class M2>
        void solveInPlace(M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::solveInPlace(QRx,beta,m2);
        }

        template <class M2>
        void solveTransposeInPlace(M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,!istrans>::solveInPlace(QRx,beta,m2);
        }

        template <class M2>
        void makeInverse(M2& minv) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            QRHelper<isvalid,istrans>::makeInverse(QRx,beta,minv); 
        }

        typename QRD<M>::qr_type SmallQRx;
        qrx_type QRx;
        beta_type beta;
        int detq;
    };
    
    template <class M>
    struct QRD_Impl<false,M>
    {
        enum { rmorcm = M::_rowmajor || M::_colmajor };
        typedef typename QRD<M>::qr_type::view_type qrx_type;

        QRD_Impl(const M& A, bool _inplace) :
            // if A is short, need to transpose
            istrans(A.colsize() < A.rowsize()),
            // inplace only if matrix is rowmajor or colmajor
            inplace(rmorcm && _inplace),
            // Aptr is the pointer to new storage if any
            Aptr( inplace ? 0 : A.rowsize()*A.rowsize() ),
            // QRx views this memory as the QR matrix
            QRx(
                inplace ? A.nonConst().ptr() : Aptr.get() , // ptr
                istrans ? A.rowsize() : A.colsize() ,  // colsize
                istrans ? A.colsize() : A.rowsize() ,  // rowsize
                inplace ? (istrans ? A.stepi() : A.stepj()) : 1 , // stepi
                ( inplace ? (istrans ? A.stepj() : A.stepi()) : 
                  int(A.rowsize()) ) // stepj
            )
            {
                if (!inplace) {
                    if (istrans) QRx = A.transpose();
                    else QRx = A;
                } else {
                    Maybe<M::_conj>::conjself(QRx);
                }
            }

        // Unlike for LU, istrans isn't known at compile time, at least
        // in the non-small Impl.
        // So we moved the following functions into the Impl.
        // The above small Impl uses the istrans parameter of QRHelper at
        // compile time, and here we check for it at run time.
        
        template <class M1, class M2>
        void solve(const M1& m1, M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans) {
                QRHelper<isvalid,true>::solve(QRx,beta,m1,m2);
            } else {
                QRHelper<isvalid,false>::solve(QRx,beta,m1,m2);
            }
        }

        template <class M1, class M2>
        void solveTranspose(const M1& m1, M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans) {
                QRHelper<isvalid,false>::solve(QRx,beta,m1,m2);
            } else {
                QRHelper<isvalid,true>::solve(QRx,beta,m1,m2);
            }
        }

        template <class M2>
        void solveInPlace(M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (pimpl->istrans) {
                QRHelper<isvalid,true>::solveInPlace(QRx,beta,m2);
            } else {
                QRHelper<isvalid,false>::solveInPlace(QRx,beta,m2);
            }
        }

        template <class M2>
        void solveTransposeInPlace(M2& m2) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (pimpl->istrans) {
                QRHelper<isvalid,false>::solveInPlace(QRx,beta,m2);
            } else {
                QRHelper<isvalid,true>::solveInPlace(QRx,beta,m2);
            }
        }


        template <class M2>
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
        {
            const bool isvalid = M::isreal || M2::iscomplex;
            if (istrans) {
                QRHelper<isvalid,true>::makeInverse(QRx,beta,minv);
            } else {
                QRHelper<isvalid,false>::makeInverse(QRx,beta,minv);
            }
        }

        const bool istrans;
        const bool inplace;
        AlignedArray<typename M::value_type> Aptr;
        qrx_type QRx;
        beta_type beta;
        int detq;
    };

    template <class M>
    QRD<M>::QRD(const M& A, bool inplace) :
        pimpl(new QRD_Impl<small,M>(A,inplace)) 
    {
        QR_Decompose(pimpl->QRx,pimpl->beta,pimpl->detq);
    }

    template <class M>
    QRD<M>::QRD(const QRD<M>& rhs) : pimpl(rhs.pimpl.release()) {}

    template <class M>
    QRD<M>::~QRD() {}

    template <class M> template <class M1, class M2>
    void QRD<M>::solve(
        const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.colsize() == colsize());
        TMVAssert(m2.colsize() == rowsize());
        pimpl->solve(m1,m2);
    }

    template <class M> template <class V1, class V2>
    void QRD<M>::solve(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V1::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVAssert(v1.size() == colsize());
        TMVAssert(v2.size() == rowsize());
        pimpl->solve(v1,v2);
    }

    template <class M> template <class M2>
    void QRD<M>::solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(m2.colsize() == colsize());
        pimpl->solveInPlace(m2);
    }

    template <class M> template <class V2>
    void QRD<M>::solveInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(v2.size() == rowsize());
        pimpl->solveInPlace(v2);
    }

    template <class M> template <class M1, class M2>
    void QRD<M>::solveTranspose(
        const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M::_rowsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.colsize() == rowsize());
        TMVAssert(m2.colsize() == colsize());
        pimpl->solveTranspose(m1,m2);
    }

    template <class M> template <class V1, class V2>
    void QRD<M>::solveTranspose(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V1::_size,M::_rowsize>::same));
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVAssert(v1.size() == rowsize());
        TMVAssert(v2.size() == colsize());
        pimpl->solveTranspose(v1,v2);
    }

    template <class M> template <class M2>
    void QRD<M>::solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
    {
        TMVStaticAssert((Sizes<M2::_colsize,M::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M::_rowsize>::same));
        TMVAssert(m2.colsize() == rowsize());
        TMVAssert(m2.colsize() == colsize());
        pimpl->solveTransposeInPlace(m2);
    }

    template <class M> template <class V2>
    void QRD<M>::solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
    {
        TMVStaticAssert((Sizes<V2::_size,M::_colsize>::same));
        TMVStaticAssert((Sizes<V2::_size,M::_rowsize>::same));
        TMVAssert(v2.size() == colsize());
        TMVAssert(v2.size() == rowsize());
        pimpl->solveTransposeInPlace(v2);
    }

    template <class M>
    typename M::value_type QRD<M>::det() const
    { return typename M::real_type(pimpl->detq) * getR().det(); }                  
    template <class M>
    typename M::real_type QRD<M>::logDet(typename M::value_type* sign) const
    {
        typename M::real_type ret = getR().logDet(sign);
        if (sign) *sign *= typename M::real_type(pimpl->detq);
        return ret;
    }                  

    template <class M>
    bool QRD<M>::isSingular() const 
    { return getR().isSingular(); }

    template <class M> template <class M2>
    void QRD<M>::makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
    {
        TMVStaticAssert((Sizes<M::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M::_rowsize,M2::_colsize>::same));
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        pimpl->makeInverse(minv);
    }

    template <class M> template <class M2>
    void QRD<M>::makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
    {
        TMVAssert(ata.rowsize() == rowsize());
        TMVAssert(ata.colsize() == rowsize());
        QR_InverseATA(pimpl->QRx,pimpl->beta,ata);
    }

    template <class M>
    bool QRD<M>::isTrans() const 
    { return QRD_Impl<small,M>::istrans; }

    template <class M>
    typename QRD<M>::getq_type QRD<M>::getQ() const 
    { return QRD<M>::getq_type(pimpl->QRx,pimpl->beta); }

    template <class M>
    typename QRD<M>::getr_type QRD<M>::getR() const 
    { return pimpl->QRx.upperTri(); }

    template <class M>
    typename QRD<M>::getqr_type QRD<M>::getQR() const 
    { return pimpl->QRx; }

    template <class M>
    typename QRD<M>::getbeta_type QRD<M>::getBeta() const 
    { return pimpl->beta; }

    template <class M>
    typename M::real_type QRD<M>::condition(RT normInf) const 
    {
        // TODO: This is a placeholder until I write the real function.
        // Make sure to do this before releasing the code!
        // See page 129 of Golub and van Loan.
        //
        // This produces the exact right answer, but it is way too slow!
        // The GvL algorithm is order n^2.  This is order n^3.
        Matrix<T> minv(rowsize(),colsize());
        makeInverse(minv);
        return normInf * minv.normInf();
    }

    template <class M, class M2>
    static bool CheckDecomp(
        const QRD<M>& qrd, const BaseMatrix_Calc<M2>& m, std::ostream* fout)
    {
        typedef typename M2::real_type RT;
        bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
        if (printmat) {
            *fout << "QR:\n";
            *fout << "M = "<<
                (qrd.isTrans()?m.transpose():m.view()) << std::endl;
            *fout << "Q = "<<qrd.getQ()<<std::endl;
            *fout << "R = "<<qrd.getR()<<std::endl;
        }
        typename M::copy_type qr = qrd.getQ()*qrd.getR();
        if (printmat) {
            *fout << "QR = "<<qr<<std::endl;
        }
        RT nm = pimpl->istrans ? Norm(qr-m.transpose()) : Norm(qr-m);
        nm /= Norm(qrd.getQ())*Norm(qrd.getR());
        RT kappa = qrd.condition(m.normInf());
        if (fout) {
            *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<" <? ";
            *fout << kappa<<"*"<<RT(m.colsize())<<"*"<<TMV_Epsilon<RT>();
            *fout << " = "<<kappa*RT(m.colsize())*TMV_Epsilon<RT>()<<std::endl;
        }
        return nm < kappa*RT(m.colsize())*TMV_Epsilon<RT>();
    }

    template <class M>
    size_t QRD<M>::colsize() const
    { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

    template <class M>
    size_t QRD<M>::rowsize() const
    { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }


} // namespace mv


#endif
