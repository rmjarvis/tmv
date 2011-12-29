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

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A (input as Q) into Q R
    // where Q is unitary and R is upper triangular.
    template <class T> 
    void QR_Decompose(const MatrixView<T>& Q, const UpperTriMatrixView<T>& R);

    // Decompose A into Q R, but don't return Q.
    // On output, R is returned as A.upperTri()
    template <class T> 
    void QR_Decompose(const MatrixView<T>& A);

    // Given that A0 = Q0 R0
    // Find R1, so that [ A0 ] = Q1 R1
    //                  [ A  ] 
    // On input R is R0; on output it is R1.
    template <class T> 
    void QR_Update(const UpperTriMatrixView<T>& R, const MatrixView<T>& A);

    // This just adds a single row
    template <class T> 
    inline void QR_Update(
        const UpperTriMatrixView<T>& R, const VectorView<T>& z)
    { QR_Update(R,RowVectorViewOf(z)); }

    // The opposite of an Update:
    // Given that [ A0 ] = Q1 R1
    //            [ A  ]
    // Find R0, so that A0 = Q0 R0
    // On input R is R1; on output it is R0.
    template <class T> 
    void QR_Downdate(const UpperTriMatrixView<T>& R, const MatrixView<T>& A);

    // Remove a single row
    template <class T> 
    inline void QR_Downdate(
        const UpperTriMatrixView<T>& R, const VectorView<T>& z)
    { QR_Downdate(R,RowVectorViewOf(z)); }

    template <class T> 
    class PackedQ;

    template <class T> 
    class QRDiv : public Divider<T> 
    {

    public :

        //
        // Constructors
        //

        QRDiv(const GenMatrix<T>& A, bool _inplace);
        ~QRDiv();

        //
        // Div, DivEq
        //

        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;

        template <class T1> 
        void doRDivEq(const MatrixView<T1>& m) const;

        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool isSingular() const;

        //
        // Access Decomposition
        //

        bool isTrans() const;
        PackedQ<T> getQ() const;
        ConstUpperTriMatrixView<T> getR() const;
        const GenMatrix<T>& getQR() const;
        const GenVector<T>& getBeta() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected :

        struct QRDiv_Impl;
        std::auto_ptr<QRDiv_Impl> pimpl;

        int colsize() const;
        int rowsize() const;

    private :

        QRDiv(const QRDiv<T>&);
        QRDiv<T>& operator=(const QRDiv<T>&);

    };

} // namespace mv


#endif
