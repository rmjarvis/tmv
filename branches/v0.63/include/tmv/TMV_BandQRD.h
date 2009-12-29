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
// This file contains the code for doing division of BandMatrices using 
// QR Decomposition.
//
// I don't implement the QRP method, since the permutations screw up the
// band structure.  It could be implemented using a similar technique 
// as I used for the BandLUDiv class, but for now if you want to use it 
// you need to copy the band matrix to a regular matrix first.
//


#ifndef TMV_BandQRD_H
#define TMV_BandQRD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    // Decompose A into Q R 
    // where Q is unitary and R is upper banded.
    // R must have R.nhi() >= A.nlo()+A.nhi() and R.nlo() >= 0
    // Also, A must have A.nrows() >= A.ncols()
    template <class T> 
    void QR_Decompose(
        const GenBandMatrix<T>& A, const MatrixView<T>& Q,
        const BandMatrixView<T>& R);

    // The same, but don't return Q
    template <class T> 
    void QR_Decompose(
        const GenBandMatrix<T>& A, const BandMatrixView<T>& R);

    template <class T> 
    class BandQRDiv : public Divider<T> 
    {

    public :

        //
        // Constructors
        //

        BandQRDiv(const GenBandMatrix<T>& A, bool _inplace);
        ~BandQRDiv();

        //
        // Div, DivEq
        //

        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;
        template <class T1> 
        void doRDivEq(const MatrixView<T1>& m) const;
        template <class T1, class T2> 
        void doLDiv(
            const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
        template <class T1, class T2> 
        void doRDiv(
            const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool isSingular() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Access Decomposition
        //

        bool isTrans() const;
        Matrix<T> getQ() const;
        ConstBandMatrixView<T> getR() const;
        const GenBandMatrix<T>& getQR() const;
        const GenVector<T>& getQBeta() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

        TMV_DEPRECATED(bool IsTrans() const)
        { return isTrans(); }
        TMV_DEPRECATED(Matrix<T> GetQ() const)
        { return getQ(); }
        TMV_DEPRECATED(ConstBandMatrixView<T> GetR() const)
        { return getR(); }
        TMV_DEPRECATED(const GenBandMatrix<T>& GetQR() const)
        { return getQR(); }
        TMV_DEPRECATED(const GenVector<T>& GetQbeta() const)
        { return getQBeta(); }

    private :

        struct BandQRDiv_Impl;
        std::auto_ptr<BandQRDiv_Impl> pimpl;

        size_t colsize() const;
        size_t rowsize() const;

    private :

        BandQRDiv(const BandQRDiv<T>&);
        BandQRDiv<T>& operator=(const BandQRDiv<T>&);

    };

} // namespace mv

#endif
