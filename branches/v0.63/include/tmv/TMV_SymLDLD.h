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
// This file contains the code for doing division of Symmetric/Hermitian
// matrices using LDLt Decomposition.  (Although it is still referred
// to as LU decomposition for the purposes of the DivideUsing method.)
//
// This Decomposition is appropriate for symmetric matrices which may
// not be positive definite (ie. they may have negative eigenvalues).
// The matrix is decomposed into P * L * D * Lt * Pt.
//
// P is a permutation.
// L is a unit diagonal lower triangle matrix.
// D is a pseudo-diagonal matrix - meaning that there may be 2x2 blocks
//   occasionally along the diagonal.
// Lt is L.adjoint for Hermitian matrices, or L.transpose for symmetric.
//
// The determinant of A is just the determinant of D.
//


#ifndef TMV_SymLDLD_H
#define TMV_SymLDLD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseSymBandMatrix.h"

namespace tmv {

    // Decompose A into P * L * D * Lt * Pt
    // where L is lower triangular with unit diagonal,
    // D is symmetric or hermitian tri-diagonal,
    // and P is a permutation.
    // The hermitian-ness of D must be the same as A.
    // L is returned as A.lowerTri(UnitDiag).
    template <class T> 
    void LDL_Decompose(
        const SymMatrixView<T>& A, const SymBandMatrixView<T>& D, int* P);

    template <class T> 
    class SymLDLDiv : public SymDivider<T> 
    {

    public :

        //
        // Constructors
        //

        SymLDLDiv(const GenSymMatrix<T>& A, bool _inplace);
        ~SymLDLDiv();

        //
        // Divider Versions of DivEq and Div
        //

        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;

        template <class T1> 
        void doRDivEq(const MatrixView<T1>& m) const;

        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;

        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;

        template <class T1> 
        void doMakeInverse(const SymMatrixView<T1>& minv) const;

        inline void makeInverse(
            const SymMatrixView<TMV_RealType(T)>& sinv) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(sinv.size() == colsize());
            doMakeInverse(sinv);
        }
        inline void makeInverse(
            const SymMatrixView<TMV_ComplexType(T)>& sinv) const
        {
            TMVAssert(sinv.size() == colsize());
            TMVAssert(sinv.isherm() == isherm());
            TMVAssert(sinv.issym() == issym());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool isSingular() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Access Decomposition
        //

        const ConstLowerTriMatrixView<T> getL() const;
        const BandMatrix<T> getD() const;
        const int* getP() const;
        const GenSymMatrix<T>& getLL() const;
        const GenVector<T>& getxD() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

        TMV_DEPRECATED(ConstLowerTriMatrixView<T> GetL() const)
        { return getL(); }
        TMV_DEPRECATED(BandMatrix<T> GetD() const)
        { return getD(); }
        TMV_DEPRECATED(const int* GetP() const)
        { return getP(); }
        TMV_DEPRECATED(const GenSymMatrix<T>& GetLL() const)
        { return getLL(); }
        TMV_DEPRECATED(const GenVector<T>& GetxD() const)
        { return getxD(); }

    private :
        struct SymLDLDiv_Impl;
        std::auto_ptr<SymLDLDiv_Impl> pimpl;

        size_t colsize() const;
        size_t rowsize() const;
        bool isherm() const;
        bool issym() const;

    private :

        SymLDLDiv(const SymLDLDiv<T>&);
        SymLDLDiv<T>& operator=(const SymLDLDiv<T>&);

    }; // SymLDLDiv

} // namespace tmv

#endif
