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
// LU Decomposition.
//
// The basics of LU decomposition for band matrices are the same as 
// for regular matrices.  However, there are a few wrinkles about doing
// it efficiently.  
//
// We leave the details to the comments in TMV_BandLUDiv.cpp, but 
// the main difference for the routines in this file is that L can
// be stored in a lower band matrix with m.nlo() subdiagonals.  
// However, U needs m.nlo() + m.nhi() superdiagonals for its storage.
//
//


#ifndef TMV_BandLUD_H
#define TMV_BandLUD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A into P L U 
    // L, U, and P must have the same size as A.
    // L should be UnitDiag.
    // U shoudl be NonUnitDiag.
    template <class T> 
    void LU_Decompose(
        const GenBandMatrix<T>& A, const LowerTriMatrixView<T>& L,
        const BandMatrixView<T>& U, int* P);

    template <class T> 
    class BandLUDiv : public Divider<T> 
    {

    public :

        //
        // Constructors
        //

        BandLUDiv(const GenBandMatrix<T>& A, bool _inplace);
        BandLUDiv(const AssignableToBandMatrix<T>& A);
        ~BandLUDiv();

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
        LowerTriMatrix<T,UnitDiag> getL() const;
        ConstBandMatrixView<T> getU() const;
        const GenBandMatrix<T>& getLU() const;
        const int* getP() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

        TMV_DEPRECATED(bool IsTrans() const)
        { return isTrans(); }
        typedef LowerTriMatrix<T,UnitDiag> UnitLowerTri;
        TMV_DEPRECATED(UnitLowerTri GetL() const)
        { return getL(); }
        TMV_DEPRECATED(ConstBandMatrixView<T> GetU() const)
        { return getU(); }
        TMV_DEPRECATED(const GenBandMatrix<T>& GetLU() const)
        { return getLU(); }
        TMV_DEPRECATED(const int* GetP() const)
        { return getP(); }

    private :

        struct BandLUDiv_Impl;
        std::auto_ptr<BandLUDiv_Impl> pimpl;

        size_t colsize() const;
        size_t rowsize() const;

    private :

        BandLUDiv(const BandLUDiv<T>&);
        BandLUDiv<T>& operator=(const BandLUDiv<T>&);

    };

} // namespace mv

#endif
