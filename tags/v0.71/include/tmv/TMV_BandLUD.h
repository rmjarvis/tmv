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
    // U should be NonUnitDiag.
    template <class T> 
    void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrixView<T> U, ptrdiff_t* P);

    // Do the decomposition in compressed form.
    template <class T> 
    void LU_Decompose(BandMatrixView<T> LUx, ptrdiff_t* P, ptrdiff_t Anhi);

    class Permutation;

    template <class T> 
    void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrixView<T> U, Permutation& P);

    template <class T> 
    void LU_Decompose(BandMatrixView<T> LUx, Permutation& P, ptrdiff_t Anhi);
    
    template <class T, int A1> 
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrix<T,A1>& L,
        BandMatrixView<T> U, Permutation& P)
    { LU_Decompose(A,L.view(),U,P); }

    template <class T, int A2> 
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrix<T,A2>& U, Permutation& P)
    { LU_Decompose(A,L,U.view(),P); }

    template <class T, int A1, int A2> 
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrix<T,A1>& L,
        BandMatrix<T,A2>& U, Permutation& P)
    { LU_Decompose(A,L.view(),U.view(),P); }

    template <class T, int A1> 
    inline void LU_Decompose(BandMatrix<T,A1>& LUx, Permutation& P, ptrdiff_t Anhi)
    { LU_Decompose(LUx.view(),P,Anhi); }
    
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
        void doLDivEq(MatrixView<T1> m) const;
        template <class T1> 
        void doRDivEq(MatrixView<T1> m) const;
        template <class T1, class T2> 
        void doLDiv(
            const GenMatrix<T1>& m, MatrixView<T2> x) const;
        template <class T1, class T2> 
        void doRDiv(
            const GenMatrix<T1>& m, MatrixView<T2> x) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(MatrixView<T1> minv) const;
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool isSingular() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Access Decomposition
        //

        bool isTrans() const;
        LowerTriMatrix<T,UnitDiag> getL() const;
        ConstBandMatrixView<T> getU() const;
        const GenBandMatrix<T>& getLU() const;
        const Permutation& getP() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

        struct BandLUDiv_Impl;
        std::auto_ptr<BandLUDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        BandLUDiv(const BandLUDiv<T>&);
        BandLUDiv<T>& operator=(const BandLUDiv<T>&);

    };

} // namespace mv

#endif
