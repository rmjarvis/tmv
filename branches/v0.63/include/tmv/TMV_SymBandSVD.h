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
#ifndef TMV_SymBandSVD_H
#define TMV_SymBandSVD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    // Find EigenValues and EigenVectors of hermitian band matrix, A.
    // For each lambda(i), A V.col(i) = lambda(i) V.col(i).
    // In other words, A * V = V * DiagMatrixViewOf(lambda)
    // Or, A = V * DiagMatrixViewOf(lambda) * V.inverse()
    // Furthermore, since A is hermitian, V.inverse() = V.adjoint().
    // On input, lambda and V must have the same size as A.
    // On output, the lambda are sorted to be increasing in value.
    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A,
        const MatrixView<T>& V, const VectorView<TMV_RealType(T)>& lambda);

    // The same, but don't return V
    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A,
        const VectorView<TMV_RealType(T)>& lambda);

    // Decompose A into U S V
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, const MatrixView<T>& U,
        const DiagMatrixView<TMV_RealType(T)>& S, const MatrixView<T>& V);

    // The same, but don't return U and/or V
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const MatrixView<T>& U, const DiagMatrixView<TMV_RealType(T)>& S);
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const DiagMatrixView<TMV_RealType(T)>& S, const MatrixView<T>& V);
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        const DiagMatrixView<TMV_RealType(T)>& S);

    // Find S such that SS = A with S positive definite.
    // A must be positive definite hermitian.
    template <class T> 
    void SquareRoot(const GenSymBandMatrix<T>& A, const SymMatrixView<T>& S);

    template <class T> 
    class HermBandSVDiv : public SymDivider<T> 
    {

    public :

        //
        // Constructors
        //

        HermBandSVDiv(const GenSymBandMatrix<T>& A);
        ~HermBandSVDiv();

        //
        // Div, DivEq
        //

        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;
        template <class T1> 
        void doRDivEq(const MatrixView<T1>& m) const;
        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m,
                    const MatrixView<T2>& x) const;
        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m, 
                    const MatrixView<T2>& x) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;
        template <class T1> 
        void doMakeInverse(const SymMatrixView<T1>& sinv) const;
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
            TMVAssert(sinv.isherm());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool isSingular() const;
        TMV_RealType(T) norm2() const;
        TMV_RealType(T) condition() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Determine which (if any) S values to zero out
        //

        void thresh(TMV_RealType(T) toler, std::ostream* debugout=0) const;
        void top(int neigen, std::ostream* debugout=0) const;
        int getKMax() const;

        //
        // Access Decomposition
        //

        ConstMatrixView<T> getU() const;
        DiagMatrix<TMV_RealType(T)> getS() const;
        Matrix<T> getV() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

        TMV_DEPRECATED(void Thresh(
                TMV_RealType(T) toler, std::ostream* debugout=0) const)
        { thresh(toler,debugout); }
        TMV_DEPRECATED(void Top(int neigen, std::ostream* debugout=0) const)
        { top(neigen,debugout); }
        TMV_DEPRECATED(int GetKMax() const)
        { return getKMax(); }
        TMV_DEPRECATED(ConstMatrixView<T> GetU() const)
        { return getU(); }
        TMV_DEPRECATED(DiagMatrix<TMV_RealType(T)> GetS() const)
        { return getS(); }
        TMV_DEPRECATED(Matrix<T> GetV() const)
        { return getV(); }

    protected:

        struct HermBandSVDiv_Impl;
        std::auto_ptr<HermBandSVDiv_Impl> pimpl;

        size_t colsize() const;
        size_t rowsize() const;

    private :

        HermBandSVDiv(const HermBandSVDiv<T>&);
        HermBandSVDiv<T>& operator=(const HermBandSVDiv<T>&);

    }; // HermBandSVDiv

    template <class T> 
    class SymBandSVDiv : public SymDivider<T> 
    {

    public :

        //
        // Constructors
        //

        SymBandSVDiv(const GenSymBandMatrix<T>& A);
        ~SymBandSVDiv();

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
        template <class T1> 
        void doMakeInverse(const SymMatrixView<T1>& sinv) const;
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
            TMVAssert(sinv.issym());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool isSingular() const;
        TMV_RealType(T) norm2() const;
        TMV_RealType(T) condition() const;

        //
        // Determine which (if any) S values to zero out
        //

        void thresh(TMV_RealType(T) toler, std::ostream* debugout=0) const;
        void top(int neigen, std::ostream* debugout=0) const;
        int getKMax() const;

        //
        // Access Decomposition
        //

        ConstMatrixView<T> getU() const;
        ConstDiagMatrixView<TMV_RealType(T)> getS() const;
        ConstMatrixView<T> getV() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

        TMV_DEPRECATED(void Thresh(
                TMV_RealType(T) toler, std::ostream* debugout=0) const)
        { thresh(toler,debugout); }
        TMV_DEPRECATED(void Top(int neigen, std::ostream* debugout=0) const)
        { top(neigen,debugout); }
        TMV_DEPRECATED(int GetKMax() const)
        { return getKMax(); }
        TMV_DEPRECATED(ConstMatrixView<T> GetU() const)
        { return getU(); }
        TMV_DEPRECATED(ConstDiagMatrixView<TMV_RealType(T)> GetS() const)
        { return getS(); }
        TMV_DEPRECATED(ConstMatrixView<T> GetV() const)
        { return getV(); }

    protected:

        struct SymBandSVDiv_Impl;
        std::auto_ptr<SymBandSVDiv_Impl> pimpl;

        size_t colsize() const;
        size_t rowsize() const;

    private :

        SymBandSVDiv(const SymBandSVDiv<T>&);
        SymBandSVDiv<T>& operator=(const SymBandSVDiv<T>&);

    }; // SymBandSVDiv

} // namespace tmv;

#endif
