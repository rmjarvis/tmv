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
        MatrixView<T> V, VectorView<TMV_RealType(T)> lambda);

    // The same, but don't return V
    template <class T> 
    void Eigen(
        const GenSymBandMatrix<T>& A,
        VectorView<TMV_RealType(T)> lambda);

    // Decompose A into U S Vt
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt);

    // The same, but don't return U and/or V
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<TMV_RealType(T)> S);
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt);
    template <class T> 
    void SV_Decompose(
        const GenSymBandMatrix<T>& A, DiagMatrixView<TMV_RealType(T)> S);

    // Find S such that SS = A with S positive definite.
    // A must be positive definite hermitian.
    template <class T> 
    void SquareRoot(const GenSymBandMatrix<T>& A, SymMatrixView<T> S);

    template <class T, int A2> 
    inline void Eigen(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> V, Vector<TMV_RealType(T),A2>& lambda)
    { Eigen(A,V,lambda.view()); }

    template <class T, int A1> 
    inline void Eigen(
        const GenSymBandMatrix<T>& A,
        Matrix<T,A1>& V, VectorView<TMV_RealType(T)> lambda)
    { Eigen(A,V.view(),lambda); }

    template <class T, int A1, int A2> 
    inline void Eigen(
        const GenSymBandMatrix<T>& A,
        Matrix<T,A1>& V, Vector<TMV_RealType(T),A2>& lambda)
    { Eigen(A,V.view(),lambda.view()); }

    template <class T, int A2> 
    inline void Eigen(
        const GenSymBandMatrix<T>& A, Vector<TMV_RealType(T),A2>& lambda)
    { Eigen(A,lambda.view()); }

    template <class T, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U,S,Vt.view()); }

    template <class T, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,U,S.view(),Vt); }

    template <class T, int A2, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, MatrixView<T> U,
        DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U,S.view(),Vt.view()); }

    template <class T, int A1> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, Matrix<T,A1>& U,
        DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt)
    { SV_Decompose(A,U.view(),S,Vt); }

    template <class T, int A1, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, Matrix<T,A1>& U,
        DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U.view(),S,Vt.view()); }

    template <class T, int A1, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, Matrix<T,A1>& U,
        DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,U.view(),S.view(),Vt); }

    template <class T, int A1, int A2, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, Matrix<T,A1>& U,
        DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U.view(),S.view(),Vt.view()); }

    template <class T, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,U,S.view()); }

    template <class T, int A1> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrixView<TMV_RealType(T)> S)
    { SV_Decompose(A,U.view(),S); }

    template <class T, int A1, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,U.view(),S.view()); }

    template <class T, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,S,Vt.view()); }

    template <class T, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,S.view(),Vt); }

    template <class T, int A2, int A3> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A,
        DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,S.view(),Vt.view()); }

    template <class T, int A2> 
    inline void SV_Decompose(
        const GenSymBandMatrix<T>& A, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,S.view()); }

    template <class T, int A2> 
    inline void SquareRoot(const GenSymBandMatrix<T>& A, HermMatrix<T,A2>& S)
    { SquareRoot(A,S.view()); }

    template <class T, int A2> 
    inline void SquareRoot(const GenSymBandMatrix<T>& A, SymMatrix<T,A2>& S)
    { SquareRoot(A,S.view()); }

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
        void doLDivEq(MatrixView<T1> m) const;
        template <class T1> 
        void doRDivEq(MatrixView<T1> m) const;
        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;
        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(MatrixView<T1> minv) const;
        template <class T1> 
        void doMakeInverse(SymMatrixView<T1> sinv) const;
        inline void makeInverse(
            SymMatrixView<TMV_RealType(T)> sinv) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(sinv.size() == colsize());
            doMakeInverse(sinv);
        }
        inline void makeInverse(
            SymMatrixView<TMV_ComplexType(T)> sinv) const
        {
            TMVAssert(sinv.size() == colsize());
            TMVAssert(sinv.isherm());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(MatrixView<T> minv) const;
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
        Matrix<T> getVt() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

        struct HermBandSVDiv_Impl;
        std::auto_ptr<HermBandSVDiv_Impl> pimpl;

        int colsize() const;
        int rowsize() const;

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
        void doLDivEq(MatrixView<T1> m) const;
        template <class T1> 
        void doRDivEq(MatrixView<T1> m) const;
        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;
        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(MatrixView<T1> minv) const;
        template <class T1> 
        void doMakeInverse(SymMatrixView<T1> sinv) const;
        inline void makeInverse(SymMatrixView<TMV_RealType(T)> sinv) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(sinv.size() == colsize());
            doMakeInverse(sinv);
        }
        inline void makeInverse(SymMatrixView<TMV_ComplexType(T)> sinv) const
        {
            TMVAssert(sinv.size() == colsize());
            TMVAssert(sinv.issym());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(MatrixView<T> minv) const;
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
        ConstMatrixView<T> getVt() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

        struct SymBandSVDiv_Impl;
        std::auto_ptr<SymBandSVDiv_Impl> pimpl;

        int colsize() const;
        int rowsize() const;

    private :

        SymBandSVDiv(const SymBandSVDiv<T>&);
        SymBandSVDiv<T>& operator=(const SymBandSVDiv<T>&);

    }; // SymBandSVDiv

} // namespace tmv;

#endif
