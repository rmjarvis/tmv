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


#ifndef TMV_SymBandMatrixArithFunc_H
#define TMV_SymBandMatrixArithFunc_H

#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_BandMatrixArithFunc.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // A = alpha * A
    template <typename T>
    inline void MultXM(const T alpha, SymBandMatrixView<T> A)
    { MultXM(alpha,A.upperBand()); }

    // B += alpha * A
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        SymBandMatrixView<T> B)
    { AddMM(alpha,A.upperBand(),B.upperBand()); }
    template <typename T, typename Ta>
    void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        BandMatrixView<T> B);
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        MatrixView<T> B)
    { AddMM(alpha,A,BandMatrixViewOf(B,A.nlo(),A.nhi())); }

    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B,
        SymBandMatrixView<T> C)
    { AddMM(alpha,A.upperBand(),beta,B.upperBand(),C.upperBand()); }
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B,
        BandMatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenBandMatrix<Tb>& B, BandMatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B, BandMatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }


    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, BandMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, BandMatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <typename T>
    class SymBandMatrixComposite : public GenSymBandMatrix<T>
    {
    public :
        inline SymBandMatrixComposite() : itsm(0) {}
        inline SymBandMatrixComposite(const SymBandMatrixComposite<T>&) :
            itsm(0) {}
        virtual inline ~SymBandMatrixComposite() {}

        // Definitions are in TMV_MultsBV.cpp
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;
        ptrdiff_t diagstep() const;

        inline SymType sym() const { return Sym; }
        inline UpLoType uplo() const { return Lower; }
        inline ConjType ct() const { return NonConj; }

    private:
        mutable AlignedArray<T> itsm1;
        mutable T* itsm;

        SymBandMatrixComposite<T>& operator=(const SymBandMatrixComposite<T>&);
    };


    template <typename T>
    class SymBandMatrixComposite<CT> :
        public BandMatrixComposite<CT>,
        virtual public AssignableToSymBandMatrix<CT>
    {
    public :
        // Need to respecify that size() is pure virtual here,
        // since GenBandMatrix defines size for compaibility
        // with AssignableToDiagMatrix, etc.
        virtual ptrdiff_t size() const = 0;

        using AssignableToSymMatrix<CT>::sym;
        using AssignableToBandMatrix<CT>::nlo;
        inline ptrdiff_t nhi() const { return nlo(); }
        inline ptrdiff_t colsize() const { return size(); }
        inline ptrdiff_t rowsize() const { return size(); }

        inline void assignToS(SymMatrixView<T> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.sym() == sym());
            this->assignTosB(SymBandMatrixViewOf(m0,nlo()));
            if (m0.size() > nlo()+1)
                BandMatrixViewOf(m0.upperTri()).diagRange(
                    nlo()+1,m0.size()).setZero();
        }
        inline void assignToS(SymMatrixView<CT> m0) const
        {
            TMVAssert(m0.size() == size());
            TMVAssert(m0.sym() == sym());
            this->assignTosB(SymBandMatrixViewOf(m0,nlo()));
            if (m0.size() > nlo()+1)
                BandMatrixViewOf(m0.upperTri()).diagRange(
                    nlo()+1,m0.size()).setZero();
        }

        inline void assignToM(MatrixView<T> m0) const
        { BandMatrixComposite<CT>::assignToM(m0); }
        inline void assignToM(MatrixView<CT> m0) const
        { BandMatrixComposite<CT>::assignToM(m0); }
        inline void assignToB(BandMatrixView<T> m0) const
        { BandMatrixComposite<CT>::assignToB(m0); }
        inline void assignToB(BandMatrixView<CT> m0) const
        { BandMatrixComposite<CT>::assignToB(m0); }
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T>
    inline void MultXM(const T alpha, SymBandMatrixView<CT> A)
    { MultXM(CT(alpha),A); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        SymBandMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        BandMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B,
        SymBandMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B,
        BandMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenBandMatrix<Tb>& B, BandMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B,
        BandMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenSymBandMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenSymBandMatrix<CT>& A,
        const CT beta, const GenSymBandMatrix<T>& B, MatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenSymBandMatrix<CT>& A,
        const CT beta, const GenSymBandMatrix<T>& B,
        BandMatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, BandMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenSymBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenSymBandMatrix<Tb>& B, BandMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenSymBandMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void MultXM(
        const CT , SymBandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& , SymBandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const CT , const GenSymBandMatrix<Tb>& , SymBandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const CT , const GenSymBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const CT , const GenSymBandMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const CT , const GenBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenBandMatrix<Ta>& ,
        const CT , const GenSymBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenSymBandMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& ,
        const GenSymBandMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const GenSymBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenSymBandMatrix<Ta>& ,
        const GenBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenBandMatrix<Ta>& ,
        const GenSymBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
