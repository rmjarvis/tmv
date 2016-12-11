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


#ifndef TMV_TriMatrixArithFunc_H
#define TMV_TriMatrixArithFunc_H

#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // A = alpha * A
    template <typename T>
    void MultXM(const T alpha, UpperTriMatrixView<T> A);
    template <typename T>
    inline void MultXM(const T alpha, LowerTriMatrixView<T> A)
    { MultXM(alpha,A.transpose()); }

    // B += alpha * A
    template <typename T, typename Ta>
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        UpperTriMatrixView<T> B);
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        LowerTriMatrixView<T> B)
    { AddMM(alpha,A.transpose(),B.transpose()); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        MatrixView<T> B)
    { AddMM(alpha,A,B.upperTri()); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        MatrixView<T> B)
    { AddMM(alpha,A.transpose(),B.lowerTri().transpose()); }
    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B,
        UpperTriMatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B,
        LowerTriMatrixView<T> C)
    { AddMM(alpha,A.transpose(),beta,B.transpose(),C.transpose()); }
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(alpha,A.transpose(),beta,B.transpose(),C.transpose()); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B.transpose(),alpha,A.transpose(),C.transpose()); }
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    void ElemMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<T> C)
    { ElemMultMM<add>(alpha,A.transpose(),B.transpose(),C.transpose()); }

    template <typename T>
    class UpperTriMatrixComposite : public GenUpperTriMatrix<T>
    {
    public:

        inline UpperTriMatrixComposite() {}
        inline UpperTriMatrixComposite(const UpperTriMatrixComposite<T>&) {}
        virtual inline ~UpperTriMatrixComposite() {}
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;
        inline ConjType ct() const { return NonConj; }

    private:
        mutable AlignedArray<T> itsm;
    };

    template <typename T>
    class LowerTriMatrixComposite : public GenLowerTriMatrix<T>
    {
    public:

        inline LowerTriMatrixComposite() {}
        inline LowerTriMatrixComposite(const LowerTriMatrixComposite<T>&) {}
        virtual inline ~LowerTriMatrixComposite() {}
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;
        inline ConjType ct() const { return NonConj; }

    private:
        mutable AlignedArray<T> itsm;
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T>
    inline void MultXM(const T alpha, UpperTriMatrixView<CT> A)
    { MultXM(CT(alpha),A); }
    template <typename T>
    inline void MultXM(const T alpha, LowerTriMatrixView<CT> A)
    { MultXM(CT(alpha),A); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        UpperTriMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        LowerTriMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B,
        UpperTriMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B,
        LowerTriMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenUpperTriMatrix<CT>& A,
        const CT beta, const GenUpperTriMatrix<T>& B,
        UpperTriMatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void MultXM(const CT , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void MultXM(const CT , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& ,
        const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& ,
        const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
