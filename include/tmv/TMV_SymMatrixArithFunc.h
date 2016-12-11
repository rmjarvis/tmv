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


#ifndef TMV_SymMatrixArithFunc_H
#define TMV_SymMatrixArithFunc_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_VectorArithFunc.h"
#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_TriMatrixArithFunc.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // A = alpha * A
    template <typename T>
    inline void MultXM(const T alpha, SymMatrixView<T> A)
    { MultXM(alpha,A.upperTri()); }

    // B += alpha * A
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        SymMatrixView<T> B)
    { AddMM(alpha,A.upperTri(),B.upperTri()); }
    template <typename T, typename Ta>
    void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A, MatrixView<T> B);

    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, SymMatrixView<T> C)
    { AddMM(alpha,A.upperTri(),beta,B.upperTri(),C.upperTri()); }
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenSymMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenSymMatrix<Tb>& B, MatrixView<T> );

    // A (+)= alpha * (x^xT) (or x* if A is Herm)
    template <bool add, typename T, typename Tx>
    void Rank1Update(
        const T alpha, const GenVector<Tx>& x, SymMatrixView<T> A);

    // B (+)= alpha * (A * AT) (or At if B is Herm)
    template <bool add, typename T, typename Ta>
    void RankKUpdate(
        const T alpha, const GenMatrix<Ta>& A, SymMatrixView<T> B);

    template <bool add, typename T, typename Ta>
    void RankKUpdate(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        SymMatrixView<T> B);
    template <bool add, typename T, typename Ta>
    void RankKUpdate(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        SymMatrixView<T> B);

    // These two don't have += forms: they must called explicitly
    // A (+)= alpha * (x^y + y^x) (or x^y* + y^x* is A is Herm)
    template <bool add, typename T, typename Tx, typename Ty>
    void Rank2Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, SymMatrixView<T> A);
    template <bool add, typename T, typename Tx, typename Ty, int A3>
    inline void Rank2Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, SymMatrix<T,A3>& A)
    { Rank2Update(alpha,x,y,A.view()); }
    template <bool add, typename T, typename Tx, typename Ty, int A3>
    inline void Rank2Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, HermMatrix<T,A3>& A)
    { Rank2Update(alpha,x,y,A.view()); }

    // C (+)= alpha * (A * BT + B*AT) (or A*Bt + B*At if C is Herm)
    template <bool add, typename T, typename Ta, typename Tb>
    void Rank2KUpdate(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb, int A3>
    inline void Rank2KUpdate(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrix<T,A3>& C)
    { Rank2KUpdate(alpha,A,B,C.view()); }
    template <bool add, typename T, typename Ta, typename Tb, int A3>
    inline void Rank2KUpdate(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, HermMatrix<T,A3>& C)
    { Rank2KUpdate(alpha,A,B,C.view()); }

    // C (+)= alpha * A * B
    // This also needs to be called explicitly.
    // This is to prevent the programmer from doing this accidentally
    // when alpha * A * B is not necessarily symmetrix/hermitian.
    template <bool add, typename T, typename Ta, typename Tb>
    void SymMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb, int A3>
    inline void SymMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrix<T,A3>& C)
    { SymMultMM(alpha,A,B,C.view()); }
    template <bool add, typename T, typename Ta, typename Tb, int A3>
    inline void SymMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, HermMatrix<T,A3>& C)
    { SymMultMM(alpha,A,B,C.view()); }

    template <typename T>
    class SymMatrixComposite : public GenSymMatrix<T>
    {
    public :
        inline SymMatrixComposite() {}
        inline SymMatrixComposite(const SymMatrixComposite<T>&) {}
        virtual inline ~SymMatrixComposite() {}

        // Definitions are in TMV_MultSV.cpp
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;

        inline SymType sym() const { return Sym; }
        inline UpLoType uplo() const { return Lower; }
        inline ConjType ct() const { return NonConj; }

    private :
        mutable AlignedArray<T> itsm;

        SymMatrixComposite<T>& operator=(const SymMatrixComposite<T>&);
    };

    template <typename T>
    class SymMatrixComposite<CT> :
        public MatrixComposite<CT>,
        virtual public AssignableToSymMatrix<CT>
    {
    public :

        inline ptrdiff_t colsize() const { return this->size(); }
        inline ptrdiff_t rowsize() const { return this->size(); }

        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        { MatrixComposite<CT>::assignToM(m0); }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        { MatrixComposite<CT>::assignToM(m0); }
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T>
    inline void MultXM(const T alpha, SymMatrixView<CT> A)
    { MultXM(CT(alpha),A); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A, SymMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, SymMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenSymMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenSymMatrix<CT>& A,
        const CT beta, const GenSymMatrix<T>& B, MatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenSymMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenSymMatrix<Ta>& A,
        const GenSymMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Tx>
    inline void Rank1Update(
        const T alpha, const GenVector<Tx>& x, SymMatrixView<CT> A)
    { Rank1Update<add>(CT(alpha),x,A); }
    template <bool add, typename T, typename Ta>
    inline void RankKUpdate(
        const T alpha, const GenMatrix<Ta>& A, SymMatrixView<CT> B)
    { RankKUpdate<add>(CT(alpha),A,B); }
    template <bool add, typename T, typename Ta>
    inline void RankKUpdate(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        SymMatrixView<CT> B)
    { RankKUpdate<add>(CT(alpha),A,B); }
    template <bool add, typename T, typename Ta>
    inline void RankKUpdate(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        SymMatrixView<CT> B)
    { RankKUpdate<add>(CT(alpha),A,B); }

    template <bool add, typename T, typename Tx, typename Ty>
    inline void Rank2Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, SymMatrixView<CT> A)
    { Rank2Update<add>(CT(alpha),x,y,A); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void Rank2KUpdate(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrixView<CT> C)
    { Rank2KUpdate<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void SymMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, SymMatrixView<CT> C)
    { SymMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenSymMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void MultXM(const CT , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymMatrix<Ta>& ,
        const CT , const GenSymMatrix<Tb>& , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymMatrix<Ta>& ,
        const CT , const GenSymMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenSymMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenSymMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenSymMatrix<Ta>& ,
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& ,
        const GenSymMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenSymMatrix<Ta>& ,
        const GenSymMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }


} // namespace tmv

#undef CT

#endif
