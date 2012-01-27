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


#ifndef TMV_TriMatrixArithFunc_H
#define TMV_TriMatrixArithFunc_H

#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenVector<Tx>& x, VectorView<T> y);
    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenVector<Tx>& x, VectorView<T> y);

    // A = alpha * A
    template <class T> 
    void MultXM(const T alpha, UpperTriMatrixView<T> A);
    template <class T> 
    inline void MultXM(const T alpha, LowerTriMatrixView<T> A)
    { MultXM(alpha,A.transpose()); }

    // B += alpha * A
    template <class T, class Ta> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        UpperTriMatrixView<T> B);
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        LowerTriMatrixView<T> B)
    { AddMM(alpha,A.transpose(),B.transpose()); }
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        MatrixView<T> B)
    { AddMM(alpha,A,B.upperTri()); }
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        MatrixView<T> B)
    { AddMM(alpha,A.transpose(),B.lowerTri().transpose()); }
    // C = alpha * A + beta * B
    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, 
        UpperTriMatrixView<T> C);
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, 
        LowerTriMatrixView<T> C)
    { AddMM(alpha,A.transpose(),beta,B.transpose(),C.transpose()); }
    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(alpha,A.transpose(),beta,B.transpose(),C.transpose()); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B.transpose(),alpha,A.transpose(),C.transpose()); }
    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<T> C);
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, class T, class Ta, class Tb> 
    void ElemMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<T> C);
    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<T> C)
    { ElemMultMM<add>(alpha,A.transpose(),B.transpose(),C.transpose()); }

    template <class T> 
    class UpperTriMatrixComposite : public GenUpperTriMatrix<T>
    {
    public:

        inline UpperTriMatrixComposite() {}
        inline UpperTriMatrixComposite(const UpperTriMatrixComposite<T>&) {}
        virtual inline ~UpperTriMatrixComposite() {}
        const T* cptr() const;
        int stepi() const;
        int stepj() const;
        inline ConjType ct() const { return NonConj; }

    private:
        mutable AlignedArray<T> itsm;
    };

    template <class T> 
    class LowerTriMatrixComposite : public GenLowerTriMatrix<T>
    {
    public:

        inline LowerTriMatrixComposite() {}
        inline LowerTriMatrixComposite(const LowerTriMatrixComposite<T>&) {}
        virtual inline ~LowerTriMatrixComposite() {}
        const T* cptr() const;
        int stepi() const;
        int stepj() const;
        inline ConjType ct() const { return NonConj; }

    private:
        mutable AlignedArray<T> itsm;
    };

    // Specialize allowed complex combinations:
    template <bool add, class T, class Ta, class Tx> 
    inline void MultMV(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }
    template <bool add, class T, class Ta, class Tx> 
    inline void MultMV(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <class T> 
    inline void MultXM(const T alpha, UpperTriMatrixView<CT> A)
    { MultXM(CT(alpha),A); }
    template <class T> 
    inline void MultXM(const T alpha, LowerTriMatrixView<CT> A)
    { MultXM(CT(alpha),A); }

    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        UpperTriMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        LowerTriMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, 
        UpperTriMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, 
        LowerTriMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const T beta, const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const T beta, const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <class T> 
    inline void AddMM(
        const CT alpha, const GenUpperTriMatrix<CT>& A,
        const CT beta, const GenUpperTriMatrix<T>& B,
        UpperTriMatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenUpperTriMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A, 
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const T alpha, const GenUpperTriMatrix<Ta>& A,
        const GenUpperTriMatrix<Tb>& B, UpperTriMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }
    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const T alpha, const GenLowerTriMatrix<Ta>& A,
        const GenLowerTriMatrix<Tb>& B, LowerTriMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMV(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMV(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T> 
    inline void MultXM(const CT , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T> 
    inline void MultXM(const CT , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta> 
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta> 
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta> 
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta> 
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenUpperTriMatrix<Ta>& ,
        const CT , const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenLowerTriMatrix<Ta>& ,
        const CT , const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , 
        const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , 
        const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenLowerTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenUpperTriMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const CT , const GenUpperTriMatrix<Ta>& , 
        const GenUpperTriMatrix<Tb>& , UpperTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const CT , const GenLowerTriMatrix<Ta>& , 
        const GenLowerTriMatrix<Tb>& , LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
