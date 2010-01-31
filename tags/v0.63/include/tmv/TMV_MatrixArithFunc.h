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


#ifndef TMV_MatrixArithFunc_H
#define TMV_MatrixArithFunc_H

#include "tmv/TMV_BaseMatrix.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, class T, class Ta, class Tx> 
    void MultMV(
        const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
        const VectorView<T>& y);

    // A *= alpha
    template <class T> 
    void MultXM(const T alpha, const MatrixView<T>& A);

    // B += alpha * A
    template <class T, class Ta> 
    void AddMM(const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B);

    // C = alpha * A + beta * B
    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C);

    // C (+)= alpha * A * B
    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        const MatrixView<T>& C);

    // A (+)= alpha * x * yT
    template <bool add, class T, class Tx, class Ty> 
    void Rank1Update(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y, 
        const MatrixView<T>& A);

    template <class T, class Ta> 
    void ElementProd(
        const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B);

    template <class T, class Ta, class Tb> 
    void AddElementProd(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, const MatrixView<T>& C);

    template <class T> 
    class MatrixComposite : public GenMatrix<T>
    {
    public:

        inline MatrixComposite() : itsm(0) {}
        inline MatrixComposite(const MatrixComposite<T>&) : itsm(0) {}
        virtual inline ~MatrixComposite() {}
        const T* cptr() const;
        int stepi() const;
        int stepj() const;
        inline ConjType ct() const { return NonConj; }
        inline bool canLinearize() const { return true; }
        size_t ls() const;

    private:
        mutable auto_array<T> itsm;
    };

    // Specialize allowed complex combinations:
    template <bool add, class T, class Ta, class Tx> 
    inline void MultMV(
        const T alpha, const GenMatrix<Ta>& A,
        const GenVector<Tx>& x, const VectorView<CT>& y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <class T> 
    inline void MultXM(const T alpha, const MatrixView<CT>& A)
    { MultXM(CT(alpha),A); }

    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A, const MatrixView<CT>& B)
    { AddMM(CT(alpha),A,B); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, const MatrixView<CT>& C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <class T> 
    inline void AddMM(
        const CT alpha, const GenMatrix<CT>& A,
        const CT beta, const GenMatrix<T>& B, const MatrixView<CT>& C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, const MatrixView<CT>& C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Tx, class Ty> 
    inline void Rank1Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, const MatrixView<CT>& A)
    { Rank1Update<add>(CT(alpha),x,y,A); }

    template <class T, class Ta> 
    inline void ElementProd(
        const T alpha, const GenMatrix<Ta>& A, const MatrixView<CT>& B)
    { ElementProd(CT(alpha),A,B); }

    template <class T, class Ta, class Tb> 
    inline void AddElementProd(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, const MatrixView<CT>& C)
    { AddElementProd(CT(alpha),A,B,C); }

    template <class T> 
    inline void AddElementProd(
        const CT alpha, const GenMatrix<CT>& A,
        const GenMatrix<T>& B, const MatrixView<CT>& C)
    { AddElementProd(alpha,B,A,C); }

    // Specialize disallowed complex combinations:
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMV(
        const CT , const GenMatrix<Ta>& ,
        const GenVector<Tb>& , const VectorView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T> 
    inline void MultXM(const CT , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta> 
    inline void AddMM(const CT , const GenMatrix<Ta>& , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , const GenMatrix<Tb>& ,
        const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void Rank1Update(
        const CT , const GenVector<Ta>& ,
        const GenVector<Tb>& , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
