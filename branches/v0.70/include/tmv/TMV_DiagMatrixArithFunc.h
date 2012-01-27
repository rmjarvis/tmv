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


#ifndef TMV_DiagMatrixArithFunc_H
#define TMV_DiagMatrixArithFunc_H

#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_VectorArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, class T, class Ta,  class Tx> 
    void MultMV(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // B += alpha * A 
    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, DiagMatrixView<T> B)
    { AddVV(alpha,A.diag(),B.diag()); }

    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    { AddVV(alpha,A.diag(),B.diag()); }

    // C = alpha * A + beta * B
    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const T beta, const GenDiagMatrix<Tb>& B, DiagMatrixView<T> C)
    { AddVV(alpha,A.diag(),beta,B.diag(),C.diag()); }

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const T beta, const GenDiagMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <class T, class Ta> 
    void MultEqMM(
        const T alpha,
        const GenDiagMatrix<Ta>& A, MatrixView<T> B);

    template <bool add, class T, class Ta, class Tb> 
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C);

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        MatrixView<T> C)
    { MultMM<add>(alpha,B,A.transpose(),C.transpose()); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<T> C)
    { MultMV<add>(alpha,A,B.diag(),C.diag()); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<T> C)
    { ElemMultVV<add>(alpha,A.diag(),B.diag(),C.diag()); }

    template <class T> 
    class DiagMatrixComposite : public GenDiagMatrix<T>
    {
    public:
        inline DiagMatrixComposite() : inst(0) {}
        inline DiagMatrixComposite(const DiagMatrixComposite<T>&) : inst(0) {}
        virtual inline ~DiagMatrixComposite() {}

    protected:
        ConstVectorView<T> cdiag() const;

    private:
        mutable std::auto_ptr<const DiagMatrix<T> > inst;
    };

    // Specialize allowed complex combinations:
    template <bool add, class T, class Ta,  class Tx> 
    inline void MultMV(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, DiagMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <class T, class Ta> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const T beta, const GenDiagMatrix<Tb>& B, DiagMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, 
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A, 
        const T beta, const GenDiagMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <class T, class Ta> 
    inline void MultEqMM(const T alpha,
                         const GenDiagMatrix<Ta>& A, MatrixView<CT> B)
    { MultEqMM(CT(alpha),A,B); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, class T, class Ta, class Tb> 
    inline void MultMV(
        const CT , const GenDiagMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenDiagMatrix<Ta>& , 
        const CT , const GenDiagMatrix<Ta>& , DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenDiagMatrix<Ta>& , 
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T, class Ta, class Tb> 
    inline void AddMM(
        const CT , const GenMatrix<Ta>& , 
        const CT , const GenDiagMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, class T, class Ta, class Tb> 
    inline void ElemMultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

}

#undef CT

#endif
