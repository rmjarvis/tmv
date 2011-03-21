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

//#define PRINTALGO_LU

#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_LUD.h"

#include "tmv/TMV_BaseMatrix_Tri.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultPV.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_AddMM.h"

namespace tmv {


    template <class T> template <class M2>
    InstLUD<T>::InstLUD(const BaseMatrix<M2>& A, bool _inplace) : 
        base(A,_inplace) {}
    template <class T> 
    InstLUD<T>::InstLUD(const InstLUD<T>& rhs) : base(rhs) {}
    template <class T> 
    InstLUD<T>::~InstLUD() {}

    template <class T> 
    void InstLUD<T>::doSolveInPlace(MatrixView<RT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstLUD<T>::doSolveInPlace(MatrixView<CT> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstLUD<T>::doSolveInPlace(
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solveInPlace(m2);  }
    template <class T> 
    void InstLUD<T>::doSolveInPlace(VectorView<RT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstLUD<T>::doSolveInPlace(VectorView<CT> v2) const 
    { base::solveInPlace(v2);  }
    template <class T> 
    void InstLUD<T>::doSolveInPlace(VectorView<CT,UNKNOWN,true> v2) const 
    { base::solveInPlace(v2);  }

    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(MatrixView<RT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(MatrixView<CT> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solveTransposeInPlace(m2); } 
    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(VectorView<RT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(VectorView<CT> v2) const 
    { base::solveTransposeInPlace(v2); } 
    template <class T> 
    void InstLUD<T>::doSolveTransposeInPlace(
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solveTransposeInPlace(v2); } 

    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<RT>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<CT>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
        MatrixView<CT> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solve(m1,m2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<RT>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<CT>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<CT,UNKNOWN,true>& v1, 
        VectorView<CT> v2) const 
    { base::solve(v1,v2); } 
    template <class T> 
    void InstLUD<T>::doSolve(
        const ConstVectorView<CT,UNKNOWN,true>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solve(v1,v2); } 

    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<RT>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<CT>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
        MatrixView<CT> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1, 
        MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const 
    { base::solveTranspose(m1,m2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<RT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<RT>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<CT>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<CT,UNKNOWN,true>& v1, 
        VectorView<CT> v2) const 
    { base::solveTranspose(v1,v2); }
    template <class T> 
    void InstLUD<T>::doSolveTranspose(
        const ConstVectorView<CT,UNKNOWN,true>& v1, 
        VectorView<CT,UNKNOWN,true> v2) const 
    { base::solveTranspose(v1,v2); }

    template <class T> 
    T InstLUD<T>::det() const
    { return base::det(); }
    template <class T> 
    typename InstLUD<T>::FT InstLUD<T>::logDet(
        typename InstLUD<T>::ZFT* sign) const
    { return base::logDet(sign); }
    template <class T> 
    bool InstLUD<T>::isSingular() const
    { return base::isSingular(); }

    template <class T> 
    void InstLUD<T>::doMakeInverse(MatrixView<RT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstLUD<T>::doMakeInverse(MatrixView<CT> minv) const 
    { base::makeInverse(minv); } 
    template <class T> 
    void InstLUD<T>::doMakeInverse(
        MatrixView<CT,UNKNOWN,UNKNOWN,true> minv) const 
    { base::makeInverse(minv); } 

    template <class T> 
    void InstLUD<T>::doMakeInverseATA(MatrixView<RT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstLUD<T>::doMakeInverseATA(MatrixView<CT> ata) const
    { base::makeInverseATA(ata); }
    template <class T> 
    void InstLUD<T>::doMakeInverseATA(
        MatrixView<CT,UNKNOWN,UNKNOWN,true> ata) const
    { base::makeInverseATA(ata); }

    template <class T> 
    typename InstLUD<T>::RT InstLUD<T>::condition(
        typename InstLUD<T>::RT normInf) const
    { return base::condition(normInf); }

    template <class T> 
    bool InstLUD<T>::preferInPlace() const
    { return base::preferInPlace(); }


#define InstFile "TMV_LUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


