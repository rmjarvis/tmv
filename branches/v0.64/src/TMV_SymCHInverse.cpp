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


//#define XDEBUG


#include "TMV_Blas.h"
#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    static void NonLapCHInverse(const SymMatrixView<T>& sinv)
    {
        TMVAssert(sinv.isherm());
        // inv = (L Lt)^-1 = Lt^-1 L^-1
        LowerTriMatrixView<T> L = sinv.lowerTri();
        L.invertSelf();
        sinv = L.adjoint() * L;
    }

#ifdef ALAP
    template <class T> 
    static inline void LapCHInverse(const SymMatrixView<T>& sinv)
    { NonLapCHInverse(sinv); }
#ifdef INST_DOUBLE
    template <> 
    void LapCHInverse(const SymMatrixView<double>& sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        LAPNAME(dpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results("dpotri");
    }
    template <> 
    void LapCHInverse(
        const SymMatrixView<std::complex<double> >& sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        LAPNAME(zpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results("zpotri");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapCHInverse(const SymMatrixView<float>& sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        LAPNAME(spotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results("spotri");
    }
    template <> 
    void LapCHInverse(
        const SymMatrixView<std::complex<float> >& sinv)
    {
        TMVAssert(sinv.isherm());
        TMVAssert(sinv.iscm());

        int n = sinv.size();
        int lda = sinv.stepj();
        LAPNAME(cpotri) (
            LAPCM (sinv.uplo() == Upper ? LAPCH_UP : LAPCH_LO),
            LAPV(n),LAPP(sinv.ptr()),LAPV(lda) LAPINFO LAP1);
        LAP_Results("cpotri");
    }
#endif // FLOAT
#endif // LAP

    template <class T, class T1> 
    void CH_Inverse(const GenSymMatrix<T>& LLx, const SymMatrixView<T1>& sinv) 
    {
        TMVAssert(LLx.size() == sinv.size());

#ifdef XDEBUG
        Matrix<T> A = LLx.lowerTri() * LLx.upperTri();
#endif

        if (sinv.size() > 0) {
            if (
#ifdef ALAP
                !isSameType(T(),T1()) || 
#endif
                !(sinv.iscm() || sinv.isrm())) {
                HermMatrix<T,Lower,ColMajor> temp(sinv.size());
                CH_Inverse(LLx,temp.view());
                sinv = temp;
#ifdef ALAP
            } else if (sinv.isrm()) {
                CH_Inverse(LLx.transpose(),sinv.transpose());
#endif
            } else {
                sinv = LLx;
#ifdef ALAP
                LapCHInverse(sinv);
#else
                NonLapCHInverse(sinv);
#endif
            }
        }

#ifdef XDEBUG
        Matrix<T1> eye = A * sinv;
        TMV_RealType(T) kappa = Norm(A) * Norm(sinv);
        if (Norm(eye-T1(1)) > 0.0001*kappa*sinv.size()) {
            cerr<<"A = "<<A<<endl;
            cerr<<"sinv = "<<sinv<<endl;
            cerr<<"A*sinv = "<<A*sinv<<endl;
            cerr<<"sinv*A = "<<sinv*A<<endl;
            cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T1(1))<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
        TMVAssert(sinv.isHermOK());
    }

#define InstFile "TMV_SymCHInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


