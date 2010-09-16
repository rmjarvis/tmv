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


#include "TMV_SymBandCHDiv.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_BandMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#else
#define TRI_DIV_BLOCKSIZE 64
#endif

    template <class T> 
    static void DoCHInverse(const SymMatrixView<T>& sinv, int nlo)
    {
        TMVAssert(sinv.isherm());
        // inv = (L Lt)^-1 = Lt^-1 L^-1
        LowerTriMatrixView<T> L = sinv.lowerTri();
        TriInverse(L.transpose(),nlo);
        sinv = L.adjoint() * L;
    }

    template <class T> 
    static void SimpleLDLt_AddXtDX(
        const SymMatrixView<T>& sinv, 
        const GenMatrix<T>& X, const GenDiagMatrix<T>& D)
    {
        const int N = D.size();
        if (N==1) {
            sinv += D(0) * X.row(0).conjugate() ^ X.row(0);
        } else {
            int No2 = N/2;
            SimpleLDLt_AddXtDX(sinv,X.rowRange(0,No2),D.subDiagMatrix(0,No2));
            SimpleLDLt_AddXtDX(sinv,X.rowRange(No2,N),D.subDiagMatrix(No2,N));
        }
    }

    template <class T> 
    static void SimpleLDLt_CombineInverse(const SymMatrixView<T>& sinv)
    {
        // This is basically the same algorithm as the LDLt_CombineInverse
        // in TMV_SymLDLInverse.cpp for the Bunch-Kaufman inverse.
        // The difference, and why this is the "simple" version, is that
        // D really is a regular DiagMatrix - not the screwy block diagonal
        // version that we had to deal with there.

        const int N = sinv.size();
        if (N > 1) {
            int No2 = N/2;
            MatrixView<T> X10 = sinv.subMatrix(No2,N,0,No2);
            LowerTriMatrixView<T> X11 = 
                sinv.lowerTri(UnitDiag).subTriMatrix(No2,N);
            DiagMatrixView<T> D(sinv.diag());

            SimpleLDLt_CombineInverse(sinv.subSymMatrix(0,No2));
            SimpleLDLt_AddXtDX(
                sinv.subSymMatrix(0,No2),X10,D.subDiagMatrix(No2,N));
            X10 = D.subDiagMatrix(No2,N) * X10;
            X10 = X11.adjoint() * X10;
            SimpleLDLt_CombineInverse(sinv.subSymMatrix(No2,N));
        }
    }

    template <class T> 
    static void DoLDL_Inverse(const SymMatrixView<T>& sinv)
    {
        LowerTriMatrixView<T> L = sinv.lowerTri(UnitDiag);
        DiagMatrixView<T> D(sinv.diag());
        TriInverse(L.transpose(),1);
        D.realPart().invertSelf();
        SimpleLDLt_CombineInverse(sinv);
    }

    template <class T, class T1> 
    void CH_Inverse(
        const GenSymBandMatrix<T>& LLx, const SymMatrixView<T1>& sinv) 
    {
        TMVAssert(LLx.size() == sinv.size());

#ifdef XDEBUG
        //cout<<"Start CH_Inverse\n";
        //cout<<"LLx = "<<TMV_Text(LLx)<<"  "<<LLx<<endl;
        //cout<<"sinv = "<<TMV_Text(sinv)<<"  "<<sinv<<endl;
        Matrix<T> A(LLx.size(),LLx.size());
        BandMatrix<T> L = LLx.lowerBand();
        A = L * L.adjoint();
#endif

        if (sinv.size() > 0) {
            if (!( sinv.iscm() || sinv.isrm() )) {
                HermMatrix<T,Lower,ColMajor> temp(sinv.size());
                CH_Inverse(LLx,temp.view());
                sinv = temp;
            } else {
                sinv = LLx;
                DoCHInverse(sinv,LLx.nlo());
            }
        }

#ifdef XDEBUG
        //cout<<"Done CH_Inverse\n";
        Matrix<T1> eye = A * sinv;
        TMV_RealType(T) kappa = Norm(A) * Norm(sinv);
        if (Norm(eye-T(1)) > 0.0001*kappa*sinv.size()) {
            cerr<<"A = "<<A<<endl;
            cerr<<"sinv = "<<sinv<<endl;
            cerr<<"A*sinv = "<<A*sinv<<endl;
            cerr<<"sinv*A = "<<sinv*A<<endl;
            cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
        TMVAssert(sinv.isHermOK());
    }

    template <class T, class T1> 
    void LDL_Inverse(
        const GenSymBandMatrix<T>& LLx, const SymMatrixView<T1>& sinv) 
    {
        TMVAssert(LLx.size() == sinv.size());
        TMVAssert(LLx.nlo() == 1);

#ifdef XDEBUG
        //cout<<"Start CH_Inverse\n";
        //cout<<"LLx = "<<TMV_Text(LLx)<<"  "<<LLx<<endl;
        //cout<<"sinv = "<<TMV_Text(sinv)<<"  "<<sinv<<endl;
        Matrix<T> A(LLx.size(),LLx.size());
        BandMatrix<T> L = LLx.lowerBand();
        L.diag().setAllTo(T(1));
        DiagMatrix<T> D(LLx.diag());
        A = L;
        A *= D;
        A *= L.adjoint();
#endif

        if (sinv.size() > 0) {
            if (!( sinv.iscm() || sinv.isrm() )) {
                HermMatrix<T,Lower,ColMajor> temp(sinv.size());
                LDL_Inverse(LLx,temp.view());
                sinv = temp;
            } else {
                sinv = LLx;
                DoLDL_Inverse(sinv);
            }
        }

#ifdef XDEBUG
        //cout<<"Done LDL_Inverse\n";
        Matrix<T1> eye = A * sinv;
        TMV_RealType(T) kappa = Norm(A) * Norm(sinv);
        if (Norm(eye-T(1)) > 0.0001*kappa*sinv.size()) {
            cerr<<"A = "<<A<<endl;
            cerr<<"sinv = "<<sinv<<endl;
            cerr<<"A*sinv = "<<A*sinv<<endl;
            cerr<<"sinv*A = "<<sinv*A<<endl;
            cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
        TMVAssert(sinv.isHermOK());
    }


#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandCHInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


