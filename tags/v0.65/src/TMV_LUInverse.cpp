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
#include "TMV_LUDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    template <class T> 
    static void NonLapLUInverse(const MatrixView<T>& minv)
    {
        // m = P L U
        // m^-1 = U^-1 L^-1 Pt
        UpperTriMatrixView<T> U = minv.upperTri();
        LowerTriMatrixView<T> L = minv.lowerTri(UnitDiag);
        U.invertSelf();
        L.invertSelf();
        //cout<<"U*Uinv = "<<U0*U<<endl;
        //cout<<"L*Linv = "<<L0*L<<endl;
        //cout<<"LUUinvLinv = "<<L0*U0*U*L<<endl;
        minv = U*L;
        //cout<<"L*U = "<<L*U<<endl;
        //cout<<"minv = "<<minv<<endl;
        //cout<<"LUminv = "<<L*U*minv<<endl;
        // Do Pt back in LU_Inverse
    }

#ifdef ALAP
    template <class T> 
    static inline void LapLUInverse(const MatrixView<T>& minv)
    { NonLapLUInverse(minv); }
#ifdef INST_DOUBLE
    template <> 
    void LapLUInverse(const MatrixView<double>& minv)
    {
        TMVAssert(minv.isSquare());
        TMVAssert(minv.iscm());

        int n = minv.colsize();
        int lda = minv.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(dgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"dgetri");
#endif
    }
    template <> 
    void LapLUInverse(const MatrixView<std::complex<double> >& minv)
    {
        TMVAssert(minv.isSquare());
        TMVAssert(minv.iscm());

        int n = minv.colsize();
        int lda = minv.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<double> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<double> > work(1);
        work.get()[0] = 0.;
        LAPNAME(zgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(zgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("zgetri");
#else
        LAP_Results(int(TMV_REAL(work[0])),n,n,lwork,"zgetri");
#endif
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapLUInverse(const MatrixView<float>& minv)
    {
        TMVAssert(minv.isSquare());
        TMVAssert(minv.iscm());

        int n = minv.colsize();
        int lda = minv.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.;
        LAPNAME(sgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(sgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"sgetri");
#endif
    }
    template <> 
    void LapLUInverse(const MatrixView<std::complex<float> >& minv)
    {
        TMVAssert(minv.isSquare());
        TMVAssert(minv.iscm());

        int n = minv.colsize();
        int lda = minv.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<float> > work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<std::complex<float> > work(1);
        work.get()[0] = 0.;
        LAPNAME(cgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(TMV_REAL(work[0]));
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(cgetri) (LAPCM LAPV(n),LAPP(minv.ptr()),LAPV(lda),
                         LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("cgetri");
#else
        LAP_Results(int(TMV_REAL(work[0])),n,n,lwork,"cgetri");
#endif
    }
#endif // FLOAT
#endif // ALAP

    template <class T, class T1> 
    void LU_Inverse(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& minv)
    {
        //std::cout<<"Start LU_Inverse"<<std::endl;
        //std::cout<<"LUx = "<<LUx<<std::endl;
        //std::cout<<"minv = "<<minv<<std::endl;
        TMVAssert(LUx.isSquare());
        TMVAssert(minv.isSquare());
        TMVAssert(minv.colsize() == LUx.colsize());
#ifdef XDEBUG
        Matrix<T> m = LUx.lowerTri(UnitDiag) * LUx.upperTri();
        m.reversePermuteRows(P);
#ifdef ALAP
        Matrix<T,ColMajor> minv2(minv.colsize(),minv.colsize());
        minv2 = LUx;
        NonLapLUInverse(minv2.view());
        minv2.reversePermuteCols(P);
#endif
#endif

        if (minv.colsize() > 0) {
            if ( !(minv.iscm()
#ifndef ALAP
                   || minv.isrm()
#endif
            )) {
                Matrix<T,ColMajor> temp(minv.colsize(),minv.colsize());
                LU_Inverse(LUx,P,temp.view());
                minv = temp;
            } else {
                minv = LUx;
#ifdef ALAP
                LapLUInverse(minv);
#else
                NonLapLUInverse(minv);
#endif
                //std::cout<<"before permute: "<<minv<<std::endl;
                minv.reversePermuteCols(P);
                //std::cout<<"after permute: "<<minv<<std::endl;
            }
        }

#ifdef XDEBUG
        TMV_RealType(T) normdiff = Norm(m*minv - T(1));
        TMV_RealType(T) kappa = Norm(m)*Norm(minv);
        if (normdiff > 0.001*kappa*minv.colsize()) {
            cerr<<"LUInverse:\n";
            cerr<<"m = "<<m<<endl;
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"P = ";
            for(size_t i=0;i<LUx.colsize();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"minv = "<<minv<<endl;
            cerr<<"m*minv = "<<m*minv<<endl;
            cerr<<"minv*m = "<<minv*m<<endl;
#ifdef ALAP
            cerr<<"Non-lap inverse = "<<minv2<<endl;
            cerr<<"m*minv2 = "<<m*minv2<<endl;
            cerr<<"minv2*m = "<<minv2*m<<endl;
#endif
            cerr<<"Norm(m*minv - 1) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;

            abort();
        }
#endif
        //std::cout<<"minv => "<<minv<<std::endl;
    }


#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_LUInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

