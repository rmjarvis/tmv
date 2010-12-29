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

    //
    // LDivEq
    //

    template <class T, class T1> 
    static void NonLapLULDivEq(
        const GenMatrix<T1>& LUx, const MatrixView<T>& m)
    {
        // Solve L U x = m:
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.colsize());
        m /= LUx.lowerTri(UnitDiag);
        m /= LUx.upperTri(NonUnitDiag);
    }


#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <class T, class T1> 
    static inline void LapLULDivEq(
        const GenMatrix<T1>& LUx, const MatrixView<T>& m)
    { NonLapLULDivEq(LUx,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapLULDivEq(
        const GenMatrix<double>& LUx, const MatrixView<double>& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.colsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());
        TMVAssert(LUx.ct()==NonConj);
        TMVAssert(m.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int lda = LUx.stepj();
        int ldb = m.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(dgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("dgetrs");
    }
    template <> 
    void LapLULDivEq(
        const GenMatrix<std::complex<double> >& LUx,
        const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.colsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());
        TMVAssert(LUx.ct()==NonConj);
        TMVAssert(m.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int lda = LUx.stepj();
        int ldb = m.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(zgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("zgetrs");
    }
#endif
#ifdef INST_FLOAT
#ifndef MKL
    template <> 
    void LapLULDivEq(const GenMatrix<float>& LUx, const MatrixView<float>& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.colsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());
        TMVAssert(LUx.ct()==NonConj);
        TMVAssert(m.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int lda = LUx.stepj();
        int ldb = m.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(sgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("sgetrs");
    }
    template <> 
    void LapLULDivEq(
        const GenMatrix<std::complex<float> >& LUx,
        const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.colsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());
        TMVAssert(LUx.ct()==NonConj);
        TMVAssert(m.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int lda = LUx.stepj();
        int ldb = m.stepj();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(cgetrs) (LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("cgetrs");
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T, class T1> 
    void LU_LDivEq(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
    {
        TMVAssert(m.colsize() == LUx.rowsize()); 
        TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef XDEBUG
        Matrix<T> m0(m);
        Matrix<T1> L = LUx.lowerTri(UnitDiag);
        Matrix<T1> U = LUx.upperTri(NonUnitDiag);
        Matrix<T1> LU = L*U;
        LU.reversePermuteRows(P);
        //cout<<"Start LU_LDivEq\n";
        //cout<<"LU = "<<LU<<endl;
        //cout<<"m0 = "<<m0<<endl;
#endif

        m.permuteRows(P); 

#ifdef ALAP
        if (m.iscm() && !m.isconj() && LUx.iscm() && !LUx.isconj())
            LapLULDivEq(LUx,m);
        else 
#endif
            NonLapLULDivEq(LUx,m); 

#ifdef XDEBUG
        //cout<<"m-> "<<m<<endl;
        Matrix<T> mm = LU*m;
        //cout<<"mm = "<<mm<<endl;
        if (!(Norm(mm-m0) <= 0.001*Norm(L)*Norm(U)*Norm(m0))) {
            cerr<<"LU_LDivEq: m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"LU = "<<LU<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"LU*m = "<<mm<<endl;
            abort();
        }
#endif
    }

    //
    // RDivEq Matrix
    //

    template <class T, class T1> 
    static void NonLapLURDivEq(
        const GenMatrix<T1>& LUx, const MatrixView<T>& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.rowsize());
        // m = m (LU)^-1 
        //   = m U^-1 L^-1
        m %= LUx.upperTri(NonUnitDiag);
        m %= LUx.lowerTri(UnitDiag);
    }

#ifdef ALAP
    template <class T, class T1> 
    static inline void LapLURDivEq(
        const GenMatrix<T1>& LUx, const MatrixView<T>& m)
    { NonLapLURDivEq(LUx,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapLURDivEq(
        const GenMatrix<double>& LUx, const MatrixView<double>& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.rowsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());
        TMVAssert(m.ct()==NonConj);
        TMVAssert(LUx.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int lda = LUx.stepj();
        int ldb = m.stepi();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(dgetrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("dgetrs");
    }
    template <> 
    void LapLURDivEq(
        const GenMatrix<std::complex<double> >& LUx,
        const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.rowsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());
        TMVAssert(LUx.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int lda = LUx.stepj();
        int ldb = m.stepi();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(zgetrs) (LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("zgetrs");
    }
#endif
#ifdef INST_FLOAT
#ifndef MKL
    template <> 
    void LapLURDivEq(
        const GenMatrix<float>& LUx, const MatrixView<float>& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.rowsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());
        TMVAssert(m.ct()==NonConj);
        TMVAssert(LUx.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int lda = LUx.stepj();
        int ldb = m.stepi();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(sgetrs) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("sgetrs");
    }
    template <> 
    void LapLURDivEq(
        const GenMatrix<std::complex<float> >& LUx,
        const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(LUx.isSquare());
        TMVAssert(LUx.rowsize() == m.rowsize());
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());
        TMVAssert(LUx.ct()==NonConj);

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int lda = LUx.stepj();
        int ldb = m.stepi();
        AlignedArray<int> ipiv(n);
#ifdef CLAP
        for(int i=0;i<n;++i) ipiv[i] = i;
#else
        for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
        LAPNAME(cgetrs) (LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
                         LAPP(LUx.cptr()),LAPV(lda),
                         LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("cgetrs");
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T, class T1> 
    void LU_RDivEq(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
        // Solve x P L U = m:
    {
#ifdef XDEBUG
        Matrix<T> m0(m);
        Matrix<T1> L = LUx.lowerTri(UnitDiag);
        Matrix<T1> U = LUx.upperTri(NonUnitDiag);
        Matrix<T1> LU = L*U;
        LU.reversePermuteRows(P);
        //cout<<"Start LU_RDivEq\n";
        //cout<<"LU = "<<LU<<endl;
        //cout<<"m0 = "<<m0<<endl;
#endif

        TMVAssert(m.rowsize() == LUx.rowsize()); 
        TMVAssert(LUx.rowsize() == LUx.colsize());

#ifdef ALAP
        if (LUx.iscm() && !LUx.isconj() && m.isrm())
            LapLURDivEq(LUx,m);
        else 
#endif
            NonLapLURDivEq(LUx,m); 

        m.reversePermuteCols(P); 

#ifdef XDEBUG
        //cout<<"m-> "<<m<<endl;
        Matrix<T> mm = m*LU;
        //cout<<"mm = "<<mm<<endl;
        if (!(Norm(mm-m0) <= 0.001*Norm(L)*Norm(U)*Norm(m0))) {
            cerr<<"LU_RDivEq: m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"LU = "<<LU<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"m*LU = "<<mm<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


