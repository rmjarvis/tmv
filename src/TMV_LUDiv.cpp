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
        const GenMatrix<T1>& LUx, MatrixView<T> m)
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
        const GenMatrix<T1>& LUx, MatrixView<T> m)
    { NonLapLULDivEq(LUx,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapLULDivEq(
        const GenMatrix<double>& LUx, MatrixView<double> m)
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
        int Lap_info=0;
        LAPNAME(dgetrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"dgetrs");
    }
    template <> 
    void LapLULDivEq(
        const GenMatrix<std::complex<double> >& LUx,
        MatrixView<std::complex<double> > m)
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
        int Lap_info=0;
        LAPNAME(zgetrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"zgetrs");
    }
#endif
#ifdef INST_FLOAT
#ifndef MKL
    template <> 
    void LapLULDivEq(const GenMatrix<float>& LUx, MatrixView<float> m)
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
        int Lap_info=0;
        LAPNAME(sgetrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"sgetrs");
    }
    template <> 
    void LapLULDivEq(
        const GenMatrix<std::complex<float> >& LUx,
        MatrixView<std::complex<float> > m)
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
        int Lap_info=0;
        LAPNAME(cgetrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"cgetrs");
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T, class T1> 
    void LU_LDivEq(
        const GenMatrix<T1>& LUx, const ptrdiff_t* P, MatrixView<T> m)
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
        const GenMatrix<T1>& LUx, MatrixView<T> m)
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
        const GenMatrix<T1>& LUx, MatrixView<T> m)
    { NonLapLURDivEq(LUx,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapLURDivEq(
        const GenMatrix<double>& LUx, MatrixView<double> m)
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
        int Lap_info=0;
        LAPNAME(dgetrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"dgetrs");
    }
    template <> 
    void LapLURDivEq(
        const GenMatrix<std::complex<double> >& LUx,
        MatrixView<std::complex<double> > m)
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
        int Lap_info=0;
        LAPNAME(zgetrs) (
            LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"zgetrs");
    }
#endif
#ifdef INST_FLOAT
#ifndef MKL
    template <> 
    void LapLURDivEq(
        const GenMatrix<float>& LUx, MatrixView<float> m)
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
        int Lap_info=0;
        LAPNAME(sgetrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"sgetrs");
    }
    template <> 
    void LapLURDivEq(
        const GenMatrix<std::complex<float> >& LUx,
        MatrixView<std::complex<float> > m)
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
        int Lap_info=0;
        LAPNAME(cgetrs) (
            LAPCM m.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.cptr()),LAPV(lda),
            LAPP(ipiv.get()),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results(Lap_info,"cgetrs");
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T, class T1> 
    void LU_RDivEq(
        const GenMatrix<T1>& LUx, const ptrdiff_t* P, MatrixView<T> m)
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


