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


#include "TMV_SymLDLDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
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

    // Note that we do not use the LAPACK division routines here.
    // LAPACK stores L and D in a complicated way which seems to be
    // intended to avoid the extra O(N) storage of xD.  
    // They store the xD vector in the subdiagonal of L, since the 
    // L matrix has 0's in these locations.
    //
    // However, this storage method makes the division _much_ slower,
    // because they don't permute all of the L matrix along the way.
    // This means that the permutations need to be done during the 
    // division routine, and they are mixed in with the multiplications,
    // which is a lot slower than doing all the permutations at the 
    // beginning and then dividing by a regular triangle matrix.
    //
    // Here are some timing measurements on my computer for 
    // 2000 x 2000 SymMatrix<double,Lower|ColMajor> decomposition and division,
    // both using my TMV code and directly calling the LAPACK routines
    // (not through the TMV library).
    // All times are in seconds - I took the fastest of 3 runs.
    //
    //                      TMV Non-Lap code     Direct LAPACK calls
    //
    // Decomposition              1.4                   1.3
    // (dsytrf)
    //
    // Divide into
    // 2000 x 2000 matrix         4.9                  54.1
    // (dsytrs)
    //
    // Divide into
    // 2000 x 200 matrix          0.54                  5.4
    // (dsytrs)
    //
    // Divide into
    // 2000 x 20 matrix           0.09                  0.17
    // (dsytrs)
    //
    // Divide into
    // 2000 x 1 vector            0.018                 0.022
    // (dsytrs)
    //
    // Find Inverse               2.9                   6.3
    // (dsytri)
    //
    // Hence our decision to forego the LAPACK routines here and
    // in calculating the inverse.  We do, however, keep the LAPACK
    // (in TMV_SymLDLDecompose.cpp) call for the decomposition, since it is 
    // comparable to my code or slightly faster.
    //
    template <class T, class T1> 
    void LDL_LDivEq(
        const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const ptrdiff_t* P, 
        MatrixView<T> m)
    {
        // Solve P L D Lt Pt x = m:
        TMVAssert(LL.size() == m.colsize());
        TMVAssert(xD.size()+1 == m.colsize());
        TMVAssert(LL.ct() == NonConj);
        TMVAssert(xD.ct() == NonConj);

        //cout<<"Start LDL_LDivEq:\n";
        //cout<<"xD = "<<xD<<endl;
#ifdef XDEBUG
        Matrix<T> m0(m);
        Matrix<T1> DD(LL.size(),LL.size(),T1(0));
        DD.diag() = LL.diag();
        DD.diag(-1) = xD;
        DD.diag(1) = LL.isherm() ? xD.conjugate() : xD.view();
        LowerTriMatrix<T1> L = LL.lowerTri(UnitDiag);
        Matrix<T1> LDL = L*DD*(LL.isherm() ? L.adjoint() : L.transpose());
        LDL.reversePermuteRows(P);
        LDL.reversePermuteCols(P);
#endif
        //cout<<"xD = "<<xD<<endl;

        m.permuteRows(P);
        m /= LL.lowerTri(UnitDiag);
        if (LL.isherm())
            PseudoDiag_LDivEq<true>(LL.diag(),xD,m);
        else
            PseudoDiag_LDivEq<false>(LL.diag(),xD,m);
        m /= LL.upperTri(UnitDiag);
        m.reversePermuteRows(P);

        //cout<<"xD = "<<xD<<endl;
#ifdef XDEBUG
        Matrix<T> mm = LDL*m;
        if (!(Norm(mm-m0) > 0.001*Norm(LDL)*Norm(m0))) {
            cerr<<"LDL_LDivEq: m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"D = "<<LL.diag()<<endl;
            cerr<<"xD = "<<xD<<endl;
            cerr<<"DD = "<<DD<<endl;
            cerr<<"LDL = "<<LDL<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"LDL*m = "<<mm<<endl;
            abort();
        }
#endif
        //cout<<"Done LDL_LDivEq: xD = "<<xD<<endl;
    }

    //
    // RDivEq Matrix
    //

    template <class T, class T1> 
    void LDL_RDivEq(
        const GenSymMatrix<T1>& LL, const GenVector<T1>& xD, const ptrdiff_t* P, 
        MatrixView<T> m)
    {
        // Solve x P L D Lt Pt = m:
        // P L Dt Lt Pt xt = mt
        TMVAssert(LL.size() == m.rowsize());
        TMVAssert(xD.size()+1 == m.rowsize());
        TMVAssert(LL.ct() == NonConj);
        TMVAssert(xD.ct() == NonConj);
        LDL_LDivEq(LL,xD,P,LL.isherm()?m.adjoint():m.transpose());
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymLDLDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


