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


#include "TMV_Blas.h"
#include "TMV_MultMM.h"
#include "TMV_MultMM_Kernel.h"
#include "TMV_Array.h"
#include "tmv/TMV_MatrixArith.h"

//#include <iostream>

namespace tmv {

    // Do the MultMM calculation for a single block of C: 
    // C += x * A * B
    template <class RT, class T, class Func>
    static void SingleBlockMultMM(
        bool add, const T& x, const GenMatrix<RT>& A,
        const GenMatrix<RT>& B, const MatrixView<T>& C,
        const int i1, const int j1, const int i2, const int j2,
        const int MB, const int NB, const int KB,
        const int size1, const int size2, const int size3,
        RT* Ap, RT* Bp, RT* Cp, bool firstA, bool firstB, 
        Func* myprod, Func* mycleanup,
        const int K, const int Kc, const int Kd)
    {
        //std::cout<<"Start SingleBlockMultMM:\n";
        //std::cout<<"i = "<<i1<<','<<i2<<"  j = "<<j1<<','<<j2<<std::endl;
        //std::cout<<"block sizes = "<<MB<<','<<NB<<','<<KB<<std::endl;
        //std::cout<<"sizes = "<<size1<<','<<size2<<','<<size3<<std::endl;
        //std::cout<<"pointers = "<<Ap<<"  "<<Bp<<"  "<<Cp<<std::endl;
        //std::cout<<"first?  "<<firstA<<"  "<<firstB<<std::endl;
        //std::cout<<"func = "<<myprod<<"  "<<mycleanup<<std::endl;
        //std::cout<<"K,Kc,kd = "<<K<<','<<Kc<<','<<Kd<<std::endl;

        //std::cout<<"Cx: "<<Cp<<"..."<<Cp+size3<<std::endl;
        VectorView<RT>(Cp,size3,1,NonConj).setZero();

        int k1 = 0;
        for (int k2=KB;k2<=K;k1=k2,k2+=KB,Ap+=size1,Bp+=size2) {
            //std::cout<<"Ax: "<<Ap<<"..."<<Ap+MB*KB<<std::endl;
            //std::cout<<"Bx: "<<Bp<<"..."<<Bp+NB*KB<<std::endl;
            MatrixView<RT> Ax(Ap,MB,KB,KB,1,RowMajor,NonConj,MB*KB);
            MatrixView<RT> Bx(Bp,KB,NB,1,KB,ColMajor,NonConj,NB*KB);
            if (firstA) Ax = A.subMatrix(i1,i2,k1,k2);
            if (firstB) Bx = B.subMatrix(k1,k2,j1,j2);
            (*myprod) (MB,NB,KB, Ap,Bp,Cp);
        }

        if (Kc) {
            //std::cout<<"Ay: "<<Ap<<"..."<<Ap+MB*Kc<<std::endl;
            //std::cout<<"By: "<<Bp<<"..."<<Bp+NB*Kc<<std::endl;
            MatrixView<RT> Ay(Ap,MB,Kc,Kd,1,RowMajor,NonConj,MB*Kc);
            MatrixView<RT> By(Bp,Kc,NB,1,Kd,ColMajor,NonConj,NB*Kc);
            if (firstA) Ay = A.subMatrix(i1,i2,k1,K);
            if (firstB) By = B.subMatrix(k1,K,j1,j2);
            (*mycleanup) (MB,NB,Kc, Ap,Bp,Cp);
        }
        //std::cout<<"Cx: "<<Cp<<"..."<<Cp+MB*NB<<std::endl;
        if (add) {
            C.subMatrix(i1,i2,j1,j2) += 
                x*ConstMatrixView<RT>(Cp,MB,NB,1,MB,ColMajor,NonConj,MB*NB);
        } else if (x != T(1)) {
            C.subMatrix(i1,i2,j1,j2) = 
                x*ConstMatrixView<RT>(Cp,MB,NB,1,MB,ColMajor,NonConj,MB*NB);
        } else {
            C.subMatrix(i1,i2,j1,j2) = 
                ConstMatrixView<RT>(Cp,MB,NB,1,MB,ColMajor,NonConj,MB*NB);
        }
    }

    template <class T>
    struct BlockHelper
    {
        typedef void KernelMultMMFunc(
            const int M, const int N, const int K,
            const T* A, const T* B, T* C);

        enum { KB = 16 };
        enum { lgKB = 4 };
        static KernelMultMMFunc* getProdFunc() 
        { return &call_multmm_16_16_16<T>; }
        static KernelMultMMFunc* getProdFuncA() 
        { return &call_multmm_16_N_16<T>; }
        static KernelMultMMFunc* getProdFuncB() 
        { return &call_multmm_M_16_16<T>; }
        static KernelMultMMFunc* getProdFuncC() 
        { return &call_multmm_M_N_K<T>; }
        static int RoundUp(int x) { return x; }
    };

#ifdef INST_DOUBLE
#ifdef __SSE2__
    template <>
    struct BlockHelper<double>
    {
        typedef double T;
        typedef void KernelMultMMFunc(
            const int M, const int N, const int K,
            const T* A, const T* B, T* C);

        enum { KB = 32 };
        enum { lgKB = 5 };
        static KernelMultMMFunc* getProdFunc() 
        { return &call_multmm_16_16_32<T>; }
        static KernelMultMMFunc* getProdFuncA() 
        { return &call_multmm_16_N_32<T>; }
        static KernelMultMMFunc* getProdFuncB() 
        { return &call_multmm_M_16_32<T>; }
        static KernelMultMMFunc* getProdFuncC() 
        { return &call_multmm_M_N_K<T>; }
        static int RoundUp(int x) 
        { return (x == 0 ? 0 : (((x-1)>>1)+1)<<1); }
    };
#endif
#endif
#ifdef INST_FLOAT
#ifdef __SSE__
    template <>
    struct BlockHelper<float>
    {
        typedef float T;
        typedef void KernelMultMMFunc(
            const int M, const int N, const int K,
            const T* A, const T* B, T* C);

        enum { KB = 64 };
        enum { lgKB = 6 };
        static KernelMultMMFunc* getProdFunc() 
        { return &call_multmm_16_16_64<T>; }
        static KernelMultMMFunc* getProdFuncA() 
        { return &call_multmm_16_N_64<T>; }
        static KernelMultMMFunc* getProdFuncB() 
        { return &call_multmm_M_16_64<T>; }
        static KernelMultMMFunc* getProdFuncC() 
        { return &call_multmm_M_N_K<T>; }
        static int RoundUp(int x) 
        { return (x == 0 ? 0 : (((x-1)>>2)+1)<<2); }
    };
#endif
#endif



    // Do the loops over the full matrix C, breaking up the calculation
    // into blocks.
    template <class RT, class T>
    static void DoBlockMultMM(
        bool add, const T& x, const GenMatrix<RT>& A,
        const GenMatrix<RT>& B, const MatrixView<T>& C)
    {
        const int M = C.colsize();
        const int N = C.rowsize();
        const int K = A.rowsize();

        const int KB = BlockHelper<RT>::KB;
        const int lgKB = BlockHelper<RT>::lgKB;

        typedef typename BlockHelper<RT>::KernelMultMMFunc Func;
        Func* prod = BlockHelper<RT>::getProdFunc();
        Func* proda = BlockHelper<RT>::getProdFuncA();
        Func* prodb = BlockHelper<RT>::getProdFuncB();
        Func* prodc = BlockHelper<RT>::getProdFuncC();

        Func* cleanup = &call_multmm_16_16_K<RT>;
        Func* cleanupabc = &call_multmm_M_N_K<RT>;

        const int Mb = (M>>4); // = M/16
        const int Nb = (N>>4); // = N/16
        const int Kb = (K>>lgKB); // = K/KB
        const int Ma = (Mb<<4); // = M/16*16
        const int Na = (Nb<<4); // = N/16*16
        const int Ka = (Kb<<lgKB); // = K/KB*KB
        const int Mc = M-Ma; // = M%16
        const int Nc = N-Na; // = N%16
        const int Kc = K-Ka; // = K%KB
        const int Kd = BlockHelper<RT>::RoundUp(Kc);
        // = Kc rounded up to multiple of 2 or 4 as required for SSE commands
        const int Ktot_d = (Kb<<lgKB) + Kd;
        TMVAssert(Ktot_d == BlockHelper<RT>::RoundUp(K));

        const int size1 = 16*KB;
        const int size1y = Mc*KB;
        const int size2 = 16*KB;
        const int size2y = Nc*KB;
        const int size3 = 16*16;

        AlignedArray<RT> C_temp(size3);
        RT*const Cp = C_temp;
        if (N >= M) {
            // Then we loop over the columns of B (in blocks).
            // We store a full copy of A in block format.
            // Each block column of B is copied one at a time into 
            // temporary storage.
            AlignedArray<RT> A_temp(M*Ktot_d);
            AlignedArray<RT> B_temp(16*Ktot_d);
            RT*const Ap0 = A_temp;
            RT*const Bp0 = B_temp;
            const int fullsize1 = 16*Ktot_d;

            for(int j=0;j<Nb;++j) {
                RT* Ap = Ap0;
                for(int i=0;i<Mb;++i) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        i*16,j*16,(i+1)*16,(j+1)*16,
                        16,16,KB,
                        size1,size2,size3,
                        Ap,Bp0,Cp, j==0,i==0, 
                        prod, cleanup,
                        K,Kc,Kd);
                    Ap += fullsize1;
                }
                if (Mc) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        Ma,j*16,M,(j+1)*16,
                        Mc,16,KB,
                        size1y,size2,size3,
                        Ap,Bp0,Cp, j==0,Mb==0, 
                        prodb, cleanupabc,
                        K,Kc,Kd);
                }
            }
            if (Nc) {
                RT* Ap = Ap0;
                for(int i=0;i<Mb;++i) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        i*16,Na,(i+1)*16,N,
                        16,Nc,KB,
                        size1,size2y,size3,
                        Ap,Bp0,Cp, Nb==0,i==0,
                        proda, cleanupabc,
                        K,Kc,Kd);
                    Ap += fullsize1;
                }
                if (Mc) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        Ma,Na,M,N,
                        Mc,Nc,KB,
                        size1y,size2y,size3,
                        Ap,Bp0,Cp, Nb==0,Mb==0,
                        prodc, cleanupabc,
                        K,Kc,Kd);
                }
            }
        } else {
            // Then we loop over the rows of A and store a full copy of B.
            AlignedArray<RT> A_temp(16*Ktot_d);
            AlignedArray<RT> B_temp(N*Ktot_d);
            RT*const Ap0 = A_temp;
            RT*const Bp0 = B_temp;
            const int fullsize2 = 16*Ktot_d;

            for(int i=0;i<Mb;++i) {
                RT* Bp = Bp0;
                for(int j=0;j<Nb;++j) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        i*16,j*16,(i+1)*16,(j+1)*16,
                        16,16,KB,
                        size1,size2,size3,
                        Ap0,Bp,Cp, j==0,i==0,
                        prod, cleanup,
                        K,Kc,Kd);
                    Bp += fullsize2;
                }
                if (Nc) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        i*16,Na,(i+1)*16,N,
                        16,Nc,KB,
                        size1,size2y,size3,
                        Ap0,Bp,Cp, Nb==0,i==0, 
                        proda, cleanupabc,
                        K,Kc,Kd);
                }
            }
            if (Mc) {
                RT* Bp = Bp0;
                for(int j=0;j<Nb;++j) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        Ma,j*16,M,(j+1)*16,
                        Mc,16,KB,
                        size1y,size2,size3,
                        Ap0,Bp,Cp, j==0,Mb==0, 
                        prodb, cleanupabc,
                        K,Kc,Kd);
                    Bp += fullsize2;
                }
                if (Nc) {
                    SingleBlockMultMM(
                        add,x,A,B,C,
                        Ma,Na,M,N,
                        Mc,Nc,KB,
                        size1y,size2y,size3,
                        Ap0,Bp,Cp, Nb==0,Mb==0, 
                        prodc, cleanupabc,
                        K,Kc,Kd);
                }
            }
        }
    }

    // Turn the various complex varieties into real versions:
#define CT std::complex<RT>
    template <class RT>
    static void DoBlockMultMM(
        bool add, const CT x, const GenMatrix<RT>& A,
        const GenMatrix<CT>& B, const MatrixView<CT>& C)
    {
        bool Bc = B.isconj();
        if (TMV_IMAG(x) == RT(0)) {
            const RT xr = TMV_REAL(x);
            DoBlockMultMM(add,xr,A,B.realPart(),C.realPart());
            DoBlockMultMM(add,Bc?-xr:xr,A,B.imagPart(),C.imagPart());
        } else {
            CT ix = std::complex<RT>(0,1) * x;
            DoBlockMultMM(add,x,A,B.realPart(),C);
            DoBlockMultMM(true,Bc?-ix:ix,A,B.imagPart(),C);
        }
    }
    template <class RT>
    static void DoBlockMultMM(
        bool add, const CT x, const GenMatrix<CT>& A,
        const GenMatrix<RT>& B, const MatrixView<CT>& C)
    {
        bool Ac = A.isconj();
        if (TMV_IMAG(x) == RT(0)) {
            const RT xr = TMV_REAL(x);
            DoBlockMultMM(add,xr,A.realPart(),B,C.realPart());
            DoBlockMultMM(add,Ac?-xr:xr,A.imagPart(),B,C.imagPart());
        } else {
            CT ix = std::complex<RT>(0,1) * x;
            DoBlockMultMM(add,x,A.realPart(),B,C);
            DoBlockMultMM(true,Ac?-ix:ix,A.imagPart(),B,C);
        }
    }
    template <class RT>
    static void DoBlockMultMM(
        bool add, const CT x, const GenMatrix<CT>& A,
        const GenMatrix<CT>& B, const MatrixView<CT>& C)
    {
        bool Ac = A.isconj();
        bool Bc = B.isconj();
        if (TMV_IMAG(x) == RT(0)) {
            const RT xr = TMV_REAL(x);
            DoBlockMultMM(
                add,xr,A.realPart(),B.realPart(),C.realPart());
            DoBlockMultMM(
                true,Ac==Bc?-xr:xr,A.imagPart(),B.imagPart(),C.realPart());
            DoBlockMultMM(
                add,Bc?-xr:xr,A.realPart(),B.imagPart(),C.imagPart());
            DoBlockMultMM(
                true,Ac?-xr:xr,A.imagPart(),B.realPart(),C.imagPart());
        } else {
            CT ix = std::complex<RT>(0,1) * x;
            DoBlockMultMM(add,x,A.realPart(),B.realPart(),C);
            DoBlockMultMM(true,Ac==Bc?-x:x,A.imagPart(),B.imagPart(),C);
            DoBlockMultMM(true,Bc?-ix:ix,A.realPart(),B.imagPart(),C);
            DoBlockMultMM(true,Ac?-ix:ix,A.imagPart(),B.realPart(),C);
        }
    }
#undef CT

    template <bool add, class T, class Ta, class Tb>
    void BlockMultMM(
        const T x, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, const MatrixView<T>& C)
    { 
        try {
            DoBlockMultMM(add,x,A,B,C); 
        } catch (std::bad_alloc) {
            // Use a simpler algorithm that doesn't allocate memory:
            if (A.iscm()) {
                if (B.iscm())
                    CCCMultMM<add>(x,A,B,C);
                else
                    CRCMultMM<add>(x,A,B,C);
            } else {
                if (B.iscm())
                    RCCMultMM<add>(x,A,B,C);
                else {
                    // There is no simple algorithm for RRC, so we 
                    // need to do some copying.  But since we have
                    // had a bad_alloc already, we only allocate 
                    // new storage in thin slices (width = 16)
                    int N=B.rowsize();
                    if (N > 16) {
                        int K=B.colsize();
                        Matrix<T,ColMajor> B1(K,16);
                        int j1=0;
                        for (int j2=16;j2<N;j1=j2,j2+=16) {
                            B1 = B.colRange(j1,j2);
                            RCCMultMM<add>(x,A,B1,C.colRange(j1,j2));
                        }
                        B1.colRange(0,N-j1) = B.colRange(j1,N);
                        RCCMultMM<add>(
                            x,A,B1.colRange(0,N-j1),C.colRange(j1,N));
                    } else {
                        Matrix<T,ColMajor> B1 = B;
                        RCCMultMM<add>(x,A,B1,C);
                    }
                }
            }
        }
    }

#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_MultMM_Block.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

