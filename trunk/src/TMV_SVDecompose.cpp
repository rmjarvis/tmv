
//#define PRINTALGO_SVD
//#define XDEBUG_SVD

#include "TMV_Blas.h"
#include "tmv/TMV_SVDecompose.h"
#include "tmv/TMV_SVDecompose_Bidiag.h"
#include "tmv/TMV_SVDecompose_QR.h"
#include "tmv/TMV_SVDecompose_DC.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ProdXV.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_SortV.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_UnpackQ.h"
#include "tmv/TMV_QRDecompose.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_SumMX.h"

namespace tmv {

    template <class Tu, class RT, class Tv> 
    static inline void DoSV_DecomposeFromBidiagonal(
        MatrixView<Tu> U, VectorView<RT>& D, VectorView<RT>& E,
        MatrixView<Tv> V, bool setUV)
    { 
#ifdef XDEBUG_SVD
        std::cout<<"Start DoSV_DecomposeFromBidiagonal\n";
        dbgcout<<"Norm(U) = "<<Norm(U)<<std::endl;
        dbgcout<<"Norm(V) = "<<Norm(V)<<std::endl;
        dbgcout<<"Norm(D) = "<<Norm(D)<<std::endl;
        dbgcout<<"Norm(E) = "<<Norm(E)<<std::endl;
#endif
        MatrixView<Tu,ColMajor> Ucm = U;
        VectorView<RT,Unit> Du = D;
        VectorView<RT,Unit> Eu = E;
#ifdef XDEBUG_SVD
        dbgcout<<"Norm(Ucm) = "<<Norm(Ucm)<<std::endl;
        dbgcout<<"Norm(Du) = "<<Norm(Du)<<std::endl;
        dbgcout<<"Norm(Eu) = "<<Norm(Eu)<<std::endl;
#endif
        if (V.iscm()) {
            MatrixView<Tv,ColMajor> Vcm = V;
#ifdef XDEBUG_SVD
            dbgcout<<"Norm(Vcm) = "<<Norm(Vcm)<<std::endl;
#endif
            InlineSV_DecomposeFromBidiagonal(Ucm,Du,Eu,Vcm,setUV);
        } else {
            MatrixView<Tv,RowMajor> Vrm = V;
#ifdef XDEBUG_SVD
            dbgcout<<"Norm(Vrm) = "<<Norm(Vrm)<<std::endl;
#endif
            InlineSV_DecomposeFromBidiagonal(Ucm,Du,Eu,Vrm,setUV);
        }
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void LapSVDecomposeFromBidiagonal(
        MatrixView<double> U, VectorView<double>& D, VectorView<double>& E,
        MatrixView<double> V, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) { 
            TMVAssert(V.rowsize() == V.colsize()); 
            TMVAssert(V.rowsize() == D.size()); 
            if (U.cptr() && setUV) TMVAssert(U.iscm() == V.iscm());
        }

        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        if (setUV) {
            char c = 'I';
            TMVAssert(U.cptr() && V.cptr());
            //std::cout<<"setUV\n";
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"V = "<<V<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
            if (U.iscm()) {
                TMVAssert(V.iscm());
                int ldu = U.stepj();
                int ldv = V.stepj();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(U.ptr()),LAPV(ldu),LAPP(V.ptr()),LAPV(ldv),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            } else {
                u = 'L';
                TMVAssert(U.isrm());
                TMVAssert(V.isrm());
                int ldu = U.stepi();
                int ldv = V.stepi();
#ifndef LAPNOWORK
                int lwork = (3*n+4)*n;
                AlignedArray<double> work(lwork);
                VectorViewOf(work.get(),lwork).setZero();
                lwork = 8*n;
                AlignedArray<int> iwork(lwork);
#endif
                LAPNAME(dbdsdc) (
                    LAPCM LAPV(u),LAPV(c),LAPV(n),
                    LAPP(D.ptr()),LAPP(E1.ptr()),
                    LAPP(V.ptr()),LAPV(ldv),LAPP(U.ptr()),LAPV(ldu),0,0
                    LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            }
        } else if (U.cptr() || V.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n,0.);
            Matrix<double,ColMajor> V1(n,n,0.);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
            //std::cout<<"U.cptr() || V.cptr()\n";
            //if (U.cptr()) std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"U1 = "<<U1<<std::endl;
            //if (V.cptr()) std::cout<<"V = "<<V<<std::endl;
            //std::cout<<"V1 = "<<V1<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (U.cptr()) U = U*U1;
            if (V.cptr()) V = V1*V;
        } else {
            //std::cout<<"!(U.cptr() || V.cptr())\n";
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
            int ldu = n;
            int ldv = n;
            char c = 'N';
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("dbdsdc");
        E = E1.subVector(0,n-1);
        //std::cout<<"Done: D => "<<D<<std::endl;
        //std::cout<<"E1 => "<<E1<<std::endl;
        //std::cout<<"E => "<<E<<std::endl;
    }
    void LapSVDecomposeFromBidiagonal(
        MatrixView<std::complex<double> > U,
        VectorView<double>& D, VectorView<double>& E,
        MatrixView<std::complex<double> > V, bool setUV)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) { 
            TMVAssert(V.rowsize() == V.colsize()); 
            TMVAssert(V.rowsize() == D.size()); 
        }

        char u = 'U';
        int n = D.size();
        Vector<double> E1(n);
        E1.subVector(0,n-1) = E;
        E1[n-1] = 0.;
        if (U.cptr() || V.cptr()) {
            char c = 'I';
            Matrix<double,ColMajor> U1(n,n,0.);
            Matrix<double,ColMajor> V1(n,n,0.);
            int ldu = U1.stepj();
            int ldv = V1.stepj();
            //std::cout<<"U.cptr() || V.cptr()\n";
            //if (U.cptr()) std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"U1 = "<<U1<<std::endl;
            //if (V.cptr()) std::cout<<"V = "<<V<<std::endl;
            //std::cout<<"V1 = "<<V1<<std::endl;
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = (3*n+4)*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                LAPP(U1.ptr()),LAPV(ldu),LAPP(V1.ptr()),LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
            if (setUV) {
                if (U.cptr()) U = U1;
                if (V.cptr()) V = V1;
            } else {
                if (U.cptr()) U = U*U1;
                if (V.cptr()) V = V1*V;
            }
        } else {
            TMVAssert(!setUV);
            int ldu = n;
            int ldv = n;
            char c = 'N';
            //std::cout<<"!(U.cptr() || V.cptr())\n";
            //std::cout<<"D = "<<D<<std::endl;
            //std::cout<<"E = "<<E<<std::endl;
            //std::cout<<"E1 = "<<E1<<std::endl;
#ifndef LAPNOWORK
            int lwork = 4*n;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
            lwork = 8*n;
            AlignedArray<int> iwork(lwork);
#endif
            LAPNAME(dbdsdc) (
                LAPCM LAPV(u),LAPV(c),LAPV(n),
                LAPP(D.ptr()),LAPP(E1.ptr()),
                0,LAPV(ldu),0,LAPV(ldv),0,0
                LAPWK(work.get()) LAPWK(iwork.get()) LAPINFO LAP1 LAP1);
        }
        LAP_Results("dbdsdc");
        E = E1.subVector(0,n-1);
        //std::cout<<"Done: D => "<<D<<std::endl;
        //std::cout<<"E1 => "<<E1<<std::endl;
        //std::cout<<"E => "<<E<<std::endl;
    }
#endif
#endif // LAP

    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool setUV)
    {
        TMVAssert(U.iscm());
        TMVAssert(V.iscm() || V.isrm());
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        if (D.size() > 0) {
            DoSV_DecomposeFromBidiagonal(U,D,E,V,setUV);
        }
    }

    template <class T, class RT>
    void InstSV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> S,
        MatrixView<T> V, T& signdet, RT& logdet, bool StoreU)
    {
#ifdef XDEBUG_SVD
        Matrix<T> m = U;
        std::cout<<"Start InstSV_Decompose\n";
        std::cout<<"U = "<<TMV_Text(U)<<std::endl;
        std::cout<<"V = "<<TMV_Text(V)<<std::endl;
        std::cout<<"StoreU = "<<StoreU<<std::endl;
        std::cout<<"ptrs = "<<U.cptr()<<" "<<V.cptr()<<std::endl;
        if (U.colsize() >= 3 && U.rowsize() >= 3)
            std::cout<<"m = "<<m.subMatrix(0,3,0,3)<<std::endl;
#endif
        if (S.size() > 0) {
            if (U.iscm()) {
                if (V.iscm() || V.isrm() || !V.cptr()) {
                    //std::cout<<"Normal Case\n";
                    MatrixView<T,ColMajor> Ucm = U;
                    if (V.iscm() || !V.cptr()) {
                        MatrixView<T,ColMajor> Vcm = V;
                        InlineSV_Decompose(Ucm,S,Vcm,signdet,logdet,StoreU); 
                    } else {
                        MatrixView<T,RowMajor> Vrm = V;
                        InlineSV_Decompose(Ucm,S,Vrm,signdet,logdet,StoreU); 
                    }
                } else {
                    //std::cout<<"Copy V\n";
                    Matrix<T,ColMajor|NoDivider> Vc = V;
                    InstSV_Decompose(U,S,Vc.xView(),signdet,logdet,StoreU);
                    InstCopy(Vc.constView().xView(),V);
                } 
            } else if (V.isrm() && V.cptr() && U.isSquare()) {
                //std::cout<<"Copy to V and transpose\n";
                TMVAssert(V.colsize() == U.rowsize());
                TMVAssert(V.rowsize() == U.colsize());
                InstCopy(U.constView(),V);
                if (StoreU) {
                    InstSV_Decompose(
                        V.transpose(),S,U.transpose(),signdet,logdet,true);
                } else {
                    MatrixView<T> Ux(0,0,0,1,1);
                    InstSV_Decompose(
                        V.transpose(),S,Ux,signdet,logdet,true);
                }
            } else if (U.isrm() && U.isSquare()) {
                //std::cout<<"Transpose Self\n";
                U.transposeSelf();
                InstSV_Decompose(
                    U.transpose(),S,V,signdet,logdet,StoreU);
                if (StoreU) U.transposeSelf();
            } else {
                //std::cout<<"Copy U\n";
                Matrix<T,ColMajor|NoDivider> Uc = U;
                InstSV_Decompose(Uc.xView(),S,V,signdet,logdet,StoreU);
                InstCopy(Uc.constView().xView(),U);
            }
        }
#ifdef XDEBUG_SVD
        std::cout<<"Done InstSV_Decompose\n";
        if (StoreU && V.cptr() && U.colsize() >= 3 && U.rowsize() >= 3) {
            Matrix<T> USV = U*S*V;
            std::cout<<"USV = "<<USV.subMatrix(0,3,0,3)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(m-USV)<<std::endl;
            if (!(Norm(m-USV) <= 1.e-4 * Norm(m))) {
                std::cout<<"U = "<<U<<std::endl;
                std::cout<<"S = "<<S<<std::endl;
                std::cout<<"V = "<<V<<std::endl;
                std::cout<<"m = "<<m<<std::endl;
                std::cout<<"USV = "<<USV<<std::endl;
                std::cout<<"Norm(m-USV) = "<<Norm(m-USV)<<std::endl;
                std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;
                abort();
            }
        }
#endif
    }


#define InstFile "TMV_SVDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


