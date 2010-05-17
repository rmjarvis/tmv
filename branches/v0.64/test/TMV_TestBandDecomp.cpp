
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestBandArith.h"

template <class T, tmv::StorageType stor> 
void TestBandDecomp()
{
    for (int mattype = 0; mattype < 7; mattype++) {
        // mattype = 0  is Square
        // mattype = 1  is NonSquare short
        // mattype = 2  is NonSquare tall
        // mattype = 3  is TriDiag
        // mattype = 4  is Singular
        // mattype = 5  is Lower BiDiag
        // mattype = 6  is Upper BiDiag

        const int N = 200;
        int M = N;
        int nlo = 12;
        int nhi = 7;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
        else if (mattype == 3) nlo = nhi = 1;
        else if (mattype == 5) { nlo = 1; nhi = 0; }
        else if (mattype == 6) { nlo = 0; nhi = 1; }

        tmv::BandMatrix<T,stor> m(M,N,nlo,nhi);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            if (i<=j+nlo && j<=i+nhi) m(i,j) = T(2+4*i-5*j);
        if (mattype != 4) {
            m(0,0) = T(14);
            if (m.nlo() >= 1) m(1,0) = T(-2);
            if (m.nhi() >= 1) m(4,5) = T(7);
            m(2,2) = T(-10);
            m.diag() *= T(30);
        }

        tmv::BandMatrix<std::complex<T>,stor> c(M,N,nlo,nhi);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            if (i<=j+nlo && j<=i+nhi) c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        if (mattype != 4) {
            c(0,0) = T(14);
            if (c.nlo() >= 1) c(1,0) = T(-2);
            if (c.nhi() >= 1) c(4,5) = T(7);
            c(2,2) = T(-10);
            c.diag() *= T(30);
        }

        T eps = EPS;
        T ceps = EPS;
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (mattype != 3) {
            T kappa = Norm(m) * Norm(m.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<std::endl;
            eps *= kappa;
            ceps *= kappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;


        // LU Decomposition
        if (mattype == 0) {
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            const int* p = m.lud().getP();
            tmv::Matrix<T> PLU = L*U;
            PLU.reversePermuteRows(p);
            if (m.lud().isTrans()) PLU.transposeSelf();
            Assert(Norm(m-PLU) < eps*Norm(m),"Band LU");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = 
                c.lud().getL();
            tmv::UpperTriMatrix<std::complex<T> > cU = c.lud().getU();
            p = c.lud().getP();
            tmv::Matrix<std::complex<T> > cPLU = cL*cU;
            cPLU.reversePermuteRows(p);
            if (c.lud().isTrans()) cPLU.transposeSelf();
            Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU");

#ifdef XTEST
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L2(M);
            tmv::BandMatrix<T,stor> U2(M,M,0,nlo+nhi);
            int p2[N];
            LU_Decompose(m,L2.view(),U2.view(),p2);
            PLU = L2*U2;
            PLU.reversePermuteRows(p2);
            Assert(Norm(m-PLU) < eps*Norm(m),"Band LU2");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL2(M);
            tmv::BandMatrix<std::complex<T>,stor> cU2(M,M,0,nlo+nhi);
            LU_Decompose(c,cL2.view(),cU2.view(),p2);
            cPLU = cL2*cU2;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU2");

            LU_Decompose(c,cL2.conjugate(),cU2.view(),p2);
            cPLU = cL2.conjugate()*cU2;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU3");

            LU_Decompose(c,cL2.view(),cU2.conjugate(),p2);
            cPLU = cL2*cU2.conjugate();
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU4");

            LU_Decompose(c,cL2.conjugate(),cU2.conjugate(),p2);
            cPLU = cL2.conjugate()*cU2.conjugate();
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU5");

            LU_Decompose(c.conjugate(),cL2.view(),cU2.view(),p2);
            cPLU = cL2*cU2;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c.conjugate()-cPLU) < ceps*Norm(c),"Band C LU6");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2.view(),p2);
            cPLU = cL2.conjugate()*cU2;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c.conjugate()-cPLU) < ceps*Norm(c),"Band C LU7");

            LU_Decompose(c.conjugate(),cL2.view(),cU2.conjugate(),p2);
            cPLU = cL2*cU2.conjugate();
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c.conjugate()-cPLU) < ceps*Norm(c),"Band C LU8");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2.conjugate(),p2);
            cPLU = cL2.conjugate()*cU2.conjugate();
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c.conjugate()-cPLU) < ceps*Norm(c),"Band C LU9");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // QR Decomposition
        if (mattype != 4) {
            tmv::Matrix<T> Q = m.qrd().getQ();
            tmv::BandMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            if (m.qrd().isTrans()) QR.transposeSelf();
            Assert(Norm(m-QR) < eps*Norm(m),"Band QR");

            tmv::Matrix<std::complex<T> > cQ = c.qrd().getQ();
            tmv::BandMatrix<std::complex<T> > cR = c.qrd().getR();
            tmv::Matrix<std::complex<T> > cQR = cQ*cR;
            if (c.qrd().isTrans()) cQR.transposeSelf();
            Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR");

#ifdef XTEST
            QR_Decompose(m,Q.view(),R.view());
            QR = Q*R;
            Assert(Norm(m-QR) < eps*Norm(m),"Band QR2");

            tmv::BandMatrix<T,stor> R2(N,N,0,nlo+nhi);
            QR_Decompose(m,R2.view());
            Assert(Norm(R-R2) < eps*Norm(R),"Band QR3");

            QR_Decompose(c,cQ.view(),cR.view());
            cQR = cQ*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR2");

            tmv::BandMatrix<std::complex<T>,stor> cR2(N,N,0,nlo+nhi);
            QR_Decompose(c,cR2.view());
            Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR3");

            QR_Decompose(c,cQ.conjugate(),cR2.view());
            cQR = cQ.conjugate()*cR2;
            Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR4");

            QR_Decompose(c,cQ.view(),cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR5");

            QR_Decompose(c,cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR6");

            QR_Decompose(c,cR2.view());
            Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR7");

            QR_Decompose(c,cR2.conjugate());
            Assert(Norm(cR-cR2.conjugate()) < ceps*Norm(cR),"Band C QR8");

            QR_Decompose(c.conjugate(),cQ.view(),cR2.view());
            cQR = cQ*cR2;
            Assert(Norm(c.conjugate()-cQR) < ceps*Norm(c),"Band C QR9");

            QR_Decompose(c.conjugate(),cR2.view());
            Assert(Norm(cR.conjugate()-cR2) < ceps*Norm(cR),"Band C QR10");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2.view());
            cQR = cQ.conjugate()*cR2;
            Assert(Norm(c.conjugate()-cQR) < ceps*Norm(c),"Band C QR11");

            QR_Decompose(c.conjugate(),cQ.view(),cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Norm(c.conjugate()-cQR) < ceps*Norm(c),"Band C QR12");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Norm(c.conjugate()-cQR) < ceps*Norm(c),"Band C QR13");

            QR_Decompose(c.conjugate(),cR2.view());
            Assert(Norm(cR.conjugate()-cR2) < ceps*Norm(cR),"Band C QR14");

            QR_Decompose(c.conjugate(),cR2.conjugate());
            Assert(Norm(cR.conjugate()-cR2.conjugate()) < ceps*Norm(cR),
                   "Band C QR15");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // SV Decomposition
        {
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            Assert(Norm(m-U*S*V) < eps*Norm(m),"SV");

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(M,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SV2");

            SV_Decompose(m,S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"SV4 S");
            Assert(Norm(U2*S2*S2*U2.transpose()-m*m.transpose()) < 
                   eps*Norm(m)*Norm(m),"SV3 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"SV5 S");
            Assert(Norm(V2.transpose()*S2*S2*V2-m.transpose()*m) < 
                   eps*Norm(m)*Norm(m),"SV5 V");

            tmv::Matrix<std::complex<T> > cU2(M,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) < eps*Norm(c),"C SV2");

            SV_Decompose(c,cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV4 S");
            Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) < 
                   ceps*Norm(c)*Norm(c),"C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV5 S");
            Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) < 
                   ceps*Norm(c)*Norm(c),"C SV5 V");

            SV_Decompose(c,cU2.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2.conjugate()*cS2*cV2) < eps*Norm(c),"C SV6");
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.conjugate());
            Assert(Norm(c-cU2*cS2*cV2.conjugate()) < eps*Norm(c),"C SV7");
            SV_Decompose(c,cU2.conjugate(),cS2.view(),cV2.conjugate());
            Assert(Norm(c-cU2.conjugate()*cS2*cV2.conjugate()) < eps*Norm(c),
                   "C SV8");

            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV9 S");
            Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) < 
                   ceps*Norm(c)*Norm(c),"C SV9 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV10 S");
            Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) < 
                   ceps*Norm(c)*Norm(c),"C SV10 V");

            SV_Decompose(c.conjugate(),cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c.conjugate()-cU2*cS2*cV2) < eps*Norm(c),"C SV11");
            SV_Decompose(c.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV12");
            SV_Decompose(c.conjugate(),cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV13 S");
            Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c.conjugate()*c.transpose()) < 
                   ceps*Norm(c)*Norm(c),"C SV13 U");
            SV_Decompose(c.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV14 S");
            Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.transpose()*c.conjugate()) < 
                   ceps*Norm(c)*Norm(c),"C SV14 V");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(c.conjugate()-cU2.conjugate()*cS2*cV2) < eps*Norm(c),
                   "C SV15");
            SV_Decompose(c.conjugate(),cU2.view(),cS2.view(),cV2.conjugate());
            Assert(Norm(c.conjugate()-cU2*cS2*cV2.conjugate()) < eps*Norm(c),
                   "C SV16");
            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view(),
                         cV2.conjugate());
            Assert(Norm(c.conjugate()-cU2.conjugate()*cS2*cV2.conjugate()) < 
                   eps*Norm(c),"C SV17");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV18 S");
            Assert(Norm(cU2.conjugate()*cS2*cS2*cU2.transpose()-
                        c.conjugate()*c.transpose())<ceps*Norm(c)*Norm(c),
                   "C SV18 U");
            SV_Decompose(c.conjugate(),cS2.view(),cV2.conjugate());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV19 S");
            Assert(Norm(cV2.transpose()*cS2*cS2*cV2.conjugate()-
                        c.transpose()*c.conjugate())<ceps*Norm(c)*Norm(c),
                   "C SV19 V");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef INST_DOUBLE
template void TestBandDecomp<double,tmv::ColMajor>();
template void TestBandDecomp<double,tmv::RowMajor>();
template void TestBandDecomp<double,tmv::DiagMajor>();
#endif
#ifdef INST_FLOAT
template void TestBandDecomp<float,tmv::ColMajor>();
template void TestBandDecomp<float,tmv::RowMajor>();
template void TestBandDecomp<float,tmv::DiagMajor>();
#endif
#ifdef INST_LONGDOUBLE
template void TestBandDecomp<long double,tmv::ColMajor>();
template void TestBandDecomp<long double,tmv::RowMajor>();
template void TestBandDecomp<long double,tmv::DiagMajor>();
#endif

