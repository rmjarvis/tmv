
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestBandArith.h"

template <class T, tmv::StorageType stor> 
void TestBandDecomp()
{
    for (int mattype = START; mattype <= 9; mattype++) {
        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<", stor = "<<
                tmv::TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Square
        // mattype = 1  is NonSquare short
        // mattype = 2  is NonSquare tall
        // mattype = 3  is TriDiag
        // mattype = 4  is Singular
        // mattype = 5  is Lower BiDiag
        // mattype = 6  is Upper BiDiag
        // mattype = 7  is Singular with seriously bad defects
        // mattype = 8  is Singular and nearly zero
        // mattype = 9  is Singular and nearly overflow

        const int N = 200;
        int M = N;
        int nlo = 12;
        int nhi = 7;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
        else if (mattype == 3) nlo = nhi = 1;
        else if (mattype == 5) { nlo = 1; nhi = 0; }
        else if (mattype == 6) { nlo = 0; nhi = 1; }
        const bool singular = mattype == 4 || mattype >= 7;
        const bool nearoverflow = mattype == 9;

        tmv::BandMatrix<T,stor> m(M,N,nlo,nhi);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            if (i<=j+nlo && j<=i+nhi) m(i,j) = T(2+4*i-5*j);
        if (!singular) {
            m(0,0) = T(14);
            if (m.nlo() >= 1) m(1,0) = T(-2);
            if (m.nhi() >= 1) m(4,5) = T(7);
            m(2,2) = T(-10);
            m.diag() *= T(30);
        }

        tmv::BandMatrix<std::complex<T>,stor> c(M,N,nlo,nhi);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            if (i<=j+nlo && j<=i+nhi) c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        if (!singular) {
            c(0,0) = T(14);
            if (c.nlo() >= 1) c(1,0) = T(-2);
            if (c.nhi() >= 1) c(4,5) = T(7);
            c(2,2) = T(-10);
            c.diag() *= T(30);
        }

        if (mattype == 7) {
            for(int i=10;i<N;i+=10) {
                m.colRange(i,N) *= T(1.e-10);
                c.colRange(i,N) *= T(1.e-10);
                m.rowRange(i+5,N) *= T(-1.e-10);
                c.rowRange(i+5,N) *= T(-1.e-10);
            }
        }

        if (mattype == 8) {
            T x = std::numeric_limits<T>::min();
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,2*x));
            for(int i=1;i<N;++i) {
                int start = i-nhi;  if (start < 0) start = 0;
                int end = i+nlo;  if (end > N) end = N;
                m.col(i,start,end) /= T(i);
                c.col(i,start,end) /= T(i);
            }
            m(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            c(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
        }

        if (mattype == 9) {
            T x = std::numeric_limits<T>::max();
            x /= N;
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,x/2));
            for(int i=1;i<N;++i) {
                int start = i-nhi;  if (start < 0) start = 0;
                int end = i+nlo;  if (end > N) end = N;
                m.col(i,start,end) /= T(i);
                c.col(i,start,end) /= T(i);
            }
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2,N) *= T(-1);
            c.diag(0,N/2,N) *= T(-1);
        }

        T eps = EPS;
        T ceps = EPS;
        T normm = Norm(m);
        T normc = Norm(c);
        if (showacc) {
            std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
            std::cout<<"norm = "<<normm<<"  "<<normc<<std::endl;
        }
        if (!singular) {
            T kappa = normm * Norm(m.inverse());
            T ckappa = normc * Norm(c.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<"  "<<ckappa<<std::endl;
            eps *= kappa;
            ceps *= ckappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;


        // LU Decomposition
        if (m.isSquare()) {
            if (showstartdone) {
                std::cout<<"LU\n";
            }
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            tmv::Permutation P = m.lud().getP();
            tmv::Matrix<T> PLU = P*L*U;
            if (m.lud().isTrans()) PLU.transposeSelf();
            Assert(Norm(m-PLU) <= eps*normm,"Band LU");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = 
                c.lud().getL();
            tmv::UpperTriMatrix<std::complex<T> > cU = c.lud().getU();
            tmv::Permutation cP = c.lud().getP();
            tmv::Matrix<std::complex<T> > cPLU = cP*cL*cU;
            if (c.lud().isTrans()) cPLU.transposeSelf();
            Assert(Norm(c-cPLU) <= ceps*normc,"Band C LU");

#ifdef XTEST
            const int Rnhi = std::min(N-1,nlo+nhi);
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L2(M);
            tmv::BandMatrix<T,stor> U2(M,M,0,Rnhi);
            tmv::Permutation P2(N);
            LU_Decompose(m,L2.view(),U2.view(),P2);
            PLU = P2*L2*U2;
            Assert(Norm(m-PLU) <= eps*normm,"Band LU2");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL2(M);
            tmv::BandMatrix<std::complex<T>,stor> cU2(M,M,0,Rnhi);
            LU_Decompose(c,cL2.view(),cU2.view(),P2);
            cPLU = P2*cL2*cU2;
            Assert(Norm(c-cPLU) <= ceps*normc,"Band C LU2");

            LU_Decompose(c,cL2.conjugate(),cU2.view(),P2);
            cPLU = P2*cL2.conjugate()*cU2;
            Assert(Norm(c-cPLU) <= ceps*normc,"Band C LU3");

            LU_Decompose(c,cL2.view(),cU2.conjugate(),P2);
            cPLU = P2*cL2*cU2.conjugate();
            Assert(Norm(c-cPLU) <= ceps*normc,"Band C LU4");

            LU_Decompose(c,cL2.conjugate(),cU2.conjugate(),P2);
            cPLU = P2*cL2.conjugate()*cU2.conjugate();
            Assert(Norm(c-cPLU) <= ceps*normc,"Band C LU5");

            LU_Decompose(c.conjugate(),cL2.view(),cU2.view(),P2);
            cPLU = P2*cL2*cU2;
            Assert(Norm(c.conjugate()-cPLU) <= ceps*normc,"Band C LU6");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2.view(),P2);
            cPLU = P2*cL2.conjugate()*cU2;
            Assert(Norm(c.conjugate()-cPLU) <= ceps*normc,"Band C LU7");

            LU_Decompose(c.conjugate(),cL2.view(),cU2.conjugate(),P2);
            cPLU = P2*cL2*cU2.conjugate();
            Assert(Norm(c.conjugate()-cPLU) <= ceps*normc,"Band C LU8");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2.conjugate(),P2);
            cPLU = P2*cL2.conjugate()*cU2.conjugate();
            Assert(Norm(c.conjugate()-cPLU) <= ceps*normc,"Band C LU9");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // QR Decomposition
        {
            if (showstartdone) {
                std::cout<<"QR\n";
            }
            tmv::Matrix<T> Q = m.qrd().getQ();
            tmv::BandMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            if (m.qrd().isTrans()) QR.transposeSelf();
            Assert(Norm(m-QR) <= eps*normm,"Band QR");

            tmv::Matrix<std::complex<T> > cQ = c.qrd().getQ();
            tmv::BandMatrix<std::complex<T> > cR = c.qrd().getR();
            tmv::Matrix<std::complex<T> > cQR = cQ*cR;
            if (c.qrd().isTrans()) cQR.transposeSelf();
            Assert(Norm(c-cQR) <= ceps*normc,"Band C QR");

#ifdef XTEST
            const int Rnhi = std::min(N-1,nlo+nhi);
            QR_Decompose(m,Q.view(),R.view());
            QR = Q*R;
            Assert(Norm(m-QR) <= eps*normm,"Band QR2");

            tmv::BandMatrix<T,stor> R2(N,N,0,Rnhi);
            QR_Decompose(m,R2.view());
            Assert(Norm(R-R2) <= eps*normm,"Band QR3");

            QR_Decompose(c,cQ.view(),cR.view());
            cQR = cQ*cR;
            Assert(Norm(c-cQR) <= ceps*normc,"Band C QR2");

            tmv::BandMatrix<std::complex<T>,stor> cR2(N,N,0,Rnhi);
            QR_Decompose(c,cR2.view());
            Assert(Norm(cR-cR2) <= ceps*normc,"Band C QR3");

            QR_Decompose(c,cQ.conjugate(),cR2.view());
            cQR = cQ.conjugate()*cR2;
            Assert(Norm(c-cQR) <= ceps*normc,"Band C QR4");

            QR_Decompose(c,cQ.view(),cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Norm(c-cQR) <= ceps*normc,"Band C QR5");

            QR_Decompose(c,cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Norm(c-cQR) <= ceps*normc,"Band C QR6");

            QR_Decompose(c,cR2.view());
            Assert(Norm(cR-cR2) <= ceps*normc,"Band C QR7");

            QR_Decompose(c,cR2.conjugate());
            Assert(Norm(cR-cR2.conjugate()) <= ceps*normc,"Band C QR8");

            QR_Decompose(c.conjugate(),cQ.view(),cR2.view());
            cQR = cQ*cR2;
            Assert(Norm(c.conjugate()-cQR) <= ceps*normc,"Band C QR9");

            QR_Decompose(c.conjugate(),cR2.view());
            Assert(Norm(cR.conjugate()-cR2) <= ceps*normc,"Band C QR10");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2.view());
            cQR = cQ.conjugate()*cR2;
            Assert(Norm(c.conjugate()-cQR) <= ceps*normc,"Band C QR11");

            QR_Decompose(c.conjugate(),cQ.view(),cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Norm(c.conjugate()-cQR) <= ceps*normc,"Band C QR12");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Norm(c.conjugate()-cQR) <= ceps*normc,"Band C QR13");

            QR_Decompose(c.conjugate(),cR2.view());
            Assert(Norm(cR.conjugate()-cR2) <= ceps*normc,"Band C QR14");

            QR_Decompose(c.conjugate(),cR2.conjugate());
            Assert(Norm(cR.conjugate()-cR2.conjugate()) <= ceps*normc,
                   "Band C QR15");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // SV Decomposition
        {
            if (showstartdone) {
                std::cout<<"SV\n";
            }
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            if (showacc) {
                std::cout<<"U = "<<U<<std::endl;
                std::cout<<"V = "<<V<<std::endl;
                std::cout<<"S = "<<S<<std::endl;
                std::cout<<"m-USV = "<<m-U*S*V<<std::endl;
                std::cout<<"Norm(m-USV) = "<<Norm(m-U*S*V)<<
                    "  cf "<<eps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.transpose()*U-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(V.transpose()*V-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(V*V.transpose()-T(1))<<
                    "  cf "<<eps<<std::endl;
            }
            Assert(Norm(m-U*S*V) <= eps*normm,"SV");

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            Assert(Norm(c-cU*cS*cV) <= ceps*normc,"C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(M,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) <= eps*normm,"SV2");

            SV_Decompose(m,S2.view());
            Assert(Norm(S2-S) <= eps*normm,"SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) <= eps*normm,"SV4 S");
            if (!nearoverflow)
                Assert(Norm(U2*S2*S2*U2.transpose()-m*m.transpose()) <= 
                       eps*normm*normm,"SV3 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) <= eps*normm,"SV5 S");
            if (!nearoverflow)
                Assert(Norm(V2.transpose()*S2*S2*V2-m.transpose()*m) <= 
                       eps*normm*normm,"SV5 V");

            tmv::Matrix<std::complex<T> > cU2(M,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) <= eps*normc,"C SV2");

            SV_Decompose(c,cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV4 S");
            if (!nearoverflow)
                Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) <= 
                       ceps*normc*normc,"C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV5 S");
            if (!nearoverflow)
                Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) <= 
                       ceps*normc*normc,"C SV5 V");

            SV_Decompose(c,cU2.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2.conjugate()*cS2*cV2) <= eps*normc,"C SV6");
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.conjugate());
            Assert(Norm(c-cU2*cS2*cV2.conjugate()) <= eps*normc,"C SV7");
            SV_Decompose(c,cU2.conjugate(),cS2.view(),cV2.conjugate());
            Assert(Norm(c-cU2.conjugate()*cS2*cV2.conjugate()) <= eps*normc,
                   "C SV8");

            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV9 S");
            if (!nearoverflow)
                Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) <= 
                       ceps*normc*normc,"C SV9 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV10 S");
            if (!nearoverflow)
                Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) <= 
                       ceps*normc*normc,"C SV10 V");

            SV_Decompose(c.conjugate(),cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c.conjugate()-cU2*cS2*cV2) <= eps*normc,"C SV11");
            SV_Decompose(c.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV12");
            SV_Decompose(c.conjugate(),cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV13 S");
            if (!nearoverflow)
                Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-
                            c.conjugate()*c.transpose()) <= 
                       ceps*normc*normc,"C SV13 U");
            SV_Decompose(c.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV14 S");
            if (!nearoverflow)
                Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-
                            c.transpose()*c.conjugate()) <= 
                       ceps*normc*normc,"C SV14 V");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view(),cV2.view());
            Assert(Norm(c.conjugate()-cU2.conjugate()*cS2*cV2) <= eps*normc,
                   "C SV15");
            SV_Decompose(c.conjugate(),cU2.view(),cS2.view(),cV2.conjugate());
            Assert(Norm(c.conjugate()-cU2*cS2*cV2.conjugate()) <= eps*normc,
                   "C SV16");
            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view(),
                         cV2.conjugate());
            Assert(Norm(c.conjugate()-cU2.conjugate()*cS2*cV2.conjugate()) <= 
                   eps*normc,"C SV17");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV18 S");
            if (!nearoverflow)
                Assert(Norm(cU2.conjugate()*cS2*cS2*cU2.transpose()-
                            c.conjugate()*c.transpose()) <= ceps*normc*normc,
                       "C SV18 U");
            SV_Decompose(c.conjugate(),cS2.view(),cV2.conjugate());
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV19 S");
            if (!nearoverflow)
                Assert(Norm(cV2.transpose()*cS2*cS2*cV2.conjugate()-
                            c.transpose()*c.conjugate()) <= ceps*normc*normc,
                       "C SV19 V");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef TEST_DOUBLE
template void TestBandDecomp<double,tmv::ColMajor>();
template void TestBandDecomp<double,tmv::RowMajor>();
template void TestBandDecomp<double,tmv::DiagMajor>();
#endif
#ifdef TEST_FLOAT
template void TestBandDecomp<float,tmv::ColMajor>();
template void TestBandDecomp<float,tmv::RowMajor>();
template void TestBandDecomp<float,tmv::DiagMajor>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDecomp<long double,tmv::ColMajor>();
template void TestBandDecomp<long double,tmv::RowMajor>();
template void TestBandDecomp<long double,tmv::DiagMajor>();
#endif

