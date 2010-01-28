
#include "../src/TMV_Blas.h"
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"

template <class T, tmv::StorageType stor> 
void TestMatrixDecomp()
{
    for (int mattype = 0; mattype < 4; mattype++) {
        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<", stor = "<<
                tmv::TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Square
        // mattype = 1  is NonSquare slightly tall
        // mattype = 2  is NonSquare very tall
        // mattype = 3  is Singular

        const int N = 200;
        int M = N;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;

        tmv::Matrix<T,stor> m(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
        if (mattype != 3) m /= T(10*N);
        m(3,3) = 14;
        m(3,4) = -2;
        m(0,2) = 7;
        m(M-1,N-4) = 23;
        m(M-1,N-2) = 13;
        m(M-1,N-1) = -10;
        if (mattype != 3) m.diag() *= T(30*N);

        tmv::Matrix<std::complex<T>,stor> c(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        if (mattype != 3) c /= T(10*N);
        c(3,3) = 14;
        c(3,4) = -2;
        c(0,2) = 7;
        c(M-1,N-4) = 23;
        c(M-1,N-2) = 13;
        c(M-1,N-1) = -10;
        if (mattype != 3) c.diag() *= T(30*N);

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

        if (showstartdone) {
            //std::cout<<"m = "<<m<<std::endl;
        }

        // LU Decomposition
        if (mattype == 0) {
            if (showstartdone) {
                std::cout<<"LU\n";
            }
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            const int* p = m.lud().getP();
            tmv::Matrix<T> PLU = L*U;
            PLU.reversePermuteRows(p);
            if (m.lud().isTrans()) PLU.transposeSelf();
            Assert(Norm(m-PLU) < eps*Norm(m),"LU"); 

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL =
                c.lud().getL();
            tmv::UpperTriMatrix<std::complex<T> > cU = c.lud().getU();
            p = c.lud().getP();
            tmv::Matrix<std::complex<T> > cPLU = cL*cU;
            cPLU.reversePermuteRows(p);
            if (c.lud().isTrans()) cPLU.transposeSelf();
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU"); 

#ifdef XTEST
            tmv::Matrix<T,stor> m2 = m;
            int p2[N];
            LU_Decompose(m2.view(),p2);
            L = m2.lowerTri(tmv::UnitDiag);
            U = m2.upperTri();
            PLU = L*U;
            PLU.reversePermuteRows(p2);
            Assert(Norm(m-PLU) < eps*Norm(m),"LU2"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            LU_Decompose(c2.view(),p2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cU = c2.upperTri();
            cPLU = cL*cU;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU2"); 

            c2.conjugate() = c;
            LU_Decompose(c2.conjugate(),p2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cU = c2.conjugate().upperTri();
            cPLU = cL*cU;
            cPLU.reversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU3"); 
#endif
            std::cout<<"."; std::cout.flush();
        }

        // QR Decomposition
        if (mattype != 4) {
            if (showstartdone) {
                std::cout<<"QR\n";
            }
            tmv::Matrix<T,stor> Q = m.qrd().getQ();
            tmv::UpperTriMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            Assert(Norm(m-QR) < eps*Norm(m),"QR"); 
            Assert(Norm(m-m.qrd().getQ()*m.qrd().getR()) < eps*Norm(m),
                   "QR - PackedQ"); 
            Assert(Norm(Q.transpose()*Q-T(1)) < T(N)*eps,"QR - QtQ"); 
            Assert(Norm(Q.transpose()*m.qrd().getQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (second)"); 
            Assert(Norm(Q / m.qrd().getQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (first)"); 

            tmv::Matrix<std::complex<T>,stor> cQ = c.qrd().getQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.qrd().getR();
            tmv::Matrix<std::complex<T> > cQR = cQ*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR"); 
            Assert(Norm(c-c.qrd().getQ()*c.qrd().getR()) < ceps*Norm(c),
                   "C QR - PackedQ"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QR - QtQ"); 
            Assert(Norm(cQ.adjoint()*c.qrd().getQ()-T(1)) < T(N)*ceps,
                   "C QR - QtQ - PackedQ (second)"); 
            Assert(Norm(cQ / c.qrd().getQ()-T(1)) < T(N)*ceps,
                   "C QR - QtQ - PackedQ (first)"); 

#ifdef XTEST
            Q = m;
            QR_Decompose(Q.view(),R.view());
            QR = Q*R;
            Assert(Norm(m-QR) < eps*Norm(m),"QR2"); 
            Assert(Norm(Q.transpose()*Q-T(1)) < T(N)*eps,"QR2 - QtQ"); 

            Q = m;
            QR_Decompose(Q.view());
            Assert(Norm(R-Q.upperTri()) < eps*Norm(R),"QR3"); 

            cQ = c;
            QR_Decompose(cQ.view(),cR.view());
            cQR = cQ*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR2"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QR2 - QtQ"); 

            cQ = c;
            QR_Decompose(cQ.view());
            Assert(Norm(cR-cQ.upperTri()) < ceps*Norm(cR),"C QR3"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR.view());
            cQR = cQ.conjugate()*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR4"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QR4 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate());
            Assert(Norm(cR-cQ.conjugate().upperTri()) < ceps*Norm(cR),"C QR5"); 

            cQ = c;
            QR_Decompose(cQ.view(),cR.conjugate());
            cQR = cQ*cR.conjugate();
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR6"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QR6 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR.conjugate());
            cQR = cQ.conjugate()*cR.conjugate();
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR7"); 
#endif
            std::cout<<"."; std::cout.flush();
        }

        // QRP Decomposition
        for (int istrict = 0; istrict <= 1; istrict++) {
            if (showstartdone) {
                std::cout<<"QRP\n";
            }
            bool strict = istrict == 1;
            tmv::QRPDiv<T>::StrictQRP = strict;
            tmv::QRPDiv<std::complex<T> >::StrictQRP = strict;
            if (istrict == 1) { m.unsetDiv(); c.unsetDiv(); }

            tmv::Matrix<T,stor> Q = m.qrpd().getQ();
            tmv::UpperTriMatrix<T> R = m.qrpd().getR();
            const int* p = m.qrpd().getP();
            tmv::Matrix<T> QRP = Q*R;
            QRP.reversePermuteCols(p);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP"); 
            QRP = m.qrpd().getQ()*R;
            QRP.reversePermuteCols(p);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP - Packed"); 
            Assert(Norm(Q.transpose()*Q-T(1)) < T(N)*eps,"QRP - QtQ"); 
            Assert(Norm(Q.transpose()*m.qrpd().getQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (first)"); 
            Assert(Norm(Q / m.qrpd().getQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (second)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP - strict"); 
            }
#endif

            tmv::Matrix<std::complex<T>,stor> cQ = c.qrpd().getQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.qrpd().getR();
            p = c.qrpd().getP();
            tmv::Matrix<std::complex<T> > cQRP = cQ*cR;
            cQRP.reversePermuteCols(p);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP - QtQ"); 
            Assert(Norm(cQ.adjoint()*c.qrpd().getQ()-T(1)) < T(N)*ceps,
                   "C QRP - QtQ - PackedQ (first)"); 
            Assert(Norm(cQ / c.qrpd().getQ()-T(1)) < T(N)*ceps,
                   "C QRP - QtQ - PackedQ (second)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP - strict"); 
            }
#endif

#ifdef XTEST
            Q = m;
            int p2[N];
            QRP_Decompose(Q.view(),R.view(),p2,strict);
            QRP = Q*R;
            QRP.reversePermuteCols(p2);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP2"); 
            Assert(Norm(Q.transpose()*Q-T(1)) < T(N)*eps,"QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP2 - strict"); 
            }
#endif

            Q = m;
            QRP_Decompose(Q.view(),strict);
            Assert(Norm(R-Q.upperTri()) < eps*Norm(R),"QRP3"); 

            cQ = c;
            QRP_Decompose(cQ.view(),cR.view(),p2,strict);
            cQRP = cQ*cR;
            cQRP.reversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP2"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP2 - strict"); 
            }
#endif
            cQ = c;
            QRP_Decompose(cQ.view(),strict);
            Assert(Norm(cR-cQ.upperTri()) < ceps*Norm(cR),"C QRP3"); 

            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),cR.view(),p2,strict);
            cQRP = cQ.conjugate()*cR;
            cQRP.reversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP4"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP4 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP4 - strict"); 
            }
#endif

            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),strict);
            Assert(Norm(cR-cQ.conjugate().upperTri()) < ceps*Norm(cR),"C QRP5"); 

            cQ = c;
            QRP_Decompose(cQ.view(),cR.conjugate(),p2,strict);
            cQRP = cQ*cR.conjugate();
            cQRP.reversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP6"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP6 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP6 - strict"); 
            }
#endif


            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),cR.conjugate(),p2,strict);
            cQRP = cQ.conjugate()*cR.conjugate();
            cQRP.reversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP7"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP7 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP7 - strict"); 
            }
#endif
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
                std::cout<<"Norm(m-USV) = "<<Norm(m-U*S*V)<<
                    "  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.transpose()*U-T(1))<<
                    "  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(V.transpose()*V-T(1))<<
                    "  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(V*V.transpose()-T(1))<<
                    "  cf "<<eps*Norm(m)<<std::endl;
            }
            Assert(Norm(m-U*S*V) < eps*Norm(m),"SV"); 
            Assert(Norm(U.transpose()*U-T(1)) < eps*Norm(m),"SV - UtU"); 
            Assert(Norm(V.transpose()*V-T(1)) < eps*Norm(m),"SV - VtV"); 
            Assert(Norm(V*V.transpose()-T(1)) < eps*Norm(m),"SV - VVt"); 

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            if (showacc) {
                std::cout<<"Norm(c-USV) = "<<Norm(c-cU*cS*cV)<<
                    "  cf "<<ceps*Norm(m)<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<
                    "  cf "<<ceps*Norm(c)<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(cV.adjoint()*cV-T(1))<<
                    "  cf "<<ceps*Norm(c)<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(cV*cV.adjoint()-T(1))<<
                    "  cf "<<ceps*Norm(c)<<std::endl;
            }
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"C SV"); 
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps*Norm(c),"C SV - UtU"); 
            Assert(Norm(cV.adjoint()*cV-T(1)) < ceps*Norm(c),"C SV - VtV"); 
            Assert(Norm(cV*cV.adjoint()-T(1)) < ceps*Norm(c),"C SV - VVt"); 

#ifdef XTEST
            tmv::Matrix<T,stor> U2 = m;
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(U2.view(),S2.view(),V2.view(),true);
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SV2"); 
            Assert(Norm(U2.transpose()*U2-T(1)) < eps*Norm(m),"SV2 - UtU"); 
            Assert(Norm(V2.transpose()*V2-T(1)) < eps*Norm(m),"SV2 - VtV"); 
            Assert(Norm(V2*V2.transpose()-T(1)) < eps*Norm(m),"SV2 - VVt"); 

            tmv::Matrix<T,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view(),false);
            Assert(Norm(S2-S) < eps*Norm(S),"SV3"); 
            U2 = m;
            SV_Decompose(U2.view(),S2.view(),true);
            Assert(Norm(S2-S) < eps*Norm(S),"SV4 S"); 
            Assert(Norm(U2*S2*S2*U2.transpose()-m*m.transpose()) <  
                   eps*Norm(m*m.transpose()),"SV4 U"); 
            m2 = m;
            SV_Decompose(m2.view(),S2.view(),V2.view(),false);
            Assert(Norm(S2-S) < eps*Norm(S),"SV5 S"); 
            Assert(Norm(V2.transpose()*S2*S2*V2-m.transpose()*m) < 
                   eps*Norm(m.transpose()*m),"SV5 V"); 

            tmv::Matrix<std::complex<T>,stor> cU2 = c;
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(cU2.view(),cS2.view(),cV2.view(),true);
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(m),"C SV2"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV2 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV2 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) < ceps*Norm(m),"C SV2 - VVt"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV3"); 
            cU2 = c;
            SV_Decompose(cU2.view(),cS2.view(),true);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV4 S"); 
            Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) <  
                   ceps*Norm(c*c.adjoint()),"C SV4 U"); 
            c2 = c;
            SV_Decompose(c2.view(),cS2.view(),cV2.view(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV5 S"); 
            Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) < 
                   ceps*Norm(c.adjoint()*c),"C SV5 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),cV2.view(),true);
            Assert(Norm(c-cU2.conjugate()*cS2*cV2) < ceps*Norm(m),"C SV6"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV6 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV6 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) < ceps*Norm(m),"C SV6 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV7"); 
            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),true);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV8 S"); 
            Assert(
                Norm(cU2.conjugate()*cS2*cS2*cU2.transpose()-c*c.adjoint()) <
                ceps*Norm(c*c.adjoint()),"C SV8 U"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),cV2.view(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV9 S"); 
            Assert(Norm(cV2.adjoint()*cS2*cS2*cV2-c.adjoint()*c) < 
                   ceps*Norm(c.adjoint()*c),"C SV9 V"); 

            cU2 = c;
            SV_Decompose(cU2.view(),cS2.view(),cV2.conjugate(),true);
            Assert(Norm(c-cU2*cS2*cV2.conjugate()) < ceps*Norm(m),"C SV10"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV10 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV10 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) < ceps*Norm(m),"C SV10 - VVt"); 
            c2 = c;
            SV_Decompose(c2.view(),cS2.view(),cV2.conjugate(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV11 S"); 
            Assert(
                Norm(cV2.transpose()*cS2*cS2*cV2.conjugate()-c.adjoint()*c) < 
                ceps*Norm(c.adjoint()*c),"C SV9 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),cV2.conjugate(),true);
            Assert(
                Norm(c-cU2.conjugate()*cS2*cV2.conjugate()) < ceps*Norm(m),
                "C SV12"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV12 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV12 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) < ceps*Norm(m),"C SV12 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),cV2.conjugate(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV13 S"); 
            Assert(
                Norm(cV2.transpose()*cS2*cS2*cV2.conjugate()-c.adjoint()*c) < 
                ceps*Norm(c.adjoint()*c),"C SV13 V"); 
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef INST_DOUBLE
template void TestMatrixDecomp<double,tmv::RowMajor>();
template void TestMatrixDecomp<double,tmv::ColMajor>();
#endif
#ifdef INST_FLOAT
template void TestMatrixDecomp<float,tmv::RowMajor>();
template void TestMatrixDecomp<float,tmv::ColMajor>();
#endif
#ifdef INST_LONGDOUBLE
template void TestMatrixDecomp<long double,tmv::RowMajor>();
template void TestMatrixDecomp<long double,tmv::ColMajor>();
#endif
#ifdef INST_INT
template void TestMatrixDecomp<int,tmv::RowMajor>();
template void TestMatrixDecomp<int,tmv::ColMajor>();
#endif
