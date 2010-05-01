
#include "../src/TMV_Blas.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

template <class T, tmv::StorageType stor> void TestMatrixDecomp()
{
    for (int mattype = 0; mattype < 4; mattype++) {
        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<", stor = "<<tmv::TMV_Text(stor)<<std::endl;
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
            T kappa = Norm(m) * Norm(m.Inverse());
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
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
            tmv::UpperTriMatrix<T> U = m.LUD().GetU();
            const int* p = m.LUD().GetP();
            tmv::Matrix<T> PLU = L*U;
            PLU.ReversePermuteRows(p);
            if (m.LUD().IsTrans()) PLU.TransposeSelf();
            Assert(Norm(m-PLU) < eps*Norm(m),"LU"); 

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
            tmv::UpperTriMatrix<std::complex<T> > cU = c.LUD().GetU();
            p = c.LUD().GetP();
            tmv::Matrix<std::complex<T> > cPLU = cL*cU;
            cPLU.ReversePermuteRows(p);
            if (c.LUD().IsTrans()) cPLU.TransposeSelf();
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU"); 

#if (XTEST & 2)
            tmv::Matrix<T,stor> m2 = m;
            int p2[N];
            LU_Decompose(m2.View(),p2);
            L = m2.LowerTri(tmv::UnitDiag);
            U = m2.UpperTri();
            PLU = L*U;
            PLU.ReversePermuteRows(p2);
            Assert(Norm(m-PLU) < eps*Norm(m),"LU2"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            LU_Decompose(c2.View(),p2);
            cL = c2.LowerTri(tmv::UnitDiag);
            cU = c2.UpperTri();
            cPLU = cL*cU;
            cPLU.ReversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU2"); 

            c2.Conjugate() = c;
            LU_Decompose(c2.Conjugate(),p2);
            cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
            cU = c2.Conjugate().UpperTri();
            cPLU = cL*cU;
            cPLU.ReversePermuteRows(p2);
            Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU3"); 
#endif
            std::cout<<"."; std::cout.flush();
        }

        // QR Decomposition
        if (mattype != 4) {
            if (showstartdone) {
                std::cout<<"QR\n";
            }
            tmv::Matrix<T,stor> Q = m.QRD().GetQ();
            tmv::UpperTriMatrix<T> R = m.QRD().GetR();
            tmv::Matrix<T> QR = Q*R;
            Assert(Norm(m-QR) < eps*Norm(m),"QR"); 
            Assert(Norm(m-m.QRD().GetQ()*m.QRD().GetR()) < eps*Norm(m),
                   "QR - PackedQ"); 
            Assert(Norm(Q.Transpose()*Q-T(1)) < T(N)*eps,"QR - QtQ"); 
            Assert(Norm(Q.Transpose()*m.QRD().GetQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (second)"); 
            Assert(Norm(Q / m.QRD().GetQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (first)"); 

            tmv::Matrix<std::complex<T>,stor> cQ = c.QRD().GetQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.QRD().GetR();
            tmv::Matrix<std::complex<T> > cQR = cQ*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR"); 
            Assert(Norm(c-c.QRD().GetQ()*c.QRD().GetR()) < ceps*Norm(c),
                   "C QR - PackedQ"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QR - QtQ"); 
            Assert(Norm(cQ.Adjoint()*c.QRD().GetQ()-T(1)) < T(N)*ceps,
                   "C QR - QtQ - PackedQ (second)"); 
            Assert(Norm(cQ / c.QRD().GetQ()-T(1)) < T(N)*ceps,
                   "C QR - QtQ - PackedQ (first)"); 

#if (XTEST & 2)
            Q = m;
            QR_Decompose(Q.View(),R.View());
            QR = Q*R;
            Assert(Norm(m-QR) < eps*Norm(m),"QR2"); 
            Assert(Norm(Q.Transpose()*Q-T(1)) < T(N)*eps,"QR2 - QtQ"); 

            Q = m;
            QR_Decompose(Q.View());
            Assert(Norm(R-Q.UpperTri()) < eps*Norm(R),"QR3"); 

            cQ = c;
            QR_Decompose(cQ.View(),cR.View());
            cQR = cQ*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR2"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QR2 - QtQ"); 

            cQ = c;
            QR_Decompose(cQ.View());
            Assert(Norm(cR-cQ.UpperTri()) < ceps*Norm(cR),"C QR3"); 

            cQ.Conjugate() = c;
            QR_Decompose(cQ.Conjugate(),cR.View());
            cQR = cQ.Conjugate()*cR;
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR4"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QR4 - QtQ"); 

            cQ.Conjugate() = c;
            QR_Decompose(cQ.Conjugate());
            Assert(Norm(cR-cQ.Conjugate().UpperTri()) < ceps*Norm(cR),"C QR5"); 

            cQ = c;
            QR_Decompose(cQ.View(),cR.Conjugate());
            cQR = cQ*cR.Conjugate();
            Assert(Norm(c-cQR) < ceps*Norm(c),"C QR6"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QR6 - QtQ"); 

            cQ.Conjugate() = c;
            QR_Decompose(cQ.Conjugate(),cR.Conjugate());
            cQR = cQ.Conjugate()*cR.Conjugate();
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
            if (istrict == 1) { m.UnSetDiv(); c.UnSetDiv(); }

            tmv::Matrix<T,stor> Q = m.QRPD().GetQ();
            tmv::UpperTriMatrix<T> R = m.QRPD().GetR();
            const int* p = m.QRPD().GetP();
            tmv::Matrix<T> QRP = Q*R;
            QRP.ReversePermuteCols(p);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP"); 
            QRP = m.QRPD().GetQ()*R;
            QRP.ReversePermuteCols(p);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP - Packed"); 
            Assert(Norm(Q.Transpose()*Q-T(1)) < T(N)*eps,"QRP - QtQ"); 
            Assert(Norm(Q.Transpose()*m.QRPD().GetQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (first)"); 
            Assert(Norm(Q / m.QRPD().GetQ()-T(1)) < T(N)*eps,
                   "QR - QtQ - PackedQ (second)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP - strict"); 
            }
#endif

            tmv::Matrix<std::complex<T>,stor> cQ = c.QRPD().GetQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.QRPD().GetR();
            p = c.QRPD().GetP();
            tmv::Matrix<std::complex<T> > cQRP = cQ*cR;
            cQRP.ReversePermuteCols(p);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP - QtQ"); 
            Assert(Norm(cQ.Adjoint()*c.QRPD().GetQ()-T(1)) < T(N)*ceps,
                   "C QRP - QtQ - PackedQ (first)"); 
            Assert(Norm(cQ / c.QRPD().GetQ()-T(1)) < T(N)*ceps,
                   "C QRP - QtQ - PackedQ (second)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP - strict"); 
            }
#endif

#if (XTEST & 2)
            Q = m;
            int p2[N];
            QRP_Decompose(Q.View(),R.View(),p2,strict);
            QRP = Q*R;
            QRP.ReversePermuteCols(p2);
            Assert(Norm(m-QRP) < eps*Norm(m),"QRP2"); 
            Assert(Norm(Q.Transpose()*Q-T(1)) < T(N)*eps,"QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP2 - strict"); 
            }
#endif

            Q = m;
            QRP_Decompose(Q.View(),strict);
            Assert(Norm(R-Q.UpperTri()) < eps*Norm(R),"QRP3"); 

            cQ = c;
            QRP_Decompose(cQ.View(),cR.View(),p2,strict);
            cQRP = cQ*cR;
            cQRP.ReversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP2"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP2 - strict"); 
            }
#endif
            cQ = c;
            QRP_Decompose(cQ.View(),strict);
            Assert(Norm(cR-cQ.UpperTri()) < ceps*Norm(cR),"C QRP3"); 

            cQ.Conjugate() = c;
            QRP_Decompose(cQ.Conjugate(),cR.View(),p2,strict);
            cQRP = cQ.Conjugate()*cR;
            cQRP.ReversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP4"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP4 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP4 - strict"); 
            }
#endif

            cQ.Conjugate() = c;
            QRP_Decompose(cQ.Conjugate(),strict);
            Assert(Norm(cR-cQ.Conjugate().UpperTri()) < ceps*Norm(cR),"C QRP5"); 

            cQ = c;
            QRP_Decompose(cQ.View(),cR.Conjugate(),p2,strict);
            cQRP = cQ*cR.Conjugate();
            cQRP.ReversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP6"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP6 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP6 - strict"); 
            }
#endif


            cQ.Conjugate() = c;
            QRP_Decompose(cQ.Conjugate(),cR.Conjugate(),p2,strict);
            cQRP = cQ.Conjugate()*cR.Conjugate();
            cQRP.ReversePermuteCols(p2);
            Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP7"); 
            Assert(Norm(cQ.Adjoint()*cQ-T(1)) < T(N)*ceps,"C QRP7 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP7 - strict"); 
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
            tmv::Matrix<T> U = m.SVD().GetU();
            tmv::DiagMatrix<T> S = m.SVD().GetS();
            tmv::Matrix<T> V = m.SVD().GetV();
            if (showacc) {
                std::cout<<"Norm(m-USV) = "<<Norm(m-U*S*V)<<"  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.Transpose()*U-T(1))<<"  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(V.Transpose()*V-T(1))<<"  cf "<<eps*Norm(m)<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(V*V.Transpose()-T(1))<<"  cf "<<eps*Norm(m)<<std::endl;
            }
            Assert(Norm(m-U*S*V) < eps*Norm(m),"SV"); 
            Assert(Norm(U.Transpose()*U-T(1)) < eps*Norm(m),"SV - UtU"); 
            Assert(Norm(V.Transpose()*V-T(1)) < eps*Norm(m),"SV - VtV"); 
            Assert(Norm(V*V.Transpose()-T(1)) < eps*Norm(m),"SV - VVt"); 

            tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
            tmv::DiagMatrix<T> cS = c.SVD().GetS();
            tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
            if (showacc) {
                std::cout<<"Norm(c-USV) = "<<Norm(c-cU*cS*cV)<<"  cf "<<ceps*Norm(m)<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(cU.Adjoint()*cU-T(1))<<"  cf "<<ceps*Norm(c)<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(cV.Adjoint()*cV-T(1))<<"  cf "<<ceps*Norm(c)<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(cV*cV.Adjoint()-T(1))<<"  cf "<<ceps*Norm(c)<<std::endl;
            }
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"C SV"); 
            Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps*Norm(c),"C SV - UtU"); 
            Assert(Norm(cV.Adjoint()*cV-T(1)) < ceps*Norm(c),"C SV - VtV"); 
            Assert(Norm(cV*cV.Adjoint()-T(1)) < ceps*Norm(c),"C SV - VVt"); 

#if (XTEST & 2)
            tmv::Matrix<T,stor> U2 = m;
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(U2.View(),S2.View(),V2.View(),true);
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SV2"); 
            Assert(Norm(U2.Transpose()*U2-T(1)) < eps*Norm(m),"SV2 - UtU"); 
            Assert(Norm(V2.Transpose()*V2-T(1)) < eps*Norm(m),"SV2 - VtV"); 
            Assert(Norm(V2*V2.Transpose()-T(1)) < eps*Norm(m),"SV2 - VVt"); 

            tmv::Matrix<T,stor> m2 = m;
            SV_Decompose(m2.View(),S2.View(),false);
            Assert(Norm(S2-S) < eps*Norm(S),"SV3"); 
            U2 = m;
            SV_Decompose(U2.View(),S2.View(),true);
            Assert(Norm(S2-S) < eps*Norm(S),"SV4 S"); 
            Assert(Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose()) <  
                   eps*Norm(m*m.Transpose()),"SV4 U"); 
            m2 = m;
            SV_Decompose(m2.View(),S2.View(),V2.View(),false);
            Assert(Norm(S2-S) < eps*Norm(S),"SV5 S"); 
            Assert(Norm(V2.Transpose()*S2*S2*V2-m.Transpose()*m) < 
                   eps*Norm(m.Transpose()*m),"SV5 V"); 

            tmv::Matrix<std::complex<T>,stor> cU2 = c;
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(cU2.View(),cS2.View(),cV2.View(),true);
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(m),"C SV2"); 
            Assert(Norm(cU2.Adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV2 - UtU"); 
            Assert(Norm(cV2.Adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV2 - VtV"); 
            Assert(Norm(cV2*cV2.Adjoint()-T(1)) < ceps*Norm(m),"C SV2 - VVt"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            SV_Decompose(c2.View(),cS2.View(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV3"); 
            cU2 = c;
            SV_Decompose(cU2.View(),cS2.View(),true);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV4 S"); 
            Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) <  
                   ceps*Norm(c*c.Adjoint()),"C SV4 U"); 
            c2 = c;
            SV_Decompose(c2.View(),cS2.View(),cV2.View(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV5 S"); 
            Assert(Norm(cV2.Adjoint()*cS2*cS2*cV2-c.Adjoint()*c) < 
                   ceps*Norm(c.Adjoint()*c),"C SV5 V"); 

            cU2.Conjugate() = c;
            SV_Decompose(cU2.Conjugate(),cS2.View(),cV2.View(),true);
            Assert(Norm(c-cU2.Conjugate()*cS2*cV2) < ceps*Norm(m),"C SV6"); 
            Assert(Norm(cU2.Adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV6 - UtU"); 
            Assert(Norm(cV2.Adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV6 - VtV"); 
            Assert(Norm(cV2*cV2.Adjoint()-T(1)) < ceps*Norm(m),"C SV6 - VVt"); 
            c2.Conjugate() = c;
            SV_Decompose(c2.Conjugate(),cS2.View(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV7"); 
            cU2.Conjugate() = c;
            SV_Decompose(cU2.Conjugate(),cS2.View(),true);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV8 S"); 
            Assert(Norm(cU2.Conjugate()*cS2*cS2*cU2.Transpose()-c*c.Adjoint()) <
                   ceps*Norm(c*c.Adjoint()),"C SV8 U"); 
            c2.Conjugate() = c;
            SV_Decompose(c2.Conjugate(),cS2.View(),cV2.View(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV9 S"); 
            Assert(Norm(cV2.Adjoint()*cS2*cS2*cV2-c.Adjoint()*c) < 
                   ceps*Norm(c.Adjoint()*c),"C SV9 V"); 

            cU2 = c;
            SV_Decompose(cU2.View(),cS2.View(),cV2.Conjugate(),true);
            Assert(Norm(c-cU2*cS2*cV2.Conjugate()) < ceps*Norm(m),"C SV10"); 
            Assert(Norm(cU2.Adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV10 - UtU"); 
            Assert(Norm(cV2.Adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV10 - VtV"); 
            Assert(Norm(cV2*cV2.Adjoint()-T(1)) < ceps*Norm(m),"C SV10 - VVt"); 
            c2 = c;
            SV_Decompose(c2.View(),cS2.View(),cV2.Conjugate(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV11 S"); 
            Assert(Norm(cV2.Transpose()*cS2*cS2*cV2.Conjugate()-c.Adjoint()*c) < 
                   ceps*Norm(c.Adjoint()*c),"C SV9 V"); 

            cU2.Conjugate() = c;
            SV_Decompose(cU2.Conjugate(),cS2.View(),cV2.Conjugate(),true);
            Assert(Norm(c-cU2.Conjugate()*cS2*cV2.Conjugate()) < ceps*Norm(m),"C SV12"); 
            Assert(Norm(cU2.Adjoint()*cU2-T(1)) < ceps*Norm(m),"C SV12 - UtU"); 
            Assert(Norm(cV2.Adjoint()*cV2-T(1)) < ceps*Norm(m),"C SV12 - VtV"); 
            Assert(Norm(cV2*cV2.Adjoint()-T(1)) < ceps*Norm(m),"C SV12 - VVt"); 
            c2.Conjugate() = c;
            SV_Decompose(c2.Conjugate(),cS2.View(),cV2.Conjugate(),false);
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV13 S"); 
            Assert(Norm(cV2.Transpose()*cS2*cS2*cV2.Conjugate()-c.Adjoint()*c) < 
                   ceps*Norm(c.Adjoint()*c),"C SV13 V"); 
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef TEST_DOUBLE
template void TestMatrixDecomp<double,tmv::RowMajor>();
template void TestMatrixDecomp<double,tmv::ColMajor>();
#endif
#ifdef TEST_FLOAT
template void TestMatrixDecomp<float,tmv::RowMajor>();
template void TestMatrixDecomp<float,tmv::ColMajor>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrixDecomp<long double,tmv::RowMajor>();
template void TestMatrixDecomp<long double,tmv::ColMajor>();
#endif
#ifdef TEST_INT
template void TestMatrixDecomp<int,tmv::RowMajor>();
template void TestMatrixDecomp<int,tmv::ColMajor>();
#endif
