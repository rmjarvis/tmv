
#define START 0

#include "../src/TMV_Blas.h"
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"

template <class T, tmv::StorageType stor> 
void TestMatrixDecomp()
{
    for (int mattype = START; mattype <= 6; mattype++) {
#if !(XTEST & 64)
        if (mattype >= 4) break;
#endif

        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<", stor = "<<
                tmv::TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Square
        // mattype = 1  is NonSquare slightly tall
        // mattype = 2  is NonSquare very tall
        // mattype = 3  is Singular
        // mattype = 4  is Singular with seriously bad defects
        // mattype = 5  is Singular and nearly zero
        // mattype = 6  is Singular and nearly overflow

        const int N = 200;
        int M = N;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
        const bool singular = mattype >= 3;
        const bool baddefect = mattype == 4;
        const bool nearunderflow = mattype == 5;
        const bool nearoverflow = mattype == 6;

        tmv::Matrix<T,stor> m(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
        if (!singular) m /= T(10*N);
        m(3,3) = 14;
        m(3,4) = -2;
        m(0,2) = 7;
        m(M-1,N-4) = 23;
        m(M-1,N-2) = 13;
        m(M-1,N-1) = -10;
        if (singular) { m.col(1).setZero(); m.row(7) = m.row(6); }
        else m.diag() *= T(30*N);

        tmv::Matrix<std::complex<T>,stor> c(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        if (!singular) c /= T(10*N);
        c(3,3) = 14;
        c(3,4) = -2;
        c(0,2) = 7;
        c(M-1,N-4) = 23;
        c(M-1,N-2) = 13;
        c(M-1,N-1) = -10;
        if (singular) { c.col(1).setZero(); c.row(7) = c.row(6); }
        else c.diag() *= T(30*N);

        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.colRange(i,N) *= T(1.e-10);
                c.colRange(i,N) *= T(1.e-10);
                m.rowRange(i+5,N) *= T(-1.e-10);
                c.rowRange(i+5,N) *= T(-1.e-10);
            }
        }

        if (nearunderflow) {
            T x = std::numeric_limits<T>::min();
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
            m(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            c(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subMatrix(N-4,N,N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
            c.subMatrix(N-4,N,N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
        }

        if (nearoverflow) {
            T x = std::numeric_limits<T>::max();
            x /= N;
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,x/2));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
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

        if (showstartdone) {
            //std::cout<<"m = "<<m<<std::endl;
        }

        // LU Decomposition
        if (m.isSquare()) do {
            if (showstartdone) {
                std::cout<<"LU\n";
            }
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            tmv::Permutation P = m.lud().getP();
            tmv::Matrix<T> PLU = P*L*U;
            if (m.lud().isTrans()) PLU.transposeSelf();
            Assert(Norm(m-PLU) <= eps*normm,"LU"); 

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL =
                c.lud().getL();
            tmv::UpperTriMatrix<std::complex<T> > cU = c.lud().getU();
            tmv::Permutation cP = c.lud().getP();
            tmv::Matrix<std::complex<T> > cPLU = cP*cL*cU;
            if (c.lud().isTrans()) cPLU.transposeSelf();
            Assert(Norm(c-cPLU) <= ceps*normc,"C LU"); 

#if (XTEST & 16)
            tmv::Matrix<T,stor> m2 = m;
            // TODO: Need to make this syntax work better with Permutation.
            int p2[N];
            LU_Decompose(m2.view(),p2);
            tmv::Permutation P2(N,p2,true);
            L = m2.lowerTri(tmv::UnitDiag);
            U = m2.upperTri();
            PLU = P2*L*U;
            Assert(Norm(m-PLU) <= eps*normm,"LU2"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            LU_Decompose(c2.view(),p2);
            tmv::Permutation cP2(N,p2,true);
            cL = c2.lowerTri(tmv::UnitDiag);
            cU = c2.upperTri();
            cPLU = cP2*cL*cU;
            Assert(Norm(c-cPLU) <= ceps*normc,"C LU2"); 

            c2.conjugate() = c;
            LU_Decompose(c2.conjugate(),p2);
            tmv::Permutation cP3(N,p2,true);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cU = c2.conjugate().upperTri();
            cPLU = cP3*cL*cU;
            Assert(Norm(c-cPLU) <= ceps*normc,"C LU3"); 
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // QR Decomposition
        do {
            // QR Decomposition always works.  It just can't be used for
            // solving equations if the matrix is singular.
            if (showstartdone) {
                std::cout<<"QR\n";
            }
            tmv::Matrix<T,stor> Q = m.qrd().getQ();
            tmv::UpperTriMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            Assert(Norm(m-QR) <= eps*normm,"QR"); 
            Assert(Norm(m-m.qrd().getQ()*m.qrd().getR()) <= eps*normm,
                   "QR - PackedQ"); 
            Assert(Norm(Q.transpose()*Q-T(1)) <= T(N)*eps,"QR - QtQ"); 
            Assert(Norm(Q.transpose()*m.qrd().getQ()-T(1)) <= T(N)*eps,
                   "QR - QtQ - PackedQ (1)"); 
            Assert(Norm(Q / m.qrd().getQ()-T(1)) <= T(N)*eps,
                   "QR - QtQ - PackedQ (2)"); 

            tmv::Matrix<std::complex<T>,stor> cQ = c.qrd().getQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.qrd().getR();
            tmv::Matrix<std::complex<T> > cQR = cQ*cR;
            if (showacc) {
                std::cout<<"Norm(c-cQR) = "<<Norm(c-cQR)<<std::endl;
                std::cout<<"cf "<<ceps*normc<<std::endl;
                std::cout<<"Norm(c-cQR) (packed) = "<<
                    Norm(c-c.qrd().getQ()*c.qrd().getR())<<std::endl;
                std::cout<<"cf "<<ceps*normc<<std::endl;
                std::cout<<"Norm(cQt*cQ-1) = "<<
                    Norm(cQ.adjoint()*cQ-T(1))<<std::endl;
                std::cout<<"cf "<<T(N)*ceps<<std::endl;
            }
            Assert(Norm(c-cQR) <= ceps*normc,"C QR"); 
            Assert(Norm(c-c.qrd().getQ()*c.qrd().getR()) <= ceps*normc,
                   "C QR - PackedQ"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QR - QtQ"); 
            Assert(Norm(cQ.adjoint()*c.qrd().getQ()-T(1)) <= T(N)*ceps,
                   "C QR - QtQ - PackedQ (1)"); 
            Assert(Norm(cQ / c.qrd().getQ()-T(1)) <= T(N)*ceps,
                   "C QR - QtQ - PackedQ (2)"); 

#if (XTEST & 16)
            Q = m;
            QR_Decompose(Q.view(),R.view());
            QR = Q*R;
            Assert(Norm(m-QR) <= eps*normm,"QR2"); 
            Assert(Norm(Q.transpose()*Q-T(1)) <= T(N)*eps,"QR2 - QtQ"); 

            Q = m;
            QR_Decompose(Q.view());
            tmv::UpperTriMatrix<T> R2 = Q.upperTri();
            if (baddefect)
                // Sometimes R,R2 come out a bit different, but RtR = mtm,
                // so I think it's ok.  It's just that rounding errors can
                // make the R's differ by more than eps*normm.
                Assert(Norm(m.adjoint()*m-R2.adjoint()*R2) <= eps*normm*normm,
                       "QR3 (RtR)");
            else
                Assert(Norm(R-R2) <= eps*normm,"QR3"); 

            cQ = c;
            QR_Decompose(cQ.view(),cR.view());
            cQR = cQ*cR;
            Assert(Norm(c-cQR) <= ceps*normc,"C QR2"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QR2 - QtQ"); 

            cQ = c;
            QR_Decompose(cQ.view());
            tmv::UpperTriMatrix<std::complex<T> > cR2 = cQ.upperTri();
            if (baddefect)
                Assert(Norm(c.adjoint()*c-cR2.adjoint()*cR2) <= 
                       ceps*normc*normc,"QR3 (RtR)");
            else
                Assert(Norm(cR-cR2) <= ceps*normc,"C QR3"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR.view());
            cQR = cQ.conjugate()*cR;
            Assert(Norm(c-cQR) <= ceps*normc,"C QR4"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QR4 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate());
            cR2 = cQ.conjugate().upperTri();
            if (baddefect)
                Assert(Norm(c.adjoint()*c-cR2.adjoint()*cR2) <= 
                       ceps*normc*normc,"QR5 (RtR)");
            else
                Assert(Norm(cR-cR2) <= ceps*normc,"C QR5"); 

            cQ = c;
            QR_Decompose(cQ.view(),cR.conjugate());
            cQR = cQ*cR.conjugate();
            Assert(Norm(c-cQR) <= ceps*normc,"C QR6"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QR6 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR.conjugate());
            cQR = cQ.conjugate()*cR.conjugate();
            Assert(Norm(c-cQR) <= ceps*normc,"C QR7"); 
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

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
            tmv::Permutation P = m.qrpd().getP();
            tmv::Matrix<T> QRP = Q*R*P;
            Assert(Norm(m-QRP) <= eps*normm,"QRP"); 
            QRP = m.qrpd().getQ()*R*P;
            Assert(Norm(m-QRP) <= eps*normm,"QRP - Packed"); 
            Assert(Norm(Q.transpose()*Q-T(1)) <= T(N)*eps,"QRP - QtQ"); 
            Assert(Norm(Q.transpose()*m.qrpd().getQ()-T(1)) <= T(N)*eps,
                   "QR - QtQ - PackedQ (1)"); 
            Assert(Norm(Q / m.qrpd().getQ()-T(1)) <= T(N)*eps,
                   "QR - QtQ - PackedQ (2)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP - strict"); 
            }
#endif

            tmv::Matrix<std::complex<T>,stor> cQ = c.qrpd().getQ();
            tmv::UpperTriMatrix<std::complex<T> > cR = c.qrpd().getR();
            tmv::Permutation cP = c.qrpd().getP();
            tmv::Matrix<std::complex<T> > cQRP = cQ*cR*cP;
            Assert(Norm(c-cQRP) <= ceps*normc,"C QRP"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QRP - QtQ"); 
            Assert(Norm(cQ.adjoint()*c.qrpd().getQ()-T(1)) <= T(N)*ceps,
                   "C QRP - QtQ - PackedQ (1)"); 
            Assert(Norm(cQ / c.qrpd().getQ()-T(1)) <= T(N)*ceps,
                   "C QRP - QtQ - PackedQ (2)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP - strict"); 
            }
#endif

#if (XTEST & 16)
            Q = m;
            int p2[N];
            QRP_Decompose(Q.view(),R.view(),p2,strict);
            tmv::Permutation P2(N,p2);
            QRP = Q*R*P2;
            Assert(Norm(m-QRP) <= eps*normm,"QRP2"); 
            Assert(Norm(Q.transpose()*Q-T(1)) <= T(N)*eps,"QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP2 - strict"); 
            }
#endif

            Q = m;
            QRP_Decompose(Q.view(),strict);
            tmv::UpperTriMatrix<T> R2 = Q.upperTri();
            Assert(Norm(R-R2) <= eps*normm,"QRP3"); 

            cQ = c;
            QRP_Decompose(cQ.view(),cR.view(),p2,strict);
            tmv::Permutation cP2(N,p2);
            cQRP = cQ*cR*cP2;
            Assert(Norm(c-cQRP) <= ceps*normc,"C QRP2"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP2 - strict"); 
            }
#endif
            cQ = c;
            QRP_Decompose(cQ.view(),strict);
            tmv::UpperTriMatrix<std::complex<T> > cR2 = cQ.upperTri();
            Assert(Norm(cR-cR2) <= ceps*normc,"C QRP3"); 

            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),cR.view(),p2,strict);
            tmv::Permutation cP3(N,p2);
            cQRP = cQ.conjugate()*cR*cP3;
            Assert(Norm(c-cQRP) <= ceps*normc,"C QRP4"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QRP4 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP4 - strict"); 
            }
#endif

            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),strict);
            cR2 = cQ.conjugate().upperTri();
            Assert(Norm(cR-cR2) <= ceps*normc,"C QRP5"); 

            cQ = c;
            QRP_Decompose(cQ.view(),cR.conjugate(),p2,strict);
            tmv::Permutation cP4(N,p2);
            cQRP = cQ*cR.conjugate()*cP4;
            Assert(Norm(c-cQRP) <= ceps*normc,"C QRP6"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QRP6 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP6 - strict"); 
            }
#endif


            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),cR.conjugate(),p2,strict);
            tmv::Permutation cP5(N,p2);
            cQRP = cQ.conjugate()*cR.conjugate()*cP5;
            Assert(Norm(c-cQRP) <= ceps*normc,"C QRP7"); 
            Assert(Norm(cQ.adjoint()*cQ-T(1)) <= T(N)*ceps,"C QRP7 - QtQ"); 
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
        do {
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
            Assert(Norm(U.transpose()*U-T(1)) <= eps,"SV - UtU"); 
            Assert(Norm(V.transpose()*V-T(1)) <= eps,"SV - VtV"); 
            Assert(Norm(V*V.transpose()-T(1)) <= eps,"SV - VVt"); 

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            if (showacc) {
                std::cout<<"Norm(c-USV) = "<<Norm(c-cU*cS*cV)<<
                    "  cf "<<ceps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<
                    "  cf "<<ceps<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(cV.adjoint()*cV-T(1))<<
                    "  cf "<<ceps<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(cV*cV.adjoint()-T(1))<<
                    "  cf "<<ceps<<std::endl;
            }
            Assert(Norm(c-cU*cS*cV) <= ceps*normc,"C SV"); 
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C SV - UtU"); 
            Assert(Norm(cV.adjoint()*cV-T(1)) <= ceps,"C SV - VtV"); 
            Assert(Norm(cV*cV.adjoint()-T(1)) <= ceps,"C SV - VVt"); 

#if (XTEST & 16)
            tmv::Matrix<T,stor> U2 = m;
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(U2.view(),S2.view(),V2.view(),true);
            Assert(Norm(m-U2*S2*V2) <= eps*normm,"SV2"); 
            Assert(Norm(U2.transpose()*U2-T(1)) <= eps,"SV2 - UtU"); 
            Assert(Norm(V2.transpose()*V2-T(1)) <= eps,"SV2 - VtV"); 
            Assert(Norm(V2*V2.transpose()-T(1)) <= eps,"SV2 - VVt"); 

            tmv::Matrix<T,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view(),false);
            Assert(Norm(S2-S) <= eps*normm,"SV3"); 
            U2 = m;
            SV_Decompose(U2.view(),S2.view(),true);
            Assert(Norm(S2-S) <= eps*normm,"SV4 S"); 
            if (!nearoverflow)
                Assert(Norm(m*m.transpose()-U2*S2*S2*U2.transpose()) <=  
                       eps*normm*normm,"SV4 U"); 
            m2 = m;
            SV_Decompose(m2.view(),S2.view(),V2.view(),false);
            Assert(Norm(S2-S) <= eps*normm,"SV5 S"); 
            if (!nearoverflow)
                Assert(Norm(m.transpose()*m-V2.transpose()*S2*S2*V2) <= 
                       eps*normm*normm,"SV5 V"); 

            tmv::Matrix<std::complex<T>,stor> cU2 = c;
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(cU2.view(),cS2.view(),cV2.view(),true);
            Assert(Norm(c-cU2*cS2*cV2) <= ceps*normm,"C SV2"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) <= ceps,"C SV2 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) <= ceps,"C SV2 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) <= ceps,"C SV2 - VVt"); 

            tmv::Matrix<std::complex<T>,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV3"); 
            cU2 = c;
            SV_Decompose(cU2.view(),cS2.view(),true);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV4 S"); 
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-cU2*cS2*cS2*cU2.adjoint()) <=  
                       ceps*normc*normc,"C SV4 U"); 
            c2 = c;
            SV_Decompose(c2.view(),cS2.view(),cV2.view(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV5 S"); 
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) <= 
                       ceps*normc*normc,"C SV5 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),cV2.view(),true);
            Assert(Norm(c-cU2.conjugate()*cS2*cV2) <= ceps*normm,"C SV6"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) <= ceps,"C SV6 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) <= ceps,"C SV6 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) <= ceps,"C SV6 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV7"); 
            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),true);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV8 S"); 
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-
                            cU2.conjugate()*cS2*cS2*cU2.transpose()) <=
                       ceps*normc*normc,"C SV8 U"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),cV2.view(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV9 S"); 
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) <= 
                       ceps*normc*normc,"C SV9 V"); 

            cU2 = c;
            SV_Decompose(cU2.view(),cS2.view(),cV2.conjugate(),true);
            Assert(Norm(c-cU2*cS2*cV2.conjugate()) <= ceps*normm,"C SV10"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) <= ceps,"C SV10 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) <= ceps,"C SV10 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) <= ceps,"C SV10 - VVt"); 
            c2 = c;
            SV_Decompose(c2.view(),cS2.view(),cV2.conjugate(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV11 S"); 
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-
                            cV2.transpose()*cS2*cS2*cV2.conjugate()) <=
                       ceps*normc*normc,"C SV9 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2.view(),cV2.conjugate(),true);
            Assert(Norm(c-cU2.conjugate()*cS2*cV2.conjugate()) <= ceps*normm,
                   "C SV12"); 
            Assert(Norm(cU2.adjoint()*cU2-T(1)) <= ceps,"C SV12 - UtU"); 
            Assert(Norm(cV2.adjoint()*cV2-T(1)) <= ceps,"C SV12 - VtV"); 
            Assert(Norm(cV2*cV2.adjoint()-T(1)) <= ceps,"C SV12 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2.view(),cV2.conjugate(),false);
            Assert(Norm(cS2-cS) <= ceps*normc,"C SV13 S"); 
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-
                            cV2.transpose()*cS2*cS2*cV2.conjugate()) <=
                       ceps*normc*normc,"C SV13 V"); 
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);
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
