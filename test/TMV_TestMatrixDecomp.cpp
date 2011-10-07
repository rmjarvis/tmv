
#define START 0

#undef NDEBUG
#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test_1.h"

template <class T, tmv::StorageType stor> 
void TestMatrixDecomp()
{
    typedef typename tmv::Traits<T>::float_type FT;
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 6; mattype++) {
#if !(XTEST & 64) || defined(LAP)
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

#if 1
        const int N = 200;
        int M = N;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
#elif 0
        const int N = 8;
        int M = N;
        if (mattype == 1) M = 10;
        else if (mattype == 2) M = 17;
#else
        const int N = 150;
        int M = N;
        if (mattype == 1) N+5;
        else if (mattype == 2) 5*N/2;
#endif

        const bool singular = mattype >= 3;
        const bool baddefect = mattype == 4;
        const bool nearunderflow = mattype == 5;
        const bool nearoverflow = mattype == 6;
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();

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

        tmv::Matrix<CT,stor> c(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            c(i,j) = CT(2+4*i-5*j,3-i);
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
            T x = Tmin;
            m.setAllTo(x);
            c.setAllTo(CT(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
            c.subMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.setAllTo(CT(x,x/2));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2,N) *= T(-1);
            c.diag(0,N/2,N) *= T(-1);
        }

        FT eps = EPS;
        FT ceps = EPS;
        FT normm = Norm(m);
        FT normc = Norm(c);
        if (showacc) {
            std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
            std::cout<<"norm = "<<normm<<"  "<<normc<<std::endl;
        }
        if (!singular) {
            FT kappa = normm * Norm(m.inverse());
            FT ckappa = normc * Norm(c.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<"  "<<ckappa<<std::endl;
            eps *= kappa;
            ceps *= ckappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= FT(10*N);
            ceps *= FT(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        if (showstartdone) {
            //std::cout<<"m = "<<m<<std::endl;
        }

        std::ostream* dbgout = showdiv ? &std::cout : 0;

        // LU Decomposition
        if (m.isSquare()) do {
            if (showstartdone) {
                std::cout<<"LU"<<std::endl;
            }
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            tmv::Permutation P = m.lud().getP();
            tmv::Matrix<T> PLU = P*L*U;
            if (m.lud().isTrans()) PLU.transposeSelf();
            if (showacc) {
                std::cout<<"Norm(m-PLU) = "<<Norm(m-PLU)<<std::endl;
                std::cout<<"cf "<<eps*normm<<std::endl;
            }
            Assert(CheckDecomp(m.lud(),m,dbgout),"LU CheckDecomp");
            Assert(Equal(m,PLU,eps*normm),"LU"); 

            tmv::LowerTriMatrix<CT,tmv::UnitDiag> cL =
                c.lud().getL();
            tmv::UpperTriMatrix<CT> cU = c.lud().getU();
            tmv::Permutation cP = c.lud().getP();
            tmv::Matrix<CT> cPLU = cP*cL*cU;
            if (c.lud().isTrans()) cPLU.transposeSelf();
            if (showacc) {
                std::cout<<"Norm(c-cPLU) = "<<Norm(c-cPLU)<<std::endl;
                std::cout<<"cf "<<ceps*normc<<std::endl;
            }
            Assert(CheckDecomp(c.lud(),c,dbgout),"C LU CheckDecomp");
            Assert(Equal(c,cPLU,ceps*normc),"C LU"); 

#if (XTEST & 16)
            tmv::Matrix<T,stor> m2 = m;
            tmv::Permutation P2(N);
            LU_Decompose(m2,P2);
            L = m2.lowerTri(tmv::UnitDiag);
            U = m2.upperTri();
            PLU = P2*L*U;
            Assert(Equal(m,PLU,eps*normm),"LU2"); 

            tmv::Matrix<CT,stor> c2 = c;
            tmv::Permutation cP2(N);
            LU_Decompose(c2,cP2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cU = c2.upperTri();
            cPLU = cP2*cL*cU;
            Assert(Equal(c,cPLU,ceps*normc),"C LU2"); 

            c2.conjugate() = c;
            tmv::Permutation cP3(N);
            LU_Decompose(c2.conjugate(),cP3);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cU = c2.conjugate().upperTri();
            cPLU = cP3*cL*cU;
            Assert(Equal(c,cPLU,ceps*normc),"C LU3"); 
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // QR Decomposition
        do {
            // QR Decomposition always works.  It just can't be used for
            // solving equations if the matrix is singular.
            if (showstartdone) {
                std::cout<<"QR"<<std::endl;
            }
            tmv::Matrix<T,stor> Q = m.qrd().getQ();
            tmv::UpperTriMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            if (showacc) {
                //std::cout<<"m = "<<m<<std::endl;
                //std::cout<<"QR = "<<QR<<std::endl;
                //std::cout<<"m-QR = "<<(m-QR)<<std::endl;
                std::cout<<"Norm(m-QR) = "<<Norm(m-QR)<<std::endl;
            }
            Assert(CheckDecomp(m.qrd(),m,dbgout),"QR CheckDecomp");
            Assert(Equal(m,QR,eps*normm),"QR"); 
            Assert(Equal(m,m.qrd().getQ()*m.qrd().getR(),eps*normm),
                   "QR - PackedQ"); 
            Assert(Equal(Q.transpose()*Q,T(1),T(N)*eps),"QR - QtQ"); 
            Assert(Equal(Q.transpose()*m.qrd().getQ(),T(1),T(N)*eps),
                   "QR - QtQ - PackedQ (1)"); 
            Assert(Equal(Q / m.qrd().getQ(),T(1),T(N)*eps),
                   "QR - QtQ - PackedQ (2)"); 

            tmv::Matrix<CT,stor> cQ = c.qrd().getQ();
            tmv::UpperTriMatrix<CT> cR = c.qrd().getR();
            tmv::Matrix<CT> cQR = cQ*cR;
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
            Assert(CheckDecomp(c.qrd(),c,dbgout),"C QR CheckDecomp");
            Assert(Equal(c,cQR,ceps*normc),"C QR"); 
            Assert(Equal(c,c.qrd().getQ()*c.qrd().getR(),ceps*normc),
                   "C QR - PackedQ"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QR - QtQ"); 
            Assert(Equal(cQ.adjoint()*c.qrd().getQ(),T(1),T(N)*ceps),
                   "C QR - QtQ - PackedQ (1)"); 
            Assert(Equal(cQ / c.qrd().getQ(),T(1),T(N)*ceps),
                   "C QR - QtQ - PackedQ (2)"); 

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            Q = m;
            QR_Decompose(Q,R);
            QR = Q*R;
            Assert(Equal(m,QR,eps*normm),"QR2"); 
            Assert(Equal(Q.transpose()*Q,T(1),T(N)*eps),"QR2 - QtQ"); 

            Q = m;
            QR_Decompose(Q);
            tmv::UpperTriMatrix<T> R2 = Q.upperTri();
            if (showacc) {
                std::cout<<"x = "<<x<<std::endl;
                std::cout<<"mtm = "<<
                    (m.adjoint()/x).calc()*(m/x).calc()<<std::endl;
                std::cout<<"RtR = "<<
                    (R2.adjoint()/x).calc()*(R2/x).calc()<<std::endl;
                std::cout<<"Norm(diff) = "<<
                    Norm((m.adjoint()/x).calc()*(m/x).calc()-
                         (R2.adjoint()/x).calc()*(R2/x).calc())<<
                    std::endl;
                std::cout<<"cf. "<< eps*(normm/x)*(normm/x)<<std::endl;
            }
            Assert(Equal((m.adjoint()/x).calc()*(m/x).calc(),
                         (R2.adjoint()/x).calc()*(R2/x).calc(),
                         eps*(normm/x)*(normm/x)),"QR3 (RtR)");

            cQ = c;
            QR_Decompose(cQ,cR);
            cQR = cQ*cR;
            Assert(Equal(c,cQR,ceps*normc),"C QR2"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QR2 - QtQ"); 

            cQ = c;
            QR_Decompose(cQ);
            tmv::UpperTriMatrix<CT> cR2 = cQ.upperTri();
            Assert(Equal((c.adjoint()/x).calc()*(c/x).calc(),
                         (cR2.adjoint()/x).calc()*(cR2/x).calc(),
                         ceps*(normc/x)*(normc/x)),"C QR3 (RtR)");

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR);
            cQR = cQ.conjugate()*cR;
            Assert(Equal(c,cQR,ceps*normc),"C QR4"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QR4 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate());
            cR2 = cQ.conjugate().upperTri();
            Assert(Equal((c.adjoint()/x).calc()*(c/x).calc(),
                         (cR2.adjoint()/x).calc()*(cR2/x).calc(),
                         ceps*(normc/x)*(normc/x)),"C QR5 (RtR)");

            cQ = c;
            QR_Decompose(cQ,cR.conjugate());
            cQR = cQ*cR.conjugate();
            Assert(Equal(c,cQR,ceps*normc),"C QR6"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QR6 - QtQ"); 

            cQ.conjugate() = c;
            QR_Decompose(cQ.conjugate(),cR.conjugate());
            cQR = cQ.conjugate()*cR.conjugate();
            Assert(Equal(c,cQR,ceps*normc),"C QR7"); 
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // TODO: Add QR_Update, QR_Downdate here to test their 
        // algorithms for large matrices.

        // QRP Decomposition
        for (int istrict = 0; istrict <= 1; istrict++) {
            if (showstartdone) {
                std::cout<<"QRP"<<std::endl;
            }
            bool strict = istrict == 1;
            tmv::UseStrictQRP(strict);
            m.unsetDiv(); c.unsetDiv(); 

            tmv::Matrix<T,stor> Q = m.qrpd().getQ();
            tmv::UpperTriMatrix<T> R = m.qrpd().getR();
            tmv::Permutation P = m.qrpd().getP();
            tmv::Matrix<T> QRP = Q*R*P;

            if (showacc) {
                std::cout<<"Norm(m-QRP) = "<<Norm(m-QRP)<<
                    "  cf "<<eps*normm<<std::endl;
                std::cout<<"Norm(QtQ-1) = "<<Norm(Q.transpose()*Q-T(1))<<
                    "  cf "<<T(N)*eps<<std::endl;
            }
            Assert(CheckDecomp(m.qrpd(),m,dbgout),"QRP CheckDecomp");
            Assert(Equal(m,QRP,eps*normm),"QRP"); 
            QRP = m.qrpd().getQ()*R*P;
            Assert(Equal(m,QRP,eps*normm),"QRP - Packed"); 
            Assert(Equal(Q.transpose()*Q,T(1),T(N)*eps),"QRP - QtQ"); 
            Assert(Equal(Q.transpose()*m.qrpd().getQ(),T(1),T(N)*eps),
                   "QR - QtQ - PackedQ (1)"); 
            Assert(Equal(Q / m.qrpd().getQ(),T(1),T(N)*eps),
                   "QR - QtQ - PackedQ (2)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP - strict"); 
            }
#endif
            tmv::Matrix<CT,stor> cQ = c.qrpd().getQ();
            tmv::UpperTriMatrix<CT> cR = c.qrpd().getR();
            tmv::Permutation cP = c.qrpd().getP();
            tmv::Matrix<CT> cQRP = cQ*cR*cP;

            if (showacc) {
                std::cout<<"Norm(c-cQRP) = "<<Norm(c-cQRP)<<
                    "  cf "<<ceps*normc<<std::endl;
                std::cout<<"Norm(QtQ-1) = "<<Norm(cQ.adjoint()*cQ-T(1))<<
                    "  cf "<<T(N)*ceps<<std::endl;
            }
            Assert(CheckDecomp(c.qrpd(),c,dbgout),"C QRP CheckDecomp");
            Assert(Equal(c,cQRP,ceps*normc),"C QRP"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QRP - QtQ"); 
            Assert(Equal(cQ.adjoint()*c.qrpd().getQ(),T(1),T(N)*ceps),
                   "C QRP - QtQ - PackedQ (1)"); 
            Assert(Equal(cQ / c.qrpd().getQ(),T(1),T(N)*ceps),
                   "C QRP - QtQ - PackedQ (2)"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP - strict"); 
            }
#endif

#if (XTEST & 16)
            Q = m;
            tmv::Permutation P2(N);
            QRP_Decompose(Q,R,P2,strict);
            QRP = Q*R*P2;
            Assert(Equal(m,QRP,eps*normm),"QRP2"); 
            Assert(Equal(Q.transpose()*Q,T(1),T(N)*eps),"QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<R.size();i++) {
                Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP2 - strict"); 
            }
#endif

            Q = m;
            QRP_Decompose(Q,strict);
            tmv::UpperTriMatrix<T> R2 = Q.upperTri();
            Assert(Equal(R,R2,eps*normm),"QRP3"); 

            cQ = c;
            tmv::Permutation cP2(N);
            QRP_Decompose(cQ,cR,cP2,strict);
            cQRP = cQ*cR*cP2;
            Assert(Equal(c,cQRP,ceps*normc),"C QRP2"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QRP2 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP2 - strict"); 
            }
#endif
            cQ = c;
            QRP_Decompose(cQ,strict);
            tmv::UpperTriMatrix<CT> cR2 = cQ.upperTri();
            Assert(Equal(cR,cR2,ceps*normc),"C QRP3"); 

            cQ.conjugate() = c;
            tmv::Permutation cP3(N);
            QRP_Decompose(cQ.conjugate(),cR,cP3,strict);
            cQRP = cQ.conjugate()*cR*cP3;
            Assert(Equal(c,cQRP,ceps*normc),"C QRP4"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QRP4 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP4 - strict"); 
            }
#endif

            cQ.conjugate() = c;
            QRP_Decompose(cQ.conjugate(),strict);
            cR2 = cQ.conjugate().upperTri();
            Assert(Equal(cR,cR2,ceps*normc),"C QRP5"); 

            cQ = c;
            tmv::Permutation cP4(N);
            QRP_Decompose(cQ,cR.conjugate(),cP4,strict);
            cQRP = cQ*cR.conjugate()*cP4;
            Assert(Equal(c,cQRP,ceps*normc),"C QRP6"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QRP6 - QtQ"); 
#ifndef LAP
            if (strict) for(size_t i=1;i<cR.size();i++) {
                Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),
                       "C QRP6 - strict"); 
            }
#endif


            cQ.conjugate() = c;
            tmv::Permutation cP5(N);
            QRP_Decompose(cQ.conjugate(),cR.conjugate(),cP5,strict);
            cQRP = cQ.conjugate()*cR.conjugate()*cP5;
            Assert(Equal(c,cQRP,ceps*normc),"C QRP7"); 
            Assert(Equal(cQ.adjoint()*cQ,T(1),T(N)*ceps),"C QRP7 - QtQ"); 
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
                std::cout<<"SV"<<std::endl;
            }
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            if (showacc) {
                //std::cout<<"U = "<<U<<std::endl;
                //std::cout<<"V = "<<V<<std::endl;
                //std::cout<<"S = "<<S<<std::endl;
                //std::cout<<"m-USV = "<<m-U*S*V<<std::endl;
                std::cout<<"Norm(m-USV) = "<<Norm(m-U*S*V)<<
                    "  cf "<<eps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.transpose()*U-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(V.transpose()*V-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(V*V.transpose()-T(1))<<
                    "  cf "<<eps<<std::endl;
            }
            Assert(CheckDecomp(m.svd(),m,dbgout),"SV CheckDecomp");
            Assert(Equal(m,U*S*V,eps*normm),"SV"); 
            Assert(Equal(U.transpose()*U,T(1),eps),"SV - UtU"); 
            Assert(Equal(V.transpose()*V,T(1),eps),"SV - VtV"); 
            Assert(Equal(V*V.transpose(),T(1),eps),"SV - VVt"); 

            tmv::Matrix<CT> cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<CT> cV = c.svd().getV();
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
            Assert(CheckDecomp(c.svd(),c,dbgout),"C SV CheckDecomp");
            Assert(Equal(c,cU*cS*cV,ceps*normc),"C SV"); 
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C SV - UtU"); 
            Assert(Equal(cV.adjoint()*cV,T(1),ceps),"C SV - VtV"); 
            Assert(Equal(cV*cV.adjoint(),T(1),ceps),"C SV - VVt"); 

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            // Direct SV_Decompose calls:
            tmv::Matrix<T,stor> U2 = m;
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T,tmv::ColMajor> V2(N,N);
            SV_Decompose(U2,S2,V2,true);
            if (showacc) {
                std::cout<<"Norm(m-USV) = "<<Norm(m-U2*S2*V2)<<std::endl;
                std::cout<<"cf "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U2*S2*V2,eps*normm),"SV2"); 
            Assert(Equal(U2.transpose()*U2,T(1),eps),"SV2 - UtU"); 
            Assert(Equal(V2.transpose()*V2,T(1),eps),"SV2 - VtV"); 
            Assert(Equal(V2*V2.transpose(),T(1),eps),"SV2 - VVt"); 

            // No U or V
            tmv::Matrix<T,stor> m2 = m;
            SV_Decompose(m2,S2,false);
            Assert(Equal(S2,S,eps*normm),"SV3"); 
            U2 = m;
            SV_Decompose(U2,S2,true);
            Assert(Equal(S2,S,eps*normm),"SV4 S"); 
            Assert(Equal(
                    (m/x).calc() * (m.transpose()/x).calc(),
                    U2 * (S2/x).calc() * (S2/x).calc() * U2.transpose(),
                    eps*(normm/x)*(normm/x)),"SV4 U"); 
            m2 = m;
            SV_Decompose(m2,S2,V2,false);
            Assert(Equal(S2,S,eps*normm),"SV5 S"); 
            Assert(Equal(
                    (m.transpose()/x).calc() * (m/x).calc(),
                    V2.transpose() * (S2/x).calc() * (S2/x).calc() * V2,
                    eps*(normm/x)*(normm/x)),"SV5 V"); 

            // Complex:
            tmv::Matrix<CT,stor> cU2 = c;
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT,tmv::ColMajor> cV2(N,N);
            SV_Decompose(cU2,cS2,cV2,true);
            Assert(Equal(c,cU2*cS2*cV2,ceps*normm),"C SV2"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV2 - UtU"); 
            Assert(Equal(cV2.adjoint()*cV2,T(1),ceps),"C SV2 - VtV"); 
            Assert(Equal(cV2*cV2.adjoint(),T(1),ceps),"C SV2 - VVt"); 

            // No U or V
            tmv::Matrix<CT,stor> c2 = c;
            SV_Decompose(c2,cS2,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV3"); 
            cU2 = c;
            SV_Decompose(cU2,cS2,true);
            Assert(Equal(cS2,cS,ceps*normc),"C SV4 S"); 
            Assert(Equal(
                    (c/x).calc() * (c.adjoint()/x).calc(),
                    cU2 * (cS2/x).calc() * (cS2/x).calc() * cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"C SV4 U"); 
            c2 = c;
            SV_Decompose(c2,cS2,cV2,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV5 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2.adjoint() * (cS2/x).calc() * (cS2/x).calc() * cV2,
                    ceps*(normc/x)*(normc/x)),"C SV5 V"); 

            // U conjugate
            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2,cV2,true);
            Assert(Equal(c,cU2.conjugate()*cS2*cV2,ceps*normm),"C SV6"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV6 - UtU"); 
            Assert(Equal(cV2.adjoint()*cV2,T(1),ceps),"C SV6 - VtV"); 
            Assert(Equal(cV2*cV2.adjoint(),T(1),ceps),"C SV6 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV7"); 
            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2,true);
            Assert(Equal(cS2,cS,ceps*normc),"C SV8 S"); 
            Assert(Equal(
                    (c/x).calc() * (c.adjoint()/x).calc(),
                    cU2.conjugate() * (cS2/x).calc() *
                    (cS2/x).calc() * cU2.transpose(),
                    ceps*(normc/x)*(normc/x)),"C SV8 U"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2,cV2,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV9 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2.adjoint() * (cS2/x).calc() * (cS2/x).calc() * cV2,
                    ceps*(normc/x)*(normc/x)),"C SV9 V"); 

            // V conjugate
            cU2 = c;
            SV_Decompose(cU2,cS2,cV2.conjugate(),true);
            Assert(Equal(c,cU2*cS2*cV2.conjugate(),ceps*normm),"C SV10"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV10 - UtU"); 
            Assert(Equal(cV2.adjoint()*cV2,T(1),ceps),"C SV10 - VtV"); 
            Assert(Equal(cV2*cV2.adjoint(),T(1),ceps),"C SV10 - VVt"); 
            c2 = c;
            SV_Decompose(c2,cS2,cV2.conjugate(),false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV11 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2.transpose() * (cS2/x).calc() *
                    (cS2/x).calc() * cV2.conjugate(),
                    ceps*(normc/x)*(normc/x)),"C SV11 V"); 

            // U,V conjugate
            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2,cV2.conjugate(),true);
            Assert(Equal(c,cU2.conjugate()*cS2*cV2.conjugate(),ceps*normm),
                   "C SV12"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV12 - UtU"); 
            Assert(Equal(cV2.adjoint()*cV2,T(1),ceps),"C SV12 - VtV"); 
            Assert(Equal(cV2*cV2.adjoint(),T(1),ceps),"C SV12 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2,cV2.conjugate(),false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV13 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2.transpose() * (cS2/x).calc() *
                    (cS2/x).calc() * cV2.conjugate(),
                    ceps*(normc/x)*(normc/x)),"C SV13 V"); 

#if XTEST & 2
            // RowMajor V
            U2 = m;
            tmv::Matrix<T,tmv::RowMajor> V2r(N,N);
            SV_Decompose(U2,S2,V2r,true);
            if (showacc) {
                std::cout<<"Norm(m-USV) = "<<Norm(m-U2*S2*V2r)<<std::endl;
                std::cout<<"cf "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U2*S2*V2r,eps*normm),"SV6"); 
            Assert(Equal(U2.transpose()*U2,T(1),eps),"SV6 - UtU"); 
            Assert(Equal(V2r.transpose()*V2r,T(1),eps),"SV6 - VtV"); 
            Assert(Equal(V2r*V2r.transpose(),T(1),eps),"SV6 - VVt"); 

            m2 = m;
            SV_Decompose(m2,S2,V2r,false);
            Assert(Equal(S2,S,eps*normm),"SV7 S"); 
            Assert(Equal(
                    (m.transpose()/x).calc() * (m/x).calc(),
                    V2r.transpose() * (S2/x).calc() * (S2/x).calc() * V2r,
                    eps*(normm/x)*(normm/x)),"SV7 V"); 

            cU2 = c;
            tmv::Matrix<CT,tmv::RowMajor> cV2r(N,N);
            SV_Decompose(cU2,cS2,cV2r,true);
            Assert(Equal(c,cU2*cS2*cV2r,ceps*normm),"C SV14"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV14 - UtU"); 
            Assert(Equal(cV2r.adjoint()*cV2r,T(1),ceps),"C SV14 - VtV"); 
            Assert(Equal(cV2r*cV2r.adjoint(),T(1),ceps),"C SV14 - VVt"); 

            c2 = c;
            SV_Decompose(c2,cS2,cV2r,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV15 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2r.adjoint() * (cS2/x).calc() * (cS2/x).calc() * cV2r,
                    ceps*(normc/x)*(normc/x)),"C SV15 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2,cV2r,true);
            Assert(Equal(c,cU2.conjugate()*cS2*cV2r,ceps*normm),"C SV16"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV16 - UtU"); 
            Assert(Equal(cV2r.adjoint()*cV2r,T(1),ceps),"C SV16 - VtV"); 
            Assert(Equal(cV2r*cV2r.adjoint(),T(1),ceps),"C SV16 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2,cV2r,false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV17 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2r.adjoint() * (cS2/x).calc() * (cS2/x).calc() * cV2r,
                    ceps*(normc/x)*(normc/x)),"C SV17 V"); 

            cU2 = c;
            SV_Decompose(cU2,cS2,cV2r.conjugate(),true);
            Assert(Equal(c,cU2*cS2*cV2r.conjugate(),ceps*normm),"C SV18"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV18 - UtU"); 
            Assert(Equal(cV2r.adjoint()*cV2r,T(1),ceps),"C SV18 - VtV"); 
            Assert(Equal(cV2r*cV2r.adjoint(),T(1),ceps),"C SV18 - VVt"); 
            c2 = c;
            SV_Decompose(c2,cS2,cV2r.conjugate(),false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV19 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2r.transpose() * (cS2/x).calc() *
                    (cS2/x).calc() * cV2r.conjugate(),
                    ceps*(normc/x)*(normc/x)),"C SV19 V"); 

            cU2.conjugate() = c;
            SV_Decompose(cU2.conjugate(),cS2,cV2r.conjugate(),true);
            Assert(Equal(c,cU2.conjugate()*cS2*cV2r.conjugate(),ceps*normm),
                   "C SV20"); 
            Assert(Equal(cU2.adjoint()*cU2,T(1),ceps),"C SV20 - UtU"); 
            Assert(Equal(cV2r.adjoint()*cV2r,T(1),ceps),"C SV20 - VtV"); 
            Assert(Equal(cV2r*cV2r.adjoint(),T(1),ceps),"C SV20 - VVt"); 
            c2.conjugate() = c;
            SV_Decompose(c2.conjugate(),cS2,cV2r.conjugate(),false);
            Assert(Equal(cS2,cS,ceps*normc),"C SV21 S"); 
            Assert(Equal(
                    (c.adjoint()/x).calc() * (c/x).calc(),
                    cV2r.transpose() * (cS2/x).calc() *
                    (cS2/x).calc() * cV2r.conjugate(),
                    ceps*(normc/x)*(normc/x)),"C SV21 V"); 
#endif
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
