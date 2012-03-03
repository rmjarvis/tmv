
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

template <class T, tmv::StorageType stor> 
void TestBandDecomp()
{
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 9; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype >= 7) break;
#endif

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
        const bool baddefect = mattype == 7;
        const bool nearunderflow = mattype == 8;
        const bool nearoverflow = mattype == 9;
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();

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

        tmv::BandMatrix<CT,stor> c(M,N,nlo,nhi);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            if (i<=j+nlo && j<=i+nhi) c(i,j) = CT(2+4*i-5*j,3-i);
        if (!singular) {
            c(0,0) = T(14);
            if (c.nlo() >= 1) c(1,0) = T(-2);
            if (c.nhi() >= 1) c(4,5) = T(7);
            c(2,2) = T(-10);
            c.diag() *= T(30);
        }

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
                int start = i-nhi;  if (start < 0) start = 0;
                int end = i+nlo;  if (end > N) end = N;
                m.col(i,start,end) /= T(i);
                c.col(i,start,end) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subBandMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
            c.subBandMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.setAllTo(CT(x,x/2));
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
        if (m.isSquare()) do {
            if (showstartdone) {
                std::cout<<"LU"<<std::endl;
            }
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::UpperTriMatrix<T> U = m.lud().getU();
            tmv::Permutation P = m.lud().getP();
            tmv::Matrix<T> PLU = P*L*U;
            if (m.lud().isTrans()) PLU.transposeSelf();
            Assert(Equal(m,PLU,eps*normm),"Band LU");

            tmv::LowerTriMatrix<CT,tmv::UnitDiag> cL = 
                c.lud().getL();
            tmv::UpperTriMatrix<CT> cU = c.lud().getU();
            tmv::Permutation cP = c.lud().getP();
            tmv::Matrix<CT> cPLU = cP*cL*cU;
            if (c.lud().isTrans()) cPLU.transposeSelf();
            Assert(Equal(c,cPLU,ceps*normc),"Band C LU");

#if (XTEST & 16)
            const int Rnhi = std::min(N-1,nlo+nhi);
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L2(M);
            tmv::BandMatrix<T,stor> U2(M,M,0,Rnhi);
            tmv::Permutation P2(N);
            LU_Decompose(m,L2,U2,P2);
            PLU = P2*L2*U2;
            Assert(Equal(m,PLU,eps*normm),"Band LU2");

            tmv::LowerTriMatrix<CT,tmv::UnitDiag> cL2(M);
            tmv::BandMatrix<CT,stor> cU2(M,M,0,Rnhi);
            LU_Decompose(c,cL2,cU2,P2);
            cPLU = P2*cL2*cU2;
            Assert(Equal(c,cPLU,ceps*normc),"Band C LU2");

            LU_Decompose(c,cL2.conjugate(),cU2,P2);
            cPLU = P2*cL2.conjugate()*cU2;
            Assert(Equal(c,cPLU,ceps*normc),"Band C LU3");

            LU_Decompose(c,cL2,cU2.conjugate(),P2);
            cPLU = P2*cL2*cU2.conjugate();
            Assert(Equal(c,cPLU,ceps*normc),"Band C LU4");

            LU_Decompose(c,cL2.conjugate(),cU2.conjugate(),P2);
            cPLU = P2*cL2.conjugate()*cU2.conjugate();
            Assert(Equal(c,cPLU,ceps*normc),"Band C LU5");

            LU_Decompose(c.conjugate(),cL2,cU2,P2);
            cPLU = P2*cL2*cU2;
            Assert(Equal(c.conjugate(),cPLU,ceps*normc),"Band C LU6");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2,P2);
            cPLU = P2*cL2.conjugate()*cU2;
            Assert(Equal(c.conjugate(),cPLU,ceps*normc),"Band C LU7");

            LU_Decompose(c.conjugate(),cL2,cU2.conjugate(),P2);
            cPLU = P2*cL2*cU2.conjugate();
            Assert(Equal(c.conjugate(),cPLU,ceps*normc),"Band C LU8");

            LU_Decompose(c.conjugate(),cL2.conjugate(),cU2.conjugate(),P2);
            cPLU = P2*cL2.conjugate()*cU2.conjugate();
            Assert(Equal(c.conjugate(),cPLU,ceps*normc),"Band C LU9");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // QR Decomposition
        do {
            if (showstartdone) {
                std::cout<<"QR"<<std::endl;
            }
            tmv::Matrix<T> Q = m.qrd().getQ();
            tmv::BandMatrix<T> R = m.qrd().getR();
            tmv::Matrix<T> QR = Q*R;
            if (m.qrd().isTrans()) QR.transposeSelf();
            Assert(Equal(m,QR,eps*normm),"Band QR");

            tmv::Matrix<CT> cQ = c.qrd().getQ();
            tmv::BandMatrix<CT> cR = c.qrd().getR();
            tmv::Matrix<CT> cQR = cQ*cR;
            if (c.qrd().isTrans()) cQR.transposeSelf();
            Assert(Equal(c,cQR,ceps*normc),"Band C QR");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            const int Rnhi = std::min(N-1,nlo+nhi);
            QR_Decompose(m,Q,R);
            QR = Q*R;
            Assert(Equal(m,QR,eps*normm),"Band QR2");

            tmv::BandMatrix<T,stor> R2(N,N,0,Rnhi);
            QR_Decompose(m,R2);
            if (showacc) {
                std::cout<<"Norm(R-R2) = "<<Norm(R-R2)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Equal(
                    tmv::Matrix<T>(m.adjoint()/x)*tmv::Matrix<T>(m/x),
                    tmv::Matrix<T>(R2.adjoint()/x)*tmv::Matrix<T>(R2/x),
                    eps*(normm/x)*(normm/x)),"Band QR3 (RtR)");

            QR_Decompose(c,cQ,cR);
            cQR = cQ*cR;
            Assert(Equal(c,cQR,ceps*normc),"Band C QR2");

            tmv::BandMatrix<CT,stor> cR2(N,N,0,Rnhi);
            QR_Decompose(c,cR2);
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    tmv::Matrix<CT>(cR2.adjoint()/x)*tmv::Matrix<CT>(cR2/x),
                    ceps*(normc/x)*(normc/x)),"Band C QR3 (RtR)");

            QR_Decompose(c,cQ.conjugate(),cR2);
            cQR = cQ.conjugate()*cR2;
            Assert(Equal(c,cQR,ceps*normc),"Band C QR4");

            QR_Decompose(c,cQ,cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Equal(c,cQR,ceps*normc),"Band C QR5");

            QR_Decompose(c,cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Equal(c,cQR,ceps*normc),"Band C QR6");

            QR_Decompose(c,cR2);
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    tmv::Matrix<CT>(cR2.adjoint()/x)*tmv::Matrix<CT>(cR2/x),
                    ceps*(normc/x)*(normc/x)),"Band C QR7 (RtR)");

            QR_Decompose(c,cR2.conjugate());
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    tmv::Matrix<CT>(cR2.transpose()/x)*tmv::Matrix<CT>(cR2.conjugate()/x),
                    ceps*(normc/x)*(normc/x)),"Band C QR8 (RtR)");

            QR_Decompose(c.conjugate(),cQ,cR2);
            cQR = cQ*cR2;
            Assert(Equal(c.conjugate(),cQR,ceps*normc),"Band C QR9");

            QR_Decompose(c.conjugate(),cR2);
            Assert(Equal(
                    tmv::Matrix<CT>(c.transpose()/x)*tmv::Matrix<CT>(c.conjugate()/x),
                    tmv::Matrix<CT>(cR2.adjoint()/x)*tmv::Matrix<CT>(cR2/x),
                    ceps*(normc/x)*(normc/x)),"Band C QR10 (RtR)");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2);
            cQR = cQ.conjugate()*cR2;
            Assert(Equal(c.conjugate(),cQR,ceps*normc),"Band C QR11");

            QR_Decompose(c.conjugate(),cQ,cR2.conjugate());
            cQR = cQ*cR2.conjugate();
            Assert(Equal(c.conjugate(),cQR,ceps*normc),"Band C QR12");

            QR_Decompose(c.conjugate(),cQ.conjugate(),cR2.conjugate());
            cQR = cQ.conjugate()*cR2.conjugate();
            Assert(Equal(c.conjugate(),cQR,ceps*normc),"Band C QR13");

            QR_Decompose(c.conjugate(),cR2.conjugate());
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    tmv::Matrix<CT>(cR2.adjoint()/x)*tmv::Matrix<CT>(cR2/x),
                    ceps*(normc/x)*(normc/x)),"Band C QR15 (RtR)");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // SV Decomposition
        do {
            if (showstartdone) {
                std::cout<<"SV"<<std::endl;
            }
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> Vt = m.svd().getVt();
            if (showacc) {
                std::cout<<"U = "<<U<<std::endl;
                std::cout<<"Vt = "<<Vt<<std::endl;
                std::cout<<"S = "<<S<<std::endl;
                std::cout<<"m-USVt = "<<m-U*S*Vt<<std::endl;
                std::cout<<"Norm(m-USVt) = "<<Norm(m-U*S*Vt)<<
                    "  cf "<<eps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.transpose()*U-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VVt-1) = "<<Norm(Vt.transpose()*Vt-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.transpose()-T(1))<<
                    "  cf "<<eps<<std::endl;
            }
            Assert(Equal(m,U*S*Vt,eps*normm),"SV");

            tmv::Matrix<CT> cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<CT> cVt = c.svd().getVt();
            if (showacc) {
                std::cout<<"cU = "<<cU<<std::endl;
                std::cout<<"cVt = "<<cVt<<std::endl;
                std::cout<<"cS = "<<cS<<std::endl;
                std::cout<<"c-cUcScVt = "<<c-cU*cS*cVt<<std::endl;
                std::cout<<"Norm(c-cUcScVt) = "<<Norm(c-cU*cS*cVt)<<
                    "  cf "<<eps*normm<<std::endl;
                std::cout<<"Norm(cUtcU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(cVcVt-1) = "<<Norm(cVt.adjoint()*cVt-T(1))<<
                    "  cf "<<eps<<std::endl;
                std::cout<<"Norm(cVtcV-1) = "<<Norm(cVt*cVt.adjoint()-T(1))<<
                    "  cf "<<eps<<std::endl;
            }
            Assert(Equal(c,cU*cS*cVt,ceps*normc),"C SV");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            tmv::Matrix<T> U2(M,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> Vt2(N,N);
            SV_Decompose(m,U2,S2,Vt2);
            Assert(Equal(m,U2*S2*Vt2,eps*normm),"SV2");

            SV_Decompose(m,S2);
            if (showacc) {
                std::cout<<"S = "<<S.diag()<<std::endl;
                std::cout<<"S2 = "<<S2.diag()<<std::endl;
                std::cout<<"S2-S = "<<S2.diag()-S.diag()<<std::endl;
                std::cout<<"Norm(S2-S) = "<<Norm(S2-S)<<
                    "  cf "<<eps*normm<<std::endl;
            }
            Assert(Equal(S2,S,eps*normm),"SV3");
            SV_Decompose(m,U2,S2);
            Assert(Equal(S2,S,eps*normm),"SV4 S");
            Assert(Equal(
                    tmv::Matrix<T>(m/x)*tmv::Matrix<T>(m.transpose()/x),
                    U2*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*U2.transpose(),
                    eps*(normm/x)*(normm/x)),"SV3 U");
            SV_Decompose(m,S2,Vt2);
            Assert(Equal(S2,S,eps*normm),"SV5 S");
            Assert(Equal(
                    tmv::Matrix<T>(m.transpose()/x)*tmv::Matrix<T>(m/x),
                    Vt2.transpose()*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*Vt2,
                    eps*(normm/x)*(normm/x)),"SV5 Vt");

            tmv::Matrix<CT> cU2(M,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT> cVt2(N,N);
            SV_Decompose(c,cU2,cS2,cVt2);
            Assert(Equal(c,cU2*cS2*cVt2,eps*normc),"C SV2");

            SV_Decompose(c,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV3");
            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"C SV4 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV5 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"C SV5 Vt");

            SV_Decompose(c,cU2.conjugate(),cS2,cVt2);
            Assert(Equal(c,cU2.conjugate()*cS2*cVt2,eps*normc),"C SV6");
            SV_Decompose(c,cU2,cS2,cVt2.conjugate());
            Assert(Equal(c,cU2*cS2*cVt2.conjugate(),eps*normc),"C SV7");
            SV_Decompose(c,cU2.conjugate(),cS2,cVt2.conjugate());
            Assert(Equal(c,cU2.conjugate()*cS2*cVt2.conjugate(),eps*normc),
                   "C SV8");

            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV9 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"C SV9 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV10 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"C SV10 Vt");

            SV_Decompose(c.conjugate(),cU2,cS2,cVt2);
            Assert(Equal(c.conjugate(),cU2*cS2*cVt2,eps*normc),"C SV11");
            SV_Decompose(c.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV12");
            SV_Decompose(c.conjugate(),cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV13 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"C SV13 U");
            SV_Decompose(c.conjugate(),cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV14 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.transpose()/x)*tmv::Matrix<CT>(c.conjugate()/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"C SV14 Vt");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2,cVt2);
            Assert(Equal(c.conjugate(),cU2.conjugate()*cS2*cVt2,eps*normc),
                   "C SV15");
            SV_Decompose(c.conjugate(),cU2,cS2,cVt2.conjugate());
            Assert(Equal(c.conjugate(),cU2*cS2*cVt2.conjugate(),eps*normc),
                   "C SV16");
            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2,
                         cVt2.conjugate());
            Assert(Equal(
                    c.conjugate(),cU2.conjugate()*cS2*cVt2.conjugate(),
                    eps*normc),"C SV17");

            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"C SV18 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2.conjugate()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.transpose(),
                    ceps*(normc/x)*(normc/x)),"C SV18 U");
            SV_Decompose(c.conjugate(),cS2,cVt2.conjugate());
            Assert(Equal(cS2,cS,ceps*normc),"C SV19 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.transpose()/x)*tmv::Matrix<CT>(c.conjugate()/x),
                    cVt2.transpose()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2.conjugate(),
                    ceps*(normc/x)*(normc/x)),"C SV19 Vt");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);
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

