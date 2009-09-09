// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestBandArith.h"

template <class T, tmv::StorageType stor> void TestBandDecomp()
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


    // LU Decomposition
    if (mattype == 0) {
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::UpperTriMatrix<T> U = m.LUD().GetU();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> PLU = L*U;
      PLU.ReversePermuteRows(p);
      if (m.LUD().IsTrans()) PLU.TransposeSelf();
      Assert(Norm(m-PLU) < eps*Norm(m),"Band LU");

      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::UpperTriMatrix<std::complex<T> > cU = c.LUD().GetU();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cPLU = cL*cU;
      cPLU.ReversePermuteRows(p);
      if (c.LUD().IsTrans()) cPLU.TransposeSelf();
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU");

#ifdef XTEST
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L2(M);
      tmv::BandMatrix<T,stor> U2(M,M,0,nlo+nhi);
      int p2[N];
      LU_Decompose(m,L2.View(),U2.View(),p2);
      PLU = L2*U2;
      PLU.ReversePermuteRows(p2);
      Assert(Norm(m-PLU) < eps*Norm(m),"Band LU2");

      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL2(M);
      tmv::BandMatrix<std::complex<T>,stor> cU2(M,M,0,nlo+nhi);
      LU_Decompose(c,cL2.View(),cU2.View(),p2);
      cPLU = cL2*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU2");

      LU_Decompose(c,cL2.Conjugate(),cU2.View(),p2);
      cPLU = cL2.Conjugate()*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU3");

      LU_Decompose(c,cL2.View(),cU2.Conjugate(),p2);
      cPLU = cL2*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU4");

      LU_Decompose(c,cL2.Conjugate(),cU2.Conjugate(),p2);
      cPLU = cL2.Conjugate()*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c-cPLU) < ceps*Norm(c),"Band C LU5");

      LU_Decompose(c.Conjugate(),cL2.View(),cU2.View(),p2);
      cPLU = cL2*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU6");

      LU_Decompose(c.Conjugate(),cL2.Conjugate(),cU2.View(),p2);
      cPLU = cL2.Conjugate()*cU2;
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU7");

      LU_Decompose(c.Conjugate(),cL2.View(),cU2.Conjugate(),p2);
      cPLU = cL2*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU8");

      LU_Decompose(c.Conjugate(),cL2.Conjugate(),cU2.Conjugate(),p2);
      cPLU = cL2.Conjugate()*cU2.Conjugate();
      cPLU.ReversePermuteRows(p2);
      Assert(Norm(c.Conjugate()-cPLU) < ceps*Norm(c),"Band C LU9");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // QR Decomposition
    if (mattype != 4) {
      tmv::Matrix<T> Q = m.QRD().GetQ();
      tmv::BandMatrix<T> R = m.QRD().GetR();
      tmv::Matrix<T> QR = Q*R;
      if (m.QRD().IsTrans()) QR.TransposeSelf();
      Assert(Norm(m-QR) < eps*Norm(m),"Band QR");

      tmv::Matrix<std::complex<T> > cQ = c.QRD().GetQ();
      tmv::BandMatrix<std::complex<T> > cR = c.QRD().GetR();
      tmv::Matrix<std::complex<T> > cQR = cQ*cR;
      if (c.QRD().IsTrans()) cQR.TransposeSelf();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR");

#ifdef XTEST
      QR_Decompose(m,Q.View(),R.View());
      QR = Q*R;
      Assert(Norm(m-QR) < eps*Norm(m),"Band QR2");

      tmv::BandMatrix<T,stor> R2(N,N,0,nlo+nhi);
      QR_Decompose(m,R2.View());
      Assert(Norm(R-R2) < eps*Norm(R),"Band QR3");

      QR_Decompose(c,cQ.View(),cR.View());
      cQR = cQ*cR;
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR2");

      tmv::BandMatrix<std::complex<T>,stor> cR2(N,N,0,nlo+nhi);
      QR_Decompose(c,cR2.View());
      Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR3");

      QR_Decompose(c,cQ.Conjugate(),cR2.View());
      cQR = cQ.Conjugate()*cR2;
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR4");

      QR_Decompose(c,cQ.View(),cR2.Conjugate());
      cQR = cQ*cR2.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR5");

      QR_Decompose(c,cQ.Conjugate(),cR2.Conjugate());
      cQR = cQ.Conjugate()*cR2.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"Band C QR6");

      QR_Decompose(c,cR2.View());
      Assert(Norm(cR-cR2) < ceps*Norm(cR),"Band C QR7");

      QR_Decompose(c,cR2.Conjugate());
      Assert(Norm(cR-cR2.Conjugate()) < ceps*Norm(cR),"Band C QR8");

      QR_Decompose(c.Conjugate(),cQ.View(),cR2.View());
      cQR = cQ*cR2;
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR9");

      QR_Decompose(c.Conjugate(),cR2.View());
      Assert(Norm(cR.Conjugate()-cR2) < ceps*Norm(cR),"Band C QR10");

      QR_Decompose(c.Conjugate(),cQ.Conjugate(),cR2.View());
      cQR = cQ.Conjugate()*cR2;
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR11");

      QR_Decompose(c.Conjugate(),cQ.View(),cR2.Conjugate());
      cQR = cQ*cR2.Conjugate();
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR12");

      QR_Decompose(c.Conjugate(),cQ.Conjugate(),cR2.Conjugate());
      cQR = cQ.Conjugate()*cR2.Conjugate();
      Assert(Norm(c.Conjugate()-cQR) < ceps*Norm(c),"Band C QR13");

      QR_Decompose(c.Conjugate(),cR2.View());
      Assert(Norm(cR.Conjugate()-cR2) < ceps*Norm(cR),"Band C QR14");

      QR_Decompose(c.Conjugate(),cR2.Conjugate());
      Assert(Norm(cR.Conjugate()-cR2.Conjugate()) < ceps*Norm(cR),"Band C QR15");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // SV Decomposition
    {
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"SV");

      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      tmv::DiagMatrix<T> cS = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"C SV");

#ifdef XTEST
      tmv::Matrix<T> U2(M,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SV2");

      SV_Decompose(m,S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV4 S");
      Assert(Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose()) < 
          eps*Norm(m)*Norm(m),"SV3 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"SV5 S");
      Assert(Norm(V2.Transpose()*S2*S2*V2-m.Transpose()*m) < 
          eps*Norm(m)*Norm(m),"SV5 V");

      tmv::Matrix<std::complex<T> > cU2(M,N);
      tmv::DiagMatrix<T> cS2(N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),cS2.View(),cV2.View());
      Assert(Norm(c-cU2*cS2*cV2) < eps*Norm(c),"C SV2");

      SV_Decompose(c,cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV3");
      SV_Decompose(c,cU2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV4 S");
      Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) < 
          ceps*Norm(c)*Norm(c),"C SV4 U");
      SV_Decompose(c,cS2.View(),cV2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV5 S");
      Assert(Norm(cV2.Adjoint()*cS2*cS2*cV2-c.Adjoint()*c) < 
          ceps*Norm(c)*Norm(c),"C SV5 V");

      SV_Decompose(c,cU2.Conjugate(),cS2.View(),cV2.View());
      Assert(Norm(c-cU2.Conjugate()*cS2*cV2) < eps*Norm(c),"C SV6");
      SV_Decompose(c,cU2.View(),cS2.View(),cV2.Conjugate());
      Assert(Norm(c-cU2*cS2*cV2.Conjugate()) < eps*Norm(c),"C SV7");
      SV_Decompose(c,cU2.Conjugate(),cS2.View(),cV2.Conjugate());
      Assert(Norm(c-cU2.Conjugate()*cS2*cV2.Conjugate()) < eps*Norm(c),"C SV8");

      SV_Decompose(c,cU2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV9 S");
      Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) < 
          ceps*Norm(c)*Norm(c),"C SV9 U");
      SV_Decompose(c,cS2.View(),cV2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV10 S");
      Assert(Norm(cV2.Adjoint()*cS2*cS2*cV2-c.Adjoint()*c) < 
          ceps*Norm(c)*Norm(c),"C SV10 V");

      SV_Decompose(c.Conjugate(),cU2.View(),cS2.View(),cV2.View());
      Assert(Norm(c.Conjugate()-cU2*cS2*cV2) < eps*Norm(c),"C SV11");
      SV_Decompose(c.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV12");
      SV_Decompose(c.Conjugate(),cU2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV13 S");
      Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c.Conjugate()*c.Transpose()) < 
          ceps*Norm(c)*Norm(c),"C SV13 U");
      SV_Decompose(c.Conjugate(),cS2.View(),cV2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV14 S");
      Assert(Norm(cV2.Adjoint()*cS2*cS2*cV2-c.Transpose()*c.Conjugate()) < 
          ceps*Norm(c)*Norm(c),"C SV14 V");

      SV_Decompose(c.Conjugate(),cU2.Conjugate(),cS2.View(),cV2.View());
      Assert(Norm(c.Conjugate()-cU2.Conjugate()*cS2*cV2) < eps*Norm(c),"C SV15");
      SV_Decompose(c.Conjugate(),cU2.View(),cS2.View(),cV2.Conjugate());
      Assert(Norm(c.Conjugate()-cU2*cS2*cV2.Conjugate()) < eps*Norm(c),"C SV16");
      SV_Decompose(c.Conjugate(),cU2.Conjugate(),cS2.View(),cV2.Conjugate());
      Assert(Norm(c.Conjugate()-cU2.Conjugate()*cS2*cV2.Conjugate()) < 
          eps*Norm(c),"C SV17");

      SV_Decompose(c.Conjugate(),cU2.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV18 S");
      Assert(Norm(cU2.Conjugate()*cS2*cS2*cU2.Transpose()-
            c.Conjugate()*c.Transpose())<ceps*Norm(c)*Norm(c),"C SV18 U");
      SV_Decompose(c.Conjugate(),cS2.View(),cV2.Conjugate());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"C SV19 S");
      Assert(Norm(cV2.Transpose()*cS2*cS2*cV2.Conjugate()-
            c.Transpose()*c.Conjugate())<ceps*Norm(c)*Norm(c),"C SV19 V");
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

