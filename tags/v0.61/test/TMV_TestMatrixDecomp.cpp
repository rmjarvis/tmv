#ifdef NDEBUG
#undef NDEBUG
#endif
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

template <class T, tmv::StorageType stor> void TestMatrixDecomp()
{
  for (int mattype = 0; mattype < 4; mattype++) {
    // mattype = 0  is Square
    // mattype = 1  is NonSquare slightly tall
    // mattype = 2  is NonSquare very tall
    // mattype = 3  is Singular

    int M = 8;
    const int N = 8;
    if (mattype == 1) M = 11;
    else if (mattype == 2) M = 45;

    tmv::Matrix<T,stor> m(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2.+4*i-5*j);
    if (mattype != 3) m /= T(10);
    m(3,3) = T(14);
    m(3,4) = T(-2);
    m(0,2) = T(7);
    m(M-1,N-4) = T(23);
    m(M-1,N-2) = T(13);
    m(M-1,N-1) = T(-10);
    if (mattype != 3) m.diag() *= T(30);

    tmv::Matrix<std::complex<T>,stor> c(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
      c(i,j) = std::complex<T>(2.+4*i-5*j,3.-i);
    if (mattype != 3) c /= T(10);
    c(3,3) = T(14);
    c(3,4) = T(-2);
    c(0,2) = T(7);
    c(M-1,N-4) = T(23);
    c(M-1,N-2) = T(13);
    c(M-1,N-1) = T(-10);
    if (mattype != 3) c.diag() *= T(30);

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 3) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(10);
      ceps *= T(10);
    }

    // LU Decomposition
    if (mattype == 0) {
      m.DivideUsing(tmv::LU);
      m.SetDiv();
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::UpperTriMatrix<T> U = m.LUD().GetU();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> PLU = L*U;
      PLU.ReversePermuteRows(p);
      if (m.LUD().IsTrans()) PLU.TransposeSelf();
      Assert(Norm(m-PLU) < eps*Norm(m),"LU");

      tmv::Matrix<T,stor> m2 = m;
      int p2[N];
      LU_Decompose(m2.View(),p2);
      L = m2.LowerTri(tmv::UnitDiag);
      U = m2.UpperTri();
      PLU = L*U;
      PLU.ReversePermuteRows(p2);
      Assert(Norm(m-PLU) < eps*Norm(m),"LU2");

      c.DivideUsing(tmv::LU);
      c.SetDiv();
      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::UpperTriMatrix<std::complex<T> > cU = c.LUD().GetU();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cPLU = cL*cU;
      cPLU.ReversePermuteRows(p);
      if (c.LUD().IsTrans()) cPLU.TransposeSelf();
      Assert(Norm(c-cPLU) < ceps*Norm(c),"C LU");

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
    }
  
    // QR Decomposition
    if (mattype != 4) {
      m.DivideUsing(tmv::QR);
      m.SetDiv();
      tmv::Matrix<T,stor> Q = m.QRD().GetQ();
      tmv::UpperTriMatrix<T> R = m.QRD().GetR();
      tmv::Matrix<T> QR = Q*R;
      Assert(Norm(m-QR) < eps*Norm(m),"QR");
      Assert(Norm(Q.Transpose()*Q-T(1)) < eps,"QR - QtQ");

      Q = m;
      QR_Decompose(Q.View(),R.View());
      QR = Q*R;
      Assert(Norm(m-QR) < eps*Norm(m),"QR2");
      Assert(Norm(Q.Transpose()*Q-T(1)) < eps,"QR2 - QtQ");

      Q = m;
      QR_Decompose(Q.View());
      Assert(Norm(R-Q.UpperTri()) < eps*Norm(R),"QR3");

      c.DivideUsing(tmv::QR);
      c.SetDiv();
      tmv::Matrix<std::complex<T>,stor> cQ = c.QRD().GetQ();
      tmv::UpperTriMatrix<std::complex<T> > cR = c.QRD().GetR();
      tmv::Matrix<std::complex<T> > cQR = cQ*cR;
      Assert(Norm(c-cQR) < ceps*Norm(c),"C QR");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QR - QtQ");

      cQ = c;
      QR_Decompose(cQ.View(),cR.View());
      cQR = cQ*cR;
      Assert(Norm(c-cQR) < ceps*Norm(c),"C QR2");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QR2 - QtQ");

      cQ = c;
      QR_Decompose(cQ.View());
      Assert(Norm(cR-cQ.UpperTri()) < ceps*Norm(cR),"C QR3");

      cQ.Conjugate() = c;
      QR_Decompose(cQ.Conjugate(),cR.View());
      cQR = cQ.Conjugate()*cR;
      Assert(Norm(c-cQR) < ceps*Norm(c),"C QR4");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QR4 - QtQ");

      cQ.Conjugate() = c;
      QR_Decompose(cQ.Conjugate());
      Assert(Norm(cR-cQ.Conjugate().UpperTri()) < ceps*Norm(cR),"C QR5");

      cQ = c;
      QR_Decompose(cQ.View(),cR.Conjugate());
      cQR = cQ*cR.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"C QR6");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QR6 - QtQ");

      cQ.Conjugate() = c;
      QR_Decompose(cQ.Conjugate(),cR.Conjugate());
      cQR = cQ.Conjugate()*cR.Conjugate();
      Assert(Norm(c-cQR) < ceps*Norm(c),"C QR7");
    }
    
    // QRP Decomposition
    for (int istrict = 0; istrict <= 1; istrict++) {
      bool strict = istrict == 1;
      m.DivideUsing(tmv::QRP);
      tmv::QRPDiv<T>::StrictQRP = strict;
      m.ReSetDiv();
      tmv::Matrix<T,stor> Q = m.QRPD().GetQ();
      tmv::UpperTriMatrix<T> R = m.QRPD().GetR();
      const int* p = m.QRPD().GetP();
      tmv::Matrix<T> QRP = Q*R;
      QRP.ReversePermuteCols(p);
      Assert(Norm(m-QRP) < eps*Norm(m),"QRP");
      Assert(Norm(Q.Transpose()*Q-T(1)) < eps,"QRP - QtQ");
      if (strict) for(size_t i=1;i<R.size();i++)
	Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP - strict");

      Q = m;
      int p2[N];
      QRP_Decompose(Q.View(),R.View(),p2,strict);
      QRP = Q*R;
      QRP.ReversePermuteCols(p2);
      Assert(Norm(m-QRP) < eps*Norm(m),"QRP2");
      Assert(Norm(Q.Transpose()*Q-T(1)) < eps,"QRP2 - QtQ");
      if (strict) for(size_t i=1;i<R.size();i++)
	Assert(std::abs(R(i,i))<=std::abs(R(i-1,i-1)),"QRP2 - strict");

      Q = m;
      QRP_Decompose(Q.View(),strict);
      Assert(Norm(R-Q.UpperTri()) < eps*Norm(R),"QRP3");

      c.DivideUsing(tmv::QRP);
      tmv::QRPDiv<std::complex<T> >::StrictQRP = strict;
      c.ReSetDiv();
      tmv::Matrix<std::complex<T>,stor> cQ = c.QRPD().GetQ();
      tmv::UpperTriMatrix<std::complex<T> > cR = c.QRPD().GetR();
      p = c.QRPD().GetP();
      tmv::Matrix<std::complex<T> > cQRP = cQ*cR;
      cQRP.ReversePermuteCols(p);
      Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QRP - QtQ");
      if (strict) for(size_t i=1;i<cR.size();i++)
	Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP - strict");

      cQ = c;
      QRP_Decompose(cQ.View(),cR.View(),p2,strict);
      cQRP = cQ*cR;
      cQRP.ReversePermuteCols(p2);
      Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP2");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QRP2 - QtQ");
      if (strict) for(size_t i=1;i<cR.size();i++)
	Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP2 - strict");

      cQ = c;
      QRP_Decompose(cQ.View(),strict);
      Assert(Norm(cR-cQ.UpperTri()) < ceps*Norm(cR),"C QRP3");

      cQ.Conjugate() = c;
      QRP_Decompose(cQ.Conjugate(),cR.View(),p2,strict);
      cQRP = cQ.Conjugate()*cR;
      cQRP.ReversePermuteCols(p2);
      Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP4");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QRP4 - QtQ");
      if (strict) for(size_t i=1;i<cR.size();i++)
	Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP4 - strict");

      cQ.Conjugate() = c;
      QRP_Decompose(cQ.Conjugate(),strict);
      Assert(Norm(cR-cQ.Conjugate().UpperTri()) < ceps*Norm(cR),"C QRP5");

      cQ = c;
      QRP_Decompose(cQ.View(),cR.Conjugate(),p2,strict);
      cQRP = cQ*cR.Conjugate();
      cQRP.ReversePermuteCols(p2);
      Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP6");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QRP6 - QtQ");
      if (strict) for(size_t i=1;i<cR.size();i++)
	Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP6 - strict");


      cQ.Conjugate() = c;
      QRP_Decompose(cQ.Conjugate(),cR.Conjugate(),p2,strict);
      cQRP = cQ.Conjugate()*cR.Conjugate();
      cQRP.ReversePermuteCols(p2);
      Assert(Norm(c-cQRP) < ceps*Norm(c),"C QRP7");
      Assert(Norm(cQ.Adjoint()*cQ-T(1)) < eps,"C QRP7 - QtQ");
      if (strict) for(size_t i=1;i<cR.size();i++)
	Assert(std::abs(cR(i,i))<=std::abs(cR(i-1,i-1)),"C QRP7 - strict");
    }

    // SV Decomposition
    {
      m.DivideUsing(tmv::SV);
      m.SaveDiv();
      m.SetDiv();
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"SV");
      Assert(Norm(U.Transpose()*U-T(1)) < eps*Norm(m),"SV - UtU");
      Assert(Norm(V.Transpose()*V-T(1)) < eps*Norm(m),"SV - VtV");
      Assert(Norm(V*V.Transpose()-T(1)) < eps*Norm(m),"SV - VVt");

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
      Assert(Norm(S2-S) < eps*Norm(S),"SV2");
      U2 = m;
      SV_Decompose(U2.View(),S2.View(),true);
      Assert(Norm(S2-S) < eps*Norm(S),"SV3 S");
      Assert(Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose()) < 
	  eps*Norm(m*m.Transpose()),"SV3 U");
      m2 = m;
      SV_Decompose(m2.View(),S2.View(),V2.View(),false);
      Assert(Norm(S2-S) < eps*Norm(S),"SV4 S");
      Assert(Norm(V2.Transpose()*S2*S2*V2-m.Transpose()*m) < 
	  eps*Norm(m.Transpose()*m),"SV4 V");

      c.DivideUsing(tmv::SV);
      c.SaveDiv();
      c.SetDiv();
      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      S = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*S*cV) < ceps*Norm(c),"C SV");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < eps*Norm(m),"C SV - UtU");
      Assert(Norm(cV.Adjoint()*cV-T(1)) < eps*Norm(m),"C SV - VtV");
      Assert(Norm(cV*cV.Adjoint()-T(1)) < eps*Norm(m),"C SV - VVt");

      tmv::Matrix<std::complex<T>,stor> cU2 = c;
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(cU2.View(),S2.View(),cV2.View(),true);
      Assert(Norm(c-cU2*S2*cV2) < eps*Norm(m),"C SV2");
      Assert(Norm(cU2.Adjoint()*cU2-T(1)) < eps*Norm(m),"C SV2 - UtU");
      Assert(Norm(cV2.Adjoint()*cV2-T(1)) < eps*Norm(m),"C SV2 - VtV");
      Assert(Norm(cV2*cV2.Adjoint()-T(1)) < eps*Norm(m),"C SV2 - VVt");

      tmv::Matrix<std::complex<T>,stor> c2 = c;
      SV_Decompose(c2.View(),S2.View(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV3");
      cU2 = c;
      SV_Decompose(cU2.View(),S2.View(),true);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV4 S");
      Assert(Norm(cU2*S2*S2*cU2.Adjoint()-c*c.Adjoint()) < 
	  ceps*Norm(c*c.Adjoint()),"C SV4 U");
      c2 = c;
      SV_Decompose(c2.View(),S2.View(),cV2.View(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV5 S");
      Assert(Norm(cV2.Adjoint()*S2*S2*cV2-c.Adjoint()*c) < 
	  ceps*Norm(c.Adjoint()*c),"C SV5 V");

      cU2.Conjugate() = c;
      SV_Decompose(cU2.Conjugate(),S2.View(),cV2.View(),true);
      Assert(Norm(c-cU2.Conjugate()*S2*cV2) < eps*Norm(m),"C SV6");
      Assert(Norm(cU2.Adjoint()*cU2-T(1)) < eps*Norm(m),"C SV6 - UtU");
      Assert(Norm(cV2.Adjoint()*cV2-T(1)) < eps*Norm(m),"C SV6 - VtV");
      Assert(Norm(cV2*cV2.Adjoint()-T(1)) < eps*Norm(m),"C SV6 - VVt");
      c2.Conjugate() = c;
      SV_Decompose(c2.Conjugate(),S2.View(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV7");
      cU2.Conjugate() = c;
      SV_Decompose(cU2.Conjugate(),S2.View(),true);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV8 S");
      Assert(Norm(cU2.Conjugate()*S2*S2*cU2.Transpose()-c*c.Adjoint()) < 
	  ceps*Norm(c*c.Adjoint()),"C SV8 U");
      c2.Conjugate() = c;
      SV_Decompose(c2.Conjugate(),S2.View(),cV2.View(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV9 S");
      Assert(Norm(cV2.Adjoint()*S2*S2*cV2-c.Adjoint()*c) < 
	  ceps*Norm(c.Adjoint()*c),"C SV9 V");
    
      cU2 = c;
      SV_Decompose(cU2.View(),S2.View(),cV2.Conjugate(),true);
      Assert(Norm(c-cU2*S2*cV2.Conjugate()) < eps*Norm(m),"C SV10");
      Assert(Norm(cU2.Adjoint()*cU2-T(1)) < eps*Norm(m),"C SV10 - UtU");
      Assert(Norm(cV2.Adjoint()*cV2-T(1)) < eps*Norm(m),"C SV10 - VtV");
      Assert(Norm(cV2*cV2.Adjoint()-T(1)) < eps*Norm(m),"C SV10 - VVt");
      c2 = c;
      SV_Decompose(c2.View(),S2.View(),cV2.Conjugate(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV11 S");
      Assert(Norm(cV2.Transpose()*S2*S2*cV2.Conjugate()-c.Adjoint()*c) < 
	  ceps*Norm(c.Adjoint()*c),"C SV9 V");

      cU2.Conjugate() = c;
      SV_Decompose(cU2.Conjugate(),S2.View(),cV2.Conjugate(),true);
      Assert(Norm(c-cU2.Conjugate()*S2*cV2.Conjugate()) < eps*Norm(m),"C SV12");
      Assert(Norm(cU2.Adjoint()*cU2-T(1)) < eps*Norm(m),"C SV12 - UtU");
      Assert(Norm(cV2.Adjoint()*cV2-T(1)) < eps*Norm(m),"C SV12 - VtV");
      Assert(Norm(cV2*cV2.Adjoint()-T(1)) < eps*Norm(m),"C SV12 - VVt");
      c2.Conjugate() = c;
      SV_Decompose(c2.Conjugate(),S2.View(),cV2.Conjugate(),false);
      Assert(Norm(S2-S) < ceps*Norm(S),"C SV13 S");
      Assert(Norm(cV2.Transpose()*S2*S2*cV2.Conjugate()-c.Adjoint()*c) < 
	  ceps*Norm(c.Adjoint()*c),"C SV13 V");
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
