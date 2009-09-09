#ifndef TMV_TEST_H
#define TMV_TEST_H

//#define LONGDOUBLE

#include <iostream>
#include <complex>
#include <vector>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::complex;
using std::abs;
using std::vector;
using std::ostream;
using std::string;

#define EPS (10*tmv::Epsilon<T>())

inline void myerror(const char* s, const char* s2="", const char* s3="")
{
  cout << "Error in test: "<< s << s2 << s3 << endl;
  exit(1);
}

inline void Assert(bool x,string s)
{
  if (!(x)) myerror(s.c_str()); 
#ifdef SHOWTESTS
  else cout<<"Passed: "<<s<<endl;
#endif
}

#ifdef TESTDIV
template <class T, class MT> void CheckDecomposition(
    const MT& m, tmv::DivType dt) 
{
  tmv::Matrix<T> DecompM(m.colsize(),m.rowsize());
  RealType(T) nm(0);
  switch (dt) {
    case tmv::LU : 
      {
#ifdef SHOWCHECK
	cout<<"M = "<<m<<endl;
	cout<<"L = "<<m.LUD().LU_GetL()<<endl;
	cout<<"U = "<<m.LUD().LU_GetU()<<endl;
	cout<<"P = "<<Matrix<double>(m.LUD().LU_GetP())<<endl;
	cout<<"Q = "<<Matrix<double>(m.LUD().LU_GetQ())<<endl;
#endif
	DecompM = m.LUD().LU_GetP()*m.LUD().LU_GetL()*m.LUD().LU_GetU()*m.LUD().LU_GetQ();
	nm = Norm(m.LUD().LU_GetL())*Norm(m.LUD().LU_GetU());
      } break;
    case tmv::QR : case tmv::QRP :
      {
#ifdef SHOWCHECK
	cout<<"M = "<<m<<endl;
	cout<<"Q = "<<m.QRD().QR_GetQ()<<endl;
	cout<<"R = "<<m.QRD().QR_GetR()<<endl;
	cout<<"P = "<<Matrix<double>(m.QRD().QR_GetP())<<endl;
#endif
	if (m.QRD().QR_IsTrans()) {
	  DecompM = Transpose(m.QRD().QR_GetQ()*m.QRD().QR_GetR()*m.QRD().QR_GetP());
	} else {
	  DecompM = m.QRD().QR_GetQ()*m.QRD().QR_GetR()*m.QRD().QR_GetP();
	}
	nm = Norm(m.QRD().QR_GetQ())*Norm(m.QRD().QR_GetR());
      } break;
    case tmv::SV :
      {
#ifdef SHOWCHECK
	cout<<"M = "<<m<<endl;
	cout<<"U = "<<m.SVD().SV_GetU()<<endl;
	cout<<"S = "<<m.SVD().SV_GetS()<<endl;
	cout<<"V = "<<m.SVD().SV_GetV()<<endl;
#endif
	DecompM = m.SVD().SV_GetU()*tmv::Diag(m.SVD().SV_GetS())*m.SVD().SV_GetV();
	nm = Norm(m.SVD().SV_GetU())*Norm(m.SVD().SV_GetS())*Norm(m.SVD().SV_GetV());
      } break;
    default : cout<<"Invalid divtype\n";
  }
#ifdef SHOWACC
  cout << "DecompM = "<<DecompM<<endl;
  cout << "Norm(M-DecompM)/norm = "<<Norm(m-DecompM)/nm;
  cout << "  (EPS = "<<EPS<<")\n";
#endif
  Assert(Norm(m-DecompM) < EPS*nm,"DecompM != M");
}
#endif // TESTDIV

#endif
