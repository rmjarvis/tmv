#include "TMV.h"

template <class SM, class V1, class V2, class T, class T2> inline void DoTestOProd(
    const SM& a, const V1& v1, const V2& v2, 
    tmv::Matrix<T>& m1, const tmv::Matrix<T2>& m2,
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start OProd"<<std::endl;
    //std::cerr<<"a = "<<Type(a)<<" = "<<a<<std::endl;
    //std::cerr<<"v1 = "<<Type(v1)<<" = "<<v1<<std::endl;
    //std::cerr<<"v2 = "<<Type(v2)<<" = "<<v2<<std::endl;
    //std::cerr<<"m1 = "<<Type(m1)<<" = "<<m1<<std::endl;
    //std::cerr<<"m2 = "<<Type(m2)<<" = "<<m2<<std::endl;
  }

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a != m1");
  Assert(Norm((v1^v2)-m2) <= EPS*Norm(m2),label+" v1^v2 != m2");
  
  a += v1^v2;
  m1 += m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += v1^v2");
  a -= v1^v2;
  m1 -= m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= v1^v2");
  a += T(7) * (v1^v2);
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += 7 * (v1^v2)");
  a -= T(7) * (v1^v2);
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= 7 * (v1^v2)");
  a += (T(7) * v1)^v2;
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += (7*v1) ^ v2)");
  a -= (T(7) * v1)^v2;
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= (7*v1) ^ v2");
  a += v1 ^ (T(7) * v2);
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += v1 ^ (7*v2)");
  a -= v1 ^ (T(7) * v2);
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= v1 ^ (7*v2)");
  a += T2(5) * (v1^v2);
  m1 += T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += 5 * (v1^v2)");
  a -= T2(5) * (v1^v2);
  m1 -= T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= 5 * (v1^v2)");
  a += (T2(5) * v1)^v2;
  m1 += T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += (5*v1) ^ v2)");
  a -= (T2(5) * v1)^v2;
  m1 -= T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= (5*v1) ^ v2");
  a += v1 ^ (T2(5) * v2);
  m1 += T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a += v1 ^ (5*v2)");
  a -= v1 ^ (T2(5) * v2);
  m1 -= T2(5) * m2;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a -= v1 ^ (5*v2)");

  if (showstartdone)
    std::cout<<"Done OProd"<<std::endl;
}

template <class SM, class V1, class V2, class T, class T2> inline void TestOProd(
    const SM& a, const V1& v1, const V2& v2, 
    tmv::Matrix<T>& m1, const tmv::Matrix<T2>& m2, std::string label)
{
  tmv::Matrix<T> m1t = m1.Transpose();
  tmv::Matrix<T> m1c = m1.Conjugate();
  tmv::Matrix<T> m1a = m1.Adjoint();
  tmv::Matrix<T2> m2t = m2.Transpose();

  DoTestOProd(a.View(),v1,v2,m1,m2,label);
  DoTestOProd(Transpose(a),v2,v1,m1t,m2t,label+" TransA");
  DoTestOProd(Conjugate(a),v1,v2,m1c,m2,label+" ConjA");
  DoTestOProd(Adjoint(a),v2,v1,m1a,m2t,label+" AdjA");
}

