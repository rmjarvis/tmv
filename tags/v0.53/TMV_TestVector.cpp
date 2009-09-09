
#include "TMV_Test.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"

using tmv::Vector;

template <class T> void TestVectorReal()
{
  const int N = 100;

  Vector<T> v(N);

  for (int i=0; i<N; ++i) v(i) = i;

  for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting Vector");

  tmv::VectorView<T> v2 = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v2(i) == T(2*i),
      "Reading Vector with stride = 2");

  for (int i=0; i<N/2; ++i) v2[i] = i + 1000;
  for (int i=0; i<N/2; ++i) Assert(v(2*i) == T(i+1000),
      "Writing Vector with stride = 2");

  Vector<T> v3 = v2;
  for (int i=0; i<N/2; ++i) Assert(v3[i] == v2[i],
      "Copying Vector with stride = 2");

  for (int i=0; i<N; ++i) v[i] = i;
  v.Swap(2,5);
  Assert(v(2) == T(5) && v(5) == T(2),"Swapping elements of Vector");
  v.Swap(2,5);
  Assert(v(2) == T(2) && v(5) == T(5),"Swapping elements of Vector");

  T sum = N*(N-1)/2;
  Assert(SumElements(v) == sum,"Vector SumElements(v)");

  v.ReverseSelf();
  for (int i=0; i<N; ++i) Assert(v(i) == T(N-i-1),"Reversing Vector");

  for (int i=0; i<N; ++i) v(i) = i+10.;
  v(23) = 10.*N;
  v(42) = 0.25;
  v(15) = -20.*N;
  size_t imax,imin;
  Assert(MaxAbsElement(v,&imax) == T(20*N),
      "MaxAbsElement of Vector did not return correct value");
  Assert(imax == 15,
      "MaxAbsElement of Vector did not return correct index");
  Assert(MinAbsElement(v,&imin) == T(0.25),
      "MinAbsElement of Vector did not return correct value");
  Assert(imin == 42,
      "MinAbsElement of Vector did not return correct index");
  Assert(MaxElement(v,&imax) == T(10*N),
      "MaxElement of Vector did not return correct value");
  Assert(imax == 23,
      "MaxElement of Vector did not return correct index");
  Assert(MinElement(v,&imin) == T(-20*N),
      "MinElement of Vector did not return correct value");
  Assert(imin == 15,
      "MinElement of Vector did not return correct index");

  Vector<T> a(N);
  Vector<T> b(N);
  for (int i=0; i<N; ++i) a(i) = T(3+i);

  b = a;
  for (int i=0; i<N; ++i) Assert(a(i) == b(i),"Vector1 = Vector2");

  Assert(a == b,"Testing Equality of Vectors");

  b(4) = 0;
  Assert(a != b,"Vector = Vector copied address, not values");

  Vector<T,tmv::FortranStyle> af(N);
  for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
  for (int i=1; i<=N; ++i) Assert(af(i) == a(i-1),"FortranStyle Vector access");
  tmv::ConstVectorView<T,tmv::FortranStyle> afcv = af.View();
  for (int i=1; i<=N; ++i) Assert(afcv(i) == a(i-1),"FortranStyle Vector CV access");
  tmv::VectorView<T,tmv::FortranStyle> afv = af.View();
  for (int i=1; i<=N; ++i) Assert(afv(i) == a(i-1),"FortranStyle Vector V access");
  Assert(a == af,"FortransStyle Vector = CStyle Vector");
  tmv::ConstVectorView<T> afc = af.View();
  Assert(afc == a,"CStyle View of FortransStyle Vector = CStyle Vector");
  Assert(tmv::ConstVectorView<T,tmv::FortranStyle>(afc) == a,
    "FortranStyle View of Vector == CStyle Vector");

  for (int i=0; i<N; ++i) b(i) = 5.+2*i;

  v = a+b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(8+3*i),"Adding Vectors");

  v = a-b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(-2-i),"Subtracting Vectors");

  a(0) = 1.;
  b(0) = 1.;
  for (int i=1; i<10; ++i) {
    a(i) = a(i-1)*2.;
    b(i) = b(i-1)/2.;
  }
  for (int i=10; i<N; ++i) a(i) = 0.;

  Assert(a*b == T(10),"Multiplying Vectors");

  Vector<T> c(5);
  c = v.SubVector(10,70,12);
  for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubVector");

  for(int i=0;i<N;++i) a(i) = i+10.;
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  T prod = 2900;
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(abs(a*b - prod) < EPS*Norm(a)*Norm(b),"Inner Product");
  Assert(abs(Norm(a+b) - normsum) < EPS*abs(Norm1(a)+Norm1(b)),"Vector Sum");
  Assert(abs(Norm(a-b) - normdiff) < EPS*abs(Norm1(a)+Norm1(b)),"Vector Diff");
}

template <class T> void TestVectorComplex()
{
  const int N = 100;

  Vector<complex<T> > v(N);
  for (int i=0; i<N; ++i) v(i) = complex<T>(i,i+1234);

  for (int i=0; i<N; ++i) Assert(v(i).real() == T(i),
      "CVector set");
  for (int i=0; i<N; ++i) Assert(v(i).imag() == T(i+1234),
      "CVector set");

  tmv::VectorView<complex<T> > v1(v.SubVector(0,N,2));
  for (int i=0; i<N/2; ++i) Assert(v1(i) == complex<T>(2*i,2*i+1234),
      "CVector stride=2");

  for (int i=0; i<N/2; ++i) v1[i] = complex<T>(i,i+1234);
  for (int i=0; i<N/2; ++i) Assert(v[2*i] == complex<T>(i,i+1234),
      "setting CVector with stride = 2");

  for (int i=0; i<N; ++i) v(i) = complex<T>(i,i+1234);

  v.Swap(2,5);
  Assert(v[2] == complex<T>(5,5+1234),"Swap in CVector");
  Assert(v[5] == complex<T>(2,2+1234),"Swap in CVector");
  v.Swap(2,5);

  Vector<complex<T> > v2 = v.Conjugate();

  for (int i=0; i<N; ++i) Assert(v2(i) == complex<T>(i,-i-1234),
      "Conjugate CVector");
  Assert(v2 == v.Conjugate(),"Conjugate == CVector");

  Assert(abs((v*v2).imag()) < EPS,"CVector * CVector");
  T norm1 = tmv::SQRT((v*v2).real());
  T norm2 = Norm(v);
  Assert(abs(norm1 - norm2) < EPS*norm1,"Norm CVector");

  Assert(v2 == v.ConjugateSelf(),"ConjugateSelf CVector");

  Vector<T> a(N);
  for(int i=0;i<N;++i) a(i) = i+10.;
  Vector<T> b(N);
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  Vector<complex<T> > ca = a;
  Assert(Norm(ca-a) < EPS*Norm(a),"Copy real V -> complex V");
  ca *= complex<T>(3,4)/T(5);
  Vector<complex<T> > cb = b*complex<T>(3,4)/T(5);

  complex<T> prod = T(29)*complex<T>(-28,96);
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(abs(ca*cb - prod) < EPS*Norm(ca)*Norm(cb),"CInner Product");
  Assert(abs(Norm(ca+cb) - normsum) < EPS*abs(Norm(ca)+Norm(cb)),"CVector Sum");
  Assert(abs(Norm(ca-cb) - normdiff) < EPS*abs(Norm(ca)+Norm(cb)),"CVector Diff");
}

template <class V1, class T> void DoTestVa(
    const V1& a, const Vector<T>& v, string label)
{
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");
  Assert(abs(Norm1(a)-Norm1(v)) <= EPS*abs(Norm1(v)),label+" Norm1");
  Assert(abs(Norm2(a)-Norm2(v)) <= EPS*abs(Norm2(v)),label+" Norm2");
  Assert(abs(NormInf(a)-NormInf(v)) <= EPS*abs(NormInf(v)),label+" NormInf");
  Assert(abs(NormSq(a)-NormSq(v)) <= EPS*abs(NormSq(v)),label+" NormSq");
}

template <class V1, class T> void DoTestV(
    const V1& a, const Vector<T>& v, string label)
{
  DoTestVa(a,v,label);
  Vector<T> v2 = v.SubVector(0,v.size(),5);
  DoTestVa(tmv::VectorView<T>(a).SubVector(0,v.size(),5),v2,label+" Step");
  DoTestVa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,v.size()-4,5),v2,label+" Step");
  Vector<T> v3 = v.Reverse();
  DoTestVa(a.Reverse(),v3,label+" Rev");
}

template <class V1, class T, class T2> void DoTestVXa(
    const V1& a, const Vector<T>& v, T2 x, string label)
{
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");
  Assert(Norm((x*a)-(x*v)) <= EPS*Norm(v)*abs(x),label+" x*a");
  Assert(Norm((a*x)-(x*v)) <= EPS*Norm(v)*abs(x),label+" a*x");
  Assert(Norm((a/x)-(v/x)) <= EPS*Norm(v)*abs(x),label+" v/x");
}

template <class V1, class T, class T2> void DoTestVX(
    const V1& a, const Vector<T>& v, T2 x, string label)
{
  DoTestVXa(a.View(),v,x,label);
  Vector<T> v2 = v.SubVector(0,v.size(),5);
  DoTestVXa(tmv::VectorView<T>(a).SubVector(0,v.size(),5),v2,x,label+" Step");
  DoTestVXa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,v.size()-4,5),v2,x,label+" Step");
  Vector<T> v3 = v.Reverse();
  DoTestVXa(a.Reverse(),v3,x,label+" Rev");
}

template <class V1, class V2, class T, class T2> void DoTestVVa(
    const V1& a, const V2& b, const Vector<T>& v1, const Vector<T2>& v2,
    string label)
{
  Assert(Norm(a-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(b-v2) <= EPS*Norm(v2),label+" b != v2");

  Assert(Norm((a+v2)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+v");
  Assert(Norm((v1+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" v+b");
  Assert(Norm((a+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+b");
  Assert(Norm((a-v2)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-v");
  Assert(Norm((v1-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" v-b");
  Assert(Norm((a-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-b");
  Assert(abs((a*v2)-(v1*v2)) <= EPS*abs(v1*v2),label+" a*v");
  Assert(abs((v1*b)-(v1*v2)) <= EPS*abs(v1*v2),label+" v*b");
  Assert(abs((a*b)-(v1*v2)) <= EPS*abs(v1*v2),label+" a*b");
}

template <class V1, class V2, class T, class T2> void DoTestVV(
    const V1& a, const V2& b, const Vector<T>& v1, const Vector<T2>& v2,
    string label)
{
  Vector<T> v1s = v1.SubVector(0,a.size(),5);
  Vector<T> v1x = v1.SubVector(0,a.size()/5);
  Vector<T> v1r = v1.Reverse();
  Vector<T2> v2s = v2.SubVector(0,a.size(),5);
  Vector<T2> v2x = v2.SubVector(0,a.size()/5);
  Vector<T2> v2r = v2.Reverse();

  DoTestVVa(a.View(),b.View(),v1,v2,label);
  DoTestVVa(tmv::VectorView<T>(a).SubVector(0,a.size(),5),
      tmv::ConstVectorView<T2>(b).SubVector(0,a.size()/5),v1s,v2x,label+" StepA");
  DoTestVVa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,a.size()-4,5),
      tmv::ConstVectorView<T2,tmv::FortranStyle>(b).SubVector(1,a.size()/5),v1s,v2x,label+" StepA");
  DoTestVVa(tmv::VectorView<T>(a).SubVector(0,a.size()/5),
      tmv::ConstVectorView<T2>(b).SubVector(0,a.size(),5),v1x,v2s,label+" StepB");
  DoTestVVa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,a.size()/5),
      tmv::ConstVectorView<T2,tmv::FortranStyle>(b).SubVector(1,a.size()-4,5),v1x,v2s,label+" StepB");
  DoTestVVa(tmv::VectorView<T>(a).SubVector(0,a.size(),5),
      tmv::ConstVectorView<T2>(b).SubVector(0,a.size(),5),v1s,v2s,label+" StepA,StepB");
  DoTestVVa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,a.size()-4,5),
      tmv::ConstVectorView<T2,tmv::FortranStyle>(b).SubVector(1,a.size()-4,5),v1s,v2s,label+" StepA,StepB");
  DoTestVVa(a.Reverse(),b.View(),v1r,v2,label+" RevA");
  DoTestVVa(a.View(),b.Reverse(),v1,v2r,label+" RevB");
  DoTestVVa(a.Reverse(),b.Reverse(),v1r,v2r,label+" RevA,RevB");
}

template <class T> void TestVectorArith()
{
  const size_t N = 100;
  Vector<T> a(N);
  for(int i=0;i<int(N);++i) a(i) = i+10.;
  Vector<T> b(N);
  for(int i=0;i<int(N);++i) b(i) = -3.*i+2.;
  Vector<T,tmv::FortranStyle> af(N);
  for(int i=1;i<=int(N);++i) af(i) = i-1+10.;
  Vector<T,tmv::FortranStyle> bf(N);
  for(int i=1;i<=int(N);++i) bf(i) = -3.*(i-1)+2.;

  Assert(a == af, "a == af");
  Assert(b == bf, "b == bf");

  Vector<complex<T> > ca = a*complex<T>(2,-1);;
  Vector<complex<T> > cb = b*complex<T>(-5,1);
  Vector<complex<T> > caf = af*complex<T>(2,-1);;
  Vector<complex<T> > cbf = bf*complex<T>(-5,1);

  Assert(ca == caf, "a == af");
  Assert(cb == cbf, "b == bf");

  T x = 12;
  complex<T> z(9,-2);
  
  Vector<T> v1 = a;
  Vector<T> v2 = b;
  Vector<complex<T> > cv1 = ca;
  Vector<complex<T> > cv2 = cb;

  DoTestV(a.View(),v1,"Vector R");
  DoTestV(af.View(),v1,"VectorF R");
  DoTestV(ca.View(),cv1,"Vector C");
  DoTestV(caf.View(),cv1,"VectorF C");

  DoTestVX(a.View(),v1,x,"Vector R,R");
  DoTestVX(af.View(),v1,x,"VectorF R,R");
  DoTestVX(a.View(),v1,z,"Vector R,C");
  DoTestVX(af.View(),v1,z,"VectorF R,C");
  DoTestVX(ca.View(),cv1,x,"Vector C,R");
  DoTestVX(caf.View(),cv1,x,"VectorF C,R");
  DoTestVX(ca.View(),cv1,z,"Vector C,C");
  DoTestVX(caf.View(),cv1,z,"VectorF C,C");

  DoTestVV(a.View(),b,v1,v2,"Vector R,R");
  DoTestVV(af.View(),b,v1,v2,"VectorF R,R");
  DoTestVV(a.View(),cb,v1,cv2,"Vector R,C");
  DoTestVV(af.View(),cb,v1,cv2,"VectorF R,C");
  DoTestVV(ca.View(),b,cv1,v2,"Vector C,R");
  DoTestVV(caf.View(),b,cv1,v2,"VectorF C,R");
  DoTestVV(ca.View(),cb,cv1,cv2,"Vector C,C");
  DoTestVV(caf.View(),cb,cv1,cv2,"VectorF C,C");
}

template <class T> void TestAllVector()
{
  TestVectorReal<T>();
  TestVectorComplex<T>();
  TestVectorArith<T>();
  cout<<"Vector<"<<tmv::Type(T())<<"> passed all tests\n";
}

template void TestAllVector<double>();
#ifndef NOFLOAT
template void TestAllVector<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllVector<long double>();
#endif
