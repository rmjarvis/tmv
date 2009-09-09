#include "TMV_Test.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
#include <fstream>

template <class T> inline void TestVectorReal()
{
  const int N = 100;

  tmv::Vector<T> v(N);

  for (int i=0; i<N; ++i) v(i) = i;

  for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting Vector");

  tmv::VectorView<T> v2 = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v2(i) == T(2*i),
      "Reading Vector with stride = 2");

  for (int i=0; i<N/2; ++i) v2[i] = i + 1000;
  for (int i=0; i<N/2; ++i) Assert(v(2*i) == T(i+1000),
      "Writing Vector with stride = 2");

  tmv::Vector<T> v3 = v2;
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

  tmv::Vector<T> a(N);
  tmv::Vector<T> b(N);
  for (int i=0; i<N; ++i) a(i) = T(3+i);

  b = a;
  for (int i=0; i<N; ++i) Assert(a(i) == b(i),"Vector1 = Vector2");

  Assert(a == b,"Testing Equality of Vectors");

  b(4) = 0;
  Assert(a != b,"Vector = Vector copied address, not values");

  tmv::Vector<T,tmv::FortranStyle> af(N);
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

  tmv::Vector<T> c(5);
  c = v.SubVector(10,70,12);
  for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubVector");

  for(int i=0;i<N;++i) a(i) = i+10.;
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  T prod = 2900;
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(std::abs(a*b - prod) < EPS*Norm(a)*Norm(b),"Inner Product");
  Assert(std::abs(Norm(a+b) - normsum) < EPS*std::abs(Norm1(a)+Norm1(b)),"Vector Sum");
  Assert(std::abs(Norm(a-b) - normdiff) < EPS*std::abs(Norm1(a)+Norm1(b)),"Vector Diff");

  const size_t NN=20;
  T w1[NN] = {3.3,1.2,5.4,-1.2,4.3,-9.4,0.,-2.,4.,-11.5,-12.,14.,33.,1.,-9.3,-3.9,4.9,10.,-31.,1.e-33};

  tmv::Vector<T> w(NN,w1);
  tmv::Vector<T> origw = w;
  size_t perm[NN];

  if (showacc)
    std::cout<<"unsorted w = "<<w<<std::endl;
  w.Sort(perm);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1) <= w(i),"Sort real Vector");
  }
  if (showacc)
    std::cout<<"sorted w = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ABS_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(std::abs(w(i-1)) <= std::abs(w(i)),"Sort real Vector abs");
  }
  if (showacc)
    std::cout<<"sorted w abs = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1) >= w(i),"Sort real Vector descend");
  }
  if (showacc)
    std::cout<<"sorted w descend = "<<w<<std::endl;

  w.Sort(0);
  w.ReversePermute(perm);
  Assert(w==origw,"Reverse permute sorted Vector = orig");
  w.Sort(0);
  origw.Permute(perm);
  Assert(w==origw,"Permute Vector = sorted Vector");

}

template <class T> inline void TestVectorComplex()
{
  const int N = 100;

  tmv::Vector<std::complex<T> > v(N);
  for (int i=0; i<N; ++i) v(i) = std::complex<T>(i,i+1234);

  for (int i=0; i<N; ++i) Assert(v(i).real() == T(i),
      "CVector set");
  for (int i=0; i<N; ++i) Assert(v(i).imag() == T(i+1234),
      "CVector set");

  tmv::VectorView<std::complex<T> > v1(v.SubVector(0,N,2));
  for (int i=0; i<N/2; ++i) Assert(v1(i) == std::complex<T>(2*i,2*i+1234),
      "CVector stride=2");

  for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(i,i+1234);
  for (int i=0; i<N/2; ++i) Assert(v[2*i] == std::complex<T>(i,i+1234),
      "setting CVector with stride = 2");

  for (int i=0; i<N; ++i) v(i) = std::complex<T>(i,i+1234);

  v.Swap(2,5);
  Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CVector");
  Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CVector");
  v.Swap(2,5);

  tmv::Vector<std::complex<T> > v2 = v.Conjugate();

  for (int i=0; i<N; ++i) Assert(v2(i) == std::complex<T>(i,-i-1234),
      "Conjugate CVector");
  Assert(v2 == v.Conjugate(),"Conjugate == CVector");

  Assert(std::abs((v*v2).imag()) < EPS,"CVector * CVector");
  T norm1 = tmv::SQRT((v*v2).real());
  T norm2 = Norm(v);
  if (showacc) {
    std::cout<<"v = "<<v<<std::endl;
    std::cout<<"v2 = "<<v2<<std::endl;
    std::cout<<"v*v2 = "<<v*v2<<std::endl;
    std::cout<<"norm1 = "<<norm1<<std::endl;
    std::cout<<"norm2 = "<<norm2<<std::endl;
  }
  Assert(std::abs(norm1 - norm2) < EPS*norm1,"Norm CVector");

  Assert(v2 == v.ConjugateSelf(),"ConjugateSelf CVector");

  tmv::Vector<T> a(N);
  for(int i=0;i<N;++i) a(i) = i+10.;
  tmv::Vector<T> b(N);
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  tmv::Vector<std::complex<T> > ca = a;
  Assert(Norm(ca-a) < EPS*Norm(a),"Copy real V -> complex V");
  ca *= std::complex<T>(3,4)/T(5);
  tmv::Vector<std::complex<T> > cb = b*std::complex<T>(3,4)/T(5);

  std::complex<T> prod = T(29)*std::complex<T>(-28,96);
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(std::abs(ca*cb - prod) < EPS*Norm(ca)*Norm(cb),"CInner Product");
  Assert(std::abs(Norm(ca+cb) - normsum) < EPS*std::abs(Norm(ca)+Norm(cb)),"CVector Sum");
  Assert(std::abs(Norm(ca-cb) - normdiff) < EPS*std::abs(Norm(ca)+Norm(cb)),"CVector Diff");

  const size_t NN=20;
  T rw1[NN] = {3.3,1.2,5.4,-1.2,4.3,-9.4,0.,-2.,4.,-11.5,-12.,14.,33.,1.,-9.3,-3.9,4.9,10.,-31.,1.e-33};
  T iw1[NN] = {1.4,9.8,-0.2,-8.6,3.0,-4.4,3.,9.,-1.9,-11.5,11.1,-140.,-23.,11.,5.2,-3.9,4.9,99.,-71.,-0.5};

  tmv::Vector<T> rw(NN,rw1);
  tmv::Vector<T> iw(NN,iw1);
  tmv::Vector<std::complex<T> > w(NN);
  w.Real() = rw;
  w.Imag() = iw;
  tmv::Vector<std::complex<T> > origw = w;
  size_t perm[NN];

  if (showacc)
    std::cout<<"unsorted w = "<<w<<std::endl;
  w.Sort(perm);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1).real() <= w(i).real(),"Sort complex Vector");
  }
  if (showacc)
    std::cout<<"sorted w = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ABS_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(std::abs(w(i-1)) <= std::abs(w(i)),"Sort complex Vector abs");
  }
  if (showacc)
    std::cout<<"sorted w abs = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND,tmv::IMAG_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(imag(w(i-1)) >= imag(w(i)),"Sort complex Vector descend");
  }
  if (showacc)
    std::cout<<"sorted w imag descend = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1).real() >= w(i).real(),"Sort complex Vector descend");
  }
  if (showacc)
    std::cout<<"sorted w imag descend = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ARG_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(arg(w(i-1)) <= arg(w(i)),"Sort complex Vector descend");
  }
  if (showacc)
    std::cout<<"sorted w arg = "<<w<<std::endl;

  w.Sort(0);
  w.ReversePermute(perm);
  Assert(w==origw,"Reverse permute sorted Vector = orig");
  w.Sort(0);
  origw.Permute(perm);
  Assert(w==origw,"Permute Vector = sorted Vector");
}

template <class V1, class T> inline void DoTestVa(
    const V1& a, const tmv::Vector<T>& v, std::string label)
{
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");
  Assert(std::abs(Norm1(a)-Norm1(v)) <= EPS*std::abs(Norm1(v)),label+" Norm1");
  Assert(std::abs(Norm2(a)-Norm2(v)) <= EPS*std::abs(Norm2(v)),label+" Norm2");
  Assert(std::abs(NormInf(a)-NormInf(v)) <= EPS*std::abs(NormInf(v)),label+" NormInf");
  Assert(std::abs(NormSq(a)-NormSq(v)) <= EPS*std::abs(NormSq(v)),label+" NormSq");
}

template <class V1, class T> inline void DoTestV(
    const V1& a, const tmv::Vector<T>& v, std::string label)
{
  DoTestVa(a,v,label);
  tmv::Vector<T> v2 = v.SubVector(0,v.size(),5);
  DoTestVa(tmv::VectorView<T>(a).SubVector(0,v.size(),5),v2,label+" Step");
  DoTestVa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,v.size()-4,5),v2,label+" Step");
  tmv::Vector<T> v3 = v.Reverse();
  DoTestVa(a.Reverse(),v3,label+" Rev");
}

template <class V1, class T, class T2> inline void DoTestVXa(
    const V1& a, const tmv::Vector<T>& v, T2 x, std::string label)
{
  Assert(Norm(a-v) <= EPS*Norm(v),label+" a != v");
  Assert(Norm((x*a)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" x*a");
  Assert(Norm((a*x)-(x*v)) <= EPS*Norm(v)*std::abs(x),label+" a*x");
  Assert(Norm((a/x)-(v/x)) <= EPS*Norm(v)*std::abs(x),label+" v/x");
}

template <class V1, class T, class T2> inline void DoTestVX(
    const V1& a, const tmv::Vector<T>& v, T2 x, std::string label)
{
  DoTestVXa(a.View(),v,x,label);
  tmv::Vector<T> v2 = v.SubVector(0,v.size(),5);
  DoTestVXa(tmv::VectorView<T>(a).SubVector(0,v.size(),5),v2,x,label+" Step");
  DoTestVXa(tmv::VectorView<T,tmv::FortranStyle>(a).SubVector(1,v.size()-4,5),v2,x,label+" Step");
  tmv::Vector<T> v3 = v.Reverse();
  DoTestVXa(a.Reverse(),v3,x,label+" Rev");
}

template <class V1, class V2, class T, class T2> inline void DoTestVVa(
    const V1& a, const V2& b, const tmv::Vector<T>& v1, const tmv::Vector<T2>& v2,
    std::string label)
{
  Assert(Norm(a-v1) <= EPS*Norm(v1),label+" a != v1");
  Assert(Norm(b-v2) <= EPS*Norm(v2),label+" b != v2");

  Assert(Norm((a+v2)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+v");
  Assert(Norm((v1+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" v+b");
  Assert(Norm((a+b)-(v1+v2)) <= EPS*Norm(v1+v2),label+" a+b");
  Assert(Norm((a-v2)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-v");
  Assert(Norm((v1-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" v-b");
  Assert(Norm((a-b)-(v1-v2)) <= EPS*Norm(v1+v2),label+" a-b");
  Assert(std::abs((a*v2)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" a*v");
  Assert(std::abs((v1*b)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" v*b");
  Assert(std::abs((a*b)-(v1*v2)) <= EPS*std::abs(v1*v2),label+" a*b");
}

template <class V1, class V2, class T, class T2> inline void DoTestVV(
    const V1& a, const V2& b, const tmv::Vector<T>& v1, const tmv::Vector<T2>& v2,
    std::string label)
{
  tmv::Vector<T> v1s = v1.SubVector(0,a.size(),5);
  tmv::Vector<T> v1x = v1.SubVector(0,a.size()/5);
  tmv::Vector<T> v1r = v1.Reverse();
  tmv::Vector<T2> v2s = v2.SubVector(0,a.size(),5);
  tmv::Vector<T2> v2x = v2.SubVector(0,a.size()/5);
  tmv::Vector<T2> v2r = v2.Reverse();

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

template <class T> inline void TestVectorArith()
{
  const size_t N = 100;
  tmv::Vector<T> a(N);
  for(int i=0;i<int(N);++i) a(i) = i+10.;
  tmv::Vector<T> b(N);
  for(int i=0;i<int(N);++i) b(i) = -3.*i+2.;
  tmv::Vector<T,tmv::FortranStyle> af(N);
  for(int i=1;i<=int(N);++i) af(i) = i-1+10.;
  tmv::Vector<T,tmv::FortranStyle> bf(N);
  for(int i=1;i<=int(N);++i) bf(i) = -3.*(i-1)+2.;

  Assert(a == af, "a == af");
  Assert(b == bf, "b == bf");

  tmv::Vector<std::complex<T> > ca = a*std::complex<T>(2,-1);;
  tmv::Vector<std::complex<T> > cb = b*std::complex<T>(-5,1);
  tmv::Vector<std::complex<T> > caf = af*std::complex<T>(2,-1);;
  tmv::Vector<std::complex<T> > cbf = bf*std::complex<T>(-5,1);

  Assert(ca == caf, "a == af");
  Assert(cb == cbf, "b == bf");

  T x = 12;
  std::complex<T> z(9,-2);
  
  tmv::Vector<T> v1 = a;
  tmv::Vector<T> v2 = b;
  tmv::Vector<std::complex<T> > cv1 = ca;
  tmv::Vector<std::complex<T> > cv2 = cb;

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

  // Test I/O:
  
  const int N = 20;
  tmv::Vector<T> v(N);
  tmv::Vector<std::complex<T> > cv(N);
  for (int i=0; i<N; ++i) {
    v(i) = T(i+34);
    cv(i) = std::complex<T>(i,N-i);
  }

  std::ofstream fout("tmvtest_vector_io.dat");
  if (!fout) throw std::runtime_error(
      "Couldn't open tmvtest_vector_io.dat for output");
  fout << v << std::endl << cv << std::endl;
  fout.close();

  tmv::Vector<T> xv1(N);
  tmv::Vector<std::complex<T> > xcv1(N);
  std::ifstream fin("tmvtest_vector_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_vector_io.dat for input");
  fin >> xv1 >> xcv1;
  fin.close();
  Assert(v == xv1,"Vector I/O check #1");
  Assert(cv == xcv1,"CVector I/O check #1");

  std::auto_ptr<tmv::Vector<T> > xv2;
  std::auto_ptr<tmv::Vector<std::complex<T> > > xcv2;
  fin.open("tmvtest_vector_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_vector_io.dat for input");
  fin >> xv2 >> xcv2;
  fin.close();
  Assert(v == *xv2,"Vector I/O check #2");
  Assert(cv == *xcv2,"CVector I/O check #2");

  std::cout<<"Vector<"<<tmv::Type(T())<<"> passed all tests\n";
}

template void TestAllVector<double>();
#ifndef NOFLOAT
template void TestAllVector<float>();
#endif
#ifdef LONGDOUBLE
template void TestAllVector<long double>();
#endif
