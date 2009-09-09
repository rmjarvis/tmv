#include "TMV_VIt.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include <fstream>

#include "TMV_TestVectorArith.h"

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
  tmv::ConstVectorView<T> afcv_c = afcv;
  Assert(afcv_c == a,"CStyle View of FortransStyle Vector = CStyle Vector");
  Assert(afcv == a,"FortranStyle View of Vector == CStyle Vector");

  for (int i=0; i<N; ++i) b(i) = 5.+2*i;

  v = a+b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(8+3*i),"Adding Vectors");

  v = a-b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(-2-i),"Subtracting Vectors");

  // b(i) = 5+2i
  // a(i) = 3+i
  // a(i)*b(i) = 15+11i+2i^2
  // Sum = 15N + 11N(N-1)/2 + 2*N*(N-1)*(2N-1)/6
  T prod = 15*N + 11*N*(N-1)/2 + 2*N*(N-1)*(2*N-1)/6;
  Assert(a*b == prod,"Multiplying Vectors");

  tmv::Vector<T> c(5);
  c = v.SubVector(10,70,12);
  for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubVector");

  if (tmv::Epsilon<T>() == T(0)) return;

  for(int i=0;i<N;++i) a(i) = i+10;
  for(int i=0;i<N;++i) b(i) = -3*i+191;

  prod = 2900;
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

  if (tmv::Epsilon<T>() == T(0)) return;

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
  for(int i=0;i<N;++i) a(i) = i+10;
  tmv::Vector<T> b(N);
  for(int i=0;i<N;++i) b(i) = -3*i+191;

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

template <class T> inline void TestVectorArith()
{
  const size_t N = 100;
  tmv::Vector<T> a(N);
  for(int i=0;i<int(N);++i) a(i) = i+10;
  tmv::Vector<T> b(N);
  for(int i=0;i<int(N);++i) b(i) = -3*i+2;

  tmv::Vector<std::complex<T> > ca = a*std::complex<T>(2,-1);;
  tmv::Vector<std::complex<T> > cb = b*std::complex<T>(-5,1);

  tmv::VectorView<T> aa = a.View();
  tmv::VectorView<T,tmv::FortranStyle> af = a.View();
  tmv::Vector<T> a10(10*N);
  tmv::VectorView<T> as = a10.SubVector(0,10*N,10);
  as = a;
  tmv::VectorView<T,tmv::FortranStyle> ag = as;

  tmv::VectorView<std::complex<T> > caa = ca.View();
  tmv::VectorView<std::complex<T>,tmv::FortranStyle> caf = ca.View();
  tmv::Vector<std::complex<T> > ca10(10*N);
  tmv::VectorView<std::complex<T> > cas = ca10.SubVector(0,10*N,10);
  cas = ca;
  tmv::VectorView<std::complex<T>,tmv::FortranStyle> cag = cas;

  tmv::VectorView<T> bb = b.View();
  tmv::VectorView<T,tmv::FortranStyle> bf = b.View();
  tmv::Vector<T> b10(10*N);
  tmv::VectorView<T> bs = b10.SubVector(0,10*N,10);
  bs = b;
  tmv::VectorView<T,tmv::FortranStyle> bg = bs;

  tmv::VectorView<std::complex<T> > cbb = cb.View();
  tmv::VectorView<std::complex<T>,tmv::FortranStyle> cbf = cb.View();
  tmv::Vector<std::complex<T> > cb10(10*N);
  tmv::VectorView<std::complex<T> > cbs = cb10.SubVector(0,10*N,10);
  cbs = cb;
  tmv::VectorView<std::complex<T>,tmv::FortranStyle> cbg = cbs;

  TestVectorArith1<T>(aa,caa,"Vector C");
  TestVectorArith1<T>(af,caf,"Vector F");
  TestVectorArith2<T>(aa,caa,b,cbb,"Vector CC");
  TestVectorArith2<T>(af,caf,b,cbb,"Vector FC");
  TestVectorArith2<T>(aa,caa,bf,cbf,"Vector CF");
  TestVectorArith2<T>(af,caf,bf,cbf,"Vector FF");

  TestVectorArith1<T>(as,cas,"Vector C Step");
  TestVectorArith2<T>(as,cas,bb,cbb,"Vector C StepA");
  TestVectorArith2<T>(aa,caa,bs,cbs,"Vector C StepB");
  TestVectorArith2<T>(as,cas,bs,cbs,"Vector C StepAB");

  TestVectorArith1<T>(ag,cag,"Vector F Step");
  TestVectorArith2<T>(ag,cag,bf,cbf,"Vector F StepA");
  TestVectorArith2<T>(af,caf,bg,cbg,"Vector F StepB");
  TestVectorArith2<T>(ag,cag,bg,cbg,"Vector F StepAB");

  TestVectorArith1<T>(aa.Reverse(),caa.Reverse(),"Vector C Rev");
  TestVectorArith2<T>(aa.Reverse(),caa.Reverse(),bb,cbb,"Vector C RevA");
  TestVectorArith2<T>(aa,caa,bb.Reverse(),cbb.Reverse(),"Vector C RevB");
  TestVectorArith2<T>(aa.Reverse(),caa.Reverse(),bb.Reverse(),cbb.Reverse(),
      "Vector C RevAB");

  TestVectorArith1<T>(af.Reverse(),caf.Reverse(),"Vector F Rev");
  TestVectorArith2<T>(af.Reverse(),caf.Reverse(),bf,cbf,"Vector F RevA");
  TestVectorArith2<T>(af,caf,bf.Reverse(),cbf.Reverse(),"Vector F RevB");
  TestVectorArith2<T>(af.Reverse(),caf.Reverse(),bf.Reverse(),cbf.Reverse(),
      "Vector F RevAB");
}

template <class T> inline void TestVectorIO()
{
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
}

template <class T> void TestAllVector()
{
  TestVectorReal<T>();
  TestVectorComplex<T>();
  TestVectorArith<T>();
  TestVectorIO<T>();

  std::cout<<"Vector<"<<tmv::Type(T())<<"> passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllVector<double>();
#endif
#ifdef INST_FLOAT
template void TestAllVector<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllVector<long double>();
#endif
#ifdef INST_INT
template void TestAllVector<int>();
#endif
