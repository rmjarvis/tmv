#include "TMV_SmallVector.h"
#include "TMV_SmallVectorArith.h"
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include <fstream>

#include "TMV_TestVectorArith.h"

#define N 100
#define NN 20

template <class T> inline void TestSmallVectorReal()
{

  tmv::SmallVector<T,N> v;

  for (int i=0; i<N; ++i) v(i) = i;

  for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting SmallVector");

  tmv::VectorView<T> v2 = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v2(i) == T(2*i),
      "Reading SmallVector with stride = 2");

  tmv::SmallVectorView<T,N/2,2> v2b = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v2b(i) == T(2*i),
      "Reading SmallVector with stride = 2 (B)");

  for (int i=0; i<N/2; ++i) v2b[i] = i + 1000;
  for (int i=0; i<N/2; ++i) Assert(v(2*i) == T(i+1000),
      "Writing SmallVector with stride = 2");

  tmv::SmallVector<T,N/2> v3 = v2;
  for (int i=0; i<N/2; ++i) Assert(v3[i] == v2[i],
      "Copying SmallVector with stride = 2");

  tmv::SmallVector<T,N/2> v3b = v2b;
  for (int i=0; i<N/2; ++i) Assert(v3[i] == v2b[i],
      "Copying SmallVector with stride = 2 (B)");

  for (int i=0; i<N; ++i) v[i] = i;
  v.Swap(2,5);
  Assert(v(2) == T(5) && v(5) == T(2),"Swapping elements of SmallVector");
  v.Swap(2,5);
  Assert(v(2) == T(2) && v(5) == T(5),"Swapping elements of SmallVector");

  T sum = N*(N-1)/2;
  Assert(SumElements(v) == sum,"SmallVector SumElements(v)");

  v.ReverseSelf();
  for (int i=0; i<N; ++i) Assert(v(i) == T(N-i-1),"Reversing SmallVector");

  for (int i=0; i<N; ++i) v(i) = i+10.;
  v(23) = 10.*N;
  v(42) = 0.25;
  v(15) = -20.*N;
  size_t imax,imin;
  Assert(MaxAbsElement(v,&imax) == T(20*N),
      "MaxAbsElement of SmallVector did not return correct value");
  Assert(imax == 15,
      "MaxAbsElement of SmallVector did not return correct index");
  Assert(MinAbsElement(v,&imin) == T(0.25),
      "MinAbsElement of SmallVector did not return correct value");
  Assert(imin == 42,
      "MinAbsElement of SmallVector did not return correct index");
  Assert(MaxElement(v,&imax) == T(10*N),
      "MaxElement of SmallVector did not return correct value");
  Assert(imax == 23,
      "MaxElement of SmallVector did not return correct index");
  Assert(MinElement(v,&imin) == T(-20*N),
      "MinElement of SmallVector did not return correct value");
  Assert(imin == 15,
      "MinElement of SmallVector did not return correct index");

  tmv::SmallVector<T,N> a;
  tmv::SmallVector<T,N> b;
  for (int i=0; i<N; ++i) a(i) = T(3+i);

  b = a;
  for (int i=0; i<N; ++i) Assert(a(i) == b(i),"SmallVector1 = SmallVector2");

  Assert(a == b,"Testing Equality of SmallVectors");

  b(4) = 0;
  Assert(a != b,"SmallVector = SmallVector copied address, not values");

  tmv::SmallVector<T,N,tmv::FortranStyle> af;
  for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
  for (int i=1; i<=N; ++i) Assert(af(i) == a(i-1),
      "FortranStyle SmallVector access");
  tmv::ConstSmallVectorView<T,N,1,false,tmv::FortranStyle> afcv = af.View();
  for (int i=1; i<=N; ++i) Assert(afcv(i) == a(i-1),
      "FortranStyle SmallVector CV access");
  tmv::SmallVectorView<T,N,1,false,tmv::FortranStyle> afv = af.View();
  for (int i=1; i<=N; ++i) Assert(afv(i) == a(i-1),
      "FortranStyle SmallVector V access");
  Assert(a == af,"FortransStyle SmallVector = CStyle SmallVector");
  tmv::ConstSmallVectorView<T,N,1,false> afc = af.View();
  Assert(afc == a,
      "CStyle View of FortransStyle SmallVector = CStyle SmallVector");
  Assert(afcv == a,
      "FortranStyle View of SmallVector == CStyle SmallVector");

  for (int i=0; i<N; ++i) b(i) = 5.+2*i;

  v = a+b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(8+3*i),"Adding SmallVectors");

  v = a-b;
  for (int i=0; i<N; ++i) Assert(v(i) == T(-2-i),"Subtracting SmallVectors");

  a(0) = 1.;
  b(0) = 1.;
  for (int i=1; i<10; ++i) {
    a(i) = a(i-1)*2.;
    b(i) = b(i-1)/2.;
  }
  for (int i=10; i<N; ++i) a(i) = 0.;

  Assert(a*b == T(10),"Multiplying SmallVectors");

  tmv::SmallVector<T,5> c = v.SubVector(10,70,12);
  for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubSmallVector");

  for(int i=0;i<N;++i) a(i) = i+10.;
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  T prod = 2900;
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(std::abs(a*b - prod) < EPS*Norm(a)*Norm(b),"Inner Product");
  Assert(std::abs(Norm(a+b) - normsum) < EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Sum");
  Assert(std::abs(Norm(a-b) - normdiff) < EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Diff");

  T w1[NN] = {3.3,1.2,5.4,-1.2,4.3,-9.4,0.,-2.,4.,-11.5,-12.,14.,33.,1.,-9.3,-3.9,4.9,10.,-31.,1.e-33};

  tmv::SmallVector<T,NN> w(w1);
  tmv::SmallVector<T,NN> origw = w;
  size_t perm[NN];

  if (showacc)
    std::cout<<"unsorted w = "<<w<<std::endl;
  w.Sort(perm);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1) <= w(i),"Sort real SmallVector");
  }
  if (showacc)
    std::cout<<"sorted w = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ABS_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(std::abs(w(i-1)) <= std::abs(w(i)),"Sort real SmallVector abs");
  }
  if (showacc)
    std::cout<<"sorted w abs = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1) >= w(i),"Sort real SmallVector descend");
  }
  if (showacc)
    std::cout<<"sorted w descend = "<<w<<std::endl;

  w.Sort(0);
  w.ReversePermute(perm);
  if (showacc)
    std::cout<<"Reverse permuted w = "<<w<<std::endl;
  Assert(w==origw,"Reverse permute sorted SmallVector = orig");
  w.Sort(0);
  origw.Permute(perm);
  if (showacc)
    std::cout<<"Sort permuted w = "<<origw<<std::endl;
  Assert(w==origw,"Permute SmallVector = sorted SmallVector");

}

template <class T> inline void TestSmallVectorComplex()
{
  tmv::SmallVector<std::complex<T>,N> v;
  for (int i=0; i<N; ++i) v(i) = std::complex<T>(i,i+1234);

  for (int i=0; i<N; ++i) Assert(v(i).real() == T(i),
      "CSmallVector set");
  for (int i=0; i<N; ++i) Assert(v(i).imag() == T(i+1234),
      "CSmallVector set");

  tmv::VectorView<std::complex<T> > v1 = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v1(i) == std::complex<T>(2*i,2*i+1234),
      "CSmallVector stride=2");

  tmv::SmallVectorView<std::complex<T>,N/2,2> v1b = v.SubVector(0,N,2);
  for (int i=0; i<N/2; ++i) Assert(v1(i) == std::complex<T>(2*i,2*i+1234),
      "CSmallVector stride=2 (B)");

  for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(i,i+9876);
  for (int i=0; i<N/2; ++i) Assert(v[2*i] == std::complex<T>(i,i+9876),
      "setting CSmallVector with stride = 2");

  for (int i=0; i<N; ++i) v(i) = std::complex<T>(i,i+1234);

  v.Swap(2,5);
  Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CSmallVector");
  Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CSmallVector");
  v.Swap(2,5);

  tmv::SmallVector<std::complex<T>,N> v2 = v.Conjugate();

  for (int i=0; i<N; ++i) Assert(v2(i) == std::complex<T>(i,-i-1234),
      "Conjugate CSmallVector");
  Assert(v2 == v.Conjugate(),"Conjugate == CSmallVector");

  Assert(std::abs((v*v2).imag()) < EPS,"CSmallVector * CSmallVector");
  T norm1 = tmv::SQRT((v*v2).real());
  T norm2 = Norm(v);
  if (showacc) {
    std::cout<<"v = "<<v<<std::endl;
    std::cout<<"v2 = "<<v2<<std::endl;
    std::cout<<"v*v2 = "<<v*v2<<std::endl;
    std::cout<<"norm1 = "<<norm1<<std::endl;
    std::cout<<"norm2 = "<<norm2<<std::endl;
  }
  Assert(std::abs(norm1 - norm2) < EPS*norm1,"Norm CSmallVector");

  Assert(v2 == v.ConjugateSelf(),"ConjugateSelf CSmallVector");

  tmv::SmallVector<T,N> a;
  for(int i=0;i<N;++i) a(i) = i+10.;
  tmv::SmallVector<T,N> b;
  for(int i=0;i<N;++i) b(i) = -3.*i+191.;

  tmv::SmallVector<std::complex<T>,N> ca = a;
  Assert(Norm(ca-a) < EPS*Norm(a),"Copy real V -> complex V");
  ca *= std::complex<T>(3,4)/T(5);
  tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(3,4)/T(5);

  std::complex<T> prod = T(29)*std::complex<T>(-28,96);
  T normsum = tmv::SQRT(T(1373700));
  T normdiff = tmv::SQRT(T(1362100));
  Assert(std::abs(ca*cb - prod) < EPS*Norm(ca)*Norm(cb),"CInner Product");
  Assert(std::abs(Norm(ca+cb) - normsum) < EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Sum");
  Assert(std::abs(Norm(ca-cb) - normdiff) < EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Diff");

  T rw1[20] = {3.3,1.2,5.4,-1.2,4.3,-9.4,0.,-2.,4.,-11.5,-12.,14.,33.,1.,-9.3,-3.9,4.9,10.,-31.,1.e-33};
  T iw1[20] = {1.4,9.8,-0.2,-8.6,3.0,-4.4,3.,9.,-1.9,-11.5,11.1,-140.,-23.,11.,5.2,-3.9,4.9,99.,-71.,-0.5};

  tmv::SmallVector<T,NN> rw(rw1);
  tmv::SmallVector<T,NN> iw(iw1);
  tmv::SmallVector<std::complex<T>,NN> w;
  w.Real() = rw;
  w.Imag() = iw;
  tmv::SmallVector<std::complex<T>,NN> origw = w;
  size_t perm[NN];

  if (showacc)
    std::cout<<"unsorted w = "<<w<<std::endl;
  w.Sort(perm);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1).real() <= w(i).real(),"Sort complex SmallVector");
  }
  if (showacc)
    std::cout<<"sorted w = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ABS_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(std::abs(w(i-1)) <= std::abs(w(i)),"Sort complex SmallVector abs");
  }
  if (showacc)
    std::cout<<"sorted w abs = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND,tmv::IMAG_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(imag(w(i-1)) >= imag(w(i)),"Sort complex SmallVector descend");
  }
  if (showacc)
    std::cout<<"sorted w imag descend = "<<w<<std::endl;

  w.Sort(0,tmv::DESCEND);
  for(size_t i=1;i<NN;++i) {
    Assert(w(i-1).real() >= w(i).real(),"Sort complex SmallVector descend");
  }
  if (showacc)
    std::cout<<"sorted w imag descend = "<<w<<std::endl;

  w.Sort(0,tmv::ASCEND,tmv::ARG_COMP);
  for(size_t i=1;i<NN;++i) {
    Assert(arg(w(i-1)) <= arg(w(i)),"Sort complex SmallVector descend");
  }
  if (showacc)
    std::cout<<"sorted w arg = "<<w<<std::endl;

  w.Sort(0);
  w.ReversePermute(perm);
  Assert(w==origw,"Reverse permute sorted SmallVector = orig");
  w.Sort(0);
  origw.Permute(perm);
  Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> inline void TestSmallVectorArith()
{
  tmv::SmallVector<T,N> a;
  for(int i=0;i<N;++i) a(i) = i+10.;
  tmv::SmallVector<T,N> b;
  for(int i=0;i<N;++i) b(i) = -3.*i+2.;

  tmv::SmallVector<std::complex<T>,N> ca = a*std::complex<T>(2,-1);;
  tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(-5,1);

  tmv::SmallVectorView<T,N> aa = a.View();
  tmv::SmallVectorView<T,N,1,false,tmv::FortranStyle> af = aa;
  tmv::SmallVector<T,10*N> a10;
  tmv::SmallVectorView<T,N,10> as = a10.SubVector(0,10*N,10);
  tmv::SmallVectorView<T,N,10,false,tmv::FortranStyle> ag = as;
  
  tmv::SmallVectorView<T,N> bb = b.View();
  tmv::SmallVectorView<T,N,1,false,tmv::FortranStyle> bf = bb;
  tmv::SmallVector<T,10*N> b10;
  tmv::SmallVectorView<T,N,10> bs = b10.SubVector(0,10*N,10);
  tmv::SmallVectorView<T,N,10,false,tmv::FortranStyle> bg = bs;

  tmv::SmallVectorView<std::complex<T>,N> caa = ca.View();
  tmv::SmallVectorView<std::complex<T>,N,1,false,tmv::FortranStyle> caf = caa;
  tmv::SmallVector<std::complex<T>,10*N> ca10;
  tmv::SmallVectorView<std::complex<T>,N,10> cas = ca10.SubVector(0,10*N,10);
  tmv::SmallVectorView<std::complex<T>,N,10,false,tmv::FortranStyle> cag = cas;
  
  tmv::SmallVectorView<std::complex<T>,N> cbb = cb.View();
  tmv::SmallVectorView<std::complex<T>,N,1,false,tmv::FortranStyle> cbf = cbb;
  tmv::SmallVector<std::complex<T>,10*N> cb10;
  tmv::SmallVectorView<std::complex<T>,N,10> cbs = cb10.SubVector(0,10*N,10);
  tmv::SmallVectorView<std::complex<T>,N,10,false,tmv::FortranStyle> cbg = cbs;

#ifndef DONTMIXWITHREG
  tmv::VectorView<T> av = a.RegView();
  tmv::VectorView<std::complex<T> > cav = ca.RegView();
  tmv::VectorView<T> bv = b.RegView();
  tmv::VectorView<std::complex<T> > cbv = cb.RegView();
#endif

  TestVectorArith1<T>(aa,caa,"SmallVector CC");
  TestVectorArith1<T>(af,caf,"SmallVector FC");
  TestVectorArith2<T>(aa,caa,bb,cbb,"SmallVector CC");
  TestVectorArith2<T>(af,caf,bb,cbb,"SmallVector FC");
  TestVectorArith2<T>(aa,caa,bf,cbf,"SmallVector CF");
  TestVectorArith2<T>(af,caf,bf,cbf,"SmallVector FF");

#ifndef DONTMIXWITHREG
  TestVectorArith2<T>(av,cav,bb,cbb,"SmallVector/Vector");
  TestVectorArith2<T>(aa,caa,bv,cbv,"Vector/SmallVector");
#endif

  TestVectorArith1<T>(as,cas,"SmallVector C Step");
  TestVectorArith2<T>(as,cas,bb,cbb,"SmallVector C StepA");
  TestVectorArith2<T>(aa,caa,bs,cbs,"SmallVector C StepB");
  TestVectorArith2<T>(as,cas,bs,cbs,"SmallVector C StepAB");

  TestVectorArith1<T>(ag,cag,"SmallVector F Step");
  TestVectorArith2<T>(ag,cag,bf,cbf,"SmallVector F StepA");
  TestVectorArith2<T>(af,caf,bg,cbg,"SmallVector F StepB");
  TestVectorArith2<T>(ag,cag,bg,cbg,"SmallVector F StepAB");

  TestVectorArith1<T>(aa.Reverse(),caa.Reverse(),"SmallVector C Rev");
  TestVectorArith2<T>(aa.Reverse(),caa.Reverse(),bb,cbb,"SmallVector C RevA");
  TestVectorArith2<T>(aa,caa,bb.Reverse(),cbb.Reverse(),"SmallVector C RevB");
  TestVectorArith2<T>(aa.Reverse(),caa.Reverse(),bb.Reverse(),cbb.Reverse(),
      "SmallVector C RevAB");

  TestVectorArith1<T>(af.Reverse(),caf.Reverse(),"SmallVector F Rev");
  TestVectorArith2<T>(af.Reverse(),caf.Reverse(),bf,cbf,"SmallVector F RevA");
  TestVectorArith2<T>(af,caf,bf.Reverse(),cbf.Reverse(),"SmallVector F RevB");
  TestVectorArith2<T>(af.Reverse(),caf.Reverse(),bf.Reverse(),cbf.Reverse(),
      "SmallVector F RevAB");
}

template <class T> inline void TestSmallVectorIO()
{
  tmv::SmallVector<T,NN> v;
  tmv::SmallVector<std::complex<T>,NN> cv;
  for (int i=0; i<NN; ++i) {
    v(i) = T(i+34);
    cv(i) = std::complex<T>(i,N-i);
  }

  std::ofstream fout("tmvtest_smallvector_io.dat");
  if (!fout) throw std::runtime_error(
      "Couldn't open tmvtest_smallvector_io.dat for output");
  fout << v << std::endl << cv << std::endl;
  fout.close();

  tmv::SmallVector<T,NN> xv1;
  tmv::SmallVector<std::complex<T>,NN> xcv1;
  std::ifstream fin("tmvtest_smallvector_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_smallvector_io.dat for input");
  fin >> xv1 >> xcv1;
  fin.close();
  Assert(v == xv1,"SmallVector I/O check #1");
  Assert(cv == xcv1,"CSmallVector I/O check #1");

#ifndef DONTMIXWITHREG
  std::auto_ptr<tmv::Vector<T> > xv2;
  std::auto_ptr<tmv::Vector<std::complex<T> > > xcv2;
  fin.open("tmvtest_smallvector_io.dat");
  if (!fin) throw std::runtime_error(
      "Couldn't open tmvtest_smallvector_io.dat for input");
  fin >> xv2 >> xcv2;
  fin.close();
  Assert(v == *xv2,"SmallVector I/O check #2");
  Assert(cv == *xcv2,"CSmallVector I/O check #2");
#endif
}

template <class T> void TestAllSmallVector()
{
  TestSmallVectorReal<T>();
  TestSmallVectorComplex<T>();
  TestSmallVectorArith<T>();
  TestSmallVectorIO<T>();

  std::cout<<"SmallVector<"<<tmv::Type(T())<<"> passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllSmallVector<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallVector<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallVector<long double>();
#endif
#ifdef INST_INT
template void TestAllSmallVector<int>();
#endif
