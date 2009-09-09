//#define SHOWTESTS

#include "TMV_Test.h"
#include "TMV.h"

using tmv::Permutation;
using tmv::Matrix;
using tmv::Vector;
using tmv::RowMajor;
using tmv::ColMajor;

void TestPermutation()
{
  Permutation p0(5);
  for(size_t i=0; i<5; ++i) Assert(p0[i] == i,"Making identity Permutation 1");
  Assert(Det(p0) == 1,"Making identity Permutation 2");
  Assert(p0.IsValid(),"Making identity Permutation 3");

  size_t pp1[] = {1,2,3,4,0};
  size_t pp2[] = {4,0,1,2,3};

  Permutation p1(5,pp1,false);
  for(int i=0; i<5; ++i) Assert(p1[i] == pp1[i],
      "Making Permutation from size_t* 1");
  Assert(p1.IsValid(),"Making Permutation from size_t* 2");
  Assert(Det(p1) == 1,"Permutation DetermineDet");

  Permutation p2(5,pp2,false);
  Permutation p1b = p2.Inverse();
  Assert(p1b.IsValid(),"Inverse Permutation 1");
  for(int i=0; i<5; ++i) Assert(p1b[i] == pp1[i],"Inverse Permutation 2");
  Assert(Det(p1b) == 1,"Inverse Permutation 3");

  Assert(p1 == p1b,"Permutation == Permutation");

  Assert(p1*p2 == p0,"Permutation * Permutation");

  p1 = p1 * p2;
  Assert(p1.IsValid() && p1 == p0,"Permutation *= Permutation");

  p1 = p1 / p2;
  Assert(p1.IsValid() && p1 == p1b,"Permutation /= Permutation");

  Assert( (p1^4) == p2,"Permutation^n");

  size_t pp3[] = {1,2,0,4,5,6,7,8,3,10,11,12,13,9};

  Permutation p3(14,pp3,false,-1);
  Assert(p3.IsValid(),"p3 not valid");
  Permutation p3inv(p3);
  p3inv.InvertSelf();
  p3inv.InvertStorage();
  Assert(p3inv.IsValid(),"p3inv not valid");
  Assert(p3 * p3inv == Permutation(14),"Permutation.InvertSelf");
  Assert(p3inv * p3 == Permutation(14),"Permutation.InvertSelf(2)");

  Assert(p3*p3 == (p3^2),"Permutation ^2 with multiple cycles");

  Assert(p3.IsValid(),"p3 not valid");
  Permutation p4 = p3*p3*p3;
  p4 = Permutation(p4 * p4 * p4);
  Assert(p4.IsValid(),"p4 not valid");

  p3 ^= 9;
  Assert(p3 == p4,"Permutation ^= with multiple cycles");

  p3 = Permutation(14,pp3,false,-1);
  p3.SwapRows(3,10);
  Assert(p3.IsValid(),"SwapRows");

  Matrix<double> m3(p3);
  Assert(m3 == p3,"Matrix == Permutation");

  Permutation p5(14);
  p5.SwapCols(2,8).SwapCols(3,5).SwapCols(3,1).SwapCols(0,11).SwapCols(5,10);
  p5.SwapCols(9,12).SwapCols(3,13).SwapCols(7,1).SwapCols(6,11).SwapCols(5,0);
  Assert(p5.IsValid(),"SwapCols");
  Matrix<double> m5(p5);

  Vector<double> v(14);
  for(int i=0;i<14;++i) v[i] = 3*i + 1000;
  Vector<double> v2 = v;
  v = p5 * v;
  v2 = m5 * v2;
  Assert(v == v2,"V = P * V");
  v = v * p5;
  v2 = v2 * m5;
  Assert(v == v2,"V = V * P");
  Vector<double> v3 = v * p5;
  Vector<double> v4 = v2 * m5;
  Assert(v3 == v4,"V1 = V2 * P");
  v3 = p5 * v;
  v4 = m5 * v2;
  Assert(v3 == v4,"V1 = P * V2");

  p5.InvertStorage();
  v = p5 * v;
  v2 = m5 * v2;
  Assert(v == v2,"V = Pt * V");
  v = v * p5;
  v2 = v2 * m5;
  Assert(v == v2,"V = V * Pt");
  v3 = v * p5;
  v4 = v2 * m5;
  Assert(v3 == v4,"V1 = V2 * Pt");
  v3 = p5 * v;
  v4 = m5 * v2;
  Assert(v3 == v4,"V1 = Pt * V2");
  p5.InvertStorage();

  // Test SortPermutation:
  v = v*p5*p5;
  Permutation psort = SortPermutation(v);
  v2 = psort * v;
  for(size_t i=1;i<v2.size();++i) 
    Assert(v2[i] >= v2[i-1],"SortPermutation(REAL,ASC)");
  psort = SortPermutation(v,tmv::DESCEND);
  v2 = psort * v;
  for(size_t i=1;i<v2.size();++i) 
    Assert(v2[i] <= v2[i-1],"SortPermutation(REAL,DESC)");
  Vector<complex<double> > vc(v.size());
  for(size_t i=0;i<vc.size();++i) vc[i] = complex<double> (3.*(i-10.),5.*(-2.*i+12.));
  vc = vc * p5*p5;
  psort = SortPermutation(vc,tmv::ASCEND,tmv::REAL_COMP);
  Vector<complex<double> > vc2 = psort*vc;
  for(size_t i=1;i<vc2.size();++i) 
    Assert(real(vc2[i]) >= real(vc2[i-1]),"SortPermutation(REAL_COMP)");
  psort = SortPermutation(vc,tmv::ASCEND,tmv::IMAG_COMP);
  vc2 = psort*vc;
  for(size_t i=1;i<vc2.size();++i) 
    Assert(imag(vc2[i]) >= imag(vc2[i-1]),"SortPermutation(IMAG_COMP)");
  psort = SortPermutation(vc,tmv::ASCEND,tmv::ABS_COMP);
  vc2 = psort*vc;
  for(size_t i=1;i<vc2.size();++i) 
    Assert(abs(vc2[i]) >= abs(vc2[i-1]),"SortPermutation(ABS_COMP,ASC)");
  psort = SortPermutation(vc,tmv::ASCEND,tmv::ARG_COMP);
  vc2 = psort*vc;
  for(size_t i=1;i<vc2.size();++i) 
    Assert(arg(vc2[i]) >= arg(vc2[i-1]),"SortPermutation(ARG_COMP,ASC)");

  Matrix<double,RowMajor> w(14,5);
  for(int i=0;i<14;++i) for(int j=0;j<5;++j) w(i,j) = 2*i + 54*j + 1000;
  Matrix<double,ColMajor> wt = w.Transpose();
  Matrix<double,RowMajor> w2 = w;
  Matrix<double,ColMajor> w2t = w2.Transpose();
  w = p5 * w;
  w2 = m5 * w2;
  Assert(w == w2,"M = P * M");
  wt = wt * p5;
  w2t = w2t * m5;
  Assert(wt == w2t,"M = M * P");
  Matrix<double> w3 = wt * p5;
  Matrix<double> w4 = w2t * m5;
  Assert(w3 == w4,"M1 = M2 * P");
  Matrix<double> w5 = p5 * w;
  Matrix<double> w6 = m5 * w2;
  Assert(w5 == w6,"M1 = P * M2");

  p5.InvertStorage();
  w = p5 * w;
  w2 = m5 * w2;
  Assert(w == w2,"M = Pt * M");
  wt = wt * p5;
  w2t = w2t * m5;
  Assert(w == w2,"M = M * Pt");
  w3 = wt * p5;
  w4 = w2t * m5;
  Assert(w3 == w4,"M1 = M2 * Pt");
  w5 = p5 * w;
  w6 = m5 * w2;
  Assert(w5 == w6,"M1 = Pt * M2");

  Permutation p6 = p5 * p3;
  Matrix<double> m6 = m5 * m3;
  Assert(p6 == m6,"P * P");
  p5.InvertStorage();
  p6 = p5 * p3;
  Assert(p6 == m6,"Pt * P");
  p3.InvertStorage();
  p6 = p5 * p3;
  Assert(p6 == m6,"Pt * Pt");
  p5.InvertStorage();
  p6 = p5 * p3;
  Assert(p6 == m6,"P * Pt");

  p6 = p6 * p5;
  m6 = m6 * m5;
  Assert(p6 == m6,"P1 = P1 * P2");
  p5.InvertStorage();
  p6 = p6 * p5;
  m6 = m6 * m5;
  Assert(p6 == m6,"P1 = P1 * P2t");
  p6 = p5 * p6;
  m6 = m5 * m6;
  Assert(p6 == m6,"P2 = P1t * P2");
  p5.InvertStorage();
  p6 = p5 * p6;
  m6 = m5 * m6;
  Assert(p6 == m6,"P2 = P1 * P2");
  cout<<"Permutation passed all tests\n";
}

