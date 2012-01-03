
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

template <class T1, class T2> 
inline bool CanLDivEq(
    const tmv::BandMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return b.size() == a.colsize() && a.nhi() == int(a.rowsize())-1; }

template <class T1, class T2> 
inline bool CanRDivEq(
    const tmv::BandMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return b.size() == a.rowsize() && a.nhi() == int(a.rowsize())-1; }

template <class T1, class T2> 
inline bool CanLDivEq(
    const tmv::BandMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return b.size() == a.colsize() && a.nlo() == int(a.colsize())-1; }

template <class T1, class T2> 
inline bool CanRDivEq(
    const tmv::BandMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return b.size() == a.rowsize() && a.nlo() == int(a.colsize())-1; }

#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestBandDiv_D2(tmv::DivType dt)
{
    const int N = 10;

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1 /= T(10);
    a1.diag().addToAll(T(10)*N);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::UpperTriMatrix<T> u(a1);
    tmv::UpperTriMatrix<std::complex<T> > cu(ca1);
    tmv::LowerTriMatrix<T> l(a1);
    tmv::LowerTriMatrix<std::complex<T> > cl(ca1);
    tmv::UpperTriMatrixView<T> uv = u.view();
    tmv::UpperTriMatrixView<std::complex<T> > cuv = cu.view();
    tmv::LowerTriMatrixView<T> lv = l.view();
    tmv::LowerTriMatrixView<std::complex<T> > clv = cl.view();

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) {
            std::cout<<"Start D2 loop: i = "<<i<<"\n";
            std::cout<<"bi = "<<tmv::TMV_Text(b[i])<<"  "<<b[i]<<std::endl;
            std::cout<<"u = "<<tmv::TMV_Text(u)<<"  "<<u<<std::endl;
            std::cout<<"l = "<<tmv::TMV_Text(l)<<"  "<<l<<std::endl;
        }

        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixDivArith1(dt,uv,bi,cuv,cbi,"Band/UpperTriMatrix");
        TestMatrixDivArith1(dt,lv,bi,clv,cbi,"Band/LowerTriMatrix");
    }
}

#ifdef TEST_DOUBLE
template void TestBandDiv_D2<double>(tmv::DivType dt);
#endif
#ifdef TEST_FLOAT
template void TestBandDiv_D2<float>(tmv::DivType dt);
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandDiv_D2<long double>(tmv::DivType dt);
#endif
