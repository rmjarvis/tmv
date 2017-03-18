#define START 0

#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

template <class M1, class M2> 
inline bool CanAddEq(
    const tmv::BaseMatrix_Diag_Mutable<M1>& a,
    const tmv::BaseMatrix_Band<M2>& b)
{ 
    return a.colsize() == b.colsize() && a.rowsize() == b.rowsize() &&
        b.nhi() == 0 && b.nlo() == 0;
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Diag<M1>& a, const tmv::BaseMatrix_Band<M2>& b,
    const tmv::BaseMatrix_Diag_Mutable<M3>& c)
{ 
    return b.colsize() == a.size() && b.rowsize() == a.size() && 
        c.size() == a.size() && b.nlo() == 0 && b.nhi() == 0;
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Diag<M2>& b,
    const tmv::BaseMatrix_Diag_Mutable<M3>& c)
{ 
    return a.colsize() == b.size() && a.rowsize() == b.size() &&
        c.size() == a.size() && a.nlo() == 0 && a.nhi() == 0;
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Band<M2>& b,
    const tmv::BaseMatrix_Diag_Mutable<M3>& c)
{ 
    return a.colsize() == c.size() && b.rowsize() == c.size() && 
        a.rowsize() == b.colsize() && 
        a.nlo() == 0 && a.nhi() == 0 && b.nlo() == 0 && b.nhi() == 0;
}

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_C2()
{
#if (XTEST & 2)
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    const int N = b[0].rowsize();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
        ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);

    tmv::DiagMatrix<T> d1(a1);
    tmv::DiagMatrix<std::complex<T> > cd1(ca1);
    tmv::DiagMatrixView<T> d1v = d1.view();
    tmv::DiagMatrixView<std::complex<T> > cd1v = cd1.view();

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"bi = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixArith4(d1v,cd1v,bi,cbi,"Diag/Band");
        TestMatrixArith5(d1v,cd1v,bi,cbi,"Diag/Band");
        TestMatrixArith6x(d1v,cd1v,bi,cbi,"Diag/Band");
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_C2<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_C2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_C2<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_C2<int>();
#endif
