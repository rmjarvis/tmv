#define START 0

#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestBandArith.h"

template <class M1, class M2> 
inline bool CanAddEq(
    const tmv::BaseMatrix_Band_Mutable<M1>& a, 
    const tmv::BaseMatrix_Tri<M2>& b)
{ 
    return a.colsize() == b.colsize() && a.rowsize() == b.rowsize() && 
        a.nhi() >= b.nhi() && a.nlo() >= b.nlo();
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Tri<M1>& a, const tmv::BaseMatrix_Band<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == c.colsize()-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == c.rowsize()-1);
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Tri<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == c.colsize()-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == c.rowsize()-1);
}

template <class M1, class M2, class M3> 
inline bool CanMult(
    const tmv::BaseMatrix_Tri<M1>& a, const tmv::BaseMatrix_Tri<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == b.colsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() &&
        (c.nlo() >= a.nlo() + b.nlo() || c.nlo() == c.colsize()-1) &&
        (c.nhi() >= a.nhi() + b.nhi() || c.nhi() == c.rowsize()-1);
}

template <class M1, class M2, class M3> 
inline bool CanElemMult(
    const tmv::BaseMatrix_Band<M1>& a, const tmv::BaseMatrix_Tri<M2>& b,
    const tmv::BaseMatrix_Band_Mutable<M3>& c)
{ 
    return a.rowsize() == c.rowsize() && a.colsize() == c.colsize() &&
        b.rowsize() == c.rowsize() && b.colsize() == c.colsize() &&
        (c.nlo() >= a.nlo() || M2::_upper) && 
        (c.nhi() >= a.nhi() || !M2::_upper);
}

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestBandMatrixArith_D1()
{
    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    const int N = b[0].rowsize();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j)
        ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor> u1(a1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::RowMajor> cu1(ca1);
    tmv::UpperTriMatrixView<T> u1v = u1.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu1v = cu1.view();
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> u1x = u1v;
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> cu1x = cu1v;

#if (XTEST & 2)
    tmv::UpperTriMatrix<T,tmv::UnitDiag|tmv::RowMajor> u2(a1);
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag|tmv::ColMajor> u3(a1);
    tmv::UpperTriMatrix<T,tmv::UnitDiag|tmv::ColMajor> u4(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor> l1(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag|tmv::RowMajor> l2(a1);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag|tmv::ColMajor> l3(a1);
    tmv::LowerTriMatrix<T,tmv::UnitDiag|tmv::ColMajor> l4(a1);

    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::RowMajor> cu2(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::ColMajor> cu3(ca1);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::ColMajor> cu4(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::RowMajor> cl1(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::RowMajor> cl2(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag|tmv::ColMajor> cl3(ca1);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag|tmv::ColMajor> cl4(ca1);

    tmv::UpperTriMatrixView<T> u2v = u2.view();
    tmv::UpperTriMatrixView<T> u3v = u3.view();
    tmv::UpperTriMatrixView<T> u4v = u4.view();
    tmv::LowerTriMatrixView<T> l1v = l1.view();
    tmv::LowerTriMatrixView<T> l2v = l2.view();
    tmv::LowerTriMatrixView<T> l3v = l3.view();
    tmv::LowerTriMatrixView<T> l4v = l4.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu2v = cu2.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu3v = cu3.view();
    tmv::UpperTriMatrixView<std::complex<T> > cu4v = cu4.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl1v = cl1.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl2v = cl2.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl3v = cl3.view();
    tmv::LowerTriMatrixView<std::complex<T> > cl4v = cl4.view();
#endif

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) {
            std::cerr<<"Start loop "<<i<<std::endl;
            std::cerr<<"bi = "<<b[i]<<std::endl;
        }
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];

        TestMatrixArith4(bi,cbi,u1v,cu1v,"Band/UpperTri");
        TestMatrixArith5(bi,cbi,u1v,cu1v,"Band/UpperTri");
        TestMatrixArith6x(bi,cbi,u1v,cu1v,"Band/UpperTri");
#if (XTEST & 2)
        TestMatrixArith4(bi,cbi,l1v,cl1v,"Band/LowerTri");
        TestMatrixArith5(bi,cbi,l1v,cl1v,"Band/LowerTri");
        TestMatrixArith6x(bi,cbi,l1v,cl1v,"Band/LowerTri");
        TestMatrixArith4(bi,cbi,u2v,cu2v,"Band/UpperTri");
        TestMatrixArith5(bi,cbi,u2v,cu2v,"Band/UpperTri");
        TestMatrixArith6x(bi,cbi,u2v,cu2v,"Band/UpperTri");
        TestMatrixArith4(bi,cbi,l2v,cl2v,"Band/LowerTri");
        TestMatrixArith5(bi,cbi,l2v,cl2v,"Band/LowerTri");
        TestMatrixArith6x(bi,cbi,l2v,cl2v,"Band/LowerTri");
        TestMatrixArith4(bi,cbi,u3v,cu3v,"Band/UpperTri");
        TestMatrixArith5(bi,cbi,u3v,cu3v,"Band/UpperTri");
        TestMatrixArith6x(bi,cbi,u3v,cu3v,"Band/UpperTri");
        TestMatrixArith4(bi,cbi,l3v,cl3v,"Band/LowerTri");
        TestMatrixArith5(bi,cbi,l3v,cl3v,"Band/LowerTri");
        TestMatrixArith6x(bi,cbi,l3v,cl3v,"Band/LowerTri");
        TestMatrixArith4(bi,cbi,u4v,cu4v,"Band/UpperTri");
        TestMatrixArith5(bi,cbi,u4v,cu4v,"Band/UpperTri");
        TestMatrixArith6x(bi,cbi,u4v,cu4v,"Band/UpperTri");
        TestMatrixArith4(bi,cbi,l4v,cl4v,"Band/LowerTri");
        TestMatrixArith5(bi,cbi,l4v,cl4v,"Band/LowerTri");
        TestMatrixArith6x(bi,cbi,l4v,cl4v,"Band/LowerTri");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestBandMatrixArith_D1<double>();
#endif
#ifdef TEST_FLOAT
template void TestBandMatrixArith_D1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestBandMatrixArith_D1<long double>();
#endif
#ifdef TEST_INT
template void TestBandMatrixArith_D1<int>();
#endif
