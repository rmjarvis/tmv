#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOADDEQ
#define NOMULTEQ
#define SYMOPROD
#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_B1()
{

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);

    const int N = s[0].size();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::Vector<T> v = a1.col(0);
    tmv::VectorView<T> vv = v.view();
    tmv::Vector<std::complex<T> > cv = ca1.col(0);
    tmv::VectorView<std::complex<T> > cvv = cv.view();

    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

#if (XTEST & 2)
    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
        std::complex<T>(1-3*i+6*j,-4+2*j);

    tmv::Vector<T> v5(5*N);
    tmv::VectorView<T> vs = v5.subVector(0,5*N,5);
    vs = v;
    tmv::Vector<std::complex<T> > cv5(5*N);
    tmv::VectorView<std::complex<T> > cvs = cv5.subVector(0,5*N,5);
    cvs = cv;

    tmv::Matrix<T,tmv::RowMajor> a3 = a2.rowRange(0,N);
    tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4);
    tmv::Matrix<T,tmv::RowMajor> a4 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca4 = a4;

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();

    tmv::Matrix<T> a1x = a1;
    tmv::Matrix<T> a3x = a3;
    tmv::Matrix<T> a4x = a4;
    tmv::Matrix<std::complex<T> > ca1x = ca1;
    tmv::Matrix<std::complex<T> > ca3x = ca3;
    tmv::Matrix<std::complex<T> > ca4x = ca4;
#endif

    symoprod = true;

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }

        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixArith4(si,csi,a1v,ca1v,"Sym/SquareM");
        TestMatrixArith5(si,csi,a1v,ca1v,"Sym/SquareM");
        TestMatrixArith6x(si,csi,a1v,ca1v,"Sym/SquareM");
        if (csi.isherm())
            TestMatrixArith7(si,csi,vv,cvv,vv,cvv.conjugate(),"Herm");
        else
            TestMatrixArith7(si,csi,vv,cvv,vv,cvv,"Sym");
#if (XTEST & 2)
        TestMatrixArith4(si,csi,a3v,ca3v,"Sym/NonSquareM");
        TestMatrixArith5(si,csi,a3v,ca3v,"Sym/NonSquareM");
        TestMatrixArith6x(si,csi,a3v,ca3v,"Sym/NonSquareM");
        TestMatrixArith4(si,csi,a4v,ca4v,"Sym/DegenerateM");
        TestMatrixArith5(si,csi,a4v,ca4v,"Sym/DegenerateM");
        TestMatrixArith6x(si,csi,a4v,ca4v,"Sym/DegenerateM");
        if (csi.isherm())
            TestMatrixArith7(si,csi,vs,cvs,vs,cvs.conjugate(),"Herm");
        else
            TestMatrixArith7(si,csi,vs,cvs,vs,cvs,"Sym");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_B1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_B1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_B1<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_B1<int>();
#endif
