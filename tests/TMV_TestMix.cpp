
#include <fstream>
#include "TMV_Test.h"
#include "TMV_Test_3.h"

bool XXDEBUG1 = false;
bool XXDEBUG2 = false;
bool XXDEBUG3 = false;
bool XXDEBUG4 = false;
bool XXDEBUG5 = false;
bool XXDEBUG6 = false;
bool XXDEBUG7 = false;
bool XXDEBUG8 = false;
bool XXDEBUG9 = false;

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool showstartdone = false;
bool donorm2 = true;
bool symoprod = false;
bool dontthrow = false;
std::string lastsuccess = "";

#include "TMV_Small.h"
#define NONSQUARE
#define INORDER
#define NOELEMMULT
#include "TMV_TestVectorArith.h"
#include "TMV_TestMatrixArith.h"

#define M 8
#define N 5

int main() try 
{
    //showacc=true;
    //showdiv=true;
    //showtests=true;
    //showstartdone=true;

    tmv::SmallVector<double,N> svn;
    for(int i=0;i<N;++i) svn(i) = i+10.;
    tmv::SmallVector<double,M> svm;
    for(int i=0;i<M;++i) svm(i) = i+10.;

    tmv::SmallVector<std::complex<double>,N> csvn =
        svn*std::complex<double>(2,-1);
    tmv::SmallVector<std::complex<double>,M> csvm =
        svm*std::complex<double>(2,-1);

    tmv::Vector<double> vn = svn;
    tmv::Vector<double> vm = svm;
    tmv::Vector<std::complex<double> > cvn = csvn;
    tmv::Vector<std::complex<double> > cvm = csvm;

    tmv::Vector<double> vx(M*N);
    tmv::Vector<double> vy(M*N);
    tmv::VectorView<double> vvn = vx.subVector(0,M*N,M);
    tmv::VectorView<double> vvm = vy.subVector(0,M*N,N);
    tmv::Vector<std::complex<double> > cvx(M*N);
    tmv::Vector<std::complex<double> > cvy(M*N);
    tmv::VectorView<std::complex<double> > cvvn = cvx.subVector(0,M*N,M);
    tmv::VectorView<std::complex<double> > cvvm = cvy.subVector(0,M*N,N);
    vvn = svn;
    cvvn = csvn;
    vvm = svm;
    cvvm = csvm;

#if 1
    TestVectorArith2(svn,csvn,svn,csvn,"SmallVector/SmallVector");
    TestVectorArith2(vn,cvn,svn,csvn,"SmallVector/Vector");
    TestVectorArith2(svn,csvn,vn,cvn,"Vector/SmallVector");
    TestVectorArith2(vn,cvvn,svn,csvn,"SmallVector/VectorView");
    TestVectorArith2(svn,csvn,vvn,cvvn,"VectorView/SmallVector");
#endif

    tmv::SmallMatrix<double,M,N> smmn;
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        smmn(i,j) = 2.9+4.3*i-5.1*j;
    }
    tmv::SmallMatrix<std::complex<double>,M,N> csmmn =
        std::complex<double>(3,2)*smmn;
    tmv::SmallMatrix<double,N,M> smnm = smmn.transpose();
    tmv::SmallMatrix<double,N,N> smnn = smmn.rowRange(0,N);
    tmv::SmallMatrix<std::complex<double>,N,M> csmnm = csmmn.transpose();
    tmv::SmallMatrix<std::complex<double>,N,N> csmnn = csmmn.rowRange(0,N);

    tmv::Matrix<double> mmn = smmn;
    tmv::Matrix<double> mnm = smnm;
    tmv::Matrix<double> mnn = smnn;
    tmv::Matrix<std::complex<double> > cmmn = csmmn;
    tmv::Matrix<std::complex<double> > cmnm = csmnm;
    tmv::Matrix<std::complex<double> > cmnn = csmnn;

    tmv::Matrix<double> mx(M*N,M*N);
    tmv::Matrix<double> my(M*N,M*N);
    tmv::Matrix<double> mz(M*N,M*N);
    tmv::MatrixView<double> mvmn = mx.subMatrix(0,M*N,0,M*N,N,M);
    tmv::MatrixView<double> mvnm = my.subMatrix(0,M*N,0,M*N,M,N);
    tmv::MatrixView<double> mvnn = mz.subMatrix(0,M*N,0,M*N,M,M);
    tmv::Matrix<std::complex<double> > cmx(M*N,M*N);
    tmv::Matrix<std::complex<double> > cmy(M*N,M*N);
    tmv::Matrix<std::complex<double> > cmz(M*N,M*N);
    tmv::MatrixView<std::complex<double> > cmvmn =
        cmx.subMatrix(0,M*N,0,M*N,N,M);
    tmv::MatrixView<std::complex<double> > cmvnm =
        cmy.subMatrix(0,M*N,0,M*N,M,N);
    tmv::MatrixView<std::complex<double> > cmvnn =
        cmz.subMatrix(0,M*N,0,M*N,M,M);
    mvmn = mmn;
    mvnm = mnm;
    mvnn = mnn;
    cmvmn = cmmn;
    cmvnm = cmnm;
    cmvnn = cmnn;

#if 1
    TestMatrixArith2a(smmn,csmmn,svn,csvn,svm,csvm,
                      "MN SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith2a(smmn,csmmn,svn,csvn,vm,cvm,
                      "MN SmallMatrix,SmallVector,Vector");
    TestMatrixArith2a(smmn,csmmn,vn,cvn,svm,csvm,
                      "MN SmallMatrix,Vector,SmallVector");
#if (XTEST & 2)
    TestMatrixArith2a(smmn,csmmn,vn,cvn,vm,cvm,
                      "MN SmallMatrix,Vector,Vector");
    TestMatrixArith2a(smmn,csmmn,svn,csvn,vvm,cvvm,
                      "MN SmallMatrix,SmallVector,VectorView");
    TestMatrixArith2a(smmn,csmmn,vvn,cvvn,svm,csvm,
                      "MN SmallMatrix,VectorView,SmallVector");
    TestMatrixArith2a(smmn,csmmn,vvn,cvvn,vvm,cvvm,
                      "MN SmallMatrix,VectorView,VectorView");
    TestMatrixArith2a(mmn,cmmn,svn,csvn,vm,cvm,
                      "MN Matrix,SmallVector,Vector");
    TestMatrixArith2a(mmn,cmmn,vn,cvn,svm,csvm,
                      "MN Matrix,Vector,SmallVector");
    TestMatrixArith2a(mmn,cmmn,svn,csvn,vvm,cvvm,
                      "MN Matrix,SmallVector,VectorView");
    TestMatrixArith2a(mmn,cmmn,vvn,cvvn,svm,csvm,
                      "MN Matrix,VectorView,SmallVector");
    TestMatrixArith2a(mvmn,cmvmn,svn,csvn,vm,cvm,
                      "MN MatrixView,SmallVector,Vector");
    TestMatrixArith2a(mvmn,cmvmn,vn,cvn,svm,csvm,
                      "MN MatrixView,Vector,SmallVector");
    TestMatrixArith2a(mvmn,cmvmn,svn,csvn,vvm,cvvm,
                      "MN MatrixView,SmallVector,VectorView");
    TestMatrixArith2a(mvmn,cmvmn,vvn,cvvn,svm,csvm,
                      "MN MatrixView,VectorView,SmallVector");

    TestMatrixArith2a(smnm,csmnm,svm,csvm,svn,csvn,
                      "NM SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith2a(smnm,csmnm,svm,csvm,vn,cvn,
                      "NM SmallMatrix,SmallVector,Vector");
    TestMatrixArith2a(smnm,csmnm,vm,cvm,svn,csvn,
                      "NM SmallMatrix,Vector,SmallVector");
    TestMatrixArith2a(smnm,csmnm,vm,cvm,vn,cvn,
                      "NM SmallMatrix,Vector,Vector");
    TestMatrixArith2a(smnm,csmnm,svm,csvm,vvn,cvvn,
                      "NM SmallMatrix,SmallVector,VectorView");
    TestMatrixArith2a(smnm,csmnm,vvm,cvvm,svn,csvn,
                      "NM SmallMatrix,VectorView,SmallVector");
    TestMatrixArith2a(smnm,csmnm,vvm,cvvm,vvn,cvvn,
                      "NM SmallMatrix,VectorView,VectorView");
    TestMatrixArith2a(mnm,cmnm,svm,csvm,vn,cvn,
                      "NM Matrix,SmallVector,Vector");
    TestMatrixArith2a(mnm,cmnm,vm,cvm,svn,csvn,
                      "NM Matrix,Vector,SmallVector");
    TestMatrixArith2a(mnm,cmnm,svm,csvm,vvn,cvvn,
                      "NM Matrix,SmallVector,VectorView");
    TestMatrixArith2a(mnm,cmnm,vvm,cvvm,svn,csvn,
                      "NM Matrix,VectorView,SmallVector");
    TestMatrixArith2a(mvnm,cmvnm,svm,csvm,vn,cvn,
                      "NM MatrixView,SmallVector,Vector");
    TestMatrixArith2a(mvnm,cmvnm,vm,cvm,svn,csvn,
                      "NM MatrixView,Vector,SmallVector");
    TestMatrixArith2a(mvnm,cmvnm,svm,csvm,vvn,cvvn,
                      "NM MatrixView,SmallVector,VectorView");
    TestMatrixArith2a(mvnm,cmvnm,vvm,cvvm,svn,csvn,
                      "NM MatrixView,VectorView,SmallVector");
#endif
#endif

#if 1
    TestMatrixArith2b(smmn,csmmn,svn,csvn,svm,csvm,
                      "MN SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith2b(smmn,csmmn,svn,csvn,vm,cvm,
                      "MN SmallMatrix,SmallVector,Vector");
    TestMatrixArith2b(smmn,csmmn,vn,cvn,svm,csvm,
                      "MN SmallMatrix,Vector,SmallVector");
#if (XTEST & 2)
    TestMatrixArith2b(smmn,csmmn,vn,cvn,vm,cvm,
                      "MN SmallMatrix,Vector,Vector");
    TestMatrixArith2b(smmn,csmmn,svn,csvn,vvm,cvvm,
                      "MN SmallMatrix,SmallVector,VectorView");
    TestMatrixArith2b(smmn,csmmn,vvn,cvvn,svm,csvm,
                      "MN SmallMatrix,VectorView,SmallVector");
    TestMatrixArith2b(smmn,csmmn,vvn,cvvn,vvm,cvvm,
                      "MN SmallMatrix,VectorView,VectorView");
    TestMatrixArith2b(mmn,cmmn,svn,csvn,vm,cvm,
                      "MN Matrix,SmallVector,Vector");
    TestMatrixArith2b(mmn,cmmn,vn,cvn,svm,csvm,
                      "MN Matrix,Vector,SmallVector");
    TestMatrixArith2b(mmn,cmmn,svn,csvn,vvm,cvvm,
                      "MN Matrix,SmallVector,VectorView");
    TestMatrixArith2b(mmn,cmmn,vvn,cvvn,svm,csvm,
                      "MN Matrix,VectorView,SmallVector");
    TestMatrixArith2b(mvmn,cmvmn,svn,csvn,vm,cvm,
                      "MN MatrixView,SmallVector,Vector");
    TestMatrixArith2b(mvmn,cmvmn,vn,cvn,svm,csvm,
                      "MN MatrixView,Vector,SmallVector");
    TestMatrixArith2b(mvmn,cmvmn,svn,csvn,vvm,cvvm,
                      "MN MatrixView,SmallVector,VectorView");
    TestMatrixArith2b(mvmn,cmvmn,vvn,cvvn,svm,csvm,
                      "MN MatrixView,VectorView,SmallVector");

    TestMatrixArith2b(smnm,csmnm,svm,csvm,svn,csvn,
                      "NM SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith2b(smnm,csmnm,svm,csvm,vn,cvn,
                      "NM SmallMatrix,SmallVector,Vector");
    TestMatrixArith2b(smnm,csmnm,vm,cvm,svn,csvn,
                      "NM SmallMatrix,Vector,SmallVector");
    TestMatrixArith2b(smnm,csmnm,vm,cvm,vn,cvn,
                      "NM SmallMatrix,Vector,Vector");
    TestMatrixArith2b(smnm,csmnm,svm,csvm,vvn,cvvn,
                      "NM SmallMatrix,SmallVector,VectorView");
    TestMatrixArith2b(smnm,csmnm,vvm,cvvm,svn,csvn,
                      "NM SmallMatrix,VectorView,SmallVector");
    TestMatrixArith2b(smnm,csmnm,vvm,cvvm,vvn,cvvn,
                      "NM SmallMatrix,VectorView,VectorView");
    TestMatrixArith2b(mnm,cmnm,svm,csvm,vn,cvn,
                      "NM Matrix,SmallVector,Vector");
    TestMatrixArith2b(mnm,cmnm,vm,cvm,svn,csvn,
                      "NM Matrix,Vector,SmallVector");
    TestMatrixArith2b(mnm,cmnm,svm,csvm,vvn,cvvn,
                      "NM Matrix,SmallVector,VectorView");
    TestMatrixArith2b(mnm,cmnm,vvm,cvvm,svn,csvn,
                      "NM Matrix,VectorView,SmallVector");
    TestMatrixArith2b(mvnm,cmvnm,svm,csvm,vn,cvn,
                      "NM MatrixView,SmallVector,Vector");
    TestMatrixArith2b(mvnm,cmvnm,vm,cvm,svn,csvn,
                      "NM MatrixView,Vector,SmallVector");
    TestMatrixArith2b(mvnm,cmvnm,svm,csvm,vvn,cvvn,
                      "NM MatrixView,SmallVector,VectorView");
    TestMatrixArith2b(mvnm,cmvnm,vvm,cvvm,svn,csvn,
                      "NM MatrixView,VectorView,SmallVector");
#endif
#endif

#if 1
    TestMatrixArith3a(smmn,csmmn,svn,csvn,svm,csvm,
                      "MN SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith3a(smmn,csmmn,svn,csvn,vm,cvm,
                      "MN SmallMatrix,SmallVector,Vector");
    TestMatrixArith3a(smmn,csmmn,vn,cvn,svm,csvm,
                      "MN SmallMatrix,Vector,SmallVector");
    TestMatrixArith3b(smmn,csmmn,svn,csvn,svm,csvm,
                      "MN SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith3b(smmn,csmmn,svn,csvn,vm,cvm,
                      "MN SmallMatrix,SmallVector,Vector");
    TestMatrixArith3b(smmn,csmmn,vn,cvn,svm,csvm,
                      "MN SmallMatrix,Vector,SmallVector");
#if (XTEST & 2)
    TestMatrixArith3a(smmn,csmmn,vn,cvn,vm,cvm,
                      "MN SmallMatrix,Vector,Vector");
    TestMatrixArith3a(smmn,csmmn,svn,csvn,vvm,cvvm,
                      "MN SmallMatrix,SmallVector,VectorView");
    TestMatrixArith3a(smmn,csmmn,vvn,cvvn,svm,csvm,
                      "MN SmallMatrix,VectorView,SmallVector");
    TestMatrixArith3a(smmn,csmmn,vvn,cvvn,vvm,cvvm,
                      "MN SmallMatrix,VectorView,VectorView");
    TestMatrixArith3a(mmn,cmmn,svn,csvn,vm,cvm,
                      "MN Matrix,SmallVector,Vector");
    TestMatrixArith3a(mmn,cmmn,vn,cvn,svm,csvm,
                      "MN Matrix,Vector,SmallVector");
    TestMatrixArith3a(mmn,cmmn,svn,csvn,vvm,cvvm,
                      "MN Matrix,SmallVector,VectorView");
    TestMatrixArith3a(mmn,cmmn,vvn,cvvn,svm,csvm,
                      "MN Matrix,VectorView,SmallVector");
    TestMatrixArith3a(mvmn,cmvmn,svn,csvn,vm,cvm,
                      "MN MatrixView,SmallVector,Vector");
    TestMatrixArith3a(mvmn,cmvmn,vn,cvn,svm,csvm,
                      "MN MatrixView,Vector,SmallVector");
    TestMatrixArith3a(mvmn,cmvmn,svn,csvn,vvm,cvvm,
                      "MN MatrixView,SmallVector,VectorView");
    TestMatrixArith3a(mvmn,cmvmn,vvn,cvvn,svm,csvm,
                      "MN MatrixView,VectorView,SmallVector");

    TestMatrixArith3a(smnm,csmnm,svm,csvm,svn,csvn,
                      "NM SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith3a(smnm,csmnm,svm,csvm,vn,cvn,
                      "NM SmallMatrix,SmallVector,Vector");
    TestMatrixArith3a(smnm,csmnm,vm,cvm,svn,csvn,
                      "NM SmallMatrix,Vector,SmallVector");
    TestMatrixArith3a(smnm,csmnm,vm,cvm,vn,cvn,
                      "NM SmallMatrix,Vector,Vector");
    TestMatrixArith3a(smnm,csmnm,svm,csvm,vvn,cvvn,
                      "NM SmallMatrix,SmallVector,VectorView");
    TestMatrixArith3a(smnm,csmnm,vvm,cvvm,svn,csvn,
                      "NM SmallMatrix,VectorView,SmallVector");
    TestMatrixArith3a(smnm,csmnm,vvm,cvvm,vvn,cvvn,
                      "NM SmallMatrix,VectorView,VectorView");
    TestMatrixArith3a(mnm,cmnm,svm,csvm,vn,cvn,
                      "NM Matrix,SmallVector,Vector");
    TestMatrixArith3a(mnm,cmnm,vm,cvm,svn,csvn,
                      "NM Matrix,Vector,SmallVector");
    TestMatrixArith3a(mnm,cmnm,svm,csvm,vvn,cvvn,
                      "NM Matrix,SmallVector,VectorView");
    TestMatrixArith3a(mnm,cmnm,vvm,cvvm,svn,csvn,
                      "NM Matrix,VectorView,SmallVector");
    TestMatrixArith3a(mvnm,cmvnm,svm,csvm,vn,cvn,
                      "NM MatrixView,SmallVector,Vector");
    TestMatrixArith3a(mvnm,cmvnm,vm,cvm,svn,csvn,
                      "NM MatrixView,Vector,SmallVector");
    TestMatrixArith3a(mvnm,cmvnm,svm,csvm,vvn,cvvn,
                      "NM MatrixView,SmallVector,VectorView");
    TestMatrixArith3a(mvnm,cmvnm,vvm,cvvm,svn,csvn,
                      "NM MatrixView,VectorView,SmallVector");

    TestMatrixArith3b(smmn,csmmn,vn,cvn,vm,cvm,
                      "MN SmallMatrix,Vector,Vector");
    TestMatrixArith3b(smmn,csmmn,svn,csvn,vvm,cvvm,
                      "MN SmallMatrix,SmallVector,VectorView");
    TestMatrixArith3b(smmn,csmmn,vvn,cvvn,svm,csvm,
                      "MN SmallMatrix,VectorView,SmallVector");
    TestMatrixArith3b(smmn,csmmn,vvn,cvvn,vvm,cvvm,
                      "MN SmallMatrix,VectorView,VectorView");
    TestMatrixArith3b(mmn,cmmn,svn,csvn,vm,cvm,
                      "MN Matrix,SmallVector,Vector");
    TestMatrixArith3b(mmn,cmmn,vn,cvn,svm,csvm,
                      "MN Matrix,Vector,SmallVector");
    TestMatrixArith3b(mmn,cmmn,svn,csvn,vvm,cvvm,
                      "MN Matrix,SmallVector,VectorView");
    TestMatrixArith3b(mmn,cmmn,vvn,cvvn,svm,csvm,
                      "MN Matrix,VectorView,SmallVector");
    TestMatrixArith3b(mvmn,cmvmn,svn,csvn,vm,cvm,
                      "MN MatrixView,SmallVector,Vector");
    TestMatrixArith3b(mvmn,cmvmn,vn,cvn,svm,csvm,
                      "MN MatrixView,Vector,SmallVector");
    TestMatrixArith3b(mvmn,cmvmn,svn,csvn,vvm,cvvm,
                      "MN MatrixView,SmallVector,VectorView");
    TestMatrixArith3b(mvmn,cmvmn,vvn,cvvn,svm,csvm,
                      "MN MatrixView,VectorView,SmallVector");

    TestMatrixArith3b(smnm,csmnm,svm,csvm,svn,csvn,
                      "NM SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith3b(smnm,csmnm,svm,csvm,vn,cvn,
                      "NM SmallMatrix,SmallVector,Vector");
    TestMatrixArith3b(smnm,csmnm,vm,cvm,svn,csvn,
                      "NM SmallMatrix,Vector,SmallVector");
    TestMatrixArith3b(smnm,csmnm,vm,cvm,vn,cvn,
                      "NM SmallMatrix,Vector,Vector");
    TestMatrixArith3b(smnm,csmnm,svm,csvm,vvn,cvvn,
                      "NM SmallMatrix,SmallVector,VectorView");
    TestMatrixArith3b(smnm,csmnm,vvm,cvvm,svn,csvn,
                      "NM SmallMatrix,VectorView,SmallVector");
    TestMatrixArith3b(smnm,csmnm,vvm,cvvm,vvn,cvvn,
                      "NM SmallMatrix,VectorView,VectorView");
    TestMatrixArith3b(mnm,cmnm,svm,csvm,vn,cvn,
                      "NM Matrix,SmallVector,Vector");
    TestMatrixArith3b(mnm,cmnm,vm,cvm,svn,csvn,
                      "NM Matrix,Vector,SmallVector");
    TestMatrixArith3b(mnm,cmnm,svm,csvm,vvn,cvvn,
                      "NM Matrix,SmallVector,VectorView");
    TestMatrixArith3b(mnm,cmnm,vvm,cvvm,svn,csvn,
                      "NM Matrix,VectorView,SmallVector");
    TestMatrixArith3b(mvnm,cmvnm,svm,csvm,vn,cvn,
                      "NM MatrixView,SmallVector,Vector");
    TestMatrixArith3b(mvnm,cmvnm,vm,cvm,svn,csvn,
                      "NM MatrixView,Vector,SmallVector");
    TestMatrixArith3b(mvnm,cmvnm,svm,csvm,vvn,cvvn,
                      "NM MatrixView,SmallVector,VectorView");
    TestMatrixArith3b(mvnm,cmvnm,vvm,cvvm,svn,csvn,
                      "NM MatrixView,VectorView,SmallVector");
#endif
#endif

#if 1
    TestMatrixArith4(mmn,cmmn,smmn,csmmn,
                     "MN Matrix,SmallMatrix");
    TestMatrixArith4(smmn,csmmn,mmn,cmmn,
                     "MN SmallMatrix,Matrix");
#if (XTEST & 2)
    TestMatrixArith4(mnm,cmnm,smnm,csmnm,
                     "NM Matrix,SmallMatrix");
    TestMatrixArith4(smnm,csmnm,mnm,cmnm,
                     "NM SmallMatrix,Matrix");
    TestMatrixArith4(mvmn,cmvmn,smmn,csmmn,
                     "MN MatrixView,SmallMatrix");
    TestMatrixArith4(smmn,csmmn,mvmn,cmvmn,
                     "MN SmallMatrix,MatrixView");
    TestMatrixArith4(mvnm,cmvnm,smnm,csmnm,
                     "NM MatrixView,SmallMatrix");
    TestMatrixArith4(smnm,csmnm,mvnm,cmvnm,
                     "NM SmallMatrix,MatrixView");
#endif
#endif

#if 1
    TestMatrixArith5(smmn,csmmn,smnn,csmnn,
                     "MN*NN SmallMatrix,SmallMatrix");
    TestMatrixArith5(smmn,csmmn,mnn,cmnn,
                     "MN*NN SmallMatrix,Matrix");
    TestMatrixArith5(mmn,cmmn,smnn,csmnn,
                     "MN*NN Matrix,SmallMatrix");
#endif

#if 1
    TestMatrixArith6(smmn,csmmn,smnn,csmnn,mmn,cmmn,
                     "MN*NN SmallMatrix,SmallMatrix,Matrix");
    TestMatrixArith6(smmn,csmmn,mnn,cmnn,mmn,cmmn,
                     "MN*NN SmallMatrix,Matrix,Matrix");
    TestMatrixArith6(mmn,cmmn,smnn,csmnn,smmn,csmmn,
                     "MN*NN Matrix,SmallMatrix,SmallMatrix");
#if (XTEST & 2)
    TestMatrixArith6(mmn,cmmn,mnn,cmnn,smmn,csmmn,
                     "MN*NN Matrix,Matrix,SmallMatrix");
    TestMatrixArith6(smnn,csmnn,smnm,csmnm,mnm,cmnm,
                     "NN*NM SmallMatrix,SmallMatrix,Matrix");
    TestMatrixArith6(mnn,cmnn,smnm,csmnm,mnm,cmnm,
                     "NN*NM Matrix,SmallMatrix,Matrix");
    TestMatrixArith6(smnn,csmnn,mnm,cmnm,smnm,csmnm,
                     "NN*NM SmallMatrix,Matrix,SmallMatrix");
    TestMatrixArith6(mnn,cmnn,mnm,cmnm,smnm,csmnm,
                     "NN*NM Matrix,Matrix,SmallMatrix");
    TestMatrixArith6(smnm,csmnm,smmn,csmmn,smnn,csmnn,
                     "NM*MN SmallMatrix,SmallMatrix,SmallMatrix");
    TestMatrixArith6(smnm,csmnm,smmn,csmmn,mnn,cmnn,
                     "NM*MN SmallMatrix,SmallMatrix,Matrix");
    TestMatrixArith6(smnm,csmnm,mmn,cmmn,smnn,csmnn,
                     "NM*MN SmallMatrix,Matrix,SmallMatrix");
    TestMatrixArith6(smnm,csmnm,mmn,cmmn,mnn,cmnn,
                     "NM*MN SmallMatrix,Matrix,Matrix");
    TestMatrixArith6(mnm,cmnm,smmn,csmmn,smnn,csmnn,
                     "NM*MN Matrix,SmallMatrix,SmallMatrix");
    TestMatrixArith6(mnm,cmnm,smmn,csmmn,mnn,cmnn,
                     "NM*MN Matrix,SmallMatrix,Matrix");
    TestMatrixArith6(mnm,cmnm,mmn,cmmn,smnn,csmnn,
                     "NM*MN Matrix,Matrix,SmallMatrix");
    TestMatrixArith6(mnm,cmnm,mmn,cmmn,mnn,cmnn,
                     "NM*MN Matrix,Matrix,Matrix");
    TestMatrixArith5(smmn,csmmn,mvnn,cmvnn,
                     "MN*NN SmallMatrix,MatrixView");
    TestMatrixArith5(mvmn,cmvmn,smnn,csmnn,
                     "MN*NN MatrixView,SmallMatrix");
    TestMatrixArith6(smmn,csmmn,smnn,csmnn,mvmn,cmvmn,
                     "MN*NN SmallMatrix,SmallMatrix,MatrixView");
    TestMatrixArith6(smmn,csmmn,mvnn,cmvnn,mvmn,cmvmn,
                     "MN*NN SmallMatrix,MatrixView,MatrixView");
    TestMatrixArith6(mvmn,cmvmn,smnn,csmnn,smmn,csmmn,
                     "MN*NN MatrixView,SmallMatrix,SmallMatrix");
    TestMatrixArith6(mvmn,cmvmn,mvnn,cmvnn,smmn,csmmn,
                     "MN*NN MatrixView,MatrixView,SmallMatrixView");
    TestMatrixArith6(smnn,csmnn,smnm,csmnm,mvnm,cmvnm,
                     "NN*NM SmallMatrix,SmallMatrix,MatrixView");
    TestMatrixArith6(mvnn,cmvnn,smnm,csmnm,mvnm,cmvnm,
                     "NN*NM MatrixView,SmallMatrix,MatrixView");
    TestMatrixArith6(smnn,csmnn,mvnm,cmvnm,smnm,csmnm,
                     "NN*NM SmallMatrix,MatrixView,SmallMatrix");
    TestMatrixArith6(mvnn,cmvnn,mvnm,cmvnm,smnm,csmnm,
                     "NN*NM MatrixView,MatrixView,SmallMatrixView");
    TestMatrixArith6(smnm,csmnm,smmn,csmmn,mvnn,cmvnn,
                     "NM*MN SmallMatrix,SmallMatrix,MatrixView");
    TestMatrixArith6(smnm,csmnm,mvmn,cmvmn,smnn,csmnn,
                     "NM*MN SmallMatrix,MatrixView,SmallMatrix");
    TestMatrixArith6(smnm,csmnm,mvmn,cmvmn,mvnn,cmvnn,
                     "NM*MN SmallMatrix,MatrixView,MatrixView");
    TestMatrixArith6(mvnm,cmvnm,smmn,csmmn,smnn,csmnn,
                     "NM*MN MatrixView,SmallMatrix,SmallMatrix");
    TestMatrixArith6(mvnm,cmvnm,smmn,csmmn,mvnn,cmvnn,
                     "NM*MN MatrixView,SmallMatrix,MatrixView");
    TestMatrixArith6(mvnm,cmvnm,mvmn,cmvmn,smnn,csmnn,
                     "NM*MN MatrixView,MatrixView,SmallMatrix");
    TestMatrixArith6(mvnm,cmvnm,mvmn,cmvmn,mvnn,cmvnn,
                     "NM*MN MatrixView,MatrixView,MatrixView");
#endif
#endif

#if 1
    TestMatrixArith7(smmn,csmmn,svm,csvm,svn,csvn,
                     "MN SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith7(smmn,csmmn,svm,csvm,vn,cvn,
                     "MN SmallMatrix,SmallVector,Vector");
    TestMatrixArith7(smmn,csmmn,vm,cvm,svn,csvn,
                     "MN SmallMatrix,Vector,SmallVector");
#if (XTEST & 2)
    TestMatrixArith7(smmn,csmmn,vm,cvm,vn,cvn,
                     "MN SmallMatrix,Vector,Vector");
    TestMatrixArith7(smmn,csmmn,svm,csvm,vvn,cvvn,
                     "MN SmallMatrix,SmallVector,VectorView");
    TestMatrixArith7(smmn,csmmn,vvm,cvvm,svn,csvn,
                     "MN SmallMatrix,VectorView,SmallVector");
    TestMatrixArith7(smmn,csmmn,vvm,cvvm,vvn,cvvn,
                     "MN SmallMatrix,VectorView,VectorView");
    TestMatrixArith7(mmn,cmmn,svm,csvm,vn,cvn,
                     "MN Matrix,SmallVector,Vector");
    TestMatrixArith7(mmn,cmmn,vm,cvm,svn,csvn,
                     "MN Matrix,Vector,SmallVector");
    TestMatrixArith7(mmn,cmmn,svm,csvm,vvn,cvvn,
                     "MN Matrix,SmallVector,VectorView");
    TestMatrixArith7(mmn,cmmn,vvm,cvvm,svn,csvn,
                     "MN Matrix,VectorView,SmallVector");
    TestMatrixArith7(mvmn,cmvmn,svm,csvm,vn,cvn,
                     "MN MatrixView,SmallVector,Vector");
    TestMatrixArith7(mvmn,cmvmn,vm,cvm,svn,csvn,
                     "MN MatrixView,Vector,SmallVector");
    TestMatrixArith7(mvmn,cmvmn,svm,csvm,vvn,cvvn,
                     "MN MatrixView,SmallVector,VectorView");
    TestMatrixArith7(mvmn,cmvmn,vvm,cvvm,svn,csvn,
                     "MN MatrixView,VectorView,SmallVector");

    TestMatrixArith7(smnm,csmnm,svn,csvn,svm,csvm,
                     "NM SmallMatrix,SmallVector,SmallVector");
    TestMatrixArith7(smnm,csmnm,svn,csvn,vm,cvm,
                     "NM SmallMatrix,SmallVector,Vector");
    TestMatrixArith7(smnm,csmnm,vn,cvn,svm,csvm,
                     "NM SmallMatrix,Vector,SmallVector");
    TestMatrixArith7(smnm,csmnm,vn,cvn,vm,cvm,
                     "NM SmallMatrix,Vector,Vector");
    TestMatrixArith7(smnm,csmnm,svn,csvn,vvm,cvvm,
                     "NM SmallMatrix,SmallVector,VectorView");
    TestMatrixArith7(smnm,csmnm,vvn,cvvn,svm,csvm,
                     "NM SmallMatrix,VectorView,SmallVector");
    TestMatrixArith7(smnm,csmnm,vvn,cvvn,vvm,cvvm,
                     "NM SmallMatrix,VectorView,VectorView");
    TestMatrixArith7(mnm,cmnm,svn,csvn,vm,cvm,
                     "NM Matrix,SmallVector,Vector");
    TestMatrixArith7(mnm,cmnm,vn,cvn,svm,csvm,
                     "NM Matrix,Vector,SmallVector");
    TestMatrixArith7(mnm,cmnm,svn,csvn,vvm,cvvm,
                     "NM Matrix,SmallVector,VectorView");
    TestMatrixArith7(mnm,cmnm,vvn,cvvn,svm,csvm,
                     "NM Matrix,VectorView,SmallVector");
    TestMatrixArith7(mvnm,cmvnm,svn,csvn,vm,cvm,
                     "NM MatrixView,SmallVector,Vector");
    TestMatrixArith7(mvnm,cmvnm,vn,cvn,svm,csvm,
                     "NM MatrixView,Vector,SmallVector");
    TestMatrixArith7(mvnm,cmvnm,svn,csvn,vvm,cvvm,
                     "NM MatrixView,SmallVector,VectorView");
    TestMatrixArith7(mvnm,cmvnm,vvn,cvvn,svm,csvm,
                     "NM MatrixView,VectorView,SmallVector");
#endif
#endif

    std::cerr<<"No mixing errors\n";
    return 0;
} 
#if 1
#ifndef NOTHROW
catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
    return 1;
}
#endif
catch (std::exception& e) {
    std::cerr<<e.what()<<std::endl;
    std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
    return 1;
} 
catch (...) {
    std::cerr<<"Unknown exception thrown\n";
    std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
    return 1;
}
#else
catch (double) {}
#endif

