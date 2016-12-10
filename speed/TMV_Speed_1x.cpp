
//#undef NDEBUG
//#define PRINTALGO_XV
#define TMV_INLINE

#include <iostream>
#include "TMV.h"

// How big do you want the vectors to be?
const int N = 1333;
//const int N = 13;

// Define the type to use:
#define TISFLOAT

// Define the operation to test:
#define TMV_OP(A,B,C) \
    C = A(0) * A + B(0) * B;
#define EIGEN_OP(A,B,C) \
    C = A(0) * A + B(0) * B;
#define EIGEN_OP_MIX_A(A,B,C) \
    C = A(0) * A.cast<CT>() + B(0) * B;
#define EIGEN_OP_MIX_B(A,B,C) \
    C = A(0) * A + B(0) * B.cast<CT>();
#define EIGEN_OP_MIX_C(A,B,C) \
    C = A(0) * A.cast<CT>() + B(0) * B.cast<CT>();

// Define which versions you want to test:
// Note: often it seems that accurate timings only happen when doing one
// version at a time.  I don't really undrestand why.
#define DOSMALL
#define DOEIGEN
#define DOEIGENSMALL

#define DOREAL
#define DOMIX
#define DOCOMPLEX

// Set this if you only want to do a single loop.
// Not so useful for timing, but useful for debugging.
//#define ONELOOP
//#define SIMPLEVALUES

// Set up the target number of operations and memory to use for testing
const int targetnflops = 10000000; // in real ops
const int targetmem = 1000000; // in bytes

// The rest of this is based on the above values.
// You shouldn't need to change anything else here.

#include <sys/time.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#if (defined(DOEIGEN) || defined(DOEIGENSMALL))
#include "Eigen/Core"
#include "Eigen/Array"
#endif

#ifdef TISFLOAT
typedef float RT;
typedef std::complex<float> CT;
const float ACC = 1.e-5;
#else
typedef double RT;
typedef std::complex<double> CT;
const double ACC = 1.e-12;
#endif

#ifdef ONELOOP
const int nloops1 = 10;
const int nloops2 = 10;
#else
const int nloops2 = targetmem / (N * sizeof(RT)) + 1;
const int nloops1 = targetnflops / (N * nloops2) + 1;
#endif

#ifdef DOEIGEN

#ifdef TISFLOAT
#define EIGENV VectorXf
#define EIGENCV VectorXcf
#define EIGENM MatrixXf
#define EIGENCM MatrixXcf
#define ACC 1.e-5
#else
#define EIGENV VectorXd
#define EIGENCV VectorXcd
#define EIGENM MatrixXd
#define EIGENCM MatrixXcd
#define ACC 1.e-12
#endif

#endif

#ifdef DOEIGENSMALL
const int nloops2x = (
    (N*sizeof(RT) < 1536 * 1024 && N!=10000) ? nloops2 : 0 );
#endif

static void ClearCache()
{
    static tmv::Vector<double> X(10000,8.);
    static tmv::Vector<double> Y(10000,8.);
    static tmv::Vector<double> Z(10000);
    Z = X + Y;
    if (Norm(Z) < 5.) exit(1);
}

#ifdef SIMPLEVALUES
#define RAND1 ( RT(1+i) )
#define RAND2 ( CT(1+i,2+i) )
#else
#define RAND1 ( RT(rand()) / RT(RAND_MAX) )
#define RAND2 ( CT(rand(),rand()) / RT(RAND_MAX) )
#endif

int main() try 
{
    srand(518423972);

#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    std::vector<tmv::Vector<RT> > A1a(nloops2,tmv::Vector<RT>(N));
    std::vector<tmv::Vector<RT> > B1a(nloops2,tmv::Vector<RT>(N));
    std::vector<tmv::Vector<RT> > C1a(nloops2,tmv::Vector<RT>(N));
    std::vector<tmv::Matrix<RT> > A1b(nloops2,tmv::Matrix<RT>(5,N));
    std::vector<tmv::Matrix<RT> > B1b(nloops2,tmv::Matrix<RT>(5,N));
    std::vector<tmv::Matrix<RT> > C1b(nloops2,tmv::Matrix<RT>(5,N));
    std::vector<tmv::Vector<CT> > A1c(nloops2,tmv::Vector<CT>(N));
    std::vector<tmv::Vector<CT> > B1c(nloops2,tmv::Vector<CT>(N));
    std::vector<tmv::Vector<CT> > C1c(nloops2,tmv::Vector<CT>(N));
    std::vector<tmv::Matrix<CT> > A1d(nloops2,tmv::Matrix<CT>(5,N));
    std::vector<tmv::Matrix<CT> > B1d(nloops2,tmv::Matrix<CT>(5,N));
    std::vector<tmv::Matrix<CT> > C1d(nloops2,tmv::Matrix<CT>(5,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) {
            A1a[k](i) = RAND1 * RT(2) - RT(1.0);
            B1a[k](i) = RAND1 * RT(2) - RT(1.0);
            C1a[k](i) = RAND1 * RT(2) - RT(1.0);
            A1b[k](0,i) = A1a[k](i);
            B1b[k](0,i) = B1a[k](i);
            C1b[k](0,i) = B1a[k](i);
            A1c[k](i) = RAND2 * RT(2) - CT(1.0,1.0);
            B1c[k](i) = RAND2 * RT(2) - CT(1.0,1.0);
            C1c[k](i) = RAND2 * RT(2) - CT(1.0,1.0);
            A1d[k](0,i) = B1c[k](i);
            B1d[k](0,i) = B1c[k](i);
            C1d[k](0,i) = B1c[k](i);
        }
    }
    std::cout<<"N = "<<N<<std::endl;
    std::cout<<"nloops = "<<nloops1<<" x "<<nloops2;
    std::cout<<" = "<<nloops1*nloops2<<std::endl;
#ifdef DOEIGENSMALL
    if (nloops2x == 0)
        std::cout<<"Vector is too big for the stack, so no \"Eigen Known\" tests.\n";
#endif
    std::cout<<"\nTime for: ";
    std::cout<<"  TMV    "<<" TMV Small"<<"   Eigen  "<<"Eigen Known"<<std::endl;

#ifdef DOSMALL
    std::vector<tmv::SmallVector<RT,N> > A2a(nloops2);
    std::vector<tmv::SmallVector<RT,N> > B2a(nloops2);
    std::vector<tmv::SmallVector<RT,N> > C2a(nloops2);
    std::vector<tmv::SmallMatrix<RT,5,N> > A2b(nloops2);
    std::vector<tmv::SmallMatrix<RT,5,N> > B2b(nloops2);
    std::vector<tmv::SmallMatrix<RT,5,N> > C2b(nloops2);
    std::vector<tmv::SmallVector<CT,N> > A2c(nloops2);
    std::vector<tmv::SmallVector<CT,N> > B2c(nloops2);
    std::vector<tmv::SmallVector<CT,N> > C2c(nloops2);
    std::vector<tmv::SmallMatrix<CT,5,N> > A2d(nloops2);
    std::vector<tmv::SmallMatrix<CT,5,N> > B2d(nloops2);
    std::vector<tmv::SmallMatrix<CT,5,N> > C2d(nloops2);

    for(int k=0; k<nloops2; ++k) {
        A2a[k] = A1a[k]; B2a[k] = B1a[k]; C2a[k] = C1a[k];
        A2b[k] = A1b[k]; B2b[k] = B1b[k]; C2b[k] = C1b[k];
        A2c[k] = A1c[k]; B2c[k] = B1c[k]; C2c[k] = C1c[k];
        A2d[k] = A1d[k]; B2d[k] = B1d[k]; C2d[k] = C1d[k];
    }
#endif

#ifdef DOEIGEN
    std::vector<Eigen::EIGENV> A4a(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> B4a(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENV> C4a(nloops2,Eigen::EIGENV(N));
    std::vector<Eigen::EIGENM> A4b(nloops2,Eigen::EIGENM(5,N));
    std::vector<Eigen::EIGENM> B4b(nloops2,Eigen::EIGENM(5,N));
    std::vector<Eigen::EIGENM> C4b(nloops2,Eigen::EIGENM(5,N));
    std::vector<Eigen::EIGENCV> A4c(nloops2,Eigen::EIGENCV(N));
    std::vector<Eigen::EIGENCV> B4c(nloops2,Eigen::EIGENCV(N));
    std::vector<Eigen::EIGENCV> C4c(nloops2,Eigen::EIGENCV(N));
    std::vector<Eigen::EIGENCM> A4d(nloops2,Eigen::EIGENCM(5,N));
    std::vector<Eigen::EIGENCM> B4d(nloops2,Eigen::EIGENCM(5,N));
    std::vector<Eigen::EIGENCM> C4d(nloops2,Eigen::EIGENCM(5,N));
    for(int k=0;k<nloops2;++k) {
        for(int i=0;i<N;++i) A4a[k](i) = A1a[k](i);
        for(int i=0;i<N;++i) B4a[k](i) = B1a[k](i);
        for(int i=0;i<N;++i) C4a[k](i) = C1a[k](i);
        for(int i=0;i<N;++i) A4b[k](0,i) = A1b[k](0,i);
        for(int i=0;i<N;++i) B4b[k](0,i) = B1b[k](0,i);
        for(int i=0;i<N;++i) C4b[k](0,i) = C1b[k](0,i);
        for(int i=0;i<N;++i) A4c[k](i) = A1c[k](i);
        for(int i=0;i<N;++i) B4c[k](i) = B1c[k](i);
        for(int i=0;i<N;++i) C4c[k](i) = C1c[k](i);
        for(int i=0;i<N;++i) A4d[k](0,i) = A1d[k](0,i);
        for(int i=0;i<N;++i) B4d[k](0,i) = B1d[k](0,i);
        for(int i=0;i<N;++i) C4d[k](0,i) = C1d[k](0,i);
    }
#endif

#ifdef DOEIGENSMALL
    std::vector<Eigen::Matrix<RT,N,1> > A5a;
    std::vector<Eigen::Matrix<RT,N,1> > B5a;
    std::vector<Eigen::Matrix<RT,N,1> > C5a;
    std::vector<Eigen::Matrix<RT,5,N> > A5b;
    std::vector<Eigen::Matrix<RT,5,N> > B5b;
    std::vector<Eigen::Matrix<RT,5,N> > C5b;
    std::vector<Eigen::Matrix<CT,N,1> > A5c;
    std::vector<Eigen::Matrix<CT,N,1> > B5c;
    std::vector<Eigen::Matrix<CT,N,1> > C5c;
    std::vector<Eigen::Matrix<CT,5,N> > A5d;
    std::vector<Eigen::Matrix<CT,5,N> > B5d;
    std::vector<Eigen::Matrix<CT,5,N> > C5d;
    if (nloops2x) {
        A5a.resize(nloops2x); B5a.resize(nloops2x); C5a.resize(nloops2x); 
        A5b.resize(nloops2x); B5b.resize(nloops2x); C5b.resize(nloops2x); 
        A5c.resize(nloops2x); B5c.resize(nloops2x); C5c.resize(nloops2x); 
        A5d.resize(nloops2x); B5d.resize(nloops2x); C5d.resize(nloops2x); 
    }
    for(int k=0;k<nloops2x;++k) {
        for(int i=0;i<N;++i) A5a[k](i) = A1a[k](i);
        for(int i=0;i<N;++i) B5a[k](i) = B1a[k](i);
        for(int i=0;i<N;++i) C5a[k](i) = C1a[k](i);
        for(int i=0;i<N;++i) A5b[k](0,i) = A1b[k](0,i);
        for(int i=0;i<N;++i) B5b[k](0,i) = B1b[k](0,i);
        for(int i=0;i<N;++i) C5b[k](0,i) = C1b[k](0,i);
        for(int i=0;i<N;++i) A5c[k](i) = A1c[k](i);
        for(int i=0;i<N;++i) B5c[k](i) = B1c[k](i);
        for(int i=0;i<N;++i) C5c[k](i) = C1c[k](i);
        for(int i=0;i<N;++i) A5d[k](0,i) = A1d[k](0,i);
        for(int i=0;i<N;++i) B5d[k](0,i) = B1d[k](0,i);
        for(int i=0;i<N;++i) C5d[k](0,i) = C1d[k](0,i);
    }
#endif

    timeval tp;

    double t1_reg=0., t1_small=0., t1_eigen=0., t1_smalleigen=0.;
    double t2_reg=0., t2_small=0., t2_eigen=0., t2_smalleigen=0.;
    double t3_reg=0., t3_small=0., t3_eigen=0., t3_smalleigen=0.;
    double t4_reg=0., t4_small=0., t4_eigen=0., t4_smalleigen=0.;
    double t5_reg=0., t5_small=0., t5_eigen=0., t5_smalleigen=0.;
    double t6_reg=0., t6_small=0., t6_eigen=0., t6_smalleigen=0.;
    double t7_reg=0., t7_small=0., t7_eigen=0., t7_smalleigen=0.;
    double t8_reg=0., t8_small=0., t8_eigen=0., t8_smalleigen=0.;
    double t9_reg=0., t9_small=0., t9_eigen=0., t9_smalleigen=0.;
    double t10_reg=0., t10_small=0., t10_eigen=0., t10_smalleigen=0.;
    double t11_reg=0., t11_small=0., t11_eigen=0., t11_smalleigen=0.;
    double t12_reg=0., t12_small=0., t12_eigen=0., t12_smalleigen=0.;
    double t13_reg=0., t13_small=0., t13_eigen=0., t13_smalleigen=0.;
    double t14_reg=0., t14_small=0., t14_eigen=0., t14_smalleigen=0.;
    double t15_reg=0., t15_small=0., t15_eigen=0., t15_smalleigen=0.;
    double t16_reg=0., t16_small=0., t16_eigen=0., t16_smalleigen=0.;
    double ta,tb;

    for (int i=0; i<nloops1; ++i) {

        ClearCache();

#ifdef DOREAL
#if 1 // aaa
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1a[k],C1a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2a[k],C2a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2a[k]-C1a[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aaa: Small rms = "<<rms<<std::endl;
                std::cout<<"C2a[0] = "<<C2a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                for (int k=0; k<nloops2; ++k) {
                    if (Norm(C2a[k]-C1a[k]) > ACC) {
                        std::cout<<"C2a["<<k<<"] = "<<C2a[k]<<std::endl;
                        std::cout<<"C1a["<<k<<"] = "<<C1a[k]<<std::endl;
                        std::cout<<"diff = "<<C2a[k]-C1a[k]<<std::endl;
                        std::cout<<"Norm(diff) = "<<Norm(C2a[k]-C1a[k])<<std::endl;
                    }
                }
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4a[k],B4a[k],C4a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aaa: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4a[0] = "<<C4a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5a[k],B5a[k],C5a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t1_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aaa: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5a[0] = "<<C5a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // baa
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1b[k].row(0),B1a[k],C1a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2b[k].row(0),B2a[k],C2a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2a[k]-C1a[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"baa: Small rms = "<<rms<<std::endl;
                std::cout<<"C2a[0] = "<<C2a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4b[k].row(0).transpose(),B4a[k],C4a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"baa: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4a[0] = "<<C4a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5b[k].row(0).transpose(),B5a[k],C5a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t2_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"baa: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5a[0] = "<<C5a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // aba
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1b[k].row(0),C1a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2b[k].row(0),C2a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2a[k]-C1a[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"aba: Small rms = "<<rms<<std::endl;
                std::cout<<"C2a[0] = "<<C2a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4a[k],B4b[k].row(0).transpose(),C4a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"aba: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4a[0] = "<<C4a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5a[k],B5b[k].row(0).transpose(),C5a[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t3_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5a[k](i)-C1a[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"aba: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5a[0] = "<<C5a[0]<<std::endl;
                std::cout<<"C1a[0] = "<<C1a[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // aab
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1a[k],C1b[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2a[k],C2b[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2b[k].row(0)-C1b[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aab: Small rms = "<<rms<<std::endl;
                std::cout<<"C2b[0] = "<<C2b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4a[k],B4a[k],C4b[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4b[k](0,i)-C1b[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aab: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4b[0] = "<<C4b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5a[k],B5a[k],C5b[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t4_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5b[k](0,i)-C1b[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aab: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5b[0] = "<<C5b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // bbb
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1b[k].row(0),B1b[k].row(0),C1b[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2b[k].row(0),B2b[k].row(0),C2b[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2b[k].row(0)-C1b[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"bbb: Small rms = "<<rms<<std::endl;
                std::cout<<"C2b[0] = "<<C2b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4b[k].row(0).transpose(),B4b[k].row(0).transpose(),C4b[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4b[k](0,i)-C1b[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"bbb: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4b[0] = "<<C4b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5b[k].row(0).transpose(),B5b[k].row(0).transpose(),C5b[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t5_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5b[k](0,i)-C1b[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1b[0].row(0)) + Norm(B1b[0].row(0));
            if (rms > ACC) {
                std::cout<<"bbb: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5b[0] = "<<C5b[0].row(0)<<std::endl;
                std::cout<<"C1b[0] = "<<C1b[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif
#endif

#ifdef DOMIX
#if 1 // aac
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1a[k],C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2a[k],C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aac: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_C(A4a[k],B4a[k],C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aac: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_C(A5a[k],B5a[k],C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t6_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aac: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // cac
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1c[k],B1a[k],C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2c[k],B2a[k],C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cac: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                for (int k=0; k<nloops2; ++k) {
                    if (Norm(C2c[k]-C1c[k]) > ACC) {
                        std::cout<<"C2c["<<k<<"] = "<<C2c[k]<<std::endl;
                        std::cout<<"C1c["<<k<<"] = "<<C1c[k]<<std::endl;
                        std::cout<<"diff = "<<C2c[k]-C1c[k]<<std::endl;
                        std::cout<<"Norm(diff) = "<<Norm(C2c[k]-C1c[k])<<std::endl;
                    }
                }
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_B(A4c[k],B4a[k],C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cac: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_B(A5c[k],B5a[k],C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t7_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cac: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // acc
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1c[k],C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2c[k],C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acc: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_A(A4a[k],B4c[k],C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acc: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_A(A5a[k],B5c[k],C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t8_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acc: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // aad
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1a[k],C1d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2a[k],C2d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2d[k].row(0)-C1d[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aad: Small rms = "<<rms<<std::endl;
                std::cout<<"C2d[0] = "<<C2d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                for (int k=0; k<nloops2; ++k) {
                    if (Norm(C2d[k].row(0)-C1d[k].row(0)) > ACC) {
                        std::cout<<"C2d["<<k<<"] = "<<C2d[k].row(0)<<std::endl;
                        std::cout<<"C1d["<<k<<"] = "<<C1d[k].row(0)<<std::endl;
                        std::cout<<"diff = "<<C2d[k].row(0)-C1d[k].row(0)<<std::endl;
                        std::cout<<"Norm(diff) = "<<Norm(C2d[k].row(0)-C1d[k].row(0))<<std::endl;
                    }
                }
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_C(A4a[k],B4a[k],C4d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aad: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4d[0] = "<<C4d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_C(A5a[k],B5a[k],C5d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t9_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"aad: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5d[0] = "<<C5d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // cad
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1c[k],B1a[k],C1d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2c[k],B2a[k],C2d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2d[k].row(0)-C1d[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cad: Small rms = "<<rms<<std::endl;
                std::cout<<"C2d[0] = "<<C2d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_B(A4c[k],B4a[k],C4d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cad: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4d[0] = "<<C4d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_B(A5c[k],B5a[k],C5d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t10_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1a[0]);
            if (rms > ACC) {
                std::cout<<"cad: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5d[0] = "<<C5d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // acd
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1a[k],B1c[k],C1d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2a[k],B2c[k],C2d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2d[k].row(0)-C1d[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acd: Small rms = "<<rms<<std::endl;
                std::cout<<"C2d[0] = "<<C2d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP_MIX_A(A4a[k],B4c[k],C4d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acd: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4d[0] = "<<C4d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP_MIX_A(A5a[k],B5c[k],C5d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t11_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1a[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"acd: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5d[0] = "<<C5d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif
#endif

#ifdef DOCOMPLEX
#if 1 // ccc
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1c[k],B1c[k],C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2c[k],B2c[k],C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccc: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4c[k],B4c[k],C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccc: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5c[k],B5c[k],C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t12_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccc: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // dcc
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1d[k].row(0),B1c[k],C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2d[k].row(0),B2c[k],C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"dcc: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4d[k].row(0).transpose(),B4c[k],C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"dcc: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5d[k].row(0).transpose(),B5c[k],C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t13_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"dcc: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // cdc
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1c[k],B1d[k].row(0),C1c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2c[k],B2d[k].row(0),C2c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2c[k]-C1c[k]);
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"cdc: Small rms = "<<rms<<std::endl;
                std::cout<<"C2c[0] = "<<C2c[0]<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4c[k],B4d[k].row(0).transpose(),C4c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"cdc: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4c[0] = "<<C4c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5c[k],B5d[k].row(0).transpose(),C5c[k]);

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t14_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5c[k](i)-C1c[k](i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"cdc: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5c[0] = "<<C5c[0].transpose()<<std::endl;
                std::cout<<"C1c[0] = "<<C1c[0]<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // ccd
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1c[k],B1c[k],C1d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2c[k],B2c[k],C2d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2d[k].row(0)-C1d[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccd: Small rms = "<<rms<<std::endl;
                std::cout<<"C2d[0] = "<<C2d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4c[k],B4c[k],C4d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccd: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4d[0] = "<<C4d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5c[k],B5c[k],C5d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t15_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1c[0]) + Norm(B1c[0]);
            if (rms > ACC) {
                std::cout<<"ccd: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5d[0] = "<<C5d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif

#if 1 // ddd
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A1d[k].row(0),B1d[k].row(0),C1d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_reg += tb-ta;

#ifdef DOSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            TMV_OP(A2d[k].row(0),B2d[k].row(0),C2d[k].row(0));

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_small += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) rms += NormSq(C2d[k].row(0)-C1d[k].row(0));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"ddd: Small rms = "<<rms<<std::endl;
                std::cout<<"C2d[0] = "<<C2d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGEN
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2; ++k)
            EIGEN_OP(A4d[k].row(0).transpose(),B4d[k].row(0).transpose(),C4d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_eigen += tb-ta;

        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C4d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"ddd: Eigen rms = "<<rms<<std::endl;
                std::cout<<"C4d[0] = "<<C4d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif

#ifdef DOEIGENSMALL
        gettimeofday(&tp,0);
        ta = tp.tv_sec + tp.tv_usec/1.e6;

        for (int k=0; k<nloops2x; ++k)
            EIGEN_OP(A5d[k].row(0).transpose(),B5d[k].row(0).transpose(),C5d[k].row(0).transpose());

        gettimeofday(&tp,0);
        tb = tp.tv_sec + tp.tv_usec/1.e6;
        t16_smalleigen += tb-ta;
        if (i == 0) {
            RT rms = 0.;
            for (int k=0; k<nloops2; ++k) 
                for (int i=0;i<N;++i) rms += tmv::TMV_NORM(C5d[k](0,i)-C1d[k](0,i));
            rms = sqrt(rms/nloops2);
            rms /= Norm(A1d[0].row(0)) + Norm(B1d[0].row(0));
            if (rms > ACC) {
                std::cout<<"ddd: SmallEigen rms = "<<rms<<std::endl;
                std::cout<<"C5d[0] = "<<C5d[0].row(0)<<std::endl;
                std::cout<<"C1d[0] = "<<C1d[0].row(0)<<std::endl;
                exit(1);
            }
        }
#endif
#endif
#endif


    }

#ifdef DOREAL
    std::cout<<"\naaa       "<<t1_reg<<"  "<<t1_small<<"  "<<t1_eigen<<"  "<<t1_smalleigen<<std::endl;
    std::cout<<"\nbaa       "<<t2_reg<<"  "<<t2_small<<"  "<<t2_eigen<<"  "<<t2_smalleigen<<std::endl;
    std::cout<<"\naba       "<<t3_reg<<"  "<<t3_small<<"  "<<t3_eigen<<"  "<<t3_smalleigen<<std::endl;
    std::cout<<"\naab       "<<t4_reg<<"  "<<t4_small<<"  "<<t4_eigen<<"  "<<t4_smalleigen<<std::endl;
    std::cout<<"\nbbb       "<<t5_reg<<"  "<<t5_small<<"  "<<t5_eigen<<"  "<<t5_smalleigen<<std::endl;
#endif
#ifdef DOMIX
    std::cout<<"\naac       "<<t6_reg<<"  "<<t6_small<<"  "<<t6_eigen<<"  "<<t6_smalleigen<<std::endl;
    std::cout<<"\ncac       "<<t7_reg<<"  "<<t7_small<<"  "<<t7_eigen<<"  "<<t7_smalleigen<<std::endl;
    std::cout<<"\nacc       "<<t8_reg<<"  "<<t8_small<<"  "<<t8_eigen<<"  "<<t8_smalleigen<<std::endl;
    std::cout<<"\naad       "<<t9_reg<<"  "<<t9_small<<"  "<<t9_eigen<<"  "<<t9_smalleigen<<std::endl;
    std::cout<<"\ncad       "<<t10_reg<<"  "<<t10_small<<"  "<<t10_eigen<<"  "<<t10_smalleigen<<std::endl;
    std::cout<<"\nacd       "<<t11_reg<<"  "<<t11_small<<"  "<<t11_eigen<<"  "<<t11_smalleigen<<std::endl;
#endif
#ifdef DOCOMPLEX
    std::cout<<"\nccc       "<<t12_reg<<"  "<<t12_small<<"  "<<t12_eigen<<"  "<<t12_smalleigen<<std::endl;
    std::cout<<"\ndcc       "<<t13_reg<<"  "<<t13_small<<"  "<<t13_eigen<<"  "<<t13_smalleigen<<std::endl;
    std::cout<<"\ncdc       "<<t14_reg<<"  "<<t14_small<<"  "<<t14_eigen<<"  "<<t14_smalleigen<<std::endl;
    std::cout<<"\nccd       "<<t15_reg<<"  "<<t15_small<<"  "<<t15_eigen<<"  "<<t15_smalleigen<<std::endl;
    std::cout<<"\nddd       "<<t16_reg<<"  "<<t16_small<<"  "<<t16_eigen<<"  "<<t16_smalleigen<<std::endl;
#endif

}
#if 0
catch (int) {}
#else
catch (tmv::Error& e) {
        std::cout<<"Caught error "<<e<<std::endl;
            return 1;
}
#endif

