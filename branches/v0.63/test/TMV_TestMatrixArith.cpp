
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> 
void TestAllMatrixArith()
{
    // First test that algorithms work for some relevant different sizes
#ifdef XTEST
#if 0
    const int NSIZE = 17;
    const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,30,31,32,63,111,128,137,337};
#else
    const int NSIZE = 4;
    const int sizear[NSIZE] = {1,9,111,337};
#endif
    for(int m1=0;m1<NSIZE;m1++) {
        for(int n1=0;n1<NSIZE;n1++) {
            int m = sizear[m1];
            int n = sizear[n1];
            if (showstartdone)
                std::cout<<"m,n = "<<m<<','<<n<<std::endl;
            tmv::Matrix<CT,tmv::ColMajor> c(m,n);
            tmv::Matrix<CT,tmv::ColMajor> c1(m,n);

            // Test various adds:
            {
                tmv::Matrix<CT,tmv::ColMajor> ac(m,n);
                tmv::Matrix<CT,tmv::ColMajor> bc(m,n);
                for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
                    ac(i,j) = CT(T(i+5),T(j-2));
                for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
                    bc(i,j) = CT(T(2*i-13),T(-3*j+8));
                tmv::Matrix<CT,tmv::RowMajor> ar = ac;
                tmv::Matrix<CT,tmv::RowMajor> br = bc;
                T eps = EPS*(T(1)+Norm(ac)+Norm(bc));

                c1 = ac + bc;
                Assert(Norm((c=ar+br)-c1) < eps,"ar+br");
                Assert(Norm((c=ar+bc)-c1) < eps,"ar+bc");
                Assert(Norm((c=ac+br)-c1) < eps,"ac+br");
                c1 = ac + tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=ac+bc.Conjugate())-c1) < eps,"ac+cbc");
                Assert(Norm((c=ar+br.Conjugate())-c1) < eps,"ar+cbr");
                Assert(Norm((c=ar+bc.Conjugate())-c1) < eps,"ar+cbc");
                Assert(Norm((c=ac+br.Conjugate())-c1) < eps,"ac+cbr");
                c1 = tmv::Matrix<CT>(ac.Conjugate()) + bc;
                Assert(Norm((c=ac.Conjugate()+bc)-c1) < eps,"cac+bc");
                Assert(Norm((c=ar.Conjugate()+br)-c1) < eps,"car+br");
                Assert(Norm((c=ar.Conjugate()+bc)-c1) < eps,"car+bc");
                Assert(Norm((c=ac.Conjugate()+br)-c1) < eps,"cac+br");
                c1 = tmv::Matrix<CT>(ac.Conjugate()) +
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=ac.Conjugate()+bc.Conjugate())-c1) < eps,
                       "cac+cbc");
                Assert(Norm((c=ar.Conjugate()+br.Conjugate())-c1) < eps,
                       "car+cbr");
                Assert(Norm((c=ar.Conjugate()+bc.Conjugate())-c1) < eps,
                       "car+cbc");
                Assert(Norm((c=ac.Conjugate()+br.Conjugate())-c1) < eps,
                       "cac+cbr");

                T x1(7);
                T x2(-3);
                c1 = x1*ac + x2*bc;
                Assert(Norm((c=x1*ar+x2*br)-c1) < eps,"xar+xbr");
                Assert(Norm((c=x1*ar+x2*bc)-c1) < eps,"xar+xbc");
                Assert(Norm((c=x1*ac+x2*br)-c1) < eps,"xac+xbr");
                c1 = x1*ac + x2*tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=x1*ac+x2*bc.Conjugate())-c1) < eps,"xac+xcbc");
                Assert(Norm((c=x1*ar+x2*br.Conjugate())-c1) < eps,"xar+xcbr");
                Assert(Norm((c=x1*ar+x2*bc.Conjugate())-c1) < eps,"xar+xcbc");
                Assert(Norm((c=x1*ac+x2*br.Conjugate())-c1) < eps,"xac+xcbr");
                c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) + x2*bc;
                Assert(Norm((c=x1*ac.Conjugate()+x2*bc)-c1) < eps,"xcac+xbc");
                Assert(Norm((c=x1*ar.Conjugate()+x2*br)-c1) < eps,"xcar+xbr");
                Assert(Norm((c=x1*ar.Conjugate()+x2*bc)-c1) < eps,"xcar+xbc");
                Assert(Norm((c=x1*ac.Conjugate()+x2*br)-c1) < eps,"xcac+xbr");
                c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) +
                    x2*tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=x1*ac.Conjugate()+x2*bc.Conjugate())-c1) < eps,
                       "xcac+xcbc");
                Assert(Norm((c=x1*ar.Conjugate()+x2*br.Conjugate())-c1) < eps,
                       "xcar+xcbr");
                Assert(Norm((c=x1*ar.Conjugate()+x2*bc.Conjugate())-c1) < eps,
                       "xcar+xcbc");
                Assert(Norm((c=x1*ac.Conjugate()+x2*br.Conjugate())-c1) < eps,
                       "xcac+xcbr");

                CT z1(7,-9);
                CT z2(-3,-8);
                c1 = z1*ac + z2*bc;
                Assert(Norm((c=z1*ar+z2*br)-c1) < eps,"zar+zbr");
                Assert(Norm((c=z1*ar+z2*bc)-c1) < eps,"zar+zbc");
                Assert(Norm((c=z1*ac+z2*br)-c1) < eps,"zac+zbr");
                c1 = z1*ac + z2*tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=z1*ac+z2*bc.Conjugate())-c1) < eps,"zac+zcbc");
                Assert(Norm((c=z1*ar+z2*br.Conjugate())-c1) < eps,"zar+zcbr");
                Assert(Norm((c=z1*ar+z2*bc.Conjugate())-c1) < eps,"zar+zcbc");
                Assert(Norm((c=z1*ac+z2*br.Conjugate())-c1) < eps,"zac+zcbr");
                c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) + z2*bc;
                Assert(Norm((c=z1*ac.Conjugate()+z2*bc)-c1) < eps,"zcac+zbc");
                Assert(Norm((c=z1*ar.Conjugate()+z2*br)-c1) < eps,"zcar+zbr");
                Assert(Norm((c=z1*ar.Conjugate()+z2*bc)-c1) < eps,"zcar+zbc");
                Assert(Norm((c=z1*ac.Conjugate()+z2*br)-c1) < eps,"zcac+zbr");
                c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) +
                    z2*tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=z1*ac.Conjugate()+z2*bc.Conjugate())-c1) < eps,
                       "zcac+zcbc");
                Assert(Norm((c=z1*ar.Conjugate()+z2*br.Conjugate())-c1) < eps,
                       "zcar+zcbr");
                Assert(Norm((c=z1*ar.Conjugate()+z2*bc.Conjugate())-c1) < eps,
                       "zcar+zcbc");
                Assert(Norm((c=z1*ac.Conjugate()+z2*br.Conjugate())-c1) < eps,
                       "zcac+zcbr");
            }

            // Now test multiplies
            for(int k1=0;k1<NSIZE;k1++) {
                int k = sizear[k1];
                tmv::Matrix<CT,tmv::ColMajor> ac(m,k);
                tmv::Matrix<CT,tmv::ColMajor> bc(k,n);
                for(int i=0;i<m;i++) for(int j=0;j<k;j++) 
                    ac(i,j) = CT(T(i+5),T(j-2));
                for(int i=0;i<k;i++) for(int j=0;j<n;j++) 
                    bc(i,j) = CT(T(2*i-13),T(-3*j+8));
                tmv::Matrix<CT,tmv::RowMajor> ar = ac;
                tmv::Matrix<CT,tmv::RowMajor> br = bc;
                T eps = T(10) * EPS * (T(1) + Norm(ac)*Norm(bc));

                c1 = ac * bc;
                Assert(Norm((c=ar*br)-c1) < eps,"ar*br");
                Assert(Norm((c=ar*bc)-c1) < eps,"ar*bc");
                Assert(Norm((c=ac*br)-c1) < eps,"ac*br");
                c1 = ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=ac*bc.Conjugate())-c1) < eps,"ac*cbc");
                Assert(Norm((c=ar*br.Conjugate())-c1) < eps,"ar*cbr");
                Assert(Norm((c=ar*bc.Conjugate())-c1) < eps,"ar*cbc");
                Assert(Norm((c=ac*br.Conjugate())-c1) < eps,"ac*cbr");
                c1 = tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm((c=ac.Conjugate()*bc)-c1) < eps,"cac*bc");
                Assert(Norm((c=ar.Conjugate()*br)-c1) < eps,"car*br");
                Assert(Norm((c=ar.Conjugate()*bc)-c1) < eps,"car*bc");
                Assert(Norm((c=ac.Conjugate()*br)-c1) < eps,"cac*br");
                c1 = tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=ac.Conjugate()*bc.Conjugate())-c1) < eps,
                       "cac*cbc");
                Assert(Norm((c=ar.Conjugate()*br.Conjugate())-c1) < eps,
                       "car*cbr");
                Assert(Norm((c=ar.Conjugate()*bc.Conjugate())-c1) < eps,
                       "car*cbc");
                Assert(Norm((c=ac.Conjugate()*br.Conjugate())-c1) < eps,
                       "cac*cbr");

                tmv::Matrix<CT> c0 = c1;
                c1 = c0 + ac * bc;
                Assert(Norm(((c=c0)+=ar*br)-c1) < eps,"+ar*br");
                Assert(Norm(((c=c0)+=ar*bc)-c1) < eps,"+ar*bc");
                Assert(Norm(((c=c0)+=ac*br)-c1) < eps,"+ac*br");
                c1 = c0 + ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=ac*bc.Conjugate())-c1) < eps,"+ac*cbc");
                Assert(Norm(((c=c0)+=ar*br.Conjugate())-c1) < eps,"+ar*cbr");
                Assert(Norm(((c=c0)+=ar*bc.Conjugate())-c1) < eps,"+ar*cbc");
                Assert(Norm(((c=c0)+=ac*br.Conjugate())-c1) < eps,"+ac*cbr");
                c1 = c0 + tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm(((c=c0)+=ac.Conjugate()*bc)-c1) < eps,"+cac*bc");
                Assert(Norm(((c=c0)+=ar.Conjugate()*br)-c1) < eps,"+car*br");
                Assert(Norm(((c=c0)+=ar.Conjugate()*bc)-c1) < eps,"+car*bc");
                Assert(Norm(((c=c0)+=ac.Conjugate()*br)-c1) < eps,"+cac*br");
                c1 = c0 + tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=ac.Conjugate()*bc.Conjugate())-c1) < eps,
                       "+cac*cbc");
                Assert(Norm(((c=c0)+=ar.Conjugate()*br.Conjugate())-c1) < eps,
                       "+car*cbr");
                Assert(Norm(((c=c0)+=ar.Conjugate()*bc.Conjugate())-c1) < eps,
                       "+car*cbc");
                Assert(Norm(((c=c0)+=ac.Conjugate()*br.Conjugate())-c1) < eps,
                       "+cac*cbr");

                T x1(7);
                eps *= x1;
                c1 = x1*ac * bc;
                Assert(Norm((c=x1*ar*br)-c1) < eps,"xar*br");
                Assert(Norm((c=x1*ar*bc)-c1) < eps,"xar*bc");
                Assert(Norm((c=x1*ac*br)-c1) < eps,"xac*br");
                c1 = x1*ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=x1*ac*bc.Conjugate())-c1) < eps,"xac*cbc");
                Assert(Norm((c=x1*ar*br.Conjugate())-c1) < eps,"xar*cbr");
                Assert(Norm((c=x1*ar*bc.Conjugate())-c1) < eps,"xar*cbc");
                Assert(Norm((c=x1*ac*br.Conjugate())-c1) < eps,"xac*cbr");
                c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm((c=x1*ac.Conjugate()*bc)-c1) < eps,"xcac*bc");
                Assert(Norm((c=x1*ar.Conjugate()*br)-c1) < eps,"xcar*br");
                Assert(Norm((c=x1*ar.Conjugate()*bc)-c1) < eps,"xcar*bc");
                Assert(Norm((c=x1*ac.Conjugate()*br)-c1) < eps,"xcac*br");
                c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=x1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
                       "xcac*cbc");
                Assert(Norm((c=x1*ar.Conjugate()*br.Conjugate())-c1) < eps,
                       "xcar*cbr");
                Assert(Norm((c=x1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
                       "xcar*cbc");
                Assert(Norm((c=x1*ac.Conjugate()*br.Conjugate())-c1) < eps,
                       "xcac*cbr");

                c1 = c0 + x1*ac * bc;
                Assert(Norm(((c=c0)+=x1*ar*br)-c1) < eps,"+xar*br");
                Assert(Norm(((c=c0)+=x1*ar*bc)-c1) < eps,"+xar*bc");
                Assert(Norm(((c=c0)+=x1*ac*br)-c1) < eps,"+xac*br");
                c1 = c0 + x1*ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=x1*ac*bc.Conjugate())-c1) < eps,
                       "+xac*cbc");
                Assert(Norm(((c=c0)+=x1*ar*br.Conjugate())-c1) < eps,
                       "+xar*cbr");
                Assert(Norm(((c=c0)+=x1*ar*bc.Conjugate())-c1) < eps,
                       "+xar*cbc");
                Assert(Norm(((c=c0)+=x1*ac*br.Conjugate())-c1) < eps,
                       "+xac*cbr");
                c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc)-c1) < eps,
                       "+xcac*bc");
                Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br)-c1) < eps,
                       "+xcar*br");
                Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc)-c1) < eps,
                       "+xcar*bc");
                Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br)-c1) < eps,
                       "+xcac*br");
                c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc.Conjugate())-c1)<eps,
                       "+xcac*cbc");
                Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br.Conjugate())-c1)<eps,
                       "+xcar*cbr");
                Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc.Conjugate())-c1)<eps,
                       "+xcar*cbc");
                Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br.Conjugate())-c1)<eps,
                       "+xcac*cbr");

                CT z1(7,-9);
                eps *= std::norm(z1)/x1;
                c1 = z1*ac * bc;
                Assert(Norm((c=z1*ar*br)-c1) < eps,"zar*br");
                Assert(Norm((c=z1*ar*bc)-c1) < eps,"zar*bc");
                Assert(Norm((c=z1*ac*br)-c1) < eps,"zac*br");
                c1 = z1*ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=z1*ac*bc.Conjugate())-c1) < eps,"zac*cbc");
                Assert(Norm((c=z1*ar*br.Conjugate())-c1) < eps,"zar*cbr");
                Assert(Norm((c=z1*ar*bc.Conjugate())-c1) < eps,"zar*cbc");
                Assert(Norm((c=z1*ac*br.Conjugate())-c1) < eps,"zac*cbr");
                c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm((c=z1*ac.Conjugate()*bc)-c1) < eps,"zcac*bc");
                Assert(Norm((c=z1*ar.Conjugate()*br)-c1) < eps,"zcar*br");
                Assert(Norm((c=z1*ar.Conjugate()*bc)-c1) < eps,"zcar*bc");
                Assert(Norm((c=z1*ac.Conjugate()*br)-c1) < eps,"zcac*br");
                c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm((c=z1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
                       "zcac*cbc");
                Assert(Norm((c=z1*ar.Conjugate()*br.Conjugate())-c1) < eps,
                       "zcar*cbr");
                Assert(Norm((c=z1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
                       "zcar*cbc");
                Assert(Norm((c=z1*ac.Conjugate()*br.Conjugate())-c1) < eps,
                       "zcac*cbr");

                c1 = c0 + z1*ac * bc;
                Assert(Norm(((c=c0)+=z1*ar*br)-c1) < eps,"+zar*br");
                Assert(Norm(((c=c0)+=z1*ar*bc)-c1) < eps,"+zar*bc");
                Assert(Norm(((c=c0)+=z1*ac*br)-c1) < eps,"+zac*br");
                c1 = c0 + z1*ac * tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=z1*ac*bc.Conjugate())-c1) < eps,
                       "+zac*cbc");
                Assert(Norm(((c=c0)+=z1*ar*br.Conjugate())-c1) < eps,
                       "+zar*cbr");
                Assert(Norm(((c=c0)+=z1*ar*bc.Conjugate())-c1) < eps,
                       "+zar*cbc");
                Assert(Norm(((c=c0)+=z1*ac*br.Conjugate())-c1) < eps,
                       "+zac*cbr");
                c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
                Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc)-c1) < eps,
                       "+zcac*bc");
                Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br)-c1) < eps,
                       "+zcar*br");
                Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc)-c1) < eps,
                       "+zcar*bc");
                Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br)-c1) < eps,
                       "+zcac*br");
                c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) *
                    tmv::Matrix<CT>(bc.Conjugate());
                Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc.Conjugate())-c1)<eps,
                       "+zcac*cbc");
                Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br.Conjugate())-c1)<eps,
                       "+zcar*cbr");
                Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc.Conjugate())-c1)<eps,
                       "+zcar*cbc");
                Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br.Conjugate())-c1)<eps,
                       "+zcac*cbr");
            }
        }
    }
#endif

    // Now we use the TestMatrixArith.h file to test lots of different
    // syntaxes for calling matrix arithmetic.  This tests the inline
    // parser more than the algorithms.
    tmv::Matrix<T,tmv::RowMajor> a1(4,4);
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
        a1(i,j) = T(2+4*i-5*j);
    }
    a1(0,0) = 14; 
    a1(1,0) = -2; 
    a1(2,0) = 7; 
    a1(3,0) = -10;
    a1(2,2) = 30;

    tmv::Matrix<CT,tmv::RowMajor> ca1 = a1;
    ca1(2,3) += CT(2,3);
    ca1(1,0) *= CT(0,2);
    ca1.col(1) *= CT(-1,3);
    ca1.row(3) += tmv::Vector<CT>(4,CT(1,9));
    tmv::MatrixView<T> a1v = a1.View();
    tmv::MatrixView<CT> ca1v = ca1.View();

    tmv::Matrix<T,tmv::ColMajor> a2 = a1.Transpose();
    a2.row(1) *= T(3);
    a2.col(2) -= tmv::Vector<T>(4,4.);
    tmv::Matrix<CT,tmv::ColMajor> ca2 = ca1;
    ca2 -= a2;
    ca2 *= CT(1,-2);
    tmv::MatrixView<T> a2v = a2.View();
    tmv::MatrixView<CT> ca2v = ca2.View();

    tmv::Matrix<T> a3x(12,16);
    for(int i=0;i<12;++i) for(int j=0;j<16;++j) a3x(i,j) = T(1-2*i+3*j);
    a3x.diag().AddToAll(30);
    tmv::Matrix<CT> ca3x = a3x*CT(1,-2);
    ca3x.diag().AddToAll(CT(-22,15));
    tmv::MatrixView<T> a3v = a3x.SubMatrix(0,12,0,16,3,4);
    tmv::MatrixView<CT> ca3v = ca3x.SubMatrix(0,12,0,16,3,4);

    tmv::Matrix<T> a1x(4,4);
    tmv::Matrix<CT> ca1x(4,4);

    TestMatrixArith123<T>(a1x,ca1x,a1v,ca1v,"Square 1");
    TestMatrixArith123<T>(a1x,ca1x,a2v,ca2v,"Square 2");
    TestMatrixArith123<T>(a1x,ca1x,a3v,ca3v,"Square 3");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a2v,ca2v,"Square 1");
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a1v,ca1v,"Square 2");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a1v,ca1v,"Square 3");
#ifdef XTEST
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a3v,ca3v,"Square 4");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a2v,ca2v,"Square 5");
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a3v,ca3v,"Square 6");
#endif

    tmv::Vector<T> v1 = a1.col(0);
    tmv::VectorView<T> v1v = v1.View();
    tmv::Vector<T> v15(20);
    tmv::VectorView<T> v1s = v15.SubVector(0,20,5);
    v1s = v1v;

    tmv::Vector<T> v2 = a1.row(2);
    tmv::VectorView<T> v2v = v2.View();
    tmv::Vector<T> v25(20);
    tmv::VectorView<T> v2s = v25.SubVector(0,20,5);
    v2s = v2v;

    tmv::Vector<CT> cv1 = ca1.col(0);
    tmv::VectorView<CT> cv1v = cv1.View();
    tmv::Vector<CT> cv15(20);
    tmv::VectorView<CT> cv1s = cv15.SubVector(0,20,5);
    cv1s = cv1v;

    tmv::Vector<CT> cv2 = ca1.row(2);
    tmv::VectorView<CT> cv2v = cv2.View();
    tmv::Vector<CT> cv25(20);
    tmv::VectorView<CT> cv2s = cv25.SubVector(0,20,5);
    cv2s = cv2v;

    TestMatrixArith7<T>(a1x,ca1x,a1v,ca1v,v1v,cv1v,v2v,cv2v,"Square 1");
    TestMatrixArith7<T>(a1x,ca1x,a1v,ca1v,v1s,cv1s,v2v,cv2v,"Square 2");
    TestMatrixArith7<T>(a1x,ca1x,a1v,ca1v,v1v,cv1v,v2s,cv2s,"Square 3");
    TestMatrixArith7<T>(a1x,ca1x,a1v,ca1v,v1s,cv1s,v2s,cv2s,"Square 4");
    TestMatrixArith7<T>(a1x,ca1x,a2v,ca2v,v1v,cv1v,v2v,cv2v,"Square 5");
    TestMatrixArith7<T>(a1x,ca1x,a2v,ca2v,v1s,cv1s,v2v,cv2v,"Square 6");
    TestMatrixArith7<T>(a1x,ca1x,a2v,ca2v,v1v,cv1v,v2s,cv2s,"Square 7");
    TestMatrixArith7<T>(a1x,ca1x,a2v,ca2v,v1s,cv1s,v2s,cv2s,"Square 8");
    TestMatrixArith7<T>(a1x,ca1x,a3v,ca3v,v1v,cv1v,v2v,cv2v,"Square 9");
    TestMatrixArith7<T>(a1x,ca1x,a3v,ca3v,v1s,cv1s,v2v,cv2v,"Square 10");
    TestMatrixArith7<T>(a1x,ca1x,a3v,ca3v,v1v,cv1v,v2s,cv2s,"Square 11");
    TestMatrixArith7<T>(a1x,ca1x,a3v,ca3v,v1s,cv1s,v2s,cv2s,"Square 12");

    tmv::Matrix<T,tmv::RowMajor> a4(7,4);
    for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = T(1-3*i+2*j);
    tmv::Matrix<T,tmv::ColMajor> a5 = a4.Transpose();
    a4.SubMatrix(2,6,0,4) += a1;
    a5.SubMatrix(0,4,1,5) -= a2;
    tmv::MatrixView<T> a4v = a4.View();
    tmv::MatrixView<T> a5v = a5.View();

    tmv::Matrix<CT,tmv::RowMajor> ca4 = a4*CT(1,2);
    tmv::Matrix<CT,tmv::ColMajor> ca5 = ca4.Adjoint();
    ca4.SubMatrix(2,6,0,4) += ca1;
    ca5.SubMatrix(0,4,1,5) -= ca2;
    ca4.col(1) *= CT(2,1);
    ca4.row(6).AddToAll(CT(-7,2));
    ca5.col(3) *= CT(-1,3);
    ca5.row(0).AddToAll(CT(1,9));
    tmv::MatrixView<CT> ca4v = ca4.View();
    tmv::MatrixView<CT> ca5v = ca5.View();

    tmv::Matrix<T> a4x(7,4);
    tmv::Matrix<CT> ca4x(7,4);
    tmv::Matrix<T> a5x(4,7);
    tmv::Matrix<CT> ca5x(4,7);
    TestMatrixArith123<T>(a4x,ca4x,a4v,ca4v,"NonSquare 1");
    TestMatrixArith123<T>(a5x,ca5x,a5v,ca5v,"NonSquare 2");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a4v,ca4v,"NonSquare 1");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a1v,ca1v,"NonSquare 2");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a5v,ca5v,"NonSquare 3");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a4v,ca4v,"NonSquare 4");
#ifdef XTEST
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a4v,ca4v,"NonSquare 5");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a4v,ca4v,"NonSquare 6");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a5v,ca5v,"NonSquare 7");
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a5v,ca5v,"NonSquare 8");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a5v,ca5v,"NonSquare 9");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a2v,ca2v,"NonSquare 10");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a3v,ca3v,"NonSquare 11");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a1v,ca1v,"NonSquare 12");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a2v,ca2v,"NonSquare 13");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a3v,ca3v,"NonSquare 14");
#endif

    tmv::Vector<T> v3 = a4.col(2);
    tmv::VectorView<T> v3v = v3.View();
    tmv::Vector<T> v35(35);
    tmv::VectorView<T> v3s = v35.SubVector(0,35,5);
    v3s = v3v;

    tmv::Vector<CT> cv3 = ca4.col(2);
    tmv::VectorView<CT> cv3v = cv3.View();
    tmv::Vector<CT> cv35(35);
    tmv::VectorView<CT> cv3s = cv35.SubVector(0,35,5);
    cv3s = cv3v;

    TestMatrixArith7<T>(a4x,ca4x,a4v,ca4v,v3v,cv3v,v2v,cv2v,"NonSquare 1");
    TestMatrixArith7<T>(a4x,ca4x,a4v,ca4v,v3s,cv3s,v2v,cv2v,"NonSquare 2");
    TestMatrixArith7<T>(a4x,ca4x,a4v,ca4v,v3v,cv3v,v2s,cv2s,"NonSquare 3");
    TestMatrixArith7<T>(a4x,ca4x,a4v,ca4v,v3s,cv3s,v2s,cv2s,"NonSquare 4");
    TestMatrixArith7<T>(a5x,ca5x,a5v,ca5v,v1v,cv1v,v3v,cv3v,"NonSquare 5");
    TestMatrixArith7<T>(a5x,ca5x,a5v,ca5v,v1s,cv1s,v3v,cv3v,"NonSquare 6");
    TestMatrixArith7<T>(a5x,ca5x,a5v,ca5v,v1v,cv1v,v3s,cv3s,"NonSquare 7");
    TestMatrixArith7<T>(a5x,ca5x,a5v,ca5v,v1s,cv1s,v3s,cv3s,"NonSquare 8");

#ifdef XTEST
    tmv::Matrix<T> a6(4,0,1);
    tmv::Matrix<T> a7(0,4,1);
    tmv::Matrix<CT> ca6 = a6;
    tmv::Matrix<CT> ca7 = a7;
    tmv::MatrixView<T> a6v = a6.View();
    tmv::MatrixView<T> a7v = a7.View();
    tmv::MatrixView<CT> ca6v = ca6.View();
    tmv::MatrixView<CT> ca7v = ca7.View();

    tmv::Matrix<T> a6x(4,0);
    tmv::Matrix<T> a7x(0,4);
    tmv::Matrix<CT> ca6x(4,0);
    tmv::Matrix<CT> ca7x(0,4);
    TestMatrixArith123<T>(a6x,ca6x,a6v,ca6v,"Degenerate 1");
    TestMatrixArith123<T>(a7x,ca7x,a7v,ca7v,"Degenerate 2");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a6v,ca6v,"Degenerate 1");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a1v,ca1v,"Degenerate 2");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a7v,ca7v,"Degenerate 3");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a6v,ca6v,"Degenerate 4");
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a6v,ca6v,"Degenerate 5");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a6v,ca6v,"Degenerate 6");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a6v,ca6v,"Degenerate 7");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a6v,ca6v,"Degenerate 8");
    TestMatrixArith456<T>(a1x,ca1x,a1v,ca1v,a7v,ca7v,"Degenerate 9");
    TestMatrixArith456<T>(a1x,ca1x,a2v,ca2v,a7v,ca7v,"Degenerate 10");
    TestMatrixArith456<T>(a1x,ca1x,a3v,ca3v,a7v,ca7v,"Degenerate 11");
    TestMatrixArith456<T>(a4x,ca4x,a4v,ca4v,a7v,ca7v,"Degenerate 12");
    TestMatrixArith456<T>(a5x,ca5x,a5v,ca5v,a7v,ca7v,"Degenerate 13");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a2v,ca2v,"Degenerate 14");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a3v,ca3v,"Degenerate 15");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a4v,ca4v,"Degenerate 16");
    TestMatrixArith456<T>(a6x,ca6x,a6v,ca6v,a5v,ca5v,"Degenerate 17");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a1v,ca1v,"Degenerate 18");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a2v,ca2v,"Degenerate 19");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a3v,ca3v,"Degenerate 20");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a4v,ca4v,"Degenerate 21");
    TestMatrixArith456<T>(a7x,ca7x,a7v,ca7v,a5v,ca5v,"Degenerate 22");

    tmv::Vector<T> v4 = a6.row(2);
    tmv::VectorView<T> v4v = v4.View();
    tmv::Vector<T> v45(0);
    tmv::VectorView<T> v4s = v45.SubVector(0,0,5);

    tmv::Vector<CT> cv4 = ca6.row(2);
    tmv::VectorView<CT> cv4v = cv4.View();
    tmv::Vector<CT> cv45(0);
    tmv::VectorView<CT> cv4s = cv45.SubVector(0,0,5);

    TestMatrixArith7<T>(a6x,ca6x,a6v,ca6v,v1v,cv1v,v4v,cv4v,"Degenerate 1");
    TestMatrixArith7<T>(a6x,ca6x,a6v,ca6v,v1s,cv1s,v4v,cv4v,"Degenerate 2");
    TestMatrixArith7<T>(a6x,ca6x,a6v,ca6v,v1v,cv1v,v4s,cv4s,"Degenerate 3");
    TestMatrixArith7<T>(a6x,ca6x,a6v,ca6v,v1s,cv1s,v4s,cv4s,"Degenerate 4");
    TestMatrixArith7<T>(a7x,ca7x,a7v,ca7v,v4v,cv4v,v2v,cv2v,"Degenerate 5");
    TestMatrixArith7<T>(a7x,ca7x,a7v,ca7v,v4s,cv4s,v2v,cv2v,"Degenerate 6");
    TestMatrixArith7<T>(a7x,ca7x,a7v,ca7v,v4v,cv4v,v2s,cv2s,"Degenerate 7");
    TestMatrixArith7<T>(a7x,ca7x,a7v,ca7v,v4s,cv4s,v2s,cv2s,"Degenerate 8");

#endif
}

#ifdef INST_DOUBLE
template void TestAllMatrixArith<double>();
#endif
#ifdef INST_FLOAT
template void TestAllMatrixArith<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllMatrixArith<long double>();
#endif
#ifdef INST_INT
template void TestAllMatrixArith<int>();
#endif
