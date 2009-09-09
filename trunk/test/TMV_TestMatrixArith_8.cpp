// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Mat.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestMatrixArith_8()
{
#if 0
  const int NSIZE = 15;
  //const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,30,31,32,63,64,65};
  const int sizear[NSIZE] = {1,2,3,4,5,63,64,65,66,70,111,128,137,256,337};
  //const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,32,128,509,1009,3009,5009};
#else
  const int NSIZE = 4;
  const int sizear[NSIZE] = {2,5,77,637};
#endif
  for(int m1=0;m1<NSIZE;m1++) {
    for(int n1=0;n1<NSIZE;n1++) {
      int m = sizear[m1];
      int n = sizear[n1];
      //if (showstartdone)
        std::cout<<"m,n = "<<m<<','<<n<<std::endl;
      tmv::Matrix<CT,tmv::ColMajor> c(m,n);
      tmv::Matrix<CT,tmv::ColMajor> c1(m,n);

#if 1
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
        Assert(Norm((c=ar+bc)-c1) < eps,"ar+bc");
#ifdef XTEST
        Assert(Norm((c=ar+br)-c1) < eps,"ar+br");
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
        c1 = tmv::Matrix<CT>(ac.Conjugate()) + tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm((c=ac.Conjugate()+bc.Conjugate())-c1) < eps,"cac+cbc");
        Assert(Norm((c=ar.Conjugate()+br.Conjugate())-c1) < eps,"car+cbr");
        Assert(Norm((c=ar.Conjugate()+bc.Conjugate())-c1) < eps,"car+cbc");
        Assert(Norm((c=ac.Conjugate()+br.Conjugate())-c1) < eps,"cac+cbr");

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
        c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) + x2*tmv::Matrix<CT>(bc.Conjugate());
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
        c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) + z2*tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm((c=z1*ac.Conjugate()+z2*bc.Conjugate())-c1) < eps,
            "zcac+zcbc");
        Assert(Norm((c=z1*ar.Conjugate()+z2*br.Conjugate())-c1) < eps,
            "zcar+zcbr");
        Assert(Norm((c=z1*ar.Conjugate()+z2*bc.Conjugate())-c1) < eps,
            "zcar+zcbc");
        Assert(Norm((c=z1*ac.Conjugate()+z2*br.Conjugate())-c1) < eps,
            "zcac+zcbr");
#endif
      }
#endif

#if 1
      // Now test multiplies
      for(int k1=0;k1<NSIZE;k1++) {
        int k = sizear[k1];
        if (showstartdone)
          std::cout<<"k = "<<k<<std::endl;
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
        Assert(Norm((c=ar*bc)-c1) < eps,"ar*bc");
#ifdef XTEST
        Assert(Norm((c=ar*br)-c1) < eps,"ar*br");
        Assert(Norm((c=ac*br)-c1) < eps,"ac*br");
#if 1
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
        c1 = tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm((c=ac.Conjugate()*bc.Conjugate())-c1) < eps,"cac*cbc");
        Assert(Norm((c=ar.Conjugate()*br.Conjugate())-c1) < eps,"car*cbr");
        Assert(Norm((c=ar.Conjugate()*bc.Conjugate())-c1) < eps,"car*cbc");
        Assert(Norm((c=ac.Conjugate()*br.Conjugate())-c1) < eps,"cac*cbr");

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
        c1 = c0 + tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm(((c=c0)+=ac.Conjugate()*bc.Conjugate())-c1) < eps,"+cac*cbc");
        Assert(Norm(((c=c0)+=ar.Conjugate()*br.Conjugate())-c1) < eps,"+car*cbr");
        Assert(Norm(((c=c0)+=ar.Conjugate()*bc.Conjugate())-c1) < eps,"+car*cbc");
        Assert(Norm(((c=c0)+=ac.Conjugate()*br.Conjugate())-c1) < eps,"+cac*cbr");

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
        c1 = x1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
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
        Assert(Norm(((c=c0)+=x1*ac*bc.Conjugate())-c1) < eps,"+xac*cbc");
        Assert(Norm(((c=c0)+=x1*ar*br.Conjugate())-c1) < eps,"+xar*cbr");
        Assert(Norm(((c=c0)+=x1*ar*bc.Conjugate())-c1) < eps,"+xar*cbc");
        Assert(Norm(((c=c0)+=x1*ac*br.Conjugate())-c1) < eps,"+xac*cbr");
        c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
        Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc)-c1) < eps,"+xcac*bc");
        Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br)-c1) < eps,"+xcar*br");
        Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc)-c1) < eps,"+xcar*bc");
        Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br)-c1) < eps,"+xcac*br");
        c1 = c0 + x1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm(((c=c0)+=x1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
            "+xcac*cbc");
        Assert(Norm(((c=c0)+=x1*ar.Conjugate()*br.Conjugate())-c1) < eps,
            "+xcar*cbr");
        Assert(Norm(((c=c0)+=x1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
            "+xcar*cbc");
        Assert(Norm(((c=c0)+=x1*ac.Conjugate()*br.Conjugate())-c1) < eps,
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
        c1 = z1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
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
        Assert(Norm(((c=c0)+=z1*ac*bc.Conjugate())-c1) < eps,"+zac*cbc");
        Assert(Norm(((c=c0)+=z1*ar*br.Conjugate())-c1) < eps,"+zar*cbr");
        Assert(Norm(((c=c0)+=z1*ar*bc.Conjugate())-c1) < eps,"+zar*cbc");
        Assert(Norm(((c=c0)+=z1*ac*br.Conjugate())-c1) < eps,"+zac*cbr");
        c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) * bc;
        Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc)-c1) < eps,"+zcac*bc");
        Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br)-c1) < eps,"+zcar*br");
        Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc)-c1) < eps,"+zcar*bc");
        Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br)-c1) < eps,"+zcac*br");
        c1 = c0 + z1*tmv::Matrix<CT>(ac.Conjugate()) * tmv::Matrix<CT>(bc.Conjugate());
        Assert(Norm(((c=c0)+=z1*ac.Conjugate()*bc.Conjugate())-c1) < eps,
            "+zcac*cbc");
        Assert(Norm(((c=c0)+=z1*ar.Conjugate()*br.Conjugate())-c1) < eps,
            "+zcar*cbr");
        Assert(Norm(((c=c0)+=z1*ar.Conjugate()*bc.Conjugate())-c1) < eps,
            "+zcar*cbc");
        Assert(Norm(((c=c0)+=z1*ac.Conjugate()*br.Conjugate())-c1) < eps,
            "+zcac*cbr");
#endif
#endif
      }
#endif
    }
  }
}

#ifdef TEST_DOUBLE
template void TestMatrixArith_8<double>();
#endif
#ifdef TEST_FLOAT
template void TestMatrixArith_8<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrixArith_8<long double>();
#endif
#ifdef TEST_INT
template void TestMatrixArith_8<int>();
#endif
