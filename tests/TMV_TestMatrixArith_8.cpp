
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>

#include "TMV_TestMatrixArith.h"
#define CT std::complex<T>

template <class T> void TestMatrixArith_8()
{
#if 1
    const int NSIZE = 4;
    const int sizear[NSIZE] = {2,5,77,637};
#else
    const int NSIZE = 15;
    //const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,30,31,32,63,64,65};
    //const int sizear[NSIZE] = {1,2,3,4,5,6,7,8,9,32,128,129,130,1009,3777};
    const int sizear[NSIZE] = {1,2,3,4,5,63,64,65,66,77,111,128,137,256,637};
#endif
    for(int m1=0;m1<NSIZE;m1++) for(int n1=0;n1<NSIZE;n1++) {
        int m = sizear[m1];
        int n = sizear[n1];
        if (showstartdone) std::cout<<"m,n = "<<m<<','<<n<<std::endl;
        else std::cout<<".";
        std::cout.flush();

#if 1
        // Test various adds:
        {
            tmv::Matrix<CT,tmv::ColMajor> c(m,n);
            tmv::Matrix<CT,tmv::ColMajor> c1(m,n);

            tmv::Matrix<CT,tmv::ColMajor> ac(m,n);
            tmv::Matrix<CT,tmv::ColMajor> bc(m,n);
            for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
                ac(i,j) = CT(T(i+5),T(j-2));
            for(int i=0;i<m;i++) for(int j=0;j<n;j++) 
                bc(i,j) = CT(T(2*i-13),T(-3*j+8));
            tmv::Matrix<CT,tmv::RowMajor> ar = ac;
            tmv::Matrix<CT,tmv::RowMajor> br = bc;
            T eps = EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+Norm(ac)+Norm(bc));

            c1 = ac + bc;
            Assert(Equal((c=ar+bc),c1,eps),"ar+bc");
#if (XTEST & 2)
            Assert(Equal((c=ar+br),c1,eps),"ar+br");
            Assert(Equal((c=ac+br),c1,eps),"ac+br");
            c1 = ac + tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=ac+bc.conjugate()),c1,eps),"ac+cbc");
            Assert(Equal((c=ar+br.conjugate()),c1,eps),"ar+cbr");
            Assert(Equal((c=ar+bc.conjugate()),c1,eps),"ar+cbc");
            Assert(Equal((c=ac+br.conjugate()),c1,eps),"ac+cbr");
            c1 = tmv::Matrix<CT>(ac.conjugate()) + bc;
            Assert(Equal((c=ac.conjugate()+bc),c1,eps),"cac+bc");
            Assert(Equal((c=ar.conjugate()+br),c1,eps),"car+br");
            Assert(Equal((c=ar.conjugate()+bc),c1,eps),"car+bc");
            Assert(Equal((c=ac.conjugate()+br),c1,eps),"cac+br");
            c1 = tmv::Matrix<CT>(ac.conjugate()) +
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=ac.conjugate()+bc.conjugate()),c1,eps),"cac+cbc");
            Assert(Equal((c=ar.conjugate()+br.conjugate()),c1,eps),"car+cbr");
            Assert(Equal((c=ar.conjugate()+bc.conjugate()),c1,eps),"car+cbc");
            Assert(Equal((c=ac.conjugate()+br.conjugate()),c1,eps),"cac+cbr");

            T x1(7);
            T x2(-3);
            c1 = x1*ac + x2*bc;
            Assert(Equal((c=x1*ar+x2*br),c1,eps),"xar+xbr");
            Assert(Equal((c=x1*ar+x2*bc),c1,eps),"xar+xbc");
            Assert(Equal((c=x1*ac+x2*br),c1,eps),"xac+xbr");
            c1 = x1*ac + x2*tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=x1*ac+x2*bc.conjugate()),c1,eps),"xac+xcbc");
            Assert(Equal((c=x1*ar+x2*br.conjugate()),c1,eps),"xar+xcbr");
            Assert(Equal((c=x1*ar+x2*bc.conjugate()),c1,eps),"xar+xcbc");
            Assert(Equal((c=x1*ac+x2*br.conjugate()),c1,eps),"xac+xcbr");
            c1 = x1*tmv::Matrix<CT>(ac.conjugate()) + x2*bc;
            Assert(Equal((c=x1*ac.conjugate()+x2*bc),c1,eps),"xcac+xbc");
            Assert(Equal((c=x1*ar.conjugate()+x2*br),c1,eps),"xcar+xbr");
            Assert(Equal((c=x1*ar.conjugate()+x2*bc),c1,eps),"xcar+xbc");
            Assert(Equal((c=x1*ac.conjugate()+x2*br),c1,eps),"xcac+xbr");
            c1 = x1 * tmv::Matrix<CT>(ac.conjugate()) + 
                x2 * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=x1*ac.conjugate()+x2*bc.conjugate()),c1,eps),
                   "xcac+xcbc");
            Assert(Equal((c=x1*ar.conjugate()+x2*br.conjugate()),c1,eps),
                   "xcar+xcbr");
            Assert(Equal((c=x1*ar.conjugate()+x2*bc.conjugate()),c1,eps),
                   "xcar+xcbc");
            Assert(Equal((c=x1*ac.conjugate()+x2*br.conjugate()),c1,eps),
                   "xcac+xcbr");

            CT z1(7,-9);
            CT z2(-3,-8);
            c1 = z1*ac + z2*bc;
            Assert(Equal((c=z1*ar+z2*br),c1,eps),"zar+zbr");
            Assert(Equal((c=z1*ar+z2*bc),c1,eps),"zar+zbc");
            Assert(Equal((c=z1*ac+z2*br),c1,eps),"zac+zbr");
            c1 = z1*ac + z2*tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=z1*ac+z2*bc.conjugate()),c1,eps),"zac+zcbc");
            Assert(Equal((c=z1*ar+z2*br.conjugate()),c1,eps),"zar+zcbr");
            Assert(Equal((c=z1*ar+z2*bc.conjugate()),c1,eps),"zar+zcbc");
            Assert(Equal((c=z1*ac+z2*br.conjugate()),c1,eps),"zac+zcbr");
            c1 = z1*tmv::Matrix<CT>(ac.conjugate()) + z2*bc;
            Assert(Equal((c=z1*ac.conjugate()+z2*bc),c1,eps),"zcac+zbc");
            Assert(Equal((c=z1*ar.conjugate()+z2*br),c1,eps),"zcar+zbr");
            Assert(Equal((c=z1*ar.conjugate()+z2*bc),c1,eps),"zcar+zbc");
            Assert(Equal((c=z1*ac.conjugate()+z2*br),c1,eps),"zcac+zbr");
            c1 = z1 * tmv::Matrix<CT>(ac.conjugate()) +
                z2 * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=z1*ac.conjugate()+z2*bc.conjugate()),c1,eps),
                   "zcac+zcbc");
            Assert(Equal((c=z1*ar.conjugate()+z2*br.conjugate()),c1,eps),
                   "zcar+zcbr");
            Assert(Equal((c=z1*ar.conjugate()+z2*bc.conjugate()),c1,eps),
                   "zcar+zcbc");
            Assert(Equal((c=z1*ac.conjugate()+z2*br.conjugate()),c1,eps),
                   "zcac+zcbr");
#endif
        }
#endif

#if 1
        // Now test multiplies
        for(int k1=0;k1<NSIZE;k1++) {
            tmv::Matrix<CT,tmv::ColMajor> c(m,n);
            tmv::Matrix<CT,tmv::ColMajor> c1(m,n);

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
            T eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+Norm(ac)*Norm(bc));

            c1 = ac * bc;
            if (showacc) {
                std::cout<<"a = "<<ar<<std::endl;
                std::cout<<"b = "<<br<<std::endl;
                std::cout<<"c (CC) = "<<c1<<std::endl;
                std::cout<<"c (RC) = "<<(c=ar*bc);
                std::cout<<Norm((c=ar*bc)-c1)<<std::endl;
                std::cout<<"c (CR) = "<<(c=ac*br);
                std::cout<<Norm((c=ac*br)-c1)<<std::endl;
                std::cout<<"c (RR) = "<<(c=ar*br);
                std::cout<<Norm((c=ar*br)-c1)<<std::endl;
                std::cout<<"cf eps = "<<eps<<std::endl;
            }
            Assert(Equal((c=ar*bc),c1,eps),"ar*bc");
#if (XTEST & 2)
            Assert(Equal((c=ar*br),c1,eps),"ar*br");
            Assert(Equal((c=ac*br),c1,eps),"ac*br");
            c1 = ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=ac*bc.conjugate()),c1,eps),"ac*cbc");
            Assert(Equal((c=ar*br.conjugate()),c1,eps),"ar*cbr");
            Assert(Equal((c=ar*bc.conjugate()),c1,eps),"ar*cbc");
            Assert(Equal((c=ac*br.conjugate()),c1,eps),"ac*cbr");
            c1 = tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal((c=ac.conjugate()*bc),c1,eps),"cac*bc");
            Assert(Equal((c=ar.conjugate()*br),c1,eps),"car*br");
            Assert(Equal((c=ar.conjugate()*bc),c1,eps),"car*bc");
            Assert(Equal((c=ac.conjugate()*br),c1,eps),"cac*br");
            c1 = tmv::Matrix<CT>(ac.conjugate()) *
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=ac.conjugate()*bc.conjugate()),c1,eps),"cac*cbc");
            Assert(Equal((c=ar.conjugate()*br.conjugate()),c1,eps),"car*cbr");
            Assert(Equal((c=ar.conjugate()*bc.conjugate()),c1,eps),"car*cbc");
            Assert(Equal((c=ac.conjugate()*br.conjugate()),c1,eps),"cac*cbr");

            tmv::Matrix<CT> c0 = c1;
            c1 = c0 + ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (Norm(c0)+Norm(ac)*Norm(bc));
            Assert(Equal(((c=c0)+=ar*br),c1,eps),"+ar*br");
            Assert(Equal(((c=c0)+=ar*bc),c1,eps),"+ar*bc");
            Assert(Equal(((c=c0)+=ac*br),c1,eps),"+ac*br");
            c1 = c0 + ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=ac*bc.conjugate()),c1,eps),"+ac*cbc");
            Assert(Equal(((c=c0)+=ar*br.conjugate()),c1,eps),"+ar*cbr");
            Assert(Equal(((c=c0)+=ar*bc.conjugate()),c1,eps),"+ar*cbc");
            Assert(Equal(((c=c0)+=ac*br.conjugate()),c1,eps),"+ac*cbr");
            c1 = c0 + tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal(((c=c0)+=ac.conjugate()*bc),c1,eps),"+cac*bc");
            Assert(Equal(((c=c0)+=ar.conjugate()*br),c1,eps),"+car*br");
            Assert(Equal(((c=c0)+=ar.conjugate()*bc),c1,eps),"+car*bc");
            Assert(Equal(((c=c0)+=ac.conjugate()*br),c1,eps),"+cac*br");
            c1 = c0 + tmv::Matrix<CT>(ac.conjugate()) *
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=ac.conjugate()*bc.conjugate()),c1,eps),
                   "+cac*cbc");
            Assert(Equal(((c=c0)+=ar.conjugate()*br.conjugate()),c1,eps),
                   "+car*cbr");
            Assert(Equal(((c=c0)+=ar.conjugate()*bc.conjugate()),c1,eps),
                   "+car*cbc");
            Assert(Equal(((c=c0)+=ac.conjugate()*br.conjugate()),c1,eps),
                   "+cac*cbr");

            T x1(7);
            c1 = x1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+x1*Norm(ac)*Norm(bc));
            Assert(Equal((c=x1*ar*br),c1,eps),"xar*br");
            Assert(Equal((c=x1*ar*bc),c1,eps),"xar*bc");
            Assert(Equal((c=x1*ac*br),c1,eps),"xac*br");
            c1 = x1*ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=x1*ac*bc.conjugate()),c1,eps),"xac*cbc");
            Assert(Equal((c=x1*ar*br.conjugate()),c1,eps),"xar*cbr");
            Assert(Equal((c=x1*ar*bc.conjugate()),c1,eps),"xar*cbc");
            Assert(Equal((c=x1*ac*br.conjugate()),c1,eps),"xac*cbr");
            c1 = x1*tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal((c=x1*ac.conjugate()*bc),c1,eps),"xcac*bc");
            Assert(Equal((c=x1*ar.conjugate()*br),c1,eps),"xcar*br");
            Assert(Equal((c=x1*ar.conjugate()*bc),c1,eps),"xcar*bc");
            Assert(Equal((c=x1*ac.conjugate()*br),c1,eps),"xcac*br");
            c1 = x1 * tmv::Matrix<CT>(ac.conjugate()) *
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=x1*ac.conjugate()*bc.conjugate()),c1,eps),
                   "xcac*cbc");
            Assert(Equal((c=x1*ar.conjugate()*br.conjugate()),c1,eps),
                   "xcar*cbr");
            Assert(Equal((c=x1*ar.conjugate()*bc.conjugate()),c1,eps),
                   "xcar*cbc");
            Assert(Equal((c=x1*ac.conjugate()*br.conjugate()),c1,eps),
                   "xcac*cbr");

            c1 = c0 + x1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (Norm(c0)+x1*Norm(ac)*Norm(bc));
            Assert(Equal(((c=c0)+=x1*ar*br),c1,eps),"+xar*br");
            Assert(Equal(((c=c0)+=x1*ar*bc),c1,eps),"+xar*bc");
            Assert(Equal(((c=c0)+=x1*ac*br),c1,eps),"+xac*br");
            c1 = c0 + x1*ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=x1*ac*bc.conjugate()),c1,eps),"+xac*cbc");
            Assert(Equal(((c=c0)+=x1*ar*br.conjugate()),c1,eps),"+xar*cbr");
            Assert(Equal(((c=c0)+=x1*ar*bc.conjugate()),c1,eps),"+xar*cbc");
            Assert(Equal(((c=c0)+=x1*ac*br.conjugate()),c1,eps),"+xac*cbr");
            c1 = c0 + x1*tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal(((c=c0)+=x1*ac.conjugate()*bc),c1,eps),"+xcac*bc");
            Assert(Equal(((c=c0)+=x1*ar.conjugate()*br),c1,eps),"+xcar*br");
            Assert(Equal(((c=c0)+=x1*ar.conjugate()*bc),c1,eps),"+xcar*bc");
            Assert(Equal(((c=c0)+=x1*ac.conjugate()*br),c1,eps),"+xcac*br");
            c1 = c0 + x1 * tmv::Matrix<CT>(ac.conjugate()) *
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=x1*ac.conjugate()*bc.conjugate()),c1,eps),
                   "+xcac*cbc");
            Assert(Equal(((c=c0)+=x1*ar.conjugate()*br.conjugate()),c1,eps),
                   "+xcar*cbr");
            Assert(Equal(((c=c0)+=x1*ar.conjugate()*bc.conjugate()),c1,eps),
                   "+xcar*cbc");
            Assert(Equal(((c=c0)+=x1*ac.conjugate()*br.conjugate()),c1,eps),
                   "+xcac*cbr");

            CT z1(7,-9);
            c1 = z1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+tmv::TMV_ABS2(z1)*Norm(ac)*Norm(bc));
            Assert(Equal((c=z1*ar*br),c1,eps),"zar*br");
            Assert(Equal((c=z1*ar*bc),c1,eps),"zar*bc");
            Assert(Equal((c=z1*ac*br),c1,eps),"zac*br");
            c1 = z1*ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=z1*ac*bc.conjugate()),c1,eps),"zac*cbc");
            Assert(Equal((c=z1*ar*br.conjugate()),c1,eps),"zar*cbr");
            Assert(Equal((c=z1*ar*bc.conjugate()),c1,eps),"zar*cbc");
            Assert(Equal((c=z1*ac*br.conjugate()),c1,eps),"zac*cbr");
            c1 = z1*tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal((c=z1*ac.conjugate()*bc),c1,eps),"zcac*bc");
            Assert(Equal((c=z1*ar.conjugate()*br),c1,eps),"zcar*br");
            Assert(Equal((c=z1*ar.conjugate()*bc),c1,eps),"zcar*bc");
            Assert(Equal((c=z1*ac.conjugate()*br),c1,eps),"zcac*br");
            c1 = z1 * tmv::Matrix<CT>(ac.conjugate()) * 
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal((c=z1*ac.conjugate()*bc.conjugate()),c1,eps),
                   "zcac*cbc");
            Assert(Equal((c=z1*ar.conjugate()*br.conjugate()),c1,eps),
                   "zcar*cbr");
            Assert(Equal((c=z1*ar.conjugate()*bc.conjugate()),c1,eps),
                   "zcar*cbc");
            Assert(Equal((c=z1*ac.conjugate()*br.conjugate()),c1,eps),
                   "zcac*cbr");

            c1 = c0 + z1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (Norm(c0)+tmv::TMV_ABS2(z1)*Norm(ac)*Norm(bc));
            Assert(Equal(((c=c0)+=z1*ar*br),c1,eps),"+zar*br");
            Assert(Equal(((c=c0)+=z1*ar*bc),c1,eps),"+zar*bc");
            Assert(Equal(((c=c0)+=z1*ac*br),c1,eps),"+zac*br");
            c1 = c0 + z1*ac * tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=z1*ac*bc.conjugate()),c1,eps),"+zac*cbc");
            Assert(Equal(((c=c0)+=z1*ar*br.conjugate()),c1,eps),"+zar*cbr");
            Assert(Equal(((c=c0)+=z1*ar*bc.conjugate()),c1,eps),"+zar*cbc");
            Assert(Equal(((c=c0)+=z1*ac*br.conjugate()),c1,eps),"+zac*cbr");
            c1 = c0 + z1*tmv::Matrix<CT>(ac.conjugate()) * bc;
            Assert(Equal(((c=c0)+=z1*ac.conjugate()*bc),c1,eps),"+zcac*bc");
            Assert(Equal(((c=c0)+=z1*ar.conjugate()*br),c1,eps),"+zcar*br");
            Assert(Equal(((c=c0)+=z1*ar.conjugate()*bc),c1,eps),"+zcar*bc");
            Assert(Equal(((c=c0)+=z1*ac.conjugate()*br),c1,eps),"+zcac*br");
            c1 = c0 + z1 * tmv::Matrix<CT>(ac.conjugate()) *
                tmv::Matrix<CT>(bc.conjugate());
            Assert(Equal(((c=c0)+=z1*ac.conjugate()*bc.conjugate()),c1,eps),
                   "+zcac*cbc");
            Assert(Equal(((c=c0)+=z1*ar.conjugate()*br.conjugate()),c1,eps),
                   "+zcar*cbr");
            Assert(Equal(((c=c0)+=z1*ar.conjugate()*bc.conjugate()),c1,eps),
                   "+zcar*cbc");
            Assert(Equal(((c=c0)+=z1*ac.conjugate()*br.conjugate()),c1,eps),
                   "+zcac*cbr");
#endif
        }
#endif

#if 1 
        // Do real version
        for(int k1=0;k1<NSIZE;k1++) {
            tmv::Matrix<T,tmv::ColMajor> c(m,n);
            tmv::Matrix<T,tmv::ColMajor> c1(m,n);

            int k = sizear[k1];
            if (showstartdone)
                std::cout<<"(real) k = "<<k<<std::endl;
            tmv::Matrix<T,tmv::ColMajor> ac(m,k);
            tmv::Matrix<T,tmv::ColMajor> bc(k,n);
            for(int i=0;i<m;i++) for(int j=0;j<k;j++) 
                ac(i,j) = T(i+j+5);
            for(int i=0;i<k;i++) for(int j=0;j<n;j++) 
                bc(i,j) = T(2*i-3*j+13);
            tmv::Matrix<T,tmv::RowMajor> ar = ac;
            tmv::Matrix<T,tmv::RowMajor> br = bc;
            T eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+Norm(ac)*Norm(bc));

            c1 = ac * bc;
            if (showacc)
            {
                std::cout<<"a = "<<ar<<std::endl;
                std::cout<<"b = "<<br<<std::endl;
                std::cout<<"c (CC) = "<<c1<<std::endl;
                std::cout<<"c (RC) = "<<(c=ar*bc);
                std::cout<<Norm((c=ar*bc)-c1)<<std::endl;
                std::cout<<"c (CR) = "<<(c=ac*br);
                std::cout<<Norm((c=ac*br)-c1)<<std::endl;
                std::cout<<"c (RR) = "<<(c=ar*br);
                std::cout<<Norm((c=ar*br)-c1)<<std::endl;
                std::cout<<"cf eps = "<<eps<<std::endl;
            }
            Assert(Equal((c=ar*bc),c1,eps),"ar*bc");
#if (XTEST & 2)
            Assert(Equal((c=ar*br),c1,eps),"ar*br");
            Assert(Equal((c=ac*br),c1,eps),"ac*br");

            tmv::Matrix<T> c0 = c1;
            c1 = c0 + ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (Norm(c0)+Norm(ac)*Norm(bc));
            Assert(Equal(((c=c0)+=ar*br),c1,eps),"+ar*br");
            Assert(Equal(((c=c0)+=ar*bc),c1,eps),"+ar*bc");
            Assert(Equal(((c=c0)+=ac*br),c1,eps),"+ac*br");

            T x1(7);
            c1 = x1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (T(1)+x1*Norm(ac)*Norm(bc));
            Assert(Equal((c=x1*ar*br),c1,eps),"xar*br");
            Assert(Equal((c=x1*ar*bc),c1,eps),"xar*bc");
            Assert(Equal((c=x1*ac*br),c1,eps),"xac*br");

            c1 = c0 + x1*ac * bc;
            eps = T(10) * EPS;
            if (!std::numeric_limits<T>::is_integer)
                eps *= (Norm(c0)+x1*Norm(ac)*Norm(bc));
            Assert(Equal(((c=c0)+=x1*ar*br),c1,eps),"+xar*br");
            Assert(Equal(((c=c0)+=x1*ar*bc),c1,eps),"+xar*bc");
            Assert(Equal(((c=c0)+=x1*ac*br),c1,eps),"+xac*br");
#endif
        }
#endif
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
