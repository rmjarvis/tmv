
#ifndef JMB_BLAS_H
#define JMB_BLAS_H

#ifndef NOBLAS
#define BLAS

#include "/usr/local/include/mkl.h"

#ifndef NOLAP
#define LAP

namespace tmv {

  const int LAP_BLOCKSIZE = 32;

  inline MKL_Complex16* LAP_Complex(complex<double>* ptr)
  { return reinterpret_cast<MKL_Complex16*>(ptr); }

  inline MKL_Complex16* LAP_Complex(const complex<double>* ptr)
  { return LAP_Complex(const_cast<complex<double>*>(ptr)); }

  inline MKL_Complex8* LAP_Complex(complex<float>* ptr)
  { return reinterpret_cast<MKL_Complex8*>(ptr); }

  inline MKL_Complex8* LAP_Complex(const complex<float>* ptr)
  { return LAP_Complex(const_cast<complex<float>*>(ptr)); }

  inline int* LAP_IPiv(const int n)
  {
    static int* ipiv = 0;
    static int currn = 0;
    if (n > currn) {
      if (ipiv) delete[] ipiv;
      currn = n;
      ipiv = new int[n];
      for(int i=0;i<n;++i) ipiv[i] = i+1;
    }
    return ipiv;
  }

  inline int* LAP_IWork(const int n)
  {
    static int* work = 0;
    static int currn = 0;
    if (n > currn) {
      if (work) delete [] work;
      work = new int[n];
      currn = n;
    }
    return work;
  }

  inline double* LAP_DWork(const int n)
  {
    static double* work = 0;
    static int currn = 0;
    if (n > currn) {
      if (work) delete [] work;
      work = new double[n];
      currn = n;
    }
    return work;
  }

  inline complex<double>* LAP_ZWork(const int n)
  {
    static complex<double>* work = 0;
    static int currn = 0;
    if (n > currn) {
      if (work) delete [] work;
      work = new complex<double>[n];
      currn = n;
    }
    return work;
  }

  inline float* LAP_SWork(const int n)
  {
    static float* work = 0;
    static int currn = 0;
    if (n > currn) {
      if (work) delete [] work;
      work = new float[n];
      currn = n;
    }
    return work;
  }

  inline complex<float>* LAP_CWork(const int n)
  {
    static complex<float>* work = 0;
    static int currn = 0;
    if (n > currn) {
      if (work) delete [] work;
      work = new complex<float>[n];
      currn = n;
    }
    return work;
  }

  inline void LAP_Results(const int info, const char* fn)
  {
    if (info < 0) {
      cout<<fn<<" returned info = "<<info<<endl;
      tmv_error("LaPack error: info < 0");
    }
  }

  inline void LAP_Results(const int info, const int lwork_opt, 
      const int m, const int n, const int lwork, const char* fn)
  {
    if (info < 0) {
      cout<<fn<<" returned info = "<<info<<endl;
      tmv_error("LaPack error: info < 0");
    }
#ifdef TMVDEBUG
    if (lwork_opt > lwork) {
      cout<<fn<<" requested more workspace than provided\n";
      cout<<"for matrix with m,n = "<<m<<','<<n<<endl;
      cout<<"Given: "<<lwork<<", requested "<<lwork_opt<<endl;
    }
#endif
  }

}

#endif // NOLAP

#endif // NOBLAS

#endif
