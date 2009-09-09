
#ifndef TMV_BLAS_H
#define TMV_BLAS_H

#ifdef MKL

#ifndef NOBLAS
#define BLAS
#endif

#ifndef NOLAP
#define LAP   
#define ALAP  // ALAP for ATLAS or LAPACK
#endif

#include <mkl.h>

namespace tmv {

#ifdef LAP
  inline MKL_Complex16* LAP_Complex(complex<double>* ptr)
  { return reinterpret_cast<MKL_Complex16*>(ptr); }
  inline MKL_Complex16* LAP_Complex(const complex<double>* ptr)
  { return LAP_Complex(const_cast<complex<double>*>(ptr)); }
  inline MKL_Complex16* LAP_Complex(double* ptr)
  { return reinterpret_cast<MKL_Complex16*>(ptr); }
  inline MKL_Complex16* LAP_Complex(const double* ptr)
  { return LAP_Complex(const_cast<double*>(ptr)); }
  inline MKL_Complex8* LAP_Complex(complex<float>* ptr)
  { return reinterpret_cast<MKL_Complex8*>(ptr); }
  inline MKL_Complex8* LAP_Complex(const complex<float>* ptr)
  { return LAP_Complex(const_cast<complex<float>*>(ptr)); }
  inline MKL_Complex8* LAP_Complex(float* ptr)
  { return reinterpret_cast<MKL_Complex8*>(ptr); }
  inline MKL_Complex8* LAP_Complex(const float* ptr)
  { return LAP_Complex(const_cast<float*>(ptr)); }
#endif
  
}

#endif // MKL

#ifdef ATLAS

#ifndef NOBLAS
#define BLAS
#endif

#ifndef NOLAP
#define ALAP  // ALAP for ATLAS or LAPACK
#define CLAP  // indicates to use the clapack names
#endif

extern "C" {
#include "cblas.h"
#include "clapack.h"
}

#endif // ATLAS

#ifdef ALAP
namespace tmv {

  const int LAP_BLOCKSIZE = 64;

  inline int* LAP_IPiv(const int n)
  {
    static int* ipiv = 0;
    static int currn = 0;
    if (n > currn) {
      if (ipiv) delete[] ipiv;
      currn = n;
      ipiv = new int[n];
#ifdef CLAP
      for(int i=0;i<n;++i) ipiv[i] = i;
#else
      for(int i=0;i<n;++i) ipiv[i] = i+1;
#endif
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
      tmv_error("LAPACK error: info < 0 in call to ",fn);
    }
  }

  inline void LAP_Results(const int info, const int lwork_opt, 
      const int m, const int n, const int lwork, const char* fn)
  {
    if (info < 0) {
      tmv_error("LAPACK error: info < 0 in call to ",fn);
    }
#ifdef TMVDEBUG
    if (lwork_opt > lwork) {
      cout<<"LAPACK function "<<fn<<" requested more workspace than provided\n";
      cout<<"for matrix with m,n = "<<m<<','<<n<<endl;
      cout<<"Given: "<<lwork<<", requested "<<lwork_opt<<endl;
    }
#endif
  }

}

#endif // LAP

#endif // TMV_BLAS_H
