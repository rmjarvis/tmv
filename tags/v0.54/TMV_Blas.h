
// This file sets up the various things needed to use BLAs and/or LAPACK
//
// There are 3 implementations which I have set up so far:
//   -- The Intel Math Kernel Library (Use -DMKL when compiling)
//   -- The AMD Core Math Library (Use -DACML when compiling)
//   -- ATLAS (Use -DATLAS when compiling)
//
// If you have a different BLAS or LAPACK library you may need
// to add your own section with the definitions specific to your
// library.  You should hopefully be able to use the existing options
// to guide you in how to set it up.  The trickiest bit is converting
// between C++'s complex class and the various constructs used by 
// different libraries for the complex data.
//
// The first part of this file does the things which are specific
// to each of the above options.
// There are several things which may be defined in these sections,
// which the rest of the library checks for when compiling.
// So here are what they mean:
//
//   BLAS -- Use Blas calls.
//   CBLAS -- Blas calls should use cblas_* calling convention.
//   BLAS_ -- Append underscore to Blas names.
//   BLASPTR -- Pass all arguments to Blas calls by pointer, not reference.
//   BLASSTRLEN -- Last argument of Blas calls need the length of char array.
//   BLASAMIN -- Has routines idamin, isamin
//
//   LAP -- Use LAPack calls (see below for more specific subsets here)
//   CLAP -- LAPack calls should use clapack_* calling convention.
//           Note: clapack_* calls have only been implemented for ALAP so far.
//   LAP_ -- Append underscore to LAPack names.
//   LAPPTR -- Pass all arguments to LAPack calls by pointer, not reference.
//   LAPNOWORK -- Skip the workspace arguments in LAPack calls.
//   LAPSTRLEN -- Last argument of LAPack calls need the length of char array.
//
// There are 4 subcategories of LAPack calls:
//   LAP = the regular LAPACK set (always includes ALAP below)
//         eg. *geqrf, *gbtrf, *orglq, etc.
//   ALAP = the minimal set of routines that ATLAS and regular LAPACK include:
//          eg. *getrf, *getri, *potrf, *trtri, etc.
//   ELAP = Extra LAPACK functions included by the standard CLAPACK
//          distribution, but not some other distributions, eg. ACML.
//          *lagtm, (z,c)rot, *lacrm, *lacpy, (z,c)sy(mv,r)
//          (also includes AELAP below)
//   AELAP = An extra routine that ATLAS has that the minimal LAPACK does not:
//           *lauum
//   XLAP = The extended "auxilliary and utility" LAPACK routines included by 
//          the INTEL MKL, but not some other LAPACK distributions
//          (specifically, not the standard CLAPACK distribution).
//          (dz,cs)sum1, (z,c)lacgv, *lan(ge,gb,gt,tr,sy,he)
//
// Each LAPack version below needs to #define the appropriate ones of these
// given which routines are included in the distribution.
//
// For the non-(CBLAS / CLAP) versions, there is generally an 
// implementation-specific structure for complex types.
// I define a conversion from a std::complex pointer to these types
// with LAPP(&z) and BLASP(&z).
// (There are also LAPP's for real types to get rid of consts if necessary.)
//
// The BLAS routines zdotc, zdotu, cdotc, and cdotu vary 
// somewhat in their implementation, since they returna a complex value.
// The CBLAS converntion is to append _sub to the names and put the 
// result as the last argument.  Others often put the result as the
// first argument.  Or they many return their style of complex value
// which needs to be converted to a regular complex number.
// All of these options are possible. 
//
// For your implementation, you should define one of the following:
//   BLASZDOTLAST
//   BLASZDOTFIRST
//   BLASZDOTRETURN
// to specify which of these conventions your implementation uses.
//
// You should also define
//   BLASZDOTSUB
// to be whatever is appended to the name of these function names, eg.
// #define BLASZDOTSUB _sub
//
// And, if you specify BLASZDOTRETURN, you need to define 
//   BLASCValue(x) 
// to be a functions which takes your implementation's complex type
// and returns a regular complex value.  (See ACML for an example.)
//
// The LAPACK info variable also has several different implementations.
// The standard method is as an int* parameter at the end of each
// LAPACK routine.  This is the default if nothing else is specified.
// However, CLAPACK and some others return info as the result of
// each function.  To specify this usage, define
//   LAPINFORETURN
// Also, I don't know if any that do this, but if any implementations 
// use a reference rather than a pointer for the info variable, you
// can define
//   LAPINFOREF
//
// The second part of the file uses all of these defines to define the
// actual syntax that is used in the library code to make the correct
// BLAS and LAPACK calls.  It alo sets up some functions I use to
// handle the LAPack temporary work space efficiently, and also
// to convert their error handling to one of my exceptions.
// You probably won't need to change anything in Part 2.

#ifndef TMV_BLAS_H
#define TMV_BLAS_H

//
//
// PART 1
//
//

// We need to make sure 


//
// *** MKL ***
//
#ifdef MKL

#ifdef NOBLAS
#define NOLAP
#else

#define BLAS

#include "mkl.h"

#define CBLAS // Use cblas_ calling convention

#ifndef NOLAP
#ifndef CLAPACK

#define LAP
#define ELAP
#define XLAP
#define LAPPTR

#include <complex>

namespace tmv {


  template <class T> inline T* MKL_ConvertP(T* x) { return x; }
  template <class T> inline T* MKL_ConvertP(const T* x) 
  { return MKL_ConvertP(const_cast<T*>(x)); }
  
  inline MKL_Complex16* MKL_ConvertP(std::complex<double>* ptr)
  { return reinterpret_cast<MKL_Complex16*>(ptr); }
  inline MKL_Complex8* MKL_ConvertP(std::complex<float>* ptr)
  { return reinterpret_cast<MKL_Complex8*>(ptr); }
  inline MKL_Complex16* MKL_ConvertP(const std::complex<double>* ptr)
  { return MKL_ConvertP(const_cast<std::complex<double>*>(ptr)); }
  inline MKL_Complex8* MKL_ConvertP(const std::complex<float>* ptr)
  { return MKL_ConvertP(const_cast<std::complex<float>*>(ptr)); }
}

#define LAPP(x) MKL_ConvertP(x)

#endif // !CLAPACK
#endif // LAP

#endif // BLAS

#endif // MKL

// 
// *** ACML ***
//
#ifdef ACML

#ifdef NOBLAS
#define NOLAP
#else

#define BLAS

#include "acml.h"
#include <complex> 
// complex needs to be included _after_ acml, otherwise there are problems.

/* Works with either of these sets of defines: */
#define BLASZDOTRETURN

/*
#define BLAS_
#define BLASPTR
#define BLASSTRLEN
#define BLASZDOTFIRST
#define BLASZDOTSUB sub  // This one is optional actually
*/

namespace tmv {

  template <class T> inline T* ACML_ConvertP(T* x) { return x; }
  template <class T> inline T* ACML_ConvertP(const T* x) 
  { return const_cast<T*>(x); }

  inline ::doublecomplex* ACML_ConvertP(std::complex<double>* ptr)
  { return reinterpret_cast< ::doublecomplex*>(ptr); }
  inline ::complex* ACML_ConvertP(std::complex<float>* ptr)
  { return reinterpret_cast< ::complex*>(ptr); }
  inline ::doublecomplex* ACML_ConvertP(const std::complex<double>* ptr)
  { return ACML_ConvertP(const_cast<std::complex<double>*>(ptr)); }
  inline ::complex* ACML_ConvertP(const std::complex<float>* ptr)
  { return ACML_ConvertP(const_cast<std::complex<float>*>(ptr)); }

  inline std::complex<double> ACML_ConvertToComplex(::doublecomplex z)
  { return std::complex<double>(z.real,z.imag); }
  inline std::complex<float> ACML_ConvertToComplex(::complex z)
  { return std::complex<float>(z.real,z.imag); }

}

#define BLASP(x) ACML_ConvertP(x)
#define BLASCValue(x) ACML_ConvertToComplex(x)

#ifndef NOLAP
#ifndef CLAPACK

#define LAP

/* Works with one or the other of these sets of defines: */
#define LAPNOWORK

/*
#define LAP_
#define LAPPTR
#define LAPSTRLEN
*/

#define LAPP(x) ACML_ConvertP(x)

#endif // !CLAPACK
#endif // LAP
#endif // BLAS
#endif // ACML


// 
// *** ATLAS ***
//
#ifdef ATLAS

#ifdef NOBLAS
#define NOLAP
#else

#define BLAS
#define CBLAS 

extern "C" {
#include "cblas.h"

#ifndef NOLAP
#ifndef CLAPACK

#define ALAP 
#define AELAP 
#define CLAP
#define LAPINFORETURN
#define LAPNOWORK

#include "clapack.h"

#endif // !CLAPACK
#endif // LAP

}

#endif // BLAS

#endif // ATLAS


//
// *** CLAPACK ***
//
#ifdef CLAPACK

#ifdef NOBLAS
#define NOLAP
#else

#ifndef NOLAP

#define LAP 
#define ELAP 
#define LAP_
#define LAPPTR

extern "C" {
#include "f2c.h"
#include "clapack.h"
#undef abs
#undef min
#undef max
}
#include <complex> 

namespace tmv {

  template <class T> inline T* CLAPACK_ConvertP(T* x) { return x; }
  template <class T> inline T* CLAPACK_ConvertP(const T* x) 
  { return const_cast<T*>(x); }

  inline ::doublecomplex* CLAPACK_ConvertP(std::complex<double>* ptr)
  { return reinterpret_cast< ::doublecomplex*>(ptr); }
  inline ::complex* CLAPACK_ConvertP(std::complex<float>* ptr)
  { return reinterpret_cast< ::complex*>(ptr); }
  inline ::doublecomplex* CLAPACK_ConvertP(const std::complex<double>* ptr)
  { return CLAPACK_ConvertP(const_cast<std::complex<double>*>(ptr)); }
  inline ::complex* CLAPACK_ConvertP(const std::complex<float>* ptr)
  { return CLAPACK_ConvertP(const_cast<std::complex<float>*>(ptr)); }

  inline ::integer* CLAPACK_ConvertP(int* ptr)
  { return reinterpret_cast< ::integer*>(ptr); }
  inline ::integer* CLAPACK_ConvertP(const int* ptr)
  { return CLAPACK_ConvertP(const_cast<int*>(ptr)); }

}

#define LAPP(x) CLAPACK_ConvertP(x)

#endif // LAP

#endif // BLAS

#endif // CLAPACK

//
//
// PART 2
//
//

// LAP always implies the ALAP subset
#ifdef LAP
#define ALAP
#endif // LAP

#ifdef ELAP
#define AELAP
#endif // ELAP

namespace tmv {
  // These are defined in TMV_Vector.cpp
  extern int Lap_info; 
  extern char Blas_ch_N;
  extern char Blas_ch_C;
  extern char Blas_ch_T;
  extern char Blas_ch_L;
  extern char Blas_ch_R;
  extern char Blas_ch_U;
}

#ifdef BLAS

#ifdef CBLAS
#define BLASNAME(x) cblas_ ## x
#elif defined BLAS_
#define BLASNAME(x) x ## _
#else
#define BLASNAME(x) x
#endif
#define BLASNAME1(x) BLASNAME(x)

#ifdef BLASPTR
#define BLASV(x) BLASP(&x)
#else
#define BLASV(x) x
#endif

#ifndef BLASP
#define BLASP(x) x
#endif

#ifdef CBLAS
#define BLASCM CblasColMajor,
#define BLASRM CblasRowMajor,
#define BLASCH_NT CblasNoTrans
#define BLASCH_CT CblasConjTrans
#define BLASCH_T CblasTrans
#define BLASCH_L CblasLeft
#define BLASCH_R CblasRight
#define BLASCH_U CblasUnit
#define BLASCH_NU CblasNonUnit
#define BLASCH_LO CblasLower
#define BLASCH_UP CblasUpper
#else
#define BLASCM
#define BLASCH_NT BLASV(Blas_ch_N)
#define BLASCH_CT BLASV(Blas_ch_C)
#define BLASCH_T BLASV(Blas_ch_T)
#define BLASCH_L BLASV(Blas_ch_L)
#define BLASCH_R BLASV(Blas_ch_R)
#define BLASCH_U BLASV(Blas_ch_U)
#define BLASCH_NU BLASV(Blas_ch_N)
#define BLASCH_LO BLASV(Blas_ch_L)
#define BLASCH_UP BLASV(Blas_ch_U)
#endif

#ifdef BLASSTRLEN
#define BLAS1 ,1
#else
#define BLAS1
#endif

#ifdef CBLAS
#ifndef BLASZDOTLAST
#define BLASZDOTLAST
#endif
#ifndef BLASZDOTSUB
#define BLASZDOTSUB _sub
#endif
#endif

#ifdef BLASZDOTLAST
#define BLASZDOT2(x) ,x
#else
#define BLASZDOT2(x)
#endif

#ifdef BLASZDOTFIRST
#define BLASZDOT1(x) x,
#else
#define BLASZDOT1(x)
#endif

#ifdef BLASZDOTRETURN
#define BLASZDOTSET(x,y) x = BLASCValue(y)
#else
#define BLASZDOTSET(x,y) y
#endif

#ifdef BLASZDOTSUB
#define BLASZDOTNAME(x) BLASNAME1(BLASZDOTNAME1(x,BLASZDOTSUB))
#define BLASZDOTNAME1(x,y) BLASZDOTNAME2(x,y)
#define BLASZDOTNAME2(x,y) x ## y
#else
#define BLASZDOTNAME(x) BLASNAME(x)
#endif

#endif // BLAS

#ifdef ALAP

#ifdef CLAP
#define LAPNAMEX(x) clapack_ ## x
#elif defined LAP_
#define LAPNAMEX(x) x ## _
#else
#define LAPNAMEX(x) x
#endif

#ifdef LAPINFORETURN
#define LAPNAME(x) tmv::Lap_info = LAPNAMEX(x)
#else
#define LAPNAME(x) LAPNAMEX(x)
#endif

#ifdef LAPPTR
#define LAPV(x) LAPP(&x)
#else
#define LAPV(x) x
#endif

#ifndef LAPP
#define LAPP(x) x
#endif

#ifdef CLAP
#define LAPCM CblasColMajor,
#define LAPCH_NT CblasNoTrans
#define LAPCH_CT CblasConjTrans
#define LAPCH_T CblasTrans
#define LAPCH_L CblasLeft
#define LAPCH_R CblasRight
#define LAPCH_U CblasUnit
#define LAPCH_NU CblasNonUnit
#define LAPCH_LO CblasLower
#define LAPCH_UP CblasUpper
#else
#define LAPCM
#define LAPCH_NT LAPV(Blas_ch_N)
#define LAPCH_CT LAPV(Blas_ch_C)
#define LAPCH_T LAPV(Blas_ch_T)
#define LAPCH_L LAPV(Blas_ch_L)
#define LAPCH_R LAPV(Blas_ch_R)
#define LAPCH_U LAPV(Blas_ch_U)
#define LAPCH_NU LAPV(Blas_ch_N)
#define LAPCH_LO LAPV(Blas_ch_L)
#define LAPCH_UP LAPV(Blas_ch_U)
#endif

#ifdef LAPSTRLEN
#define LAP1 ,1
#else
#define LAP1
#endif

#ifdef LAPNOWORK
#define LAPWK(x)
#define LAPVWK(x)
#else
#define LAPWK(x) ,LAPP(x)
#define LAPVWK(x) ,LAPV(x)
#endif

#ifdef LAPINFORETURN
#define LAPINFO
#elif defined LAPINFOREF
#define LAPINFO ,tmv::Lap_info
#else
#define LAPINFO ,LAPP(&tmv::Lap_info)
#endif

#ifdef CLAP
#define LAPMINUS1
#define LAPPLUS1
#else
#define LAPMINUS1 -1
#define LAPPLUS1 +1
#endif

#include "TMV_Base.h"

namespace tmv {

  inline int* LAP_IPiv(const int n)
  {
    static tmv::auto_array<int> ipiv(0);
    static int currn = 0;
    if (n > currn) {
      ipiv.reset(new int[n]);
      currn = n;
#ifdef CLAP
      for(int i=0;i<n;++i) ipiv.get()[i] = i;
#else
      for(int i=0;i<n;++i) ipiv.get()[i] = i+1;
#endif
    }
    return ipiv.get();
  }

#ifndef LAPNOWORK

  const int LAP_BLOCKSIZE = 64;

  inline int* LAP_IWork(const int n)
  {
    static tmv::auto_array<int> work(0);
    static int currn = 0;
    if (n > currn) {
      work.reset(new int[n]);
      currn = n;
    }
    return work.get();
  }

  inline double* LAP_DWork(const int n)
  {
    static tmv::auto_array<double> work(0);
    static int currn = 0;
    if (n > currn) {
      work.reset(new double[n]);
      currn = n;
    }
    return work.get();
  }

  inline std::complex<double>* LAP_ZWork(const int n)
  {
    static tmv::auto_array<std::complex<double> > work(0);
    static int currn = 0;
    if (n > currn) {
      work.reset(new std::complex<double>[n]);
      currn = n;
    }
    return work.get();
  }

  inline float* LAP_SWork(const int n)
  {
    static tmv::auto_array<float> work(0);
    static int currn = 0;
    if (n > currn) {
      work.reset(new float[n]);
      currn = n;
    }
    return work.get();
  }

  inline std::complex<float>* LAP_CWork(const int n)
  {
    static tmv::auto_array<std::complex<float> > work(0);
    static int currn = 0;
    if (n > currn) {
      work.reset(new std::complex<float>[n]);
      currn = n;
    }
    return work.get();
  }

#endif // !LAPNOWORK

  inline void LAP_Results(const char* fn)
  {
    if (Lap_info < 0) {
      throw Error("info < 0 returned by LAPACK function ",fn);
    }
  }

  inline void LAP_Results(
#ifdef TMVDEBUG
      const int lwork_opt, const int m, const int n, const int lwork, 
#else
      const int , const int , const int, const int, 
#endif
      const char* fn)
  {
    if (Lap_info < 0) {
      throw Error("info < 0 returned by LAPACK function ",fn);
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
#endif // ALAP

#endif // TMV_BLAS_H
