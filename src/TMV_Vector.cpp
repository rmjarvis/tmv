///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMV_Blas.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include "TMV_ConvertIndex.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <sstream>

namespace tmv {

#define RT TMV_RealType(T)

    // First two things from TMV_Base.h
    bool TMV_FALSE = false;

    // And also some things from TMV_Blas.h
    void LAP_Results(int Lap_info, const char* fn)
    {
        if (Lap_info < 0) {
#ifdef NOTHROW
            std::cerr<<"info < 0 returned by LAPACK function "<<fn<<std::endl;
            exit(1);
#else
            throw Error(std::string("info < 0 returned by LAPACK function ")+fn);
#endif
        }
    }

    void LAP_Results(
        int Lap_info, int lwork_opt, int m, int n, int lwork, const char* fn)
    {
        LAP_Results(Lap_info,fn);
        if (lwork_opt > lwork) { 
            std::ostringstream s;
            s << "LAPACK function " << fn << 
                " requested more workspace than provided";
            TMV_Warning(s.str());
            s.str(std::string());
            s<<"for matrix with m,n = "<<m<<','<<n<<std::endl;
            TMV_Warning(s.str());
            s.str(std::string());
            s<<"Given: "<<lwork<<", requested "<<lwork_opt<<std::endl;
            TMV_Warning(s.str());
        }
    }

    //
    // Access
    //

    template <class T> 
    T GenVector<T>::cref(ptrdiff_t i) const
    {
        const T* vi = cptr() + i*step();
        return (ct() == Conj) ? TMV_CONJ(*vi) : *vi;
    }

    template <class T, int A>
    typename VectorView<T,A>::reference VectorView<T,A>::ref(ptrdiff_t i)
    {
        T* vi = ptr() + i*step();
        return RefHelper<T>::makeRef(vi,ct());
    }

    //
    // hasSubVector
    //
    template <class T> 
    bool GenVector<T>::hasSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
    {
        if (i1==i2) return true;  // no elements
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") cannot be 0\n";
        }
        if (i1 < 0 || i1 >= size()) {
            ok = false;
            std::cerr<<"first element ("<<i1<<") must be in 0 -- "<<
                size()-1<<std::endl;
        }
        if (i2-istep < 0 || i2-istep >= size()) {
            ok = false;
            std::cerr<<"last element ("<<i2-istep<<") must be in 0 -- "<<
                size()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"range ("<<i2-i1<<") must be multiple of istep ("<<
                istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstVectorView<T,FortranStyle>::hasSubVector(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
    {
        if (i1==i2) return true;  // no elements
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") cannot be 0\n";
        }
        if (i1 < 1 || i1 > this->size()) {
            ok = false;
            std::cerr<<"first element ("<<i1<<") must be in 1 -- "<<
                this->size()<<std::endl;
        }
        if (i2 < 1 || i2 > this->size()) {
            ok = false;
            std::cerr<<"last element ("<<i2<<") must be in 1 -- "<<
                this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"range ("<<i2-i1<<") must be multiple of istep ("<<
                istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        return ok;
    }

    //
    // norm2
    //
    template <class T> 
    RT GenVector<T>::normSq(const RT scale) const
    {
        if (size() == 0) return RT(0);

        const ptrdiff_t s = step();
        if (s == 1) {
            if (isComplex(T())) return flatten().normSq(scale);
            else {
                RT sum(0);
                const T* p = cptr();
                if (scale == RT(1)) 
                    for(ptrdiff_t i=size();i>0;--i,++p) sum += TMV_NORM(*p);
                else
                    for(ptrdiff_t i=size();i>0;--i,++p) sum += TMV_NORM(scale* (*p));
                return sum;
            }
        } else if (s < 0) {
            return reverse().normSq(scale);
        } else if (s == 0) {
            return RT(size()) * TMV_NORM(scale*(*cptr()));
        } else {
            RT sum(0);
            const T* p = cptr();
            if (scale == RT(1))
                for(ptrdiff_t i=size();i>0;--i,p+=s) sum += TMV_NORM(*p);
            else
                for(ptrdiff_t i=size();i>0;--i,p+=s) sum += TMV_NORM(scale* (*p));
            return sum;
        }
    }

#ifdef INST_INT
    static int DoNorm2(const GenVector<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenVector<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    static RT DoNorm2(const GenVector<T>& v)
    { 
        const RT eps = TMV_Epsilon<T>();

        RT vmax = v.maxAbs2Element();
        if (vmax == RT(0)) return RT(0);
        else if (TMV_Underflow(vmax * vmax)) {
            // Then we need to rescale, since underflow will cause 
            // rounding errors
            // Epsilon is a pure power of 2, so this scaling doesn't 
            // introduce any no rounding errors.
            const RT inveps = RT(1)/eps;
            RT scale = inveps;
            vmax *= scale;
            const RT eps2 = eps*eps;
            while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
            return TMV_SQRT(v.normSq(scale))/scale;
        } else if (RT(1) / (vmax) == RT(0)) {
            // Then vmax is already inf, so no hope of making more accurate.
            // (And specifically, next section would cause infinite loop.)
            return vmax;
        } else if (RT(1) / (vmax*vmax) == RT(0)) {
            // Then we have overflow, so we need to rescale:
            const RT inveps = RT(1)/eps;
            RT scale = eps;
            vmax *= scale;
            while (vmax > inveps) { scale *= eps; vmax *= eps; }
            return TMV_SQRT(v.normSq(scale))/scale;
        } 
        return TMV_SQRT(v.normSq());

    }
#ifdef BLAS
#ifndef BLASNORETURN
#ifdef INST_DOUBLE
    template <> 
    double DoNorm2(const GenVector<double>& v)
    {
        int n=v.size();
        int s=v.step();
        const double* vp = v.cptr();
        return BLASNAME(dnrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
    template <> 
    double DoNorm2(const GenVector<std::complex<double> >& v)
    { 
        int n=v.size();
        int s=v.step();
        const std::complex<double>* vp = v.cptr();
        return BLASNAME(dznrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    float DoNorm2(const GenVector<float>& v)
    {
        int n=v.size();
        int s=v.step();
        const float* vp = v.cptr();
        return BLASNAME(snrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
    template <> 
    float DoNorm2(const GenVector<std::complex<float> >& v)
    { 
        int n=v.size();
        int s=v.step();
        const std::complex<float>* vp = v.cptr();
        return BLASNAME(scnrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#endif // BLASNORETURN
#endif // BLAS
    template <class T> 
    RT GenVector<T>::norm2() const
    { 
        if (size() == 0) return RT(0);
        else if (step() < 0) return DoNorm2(reverse());
        else  return DoNorm2(*this);
    }

    //
    // sumElements
    //
    template <class T> 
    T GenVector<T>::sumElements() const
    {
        if (size() == 0) return T(0);

        const ptrdiff_t s = step();
        if (s == 1) {
            T sum(0);
            const T* p = cptr();
            for(ptrdiff_t i=size();i>0;--i,++p) sum += *p;
            return this->isconj() ? TMV_CONJ(sum) : sum;
        } else if (s < 0) {
            return reverse().sumElements();
        } else if (s == 0) {
            return RT(size()) * (*cptr());
        } else {
            T sum(0);
            const T* p = cptr();
            for(ptrdiff_t i=size();i>0;--i,p+=s) sum += *p;
            return this->isconj() ? TMV_CONJ(sum) : sum;
        }
    }

    template <class T> 
    static RT DoSumAbsElements(const GenVector<T>& v)
    {
        TMVAssert(v.step() > 0);
        const ptrdiff_t s = v.step();
        RT sum(0);
        const T* p = v.cptr();
        if (s == 1) {
            for(ptrdiff_t i=v.size();i>0;--i,++p) sum += TMV_ABS(*p);
        } else {
            for(ptrdiff_t i=v.size();i>0;--i,p+=s) sum += TMV_ABS(*p);
        }
        return sum;
    }

#ifdef BLAS
#ifndef BLASNORETURN
#ifdef INST_DOUBLE
    template <> 
    double DoSumAbsElements(const GenVector<double>& v)
    { 
        int n=v.size();
        int s=v.step();
        return BLASNAME(dasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    float DoSumAbsElements(const GenVector<float>& v)
    { 
        int n=v.size();
        int s=v.step();
        return BLASNAME(sasum) (BLASV(n),BLASP(v.cptr()),BLASV(s));
    }
#endif
#ifdef XLAP
#ifdef INST_DOUBLE
    template <> 
    double DoSumAbsElements(
        const GenVector<std::complex<double> >& v)
    { 
        int n = v.size();
        int s = v.step();
        return LAPNAME(dzsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
    }
#endif
#ifdef INST_FLOAT
    template <> 
    float DoSumAbsElements(
        const GenVector<std::complex<float> >& v)
    { 
        int n = v.size();
        int s = v.step();
        return LAPNAME(scsum1) (LAPV(n),LAPP(v.cptr()),LAPV(s)); 
    }
#endif
#endif // XLAP
#endif // BLASNORETURN
#endif // BLAS

    template <class T> 
    static T DoSumAbs2Elements(const GenVector<T>& v)
    { return DoSumAbsElements(v); }

    template <class T> 
    static T DoSumAbs2Elements(const GenVector<std::complex<T> >& v)
    {
        TMVAssert(v.step() > 0);
        if (v.step() == 1) return DoSumAbsElements(v.flatten());
        else {
            const ptrdiff_t s = v.step();
            T sum(0);
            const std::complex<T>* p = v.cptr();
            for(ptrdiff_t i=v.size();i>0;--i,p+=s) sum += TMV_ABS2(*p);
            return sum;
        }
    }

#ifdef INST_INT
    static int DoSumAbsElements(const GenVector<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenVector<T>::sumAbsElements() const
    {
        if (size() == 0) return RT(0);
        else if (step() > 0) return DoSumAbsElements(*this); 
        else if (step() < 0) return DoSumAbsElements(reverse()); 
        // If s == 0, the BLAS standard is to return 0 for the sum,
        // rather than the correct value.  Weird.
        // The non-blas version works correctly, but just do it here anyway.
        else return RT(size() * TMV_ABS(*cptr()));
    }

    template <class T> 
    RT GenVector<T>::sumAbs2Elements() const
    {
        if (size() == 0) return RT(0);
        else if (step() > 0) return DoSumAbs2Elements(*this); 
        else if (step() < 0) return DoSumAbs2Elements(reverse()); 
        // If s == 0, the BLAS standard is to return 0 for the sum,
        // rather than the correct value.  Weird.
        // The non-blas version works correctly, but just do it here anyway.
        else return RT(size() * TMV_ABS2(*cptr()));
    }

    //
    // Find Min/Max Element
    //

    template <class T> 
    static T FindMinElement(const GenVector<T>& v, ptrdiff_t& imin)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        T min = *p;
        imin = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                if (TMV_REAL(*p) < TMV_REAL(min)) {
                    min = *p;
                    imin = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                if (TMV_REAL(*p) < TMV_REAL(min)) {
                    min = *p;
                    imin = i;
                }
            }
        }
        return v.isconj() ? TMV_CONJ(min) : min;
    }
    template <class T> 
    static T FindMaxElement(const GenVector<T>& v, ptrdiff_t& imax)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        T max = *p;
        imax = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                if (TMV_REAL(*p) > TMV_REAL(max)) {
                    max = *p;
                    imax = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                if (TMV_REAL(*p) > TMV_REAL(max)) {
                    max = *p;
                    imax = i;
                }
            }
        }
        return v.isconj() ? TMV_CONJ(max) : max;
    }
    template <class T> 
    static RT FindMaxAbsElement(const GenVector<T>& v, ptrdiff_t& imax)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        RT max = TMV_ABS(*p);
        imax = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                RT absval = TMV_ABS(*p);
                if (absval > max) { 
                    max = absval; 
                    imax = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                RT absval = TMV_ABS(*p);
                if (absval > max) { 
                    max = absval; 
                    imax = i;
                }
            }
        }
        return max;
    }
    template <class T> 
    static RT FindMinAbsElement(const GenVector<T>& v, ptrdiff_t& imin)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        RT min = TMV_ABS(*p);
        imin = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                RT absval = TMV_ABS(*p);
                if (absval < min) {
                    min = absval; 
                    imin = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                RT absval = TMV_ABS(*p);
                if (absval < min) {
                    min = absval; 
                    imin = i;
                }
            }
        }
        return min;
    }

    template <class T> 
    static RT FindMaxAbs2Element(const GenVector<T>& v, ptrdiff_t& imax)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        RT max = TMV_ABS2(*p);
        imax = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                RT absval = TMV_ABS2(*p);
                if (absval > max) { 
                    max = absval; 
                    imax = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                RT absval = TMV_ABS2(*p);
                if (absval > max) { 
                    max = absval; 
                    imax = i;
                }
            }
        }
        return max;
    }
    template <class T> 
    static RT FindMinAbs2Element(const GenVector<T>& v, ptrdiff_t& imin)
    {
        TMVAssert(v.size() > 0);
        TMVAssert(v.step() > 0);

        const T* p = v.cptr();
        const ptrdiff_t s = v.step();
        RT min = TMV_ABS2(*p);
        imin = 0;
        ptrdiff_t i=1;
        if (s == 1) {
            ++p;
            for(ptrdiff_t k=v.size()-1;k>0; --k,++p,++i) {
                RT absval = TMV_ABS2(*p);
                if (absval < min) {
                    min = absval; 
                    imin = i;
                }
            }
        } else {
            p += s;
            for(ptrdiff_t k=v.size()-1;k>0; --k,p+=s,++i) {
                RT absval = TMV_ABS2(*p);
                if (absval < min) {
                    min = absval; 
                    imin = i;
                }
            }
        }
        return min;
    }
#ifdef BLAS
    // These return values seem to work, so I don't guard this segment 
    // with BLASNORETURN
#ifdef INST_DOUBLE
    static double FindMaxAbsElement(const GenVector<double>& v, ptrdiff_t& imax)
    {
        int n=v.size();
        int s=v.step();
        imax = BLASNAME(idamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imax;
#endif
        // If the input vector has nan's (especially all nan's), then
        // the output imax might not be in the allowed bounds.
        // So make sure to reset it to 0, so we don't get a segfault.
        if (imax < 0 || imax >= v.size()) imax = 0;
        TMVAssert(imax < v.size());
        return TMV_ABS(v[imax]);
    }
    static double FindMaxAbs2Element(const GenVector<double>& v, ptrdiff_t& imax)
    { return FindMaxAbsElement(v,imax); }
    static double FindMaxAbs2Element(
        const GenVector<std::complex<double> >& v, ptrdiff_t& imax)
    {
        int n=v.size();
        int s=v.step();
        imax = BLASNAME(izamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imax;
#endif
        if (imax < 0 || imax >= v.size()) imax = 0;
        TMVAssert(imax < v.size());
        return TMV_ABS2(v[imax]);
    }
#ifdef BLASIDAMIN
    static double FindMinAbsElement(const GenVector<double>& v, ptrdiff_t& imin)
    {
        int n=v.size();
        int s=v.step();
        imin = BLASNAME(idamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imin;
#endif
        if (imin < 0 || imin >= v.size()) imin = 0;
        TMVAssert(imin < v.size());
        return TMV_ABS(v[imin]);
    }
    static double FindMinAbs2Element(const GenVector<double>& v, ptrdiff_t& imin)
    { return FindMinAbsElement(v,imin); }
    static double FindMinAbs2Element(
        const GenVector<std::complex<double> >& v, ptrdiff_t& imin)
    {
        int n=v.size();
        int s=v.step();
        imin = BLASNAME(izamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imin;
#endif
        if (imin < 0 || imin >= v.size()) imin = 0;
        TMVAssert(imin < v.size());
        return TMV_ABS2(v[imin]);
    }
#endif // BLASIDAMIN
#endif
#ifdef INST_FLOAT
    static float FindMaxAbsElement(const GenVector<float>& v, ptrdiff_t& imax)
    {
        int n=v.size();
        int s=v.step();
        imax = BLASNAME(isamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imax;
#endif
        if (imax < 0 || imax >= v.size()) imax = 0;
        TMVAssert(imax < v.size());
        return TMV_ABS(v[imax]);
    }
    static float FindMaxAbs2Element(const GenVector<float>& v, ptrdiff_t& imax)
    { return FindMaxAbsElement(v,imax); }
    static float FindMaxAbs2Element(
        const GenVector<std::complex<float> >& v, ptrdiff_t& imax)
    {
        int n=v.size();
        int s=v.step();
        imax = BLASNAME(icamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imax;
#endif
        if (imax < 0 || imax >= v.size()) imax = 0;
        TMVAssert(imax < v.size());
        return TMV_ABS2(v[imax]);
    }
#ifdef BLASIDAMIN
    static float FindMinAbsElement(const GenVector<float>& v, ptrdiff_t& imin)
    {
        int n=v.size();
        int s=v.step();
        imin = BLASNAME(isamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imin;
#endif
        if (imin < 0 || imin >= v.size()) imin = 0;
        TMVAssert(imin < v.size());
        return TMV_ABS(v[imin]);
    }
    static float FindMinAbs2Element(const GenVector<float>& v, ptrdiff_t& imin)
    { return FindMinAbsElement(v,imin); }
    static float FindMinAbs2Element(
        const GenVector<std::complex<float> >& v, ptrdiff_t& imin)
    {
        int n=v.size();
        int s=v.step();
        imin = BLASNAME(icamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
        --imin;
#endif
        if (imin < 0 || imin >= v.size()) imin = 0;
        TMVAssert(imin < v.size());
        return TMV_ABS2(v[imin]);
    }
#endif // BLASAMIN
#endif // FLOAT
#endif // BLAS

#ifdef INST_INT
    static int FindMinAbsElement(
        const GenVector<std::complex<int> >& , ptrdiff_t& imax)
    { TMVAssert(TMV_FALSE); imax=0; return 0; }
    static int FindMaxAbsElement(
        const GenVector<std::complex<int> >& , ptrdiff_t& imin)
    { TMVAssert(TMV_FALSE); imin=0; return 0; }
#endif

    template <class T> 
    T GenVector<T>::doMinElement(ptrdiff_t* iminout) const
    {
        if (size() == 0) {
            if (iminout) *iminout = -1;
            return T(0);
        }
        if (step() > 0) {
            ptrdiff_t imin;
            T min = FindMinElement(*this,imin);
            TMVAssert(imin < size());
            if (iminout) *iminout = imin;
            return min;
        } else if (step() == 0) {
            if (iminout) *iminout = 0;
            return *cptr();
        } else {
            T min = reverse().doMinElement(iminout);
            if (iminout) *iminout = size()-1-(*iminout);
            return min;
        }
    }
    template <class T> 
    T GenVector<T>::doMaxElement(ptrdiff_t* imaxout) const
    {
        if (size() == 0) {
            if (imaxout) *imaxout = -1;
            return T(0);
        }
        if (step() > 0) {
            ptrdiff_t imax;
            T max = FindMaxElement(*this,imax);
            TMVAssert(imax < size());
            if (imaxout) *imaxout = imax;
            return max;
        } else if (step() == 0) {
            if (imaxout) *imaxout = 0;
            return *cptr();
        } else {
            T max = reverse().doMaxElement(imaxout);
            if (imaxout) *imaxout = size()-1-(*imaxout);
            return max;
        }
    }

    template <class T> 
    RT GenVector<T>::doMinAbsElement(ptrdiff_t* iminout) const
    {
        if (size() == 0) {
            if (iminout) *iminout = -1;
            return RT(0);
        }
        if (step() > 0) {
            ptrdiff_t imin;
            RT min = FindMinAbsElement(*this,imin);
            TMVAssert(imin < size());
            if (iminout) *iminout = imin;
            return min;
        } else if (step() == 0) {
            if (iminout) *iminout = 0;
            return RT(TMV_ABS(*cptr()));
        } else {
            RT min = reverse().doMinAbsElement(iminout);
            if (iminout) *iminout = size()-1-(*iminout);
            return min;
        }
    }
    template <class T> 
    RT GenVector<T>::doMaxAbsElement(ptrdiff_t* imaxout) const
    {
        if (size() == 0) {
            if (imaxout) *imaxout = -1;
            return RT(0);
        }
        if (step() > 0) {
            ptrdiff_t imax;
            RT max = FindMaxAbsElement(*this,imax);
            TMVAssert(imax < size());
            if (imaxout) *imaxout = imax;
            return max;
        } else if (step() == 0) {
            if (imaxout) *imaxout = 0;
            return RT(TMV_ABS(*cptr()));
        } else {
            RT max = reverse().doMaxAbsElement(imaxout);
            if (imaxout) *imaxout = size()-1-(*imaxout);
            return max;
        }
    }

    template <class T> 
    RT GenVector<T>::doMinAbs2Element(ptrdiff_t* iminout) const
    {
        if (size() == 0) {
            if (iminout) *iminout = -1;
            return RT(0);
        }
        if (step() > 0) {
            ptrdiff_t imin;
            RT min = FindMinAbs2Element(*this,imin);
            TMVAssert(imin < size());
            if (iminout) *iminout = imin;
            return min;
        } else if (step() == 0) {
            if (iminout) *iminout = 0;
            return TMV_ABS2(*cptr());
        } else {
            RT min = reverse().doMinAbs2Element(iminout);
            if (iminout) *iminout = size()-1-(*iminout);
            return min;
        }
    }
    template <class T> 
    RT GenVector<T>::doMaxAbs2Element(ptrdiff_t* imaxout) const
    {
        if (size() == 0) {
            if (imaxout) *imaxout = -1;
            return RT(0);
        }
        if (step() > 0) {
            ptrdiff_t imax;
            RT max = FindMaxAbs2Element(*this,imax);
            TMVAssert(imax < size());
            if (imaxout) *imaxout = imax;
            return max;
        } else if (step() == 0) {
            if (imaxout) *imaxout = 0;
            return TMV_ABS2(*cptr());
        } else {
            RT max = reverse().doMaxAbs2Element(imaxout);
            if (imaxout) *imaxout = size()-1-(*imaxout);
            return max;
        }
    }

    //
    // Other Modifying Functions:
    //   clip
    //   setAllTo
    //   addToAll
    //   conjugateSelf
    //   reverseSelf
    //   sort
    //

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::setZero() 
    {
        if (step() == 1) std::fill_n(ptr(),size(),T(0));
        else setAllTo(T(0));
        return *this;
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::setZero() 
    {
        std::fill_n(ptr(),size(),T(0));
        return *this;
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::clip(RT thresh) 
    {
        const ptrdiff_t s = step();
        if (s < 0) {
            reverse().clip(thresh);
        } else if (s == 0) {
            if (TMV_ABS(*ptr()) < thresh) *ptr() = T(0);
        } else {
            T* p = ptr();

            if (s == 1) {
                for(ptrdiff_t i=size();i>0;--i,++p) 
                    if (TMV_ABS(*p) < thresh) *p = T(0);
            } else {
                for(ptrdiff_t i=size();i>0;--i,p+=s) 
                    if (TMV_ABS(*p) < thresh) *p = T(0);
            }
        }
        return *this; 
    }

#ifdef INST_INT
    template <>
    VectorView<std::complex<int>,CStyle>&
    VectorView<std::complex<int>,CStyle>::clip(int ) 
    { TMVAssert(TMV_FALSE); return *this; }
#endif

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::clip(RT thresh)
    { view().clip(thresh); return *this; }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::setAllTo(const T& x) 
    {
        const ptrdiff_t s = step();
        if (s < 0) {
            reverse().setAllTo(x);
        } else if (s == 0) {
#ifdef TMVFLDEBUG
            TMVAssert(ptr() >= _first);
            TMVAssert(ptr() < _last);
#endif
            *ptr() = x;
        } else if (s == 1) {
            std::fill_n(ptr(),size(),x);
        } else {
            T* p = ptr();
            if (this->isconj()) {
                for(ptrdiff_t i=size();i>0;--i,p+=s) {
#ifdef TMVFLDEBUG
                    TMVAssert(p >= _first);
                    TMVAssert(p < _last);
#endif
                    *p = TMV_CONJ(x); 
                }
            } else {
                for(ptrdiff_t i=size();i>0;--i,p+=s) {
#ifdef TMVFLDEBUG
                    TMVAssert(p >= _first);
                    TMVAssert(p < _last);
#endif
                    *p = x; 
                }
            }
        }
        return *this; 
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::setAllTo(const T& x) 
    {
        std::fill_n(ptr(),size(),x);
        return *this;
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::addToAll(const T& x) 
    {
        const ptrdiff_t s = step();
        if (s < 0) {
            reverse().addToAll(x);
        } else if (s == 0) {
            *ptr() += x;
        } else {
            T* p = ptr();

            if (this->isconj()) {
                if (s == 1) {
                    for(ptrdiff_t i=size();i>0;--i,++p) *p += TMV_CONJ(x); 
                } else {
                    for(ptrdiff_t i=size();i>0;--i,p+=s) *p += TMV_CONJ(x); 
                }
            } else {
                if (s == 1) {
                    for(ptrdiff_t i=size();i>0;--i,++p) *p += x; 
                } else {
                    for(ptrdiff_t i=size();i>0;--i,p+=s) *p += x; 
                }
            }
        }
        return *this; 
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::addToAll(const T& x)
    {
        T* p = ptr();
        for(ptrdiff_t i=size();i>0;--i,++p) *p += x;
        return *this;
    }

    template <class T> 
    static void NonLapConjugate(VectorView<std::complex<T>,CStyle> v)
    {
        TMVAssert(v.step() > 0);

        T* p = v.realPart().ptr();
        p++;

        if (v.step() == 1) {
            for(ptrdiff_t i=v.size();i>0;--i,p+=2) *p = -(*p);
        } else if (v.step() == 0) {
            *p = -*p;
        } else {
            const ptrdiff_t s = 2*v.step();
            for(ptrdiff_t i=v.size();i>0;--i,p+=s) *p = -(*p);
        }
    }
    template <class T> 
    static inline void NonLapConjugate(VectorView<T,CStyle> ) 
    {}
#ifdef ELAP
    template <class T> 
    static inline void LapConjugate(VectorView<T,CStyle> v)
    { NonLapConjugate(v); }
#ifdef INST_DOUBLE
    template <> 
    void LapConjugate(VectorView<std::complex<double>,CStyle> v)
    { 
        int n = v.size();
        int s = v.step();
        LAPNAME(zlacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapConjugate(VectorView<std::complex<float>,CStyle> v)
    {
        int n = v.size();
        int s = v.step();
        LAPNAME(clacgv) (LAPV(n),LAPP(v.ptr()),LAPV(s)); 
    }
#endif
#endif // ELAP
    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::conjugateSelf() 
    {
        if (step() < 0) {
            reverse().conjugateSelf();
        } else {
#ifdef ELAP
            LapConjugate(*this);
#else
            NonLapConjugate(*this);
#endif
        }
        return *this; 
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::conjugateSelf()
    { 
        view().conjugateSelf(); 
        return *this; 
    }


    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::DoBasis(ptrdiff_t i) 
    { 
        TMVAssert(A==CStyle);
        TMVAssert(i < size());
        setZero();
        ref(i) = T(1);
        return *this; 
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::DoBasis(ptrdiff_t i)
    {
        if (A == CStyle) {
            TMVAssert(i<size()); 
        } else { 
            TMVAssert(i>0 && i<=size()); 
        }

        const ptrdiff_t ix = (A==CStyle ? i : i-1);
        setZero();
        ref(ix) = T(1);
        return *this;
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::DoSwap(ptrdiff_t i1, ptrdiff_t i2) 
    {
        TMVAssert(i1 < size());
        TMVAssert(i2 < size());
        TMVAssert(A==CStyle);
        if (i1 != i2) {
            const ptrdiff_t s = step();
            if (s == 1) TMV_SWAP(*(ptr()+i1),*(ptr()+i2));
            else TMV_SWAP(*(ptr()+i1*s),*(ptr()+i2*s));
        }
        return *this;
    }

    template <class T, int A>
    Vector<T,A>& Vector<T,A>::DoSwap(ptrdiff_t i1, ptrdiff_t i2)
    {
        TMVAssert(i1 < size());
        TMVAssert(i2 < size());
        if (i1 != i2)  {
            if (A==CStyle) TMV_SWAP(*(ptr()+i1),*(ptr()+i2));
            else TMV_SWAP(*(ptr()+i1-1),*(ptr()+i2-1));
        }
        return *this;
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::DoPermute(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
    { 
        TMVAssert(i2 <= size());
        TMVAssert(i1 <= i2);
        TMVAssert(A==CStyle);
        for(ptrdiff_t i=i1;i<i2;++i) DoSwap(i,p[i]);
        return *this; 
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::DoReversePermute(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
    { 
        TMVAssert(i2 <= size());
        TMVAssert(i1 <= i2);
        TMVAssert(A==CStyle);
        for(ptrdiff_t i=i2;i>i1;) { --i; DoSwap(i,p[i]); }
        return *this; 
    }

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::reverseSelf() 
    {
        const ptrdiff_t s = step();
        if (s < 0) reverse().reverseSelf();
        else if (s > 0) {
            T* p1 = ptr();
            if (s == 1) 
                for(T* p2=p1+size()-1;p2>p1;++p1,--p2) TMV_SWAP(*p1,*p2);
            else
                for(T* p2=p1+s*(size()-1);p2>p1;p1+=s,p2-=s) TMV_SWAP(*p1,*p2);
        }
        return *this;
    }

    template <class T> class VTIndex
    {

    private :

        RT itsvalue;
        ptrdiff_t itsi;

    public :

        VTIndex() : itsvalue(RT(0)), itsi(0) {}

        VTIndex(T val, ptrdiff_t i, ADType ad, CompType comp) : 
            itsvalue(0), itsi(i)
        {
            bool neg = ad==Descend;
            switch(comp) {
              case RealComp : 
                   itsvalue = RT(neg ? -TMV_REAL(val) : TMV_REAL(val));
                   break;
              case AbsComp : 
                   itsvalue = RT(neg ? -TMV_ABS(val) : TMV_ABS(val));
                   break;
              case ImagComp :
                   itsvalue = RT(neg ? -TMV_IMAG(val) : TMV_IMAG(val));
                   break;
              case ArgComp : 
                   itsvalue = RT(neg ? -TMV_ARG(val) : TMV_ARG(val));
                   break;
              default : 
                   TMVAssert2(TMV_FALSE);
            }
        }

        // Use default copy, op=, destructor

        ptrdiff_t getI() const { return itsi; }
        RT getVal() const { return itsvalue; }
        bool operator<(const VTIndex& rhs) const
        { return itsvalue < rhs.itsvalue; }
        operator ptrdiff_t() const { return itsi; }

    };

    template <class T> 
    static inline std::ostream& operator<<(std::ostream& os, VTIndex<T>& i)
    { os << i.getVal(); return os; }

    template <class T> class Compare
    {
    public:

        Compare(ADType _ad, CompType _comp) : ad(_ad), comp(_comp) {}

        bool operator()(const T& x, const T& y) const 
        {
            if (ad == Ascend) {
                switch(comp) {
                  case RealComp : 
                       return TMV_REAL(x) < TMV_REAL(y);
                  case AbsComp : 
                       return TMV_ABS(x) < TMV_ABS(y);
                  case ImagComp : 
                       return TMV_IMAG(x) < TMV_IMAG(y);
                  case ArgComp : 
                       return TMV_ARG(x) < TMV_ARG(y);
                  default : 
                       TMVAssert2(TMV_FALSE); 
                       return false;
                }
            } else {
                switch(comp) {
                  case RealComp : 
                       return TMV_REAL(x) > TMV_REAL(y);
                  case AbsComp : 
                       return TMV_ABS(x) > TMV_ABS(y);
                  case ImagComp : 
                       return TMV_IMAG(x) > TMV_IMAG(y);
                  case ArgComp : 
                       return TMV_ARG(x) > TMV_ARG(y);
                  default : 
                       TMVAssert2(TMV_FALSE); 
                       return false;
                }
            } 
        }

    private:

        const ADType ad;
        const CompType comp;
    };

    template <class T, int A>
    VectorView<T,A>& VectorView<T,A>::sort(ptrdiff_t* p, ADType ad, CompType comp) 
    {
        if (tmv::Traits<T>::iscomplex) {
            if (std::numeric_limits<T>::is_integer) {
                TMVAssert(comp != AbsComp && comp != ArgComp);
            }
        } else {
            TMVAssert(comp != ImagComp && comp != ArgComp);
        }
        if (p) {
            std::vector<VTIndex<T> > newindex(size());
            const ptrdiff_t N = size();
            for(ptrdiff_t i=0;i<N;++i) newindex[i] = VTIndex<T>(ref(i),i,ad,comp);
            std::sort(newindex.begin(),newindex.end());
            ConvertIndexToPermute(size(),newindex,p);
            permute(p);
        } else {
            // Swap ad as necessary according the the conj status of the vector:
            if (this->isconj() && (comp==ImagComp || comp==ArgComp)) {
                if (ad == Ascend) ad = Descend;
                else ad = Ascend;
            }
            const Compare<T> cc(ad,comp);
            std::sort(ptr(),ptr()+size(),cc);
        }
        return *this;
    }

    //
    // Special Constructors
    //

    template <class T, int A>
    Vector<T,A> DoBasisVector(ptrdiff_t n, ptrdiff_t i)
    {
        if (A == CStyle) { TMVAssert(i<n); }
        else { TMVAssert(i>0 && i<=n); }
        Vector<T,A> temp(n,T(0));
        temp(i) = T(1);
        return temp;
    }

    //
    // Copy Vectors
    //

    template <class T>  
    void DoCopySameType(const GenVector<T>& v1, VectorView<T> v2)
    {
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v2.size()>0);
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        TMVAssert(v2.step() != -1);
        TMVAssert(v1.step() != -1 || v2.step() == 1);
        TMVAssert(v2.step() > 0 || v1.step() == 1);
        TMVAssert(!v2.isSameAs(v1));
        const T* v1ptr = v1.cptr();
        T* v2ptr = v2.ptr();
        const ptrdiff_t step1 = v1.step();
        const ptrdiff_t step2 = v2.step();

        if (step1 == 1 && step2 == 1) {
            std::copy(v1ptr,v1ptr+v2.size(),v2ptr);
        } else {
            for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = *v1ptr;
            }
        }
    }
#ifdef BLAS
#ifdef INST_DOUBLE
    template <> 
    void DoCopySameType(
        const GenVector<double>& v1, VectorView<double> v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const double* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(dcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    template <> 
    void DoCopySameType(
        const GenVector<std::complex<double> >& v1,
        VectorView<std::complex<double> > v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<double>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(zcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void DoCopySameType(
        const GenVector<float>& v1, VectorView<float> v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const float* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(scopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    template <> 
    void DoCopySameType(
        const GenVector<std::complex<float> >& v1,
        VectorView<std::complex<float> > v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<float>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(ccopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
#endif // FLOAT
#endif // BLAS

    //
    // Swap Vectors
    //

    template <class T> 
    static void conjswap(T& x, T& y)
    {
        T temp = x;
        x = TMV_CONJ(y);
        y = TMV_CONJ(temp);
    }

    template <class T> 
    static void DoSwap(VectorView<T> v1, VectorView<T> v2)
    {
        T* v1ptr = v1.ptr();
        T* v2ptr = v2.ptr();
        const ptrdiff_t step1 = v1.step();
        const ptrdiff_t step2 = v2.step();

        if (v1.isconj())
            if (step1 == 1 && step2 == 1)
                for(ptrdiff_t i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
                    conjswap(*v1ptr,*v2ptr); 
            else
                for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
                    conjswap(*v1ptr,*v2ptr); 
        else
            if (step1 == 1 && step2 == 1)
                for(ptrdiff_t i=v2.size();i>0;--i,++v1ptr,++v2ptr) 
                    TMV_SWAP(*v1ptr,*v2ptr); 
            else
                for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) 
                    TMV_SWAP(*v1ptr,*v2ptr); 
    }
#ifdef BLAS
#ifdef INST_DOUBLE
    template <> 
    void DoSwap(VectorView<double> v1, VectorView<double> v2)
    { 
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        double* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(dswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    template <> 
    void DoSwap(
        VectorView<std::complex<double> > v1, 
        VectorView<std::complex<double> > v2)
    { 
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(zswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void DoSwap(VectorView<float> v1, VectorView<float> v2)
    { 
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        float* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(sswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    template <> 
    void DoSwap(
        VectorView<std::complex<float> > v1, 
        VectorView<std::complex<float> > v2)
    { 
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(cswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
#endif
#endif // BLAS

    template <class T> 
    void Swap(VectorView<T> v1, VectorView<T> v2)
    { 
        TMVAssert2(v1.size() == v2.size());
        TMVAssert(v1.step() != 0);
        TMVAssert(v2.step() != 0);
        TMVAssert(!v2.isSameAs(v1) || v1.isconj() == v2.isconj());

        if (v1.size() > 0 && !v2.isSameAs(v1)) {
            if (shouldReverse(v1.step(),v2.step()))
                Swap(v1.reverse(),v2.reverse());
            else if (v2.isconj()) DoSwap(v1.conjugate(),v2.conjugate());
            else DoSwap(v1,v2);
        }
    }

    //
    // op ==
    //
    template <class T1, class T2> 
    bool operator==(const GenVector<T1>& v1, const GenVector<T2>& v2)
    {
        if (v1.size() != v2.size()) return false;
        else if (v2.isSameAs(v1)) return true;
        else {
            const T1* v1ptr = v1.cptr();
            const T2* v2ptr = v2.cptr();
            const ptrdiff_t step1 = v1.step();
            const ptrdiff_t step2 = v2.step();

            if (v1.isconj() == v2.isconj()) {
                if (step1 == 1 && step2 == 1) {
                    for(ptrdiff_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
                        if ( *v1ptr != *v2ptr ) return false;
                } else {
                    for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
                        if ( *v1ptr != *v2ptr ) return false;
                }
            } else {
                if (step1 == 1 && step2 == 1) {
                    for(ptrdiff_t i=v2.size();i>0;--i,++v1ptr,++v2ptr)
                        if ( *v1ptr != TMV_CONJ(*v2ptr) ) return false;
                } else {
                    for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2)
                        if ( *v1ptr != TMV_CONJ(*v2ptr) ) return false;
                }
            }
            return true;
        }
    }

    //
    // I/O
    //

    template <class T> 
    void GenVector<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = size();
        writer.begin();
        writer.writeCode("V");
        writer.writeSize(N);
        writer.writeLParen();
        for(ptrdiff_t i=0;i<N;++i) {
            if (i > 0) writer.writeSpace();
            writer.writeValue(cref(i));
        }
        writer.writeRParen();
        writer.end();
    }

#ifndef NOTHROW
    template <class T> 
    class VectorReadError : public ReadError
    {
    public :
        Vector<T> v;
        ptrdiff_t i;
        std::string exp,got;
        ptrdiff_t s;
        bool is, iseof, isbad;

        VectorReadError(std::istream& _is) throw():
            ReadError("Vector"),
            i(0), s(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        VectorReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw():
            ReadError("Vector"),
            i(0), exp(_e), got(_g), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        VectorReadError(
            ptrdiff_t _i, const GenVector<T>& _v, std::istream& _is) throw():
            ReadError("Vector"),
            v(_v), i(_i), s(v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        VectorReadError(
            ptrdiff_t _i, const GenVector<T>& _v, std::istream& _is,
            const std::string& _e, const std::string& _g) throw():
            ReadError("Vector"),
            v(_v), i(_i), exp(_e), got(_g), s(v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        VectorReadError(
            const GenVector<T>& _v, std::istream& _is, ptrdiff_t _s) throw():
            ReadError("Vector"),
            v(_v), i(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        VectorReadError(const VectorReadError<T>& rhs) throw():
            ReadError("Vector"),
            v(rhs.v), i(rhs.i), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~VectorReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Vector\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (s != v.size()) {
                os<<"Wrong size: expected "<<v.size()<<", got "<<s<<".\n";
            }
            if (!is) {
                if (iseof) {
                    os<<"Input stream reached end-of-file prematurely.\n";
                } else if (isbad) {
                    os<<"Input stream is corrupted.\n";
                } else {
                    os<<"Input stream cannot read next character.\n";
                }
            }
            if (v.size() > 0) {
                os<<"The portion of the Vector which was successfully read is: \n";
                os<<"(";
                for(ptrdiff_t ii=0;ii<i;++ii) os<<' '<<v.cref(ii)<<' ';
                os<<")\n";
            }
        }
    };
#endif

    template <class T>
    static void FinishRead(const TMV_Reader& reader, VectorView<T> v)
    {
        const ptrdiff_t n = v.size();
        std::string exp, got;
        T temp;
        if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1); 
#else
            throw VectorReadError<T>(0,v,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<n;++i) {
            if (i>0) {
                if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1); 
#else
                    throw VectorReadError<T>(i,v,reader.getis(),exp,got);
#endif
                }
            }
            if (!reader.readValue(temp)) {
#ifdef NOTHROW
                std::cerr<<"Vector Read Error: reading value\n";
                exit(1); 
#else
                throw VectorReadError<T>(i,v,reader.getis());
#endif
            }
            v.ref(i) = temp;
        }
        if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1); 
#else
            throw VectorReadError<T>(n,v,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void Vector<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("V",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1); 
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t n=size();
        if (!reader.readSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: reading size\n";
            exit(1); 
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (n != size()) resize(n);
        VectorView<T> vv = view();
        FinishRead(reader,vv);
    }

    template <class T, int A>
    void VectorView<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("V",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1); 
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t n=size();
        if (!reader.readSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: reading size\n";
            exit(1); 
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (n != size()) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: wrong size\n";
            exit(1); 
#else
            throw VectorReadError<T>(*this,reader.getis(),n);
#endif
        }
        FinishRead(reader,*this);
    }


#ifdef BLAS
#define INST_SKIP_BLAS
#endif

#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


