///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// This file defines two classes that are used by TMV for allocating
// memory.
//
// First, AlignedArray works basically like a regular new v[]
// allocation, except that if we are doing SSE or SSE2, then we enforce
// that the allocation is aligned on a 128 byte boundary.
//
// Second, StackArray is used for small arrays and emulates a normal C array
// on the stack: T v[N].  However, when N is large, it uses the 
// heap instead to avoid stack overflows.

#ifndef Array_H
#define Array_H

const int TMV_MaxStack = 1024; // bytes

namespace tmv
{
    // There doesn't seem to be any portable C++ function that guarantees
    // that the memory allocated will be aligned as necessary for 
    // SSE functions.
    // Sometimes posix_memalign or memalign does this job.
    // But it doesn't seem to be standard, since some systems don't have it.
    // So we make this simple class that simply loads a bit more memory than
    // necessary and then finds the starting point that is 16 byte aligned.

    // First the regular non-SSE version, where we don't need aligment.
    template <class T>
    struct AlignedMemory
    {
        T* p;

        AlignedMemory() : p(0) {}
        inline void allocate(const size_t n) 
        { p = new T[n]; }
        inline void deallocate()
        { if (p) delete [] p; p=0; }
        inline void swapWith(AlignedMemory<T>& rhs)
        { T* temp = p; p = rhs.p; rhs.p = temp; }
        inline T* getP() { return p; }
        inline const T* getP() const { return p; }
    };

    // Now specialize float and double if SSE commands are enabled
    // Thie solution is adapted from the web page:
    // http://stackoverflow.com/questions/227897/
    //   solve-the-memory-alignment-in-c-interview-question-that-stumped-me
    // See that page for more discussion of this kind of problem.
#ifdef __SSE__
    template <>
    struct AlignedMemory<float>
    {
        float* p;

#if 0
        AlignedMemory() : mem(0), p(0) {}
        float* mem;
        inline void allocate(const size_t n) 
        { 
            mem = new float[n+3];
            //p = mem + ((4 - (size_t(mem) & 3)) & 3);
            p = mem + ((4 - (((size_t)(mem) & 0xf)>>2)) & 3);
            TMVAssert( ((size_t)(p) & 0xf) == 0);
            TMVAssert(p <= mem + 3);
            TMVAssert(p >= mem);
        }
        inline void deallocate()
        { if (mem) delete [] mem; p=0; mem=0; }
        inline void swapWith(AlignedMemory<float>& rhs)
        {
            float* temp = p; p = rhs.p; rhs.p = temp; 
            temp = mem; mem = rhs.mem; rhs.mem = temp; 
        }
#else
        AlignedMemory() : p(0) {}
        inline void allocate(const size_t n) 
        { 
            p = reinterpret_cast<float*>(new __m128[(n+3)>>2]);
            TMVAssert( ((size_t)(p) & 0xf) == 0);
        }
        inline void deallocate()
        { if (p) delete [] reinterpret_cast<__m128*>(p); p=0; }
        inline void swapWith(AlignedMemory<float>& rhs)
        { float* temp = p; p = rhs.p; rhs.p = temp; }
#endif
        inline float* getP() { return p; }
        inline const float* getP() const { return p; }
    };
#endif
#ifdef __SSE2__
    template <>
    struct AlignedMemory<double>
    {
        double* p;

#if 0
        double* mem;
        AlignedMemory() : mem(0), p(0) {}
        inline void allocate(const size_t n)
        {
            mem = new double[n+1]; 
            //p = mem + (size_t(mem) & 1);
            p = mem + (((size_t)(mem) & 0xf)>>3);
            TMVAssert( ((size_t)(p) & 0xf) == 0);
            TMVAssert(p <= mem + 1);
            TMVAssert(p >= mem);
        }
        inline void deallocate()
        { if (mem) delete [] mem; p=0; mem=0; }
        inline void swapWith(AlignedMemory<double>& rhs)
        {
            double* temp = p; p = rhs.p; rhs.p = temp; 
            temp = mem; mem = rhs.mem; rhs.mem = temp; 
        }
#else
        AlignedMemory() : p(0) {}
        inline void allocate(const size_t n) 
        { 
            p = reinterpret_cast<double*>(new __m128d[(n+1)>>1]);
            TMVAssert( ((size_t)(p) & 0xf) == 0);
        }
        inline void deallocate()
        { if (p) delete [] reinterpret_cast<__m128d*>(p); p=0; }
        inline void swapWith(AlignedMemory<double>& rhs)
        { double* temp = p; p = rhs.p; rhs.p = temp; }
#endif
        inline double* getP() { return p; }
        inline const double* getP() const { return p; }
    };
#endif

    // Now the actual class that makes AlignedMemory work like a normal 
    // pointer.
    template <class T>
    class AlignedArray
    {
    public :

        inline AlignedArray() {}
        inline AlignedArray(const size_t n) { p.allocate(n); }
        inline ~AlignedArray() { p.deallocate(); }

        inline T& operator*() { return *getP(); }
        inline T* operator->() { return getP(); }
        inline operator T*() { return getP(); }

        inline const T& operator*() const { return *getP(); }
        inline const T* operator->() const { return getP(); }
        inline operator const T*() const { return getP(); }

        inline void swapWith(AlignedArray<T>& rhs) { p.swapWith(rhs.p); }
        inline void resize(const size_t n) { p.deallocate(); p.allocate(n); }

        inline T* getP() { return p.getP(); }
        inline const T* getP() const  { return p.getP(); }

    private :

        AlignedMemory<T> p;

        inline AlignedArray& operator=(AlignedArray& p2);
        inline AlignedArray(const AlignedArray& p2);
    };

    template <class RT>
    class AlignedArray<std::complex<RT> >
    {
    public :
        typedef std::complex<RT> T;

        inline AlignedArray() {}
        inline AlignedArray(const size_t n) { p.allocate(n<<1); }
        inline ~AlignedArray() { p.deallocate(); }

        inline T& operator*() { return *getP(); }
        inline T* operator->() { return getP(); }
        inline operator T*() { return getP(); }

        inline const T& operator*() const { return *getP(); }
        inline const T* operator->() const { return getP(); }
        inline operator const T*() const { return getP(); }

        inline void swapWith(AlignedArray<T>& rhs) { p.swapWith(rhs.p); }
        inline void resize(const size_t n) { p.deallocate(); p.allocate(n<<1); }

        inline T* getP() { return reinterpret_cast<T*>(p.getP()); }
        inline const T* getP() const 
        { return reinterpret_cast<const T*>(p.getP()); }

    private :

        AlignedMemory<RT> p;

        inline AlignedArray& operator=(AlignedArray& p2);
        inline AlignedArray(const AlignedArray& p2);
    };


    // This is a helper class that has bigN as a template parameter
    // to decide whether to use the stack or not.
    // smallN is for N < 4 or N < 2 where the SSE alignment isn't necessary,
    // to make sure we don't gratuitously use extra memory when we have
    // a lot of SmallVector<float,2>'s or something like that.
    template <class T, int N, bool bigN, bool smallN> 
    class StackArray2;

    template <class T, int N>
    struct StackArray2<T,N,false,false>
    { 
        T p[N]; 
        inline T* getP() { return p; }
        inline const T* getP() const { return p; }
    };

#ifdef __SSE__
    template <int N>
    struct StackArray2<float,N,false,false>
    {
        union { float xf[N]; __m128 xm; } xp;
        inline float* getP() { return xp.xf; }
        inline const float* getP() const { return xp.xf; }
    };
    template <int N>
    struct StackArray2<float,N,false,true>
    {
        float p[N];
        inline float* getP() { return p; }
        inline const float* getP() const { return p; }
    };
#endif
#ifdef __SSE2__
    template <int N>
    struct StackArray2<double,N,false,false>
    { 
        union { double xd[N]; __m128d xm; } xp;
        inline double* getP() { return xp.xd; }
        inline const double* getP() const { return xp.xd; }
    };
    template <int N>
    struct StackArray2<double,N,false,true>
    {
        double p[N];
        inline double* getP() { return p; }
        inline const double* getP() const { return p; }
    };
#endif

    template <class T, int N>
    struct StackArray2<T,N,true,false>
    {
        AlignedArray<T> p;
        inline StackArray2() : p(N) {}
        inline T* getP() { return p.getP(); }
        inline const T* getP() const { return p.getP(); }
    };

    // Now the real class that we use: StackArray<T,N>
    template <class T, int N>
    class StackArray
    {
    public :
        inline StackArray() {}
        inline ~StackArray() {}

        inline T& operator*() { return *getP(); }
        inline T* operator->() { return getP(); }
        inline operator T*() { return getP(); }

        inline const T& operator*() const { return *getP(); }
        inline const T* operator->() const { return getP(); }
        inline operator const T*() const { return getP(); }

        inline T* getP() { return p.getP(); }
        inline const T* getP() const { return p.getP(); }

    private :
        enum { bigN = N*int(sizeof(T)) > TMV_MaxStack };
        enum { smallN = (
#ifdef __SSE__
                Traits2<T,float>::sametype ? ( N < 4 ) :
#endif
#ifdef __SSE2__
                Traits2<T,double>::sametype ? ( N < 2 ) :
#endif
                false ) };

        StackArray2<T,N,bigN,smallN> p;

        inline StackArray& operator=(StackArray& p2);
        inline StackArray(const StackArray& p2);
    };

    template <class RT, int N>
    class StackArray<std::complex<RT>,N>
    {
    public :
        typedef std::complex<RT> T;

        inline StackArray() {}
        inline ~StackArray() {}

        inline T& operator*() { return *getP(); }
        inline T* operator->() { return getP(); }
        inline operator T*() { return getP(); }

        inline const T& operator*() const { return *getP(); }
        inline const T* operator->() const { return getP(); }
        inline operator const T*() const { return getP(); }

        inline T* getP() { return reinterpret_cast<T*>(p.getP()); }
        inline const T* getP() const 
        { return reinterpret_cast<const T*>(p.getP()); }

    private :

        enum { bigN = N*int(sizeof(T)) > TMV_MaxStack };
        enum { smallN = (
#ifdef __SSE__
                Traits2<T,float>::sametype ? ( N < 2 ) :
#endif
#ifdef __SSE2__
                Traits2<T,double>::sametype ? ( N < 1 ) :
#endif
                false ) };

        StackArray2<T,(N<<1),bigN,smallN> p;

        inline StackArray& operator=(StackArray& p2);
        inline StackArray(const StackArray& p2);
    };

} // namespace tmv

#endif
