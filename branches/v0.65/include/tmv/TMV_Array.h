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

#include <complex>

namespace tmv
{
    // There doesn't seem to be any portable C++ function that guarantees
    // that the memory allocated will be aligned as necessary for 
    // SSE functions.
    // Sometimes posix_memalign or memalign does this job.
    // But it doesn't seem to be standard, since some systems don't have it.
    // Other systems have _mm_malloc (ICPC) or _aligned_malloc (WIN).
    // So probably I should have ifdefs that check to see if one of these
    // is available and use that instead.
    //
    // But for now we have this simple class that simply loads a bit more 
    // memory than necessary and then finds the starting point that is 
    // 16 byte aligned.

    // First the regular non-SSE version, where we don't need aligment.
    template <class T>
    class AlignedMemory
    {
    public:
        AlignedMemory() : p(0) 
        { 
            //std::cout<<this<<" X constructor: p = "<<p<<std::endl; 
        }
        inline void allocate(const size_t n) 
        { 
            p = new T[n]; 
            //std::cout<<this<<" X allocate: n = "<<n<<"  p = "<<p<<std::endl; 
        }
        inline void deallocate()
        { 
            if (p) delete [] p; p=0; 
            //std::cout<<this<<" X deallocate: p = "<<p<<std::endl; 
        }
        inline void swapWith(AlignedMemory<T>& rhs)
        { T* temp = p; p = rhs.p; rhs.p = temp; }
        inline T* get() { return p; }
        inline const T* get() const { return p; }
    private:
        T* p;
    };

    // Now specialize float and double
    // We do this regardless of whether __SSE__ or __SSE2__ is defined,
    // since these might be allocated in a unit that doesn't defined them
    // and then have get() called in a unit that does.  This leads to problems!
    // So we always make the float and double allocations aligned at
    // 16 byte boundaries.
    template <>
    class AlignedMemory<float>
    {
    public:
        AlignedMemory() : p(0)
        {
            //std::cout<<this<<" F constructor: p = "<<(void*)p<<std::endl; 
        }
        inline void allocate(const size_t n) 
        { 
            p = new char[(n<<2)+15];
            //std::cout<<this<<" F allocate: p = "<<(void*)p<<std::endl;
            TMVAssert((void*)(p+(n<<2)+15) >= (void*)(get()+n));
        }
        inline void deallocate()
        {
            //std::cout<<this<<" F deallocate: p = "<<(void*)p<<", p = "<<((size_t)p)<<std::endl; 
            if (p) delete [] p; p=0;
        }
        inline void swapWith(AlignedMemory<float>& rhs)
        { char* temp = p; p = rhs.p; rhs.p = temp; }
        inline float* get() 
        {
            //std::cout<<this<<" F get: p = "<<(void*)p<<std::endl;
            float* pf = reinterpret_cast<float*>(
                p + ((0x10-((size_t)(p) & 0xf)) & ~0x10));
            //std::cout<<this<<" pf = "<<(void*)pf<<std::endl;
            TMVAssert( ((size_t)(pf) & 0xf) == 0);
            return pf;
        }
        inline const float* get() const 
        {
            const float* pf = reinterpret_cast<const float*>(
                p + ((0x10-((size_t)(p) & 0xf)) & ~0x10));
            TMVAssert( ((size_t)(pf) & 0xf) == 0);
            return pf;
        }
    private:
        char* p;
    };
    template <>
    class AlignedMemory<double>
    {
    public:
        AlignedMemory() : p(0)
        {
            //std::cout<<this<<" D constructor: p = "<<(void*)p<<std::endl; 
        }
        inline void allocate(const size_t n) 
        { 
            p = new char[(n<<3)+15];
            //std::cout<<this<<" D allocate: p = "<<(void*)p<<std::endl;
            TMVAssert((void*)(p+(n<<3)+15) >= (void*)(get()+n));
        }
        inline void deallocate()
        {
            //std::cout<<this<<" D deallocate: p = "<<(void*)p<<std::endl; 
            if (p) delete [] p; p=0;
        }
        inline void swapWith(AlignedMemory<double>& rhs)
        { char* temp = p; p = rhs.p; rhs.p = temp; }
        inline double* get() 
        {
            //std::cout<<this<<" D get: p = "<<(void*)p<<std::endl;
            double* pd = reinterpret_cast<double*>(
                p + ((0x10-((size_t)(p) & 0xf)) & ~0x10));
            //std::cout<<"pd = "<<(void*)pd<<std::endl;
            TMVAssert( ((size_t)(pd) & 0xf) == 0);
            return pd;
        }
        inline const double* get() const 
        {
            const double* pd = reinterpret_cast<const double*>(
                p + ((0x10-((size_t)(p) & 0xf)) & ~0x10));
            TMVAssert( ((size_t)(pd) & 0xf) == 0);
            return pd;
        }

    private :
        char* p;

    };

    // Now the actual class that makes AlignedMemory work like a normal 
    // pointer.
    template <class T>
    class AlignedArray
    {
    public :

        inline AlignedArray() 
        {
            //std::cout<<"AA Default constructor: "<<&p<<std::endl;
        }
        inline AlignedArray(const size_t n) 
        {
            //std::cout<<"AA Sized constructor: n = "<<n<<"  "<<&p<<std::endl;
            p.allocate(n); 
        }
        inline ~AlignedArray() 
        { 
            //std::cout<<"AA Destructor: "<<&p<<std::endl;
            p.deallocate(); 
        }

        inline T& operator*() { return *get(); }
        inline T* operator->() { return get(); }
        inline operator T*() { return get(); }

        inline const T& operator*() const { return *get(); }
        inline const T* operator->() const { return get(); }
        inline operator const T*() const { return get(); }

        inline void swapWith(AlignedArray<T>& rhs) { p.swapWith(rhs.p); }
        inline void resize(const size_t n) { p.deallocate(); p.allocate(n); }

        inline T* get() { return p.get(); }
        inline const T* get() const  { return p.get(); }

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

        inline AlignedArray()
        {
            //std::cout<<"CAA Default constructor: "<<&p<<std::endl;
        }
        inline AlignedArray(const size_t n) 
        { 
            //std::cout<<"CAA Sized constructor: n = "<<n<<"  "<<&p<<std::endl;
            p.allocate(n<<1); 
        }
        inline ~AlignedArray() 
        {
            //std::cout<<"CAA Destructor: "<<&p<<std::endl;
            p.deallocate(); 
        }

        inline T& operator*() { return *get(); }
        inline T* operator->() { return get(); }
        inline operator T*() { return get(); }

        inline const T& operator*() const { return *get(); }
        inline const T* operator->() const { return get(); }
        inline operator const T*() const { return get(); }

        inline void swapWith(AlignedArray<T>& rhs) { p.swapWith(rhs.p); }
        inline void resize(const size_t n) { p.deallocate(); p.allocate(n<<1); }

        inline T* get() { return reinterpret_cast<T*>(p.get()); }
        inline const T* get() const 
        { return reinterpret_cast<const T*>(p.get()); }

    private :

        AlignedMemory<RT> p;

        inline AlignedArray& operator=(AlignedArray& p2);
        inline AlignedArray(const AlignedArray& p2);
    };



    // The rest of this file implements memory that is allocated on 
    // the stack rather than on the heap.  Here things are a bit easier,
    // since this will always align __m128 objects on 16 byte boundaries.
    // So the way to get float or double aligned is simply to union
    // the array with one of these.
    
    // Another wrinkle though is that we don't actually want to use the
    // stack for very large arrays, since it will crash.  And very small
    // arrays don't need the alignement.
    // So here is a helper class that has bigN as a template parameter
    // to decide whether to use the stack or not.
    // And smallN is for N < 4 or N < 2 where the SSE alignment isn't 
    // necessary, to make sure we don't gratuitously use extra memory when 
    // we have a lot of SmallVector<float,2>'s or something like that.
    template <class T, int N, bool bigN, bool smallN> 
    class StackArray2;

    template <class T, int N>
    class StackArray2<T,N,false,false>
    { 
    public:
        inline T* get() { return p; }
        inline const T* get() const { return p; }
    private:
        T p[N]; 
    };

#ifdef __SSE__
    template <int N>
    class StackArray2<float,N,false,false>
    {
    public:
        inline float* get() { return xp.xf; }
        inline const float* get() const { return xp.xf; }
    private:
        union { float xf[N]; __m128 xm; } xp;
    };
    template <int N>
    class StackArray2<float,N,false,true>
    {
    public:
        inline float* get() { return p; }
        inline const float* get() const { return p; }
    private:
        float p[N];
    };
#endif
#ifdef __SSE2__
    template <int N>
    class StackArray2<double,N,false,false>
    { 
    public:
        inline double* get() { return xp.xd; }
        inline const double* get() const { return xp.xd; }
    private:
        union { double xd[N]; __m128d xm; } xp;
    };
    template <int N>
    class StackArray2<double,N,false,true>
    {
    public:
        inline double* get() { return p; }
        inline const double* get() const { return p; }
    private:
        double p[N];
    };
#endif

    template <class T, int N>
    class StackArray2<T,N,true,false>
    {
    public:
        inline StackArray2() : p(N) {}
        inline T* get() { return p.get(); }
        inline const T* get() const { return p.get(); }
    private:
        AlignedArray<T> p;
    };

    // Now the real class that we use: StackArray<T,N>
    template <class T, int N>
    class StackArray
    {
    public :
        inline StackArray() {}
        inline ~StackArray() {}

        inline T& operator*() { return *get(); }
        inline T* operator->() { return get(); }
        inline operator T*() { return get(); }

        inline const T& operator*() const { return *get(); }
        inline const T* operator->() const { return get(); }
        inline operator const T*() const { return get(); }

        inline T* get() { return p.get(); }
        inline const T* get() const { return p.get(); }

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

        inline T& operator*() { return *get(); }
        inline T* operator->() { return get(); }
        inline operator T*() { return get(); }

        inline const T& operator*() const { return *get(); }
        inline const T* operator->() const { return get(); }
        inline operator const T*() const { return get(); }

        inline T* get() { return reinterpret_cast<T*>(p.get()); }
        inline const T* get() const 
        { return reinterpret_cast<const T*>(p.get()); }

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
