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


// AlignedArray works basically like a regular new v[]
// allocation, except that if we are doing SSE or SSE2, then we enforce
// that the allocation is aligned on a 128 byte boundary.
//

#ifndef Array_H
#define Array_H

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
    // TODO: There are non-portable things like memalign, posix_memalign,
    // etc. that can be used in some cases.  I should use these
    // when available, since they are probably more efficient that 
    // the hack I do here.
#ifdef __SSE__
    template <>
    struct AlignedMemory<float>
    {
        float* mem;
        float* p;

        AlignedMemory() : mem(0), p(0) {}
        inline void allocate(const size_t n) 
        { 
            mem = new float[n+3];
            p = mem + ((0x10 - int((size_t)(mem) & 0xf))>>2);
            TMVAssert( ((size_t)(p) & 0xf) == 0 );
        }
        inline void deallocate()
        { if (mem) delete [] mem; mem=0; p=0; }
        inline void swapWith(AlignedMemory<float>& rhs)
        { 
            float* temp = p; p = rhs.p; rhs.p = temp; 
            temp = mem; mem = rhs.mem; rhs.mem = temp; 
        }
        inline float* getP() { return p; }
        inline const float* getP() const { return p; }
    };
#endif
#ifdef __SSE2__
    template <>
    struct AlignedMemory<double>
    {
        double* mem;
        double* p;

        AlignedMemory() : mem(0), p(0) {}
        inline void allocate(const size_t n) 
        { 
            mem = new double[n+1];
            p = mem + ((0x10 - int((size_t)(mem) & 0xf))>>3);
            TMVAssert( ((size_t)(p) & 0xf) == 0 );
        }
        inline void deallocate()
        { if (p) delete [] mem; mem=0; p=0; }
        inline void swapWith(AlignedMemory<double>& rhs)
        {
            double* temp = p; p = rhs.p; rhs.p = temp; 
            temp = mem; mem = rhs.mem; rhs.mem = temp; 
        }
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

} // namespace tmv

#endif
