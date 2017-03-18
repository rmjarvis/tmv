

#ifndef TMV_Prefetch_H
#define TMV_Prefetch_H

namespace tmv {

    // We have 4 functions that (might - depending on the compiler) 
    // prefetch data.  Two are for reading and two are for writing.
    // In each case the Multi version indicates that we are going to 
    // use the data multiple times, so try to keep it in cache.

    // The only supported prefetch so far is the GNU builtin prefetch.
    // TODO: Find out what the correlate on other compilers is.
    // I'm sure icpc must have something for this.

#ifdef __GNUC__
    // TODO: I'm not sure if the __GNUC__ guard is sufficient.  I might
    // need to check which version we are using.  I should find out
    // when gcc added this feature.
    TMV_INLINE void Prefetch_Read(const void* p) 
    { __builtin_prefetch(p,0,0); }
    TMV_INLINE void Prefetch_MultiRead(const void* p) 
    { __builtin_prefetch(p,0,3); }
    TMV_INLINE void Prefetch_Write(void* p) 
    { __builtin_prefetch(p,1,0); }
    TMV_INLINE void Prefetch_MultiWrite(void* p) 
    { __builtin_prefetch(p,1,3); }
#define TMV_PREFETCH
#else
    TMV_INLINE void Prefetch_Read(const void* ) {}
    TMV_INLINE void Prefetch_MultiRead(const void* ) {}
    TMV_INLINE void Prefetch_Write(void* ) {}
    TMV_INLINE void Prefetch_MultiWrite(void* ) {}
#endif

} // namespace tmv

#endif
