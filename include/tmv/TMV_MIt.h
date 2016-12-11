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


//-----------------------------------------------------------------------------
//
// This file defines the iterator classes for matrices to allow the
// m << 1,2,3
//      4,5,6
// kind of assignment to work regardless of the StorageType of the matrix.
//

#ifndef TMV_MIt_H
#define TMV_MIt_H

namespace tmv {

    template <class M>
    class RMIt
    {
    public :

        RMIt() : m(0), i(0), j(0) {}
        RMIt(M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        RMIt(const RMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        RMIt<M>& operator=(const RMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~RMIt() {}

        M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const RMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const RMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::reference operator*() const
        {
            TMVAssert(m);
            return m->ref(i,j);
        }

        inline RMIt<M>& operator++()
        {
            ++j;
            if (j == m->rowend(i)) { ++i; j = m->rowstart(i); }
            return *this;
        }
        inline RMIt<M> operator++(int)
        { RMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef typename M::value_type*         pointer;
        typedef typename M::reference           reference;

    private :

        M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    class CRMIt
    {
    public :

        CRMIt() : m(0), i(0), j(0) {}
        CRMIt(const M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        CRMIt(const CRMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        CRMIt<M>& operator=(const CRMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~CRMIt() {}

        const M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const CRMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const CRMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::value_type operator*() const
        {
            TMVAssert(m);
            return m->cref(i,j);
        }

        inline CRMIt<M>& operator++()
        {
            ++j;
            if (j == m->rowend(i)) { ++i; j = m->rowstart(i); }
            return *this;
        }
        inline CRMIt<M> operator++(int)
        { CRMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const typename M::value_type*   pointer;
        typedef const typename M::value_type&   reference;

    private :

        const M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    class CMIt
    {
    public :

        CMIt() : m(0), i(0), j(0) {}
        CMIt(M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        CMIt(const CMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        CMIt<M>& operator=(const CMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~CMIt() {}

        M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const CMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const CMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::reference operator*() const
        {
            TMVAssert(m);
            return m->ref(i,j);
        }

        inline CMIt<M>& operator++()
        {
            ++i;
            if (i == m->colend(j)) { ++j; i = m->colstart(j); }
            return *this;
        }
        inline CMIt<M> operator++(int)
        { CMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef typename M::value_type*         pointer;
        typedef typename M::reference           reference;

    private :

        M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    class CCMIt
    {
    public :

        CCMIt() : m(0), i(0), j(0) {}
        CCMIt(const M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        CCMIt(const CCMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        CCMIt<M>& operator=(const CCMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~CCMIt() {}

        const M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const CCMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const CCMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::value_type operator*() const
        {
            TMVAssert(m);
            return m->cref(i,j);
        }

        inline CCMIt<M>& operator++()
        {
            ++i;
            if (i == m->colend(j)) { ++j; i = m->colstart(j); }
            return *this;
        }
        inline CCMIt<M> operator++(int)
        { CCMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const typename M::value_type*   pointer;
        typedef const typename M::value_type&   reference;

    private :

        const M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    class DMIt
    {
    public :

        DMIt() : m(0), i(0), j(0) {}
        DMIt(M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        DMIt(const DMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        DMIt<M>& operator=(const DMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~DMIt() {}

        M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const DMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const DMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::reference operator*() const
        {
            TMVAssert(m);
            return m->ref(i,j);
        }

        inline DMIt<M>& operator++()
        {
            ++i; ++j;
            if (i == m->colsize() || j == m->rowsize()) {
                if (i > j) { i -= j+1; j = 0; }
                else { j -= i-1; i = 0; }
            }
            return *this;
        }
        inline DMIt<M> operator++(int)
        { DMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef typename M::value_type*         pointer;
        typedef typename M::reference           reference;

    private :

        M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    class CDMIt
    {
    public :

        CDMIt() : m(0), i(0), j(0) {}
        CDMIt(const M* _m, ptrdiff_t _i, ptrdiff_t _j) : m(_m), i(_i), j(_j) {}
        CDMIt(const CDMIt<M>& rhs) : m(rhs.m), i(rhs.i), j(rhs.j) {}
        CDMIt<M>& operator=(const CDMIt<M>& rhs)
        {
            TMVAssert(m == rhs.m);
            i = rhs.i;
            j = rhs.j;
            return *this;
        }
        ~CDMIt() {}

        const M* getM() const { return m; }
        ptrdiff_t getI() const { return i; }
        ptrdiff_t getJ() const { return j; }

        inline bool operator==(const CDMIt<M>& rhs) const
        { return m==rhs.m && i==rhs.i && j==rhs.j; }
        inline bool operator!=(const CDMIt<M>& rhs) const
        { return !operator==(rhs); }

        inline typename M::value_type operator*() const
        {
            TMVAssert(m);
            return m->cref(i,j);
        }

        inline CDMIt<M>& operator++()
        {
            ++i; ++j;
            if (i == m->colsize() || j == m->rowsize()) {
                if (i > j) { i -= j+1; j = 0; }
                else { j -= i-1; i = 0; }
            }
            return *this;
        }
        inline CDMIt<M> operator++(int)
        { CDMIt<M> p2 = *this; operator++(); return p2; }

        typedef std::forward_iterator_tag       iterator_category;
        typedef typename M::value_type          value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const typename M::value_type*   pointer;
        typedef const typename M::value_type&   reference;

    private :

        const M* m;
        ptrdiff_t i;
        ptrdiff_t j;

    };

    template <class M>
    inline std::string TMV_Text(const RMIt<M>& it)
    { return std::string("RMIt<") + TMV_Text(*it.getM()) + ">"; }

    template <class M>
    inline std::string TMV_Text(const CRMIt<M>& it)
    { return std::string("CRMIt<") + TMV_Text(*it.getM()) + ">"; }

    template <class M>
    inline std::string TMV_Text(const CMIt<M>& it)
    { return std::string("CMIt<") + TMV_Text(*it.getM()) + ">"; }

    template <class M>
    inline std::string TMV_Text(const CCMIt<M>& it)
    { return std::string("CCMIt<") + TMV_Text(*it.getM()) + ">"; }

    template <class M>
    inline std::string TMV_Text(const DMIt<M>& it)
    { return std::string("DMIt<") + TMV_Text(*it.getM()) + ">"; }

    template <class M>
    inline std::string TMV_Text(const CDMIt<M>& it)
    { return std::string("CDMIt<") + TMV_Text(*it.getM()) + ">"; }

} // namespace tmv

#endif
