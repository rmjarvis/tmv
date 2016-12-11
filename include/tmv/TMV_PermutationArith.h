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


#ifndef TMV_PermutationArith_H
#define TMV_PermutationArith_H

#include "TMV_Permutation.h"
#include "TMV_MatrixArith.h"

namespace tmv {

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

    template <typename T>
    class ProdPV : public VectorComposite<T>
    {
    public:
        inline ProdPV(const Permutation& _p, const GenVector<T>& _v) :
            p(_p), v(_v)
        { TMVAssert(v.size()==p.rowsize()); }
        inline ptrdiff_t size() const { return p.colsize(); }
        inline const Permutation& getP() const { return p; }
        inline const GenVector<T>& getV() const { return v; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            p.apply(v0=v);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            p.apply(v0=v);
        }
    private:
        const Permutation& p;
        const GenVector<T>& v;
    };

    template <typename T>
    class ProdVP : public VectorComposite<T>
    {
    public:
        inline ProdVP(const GenVector<T>& _v, const Permutation& _p) :
            v(_v), p(_p)
        { TMVAssert(v.size()==p.colsize()); }
        inline ptrdiff_t size() const { return p.rowsize(); }
        inline const GenVector<T>& getV() const { return v; }
        inline const Permutation& getP() const { return p; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            p.inverse().apply(v0=v);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            p.inverse().apply(v0=v);
        }
    private:
        const GenVector<T>& v;
        const Permutation& p;
    };

    template <typename T>
    inline ProdVP<T> operator*(const GenVector<T>& v, const Permutation& p)
    { return ProdVP<T>(v,p); }

    template <typename T>
    inline ProdPV<T> operator*(const Permutation& p, const GenVector<T>& v)
    { return ProdPV<T>(p,v); }

    template <typename T>
    inline VectorView<T> operator*=(VectorView<T> v, const Permutation& p)
    {
        TMVAssert(v.size() == p.colsize());
        p.inverse().apply(v);
        return v;
    }

    template <typename T>
    class ProdPM : public MatrixComposite<T>
    {
    public:
        inline ProdPM(const Permutation& _p, const BaseMatrix<T>& _m) :
            p(_p), m(_m)
        { TMVAssert(m.colsize()==p.rowsize()); }
        inline ptrdiff_t colsize() const { return p.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline const Permutation& getP() const { return p; }
        inline const BaseMatrix<T>& getM() const { return m; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.applyOnLeft(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            p.applyOnLeft(m0=m);
        }
    private:
        const Permutation& p;
        const BaseMatrix<T>& m;
    };

    template <typename T>
    class ProdMP : public MatrixComposite<T>
    {
    public:
        inline ProdMP(const BaseMatrix<T>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert(m.rowsize()==p.colsize()); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return p.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const BaseMatrix<T>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.applyOnRight(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            p.applyOnRight(m0=m);
        }
    private:
        const BaseMatrix<T>& m;
        const Permutation& p;
    };

    template <typename T>
    inline ProdMP<T> operator*(const BaseMatrix<T>& m, const Permutation& p)
    { return ProdMP<T>(m,p); }

    template <typename T>
    inline ProdPM<T> operator*(const Permutation& p, const BaseMatrix<T>& m)
    { return ProdPM<T>(p,m); }

    template <typename T>
    inline MatrixView<T> operator*=(MatrixView<T> m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.applyOnRight(m);
        return m;
    }

    template <typename T>
    class QuotVP : public VectorComposite<T>
    {
    public:
        inline QuotVP(const GenVector<T>& _v, const Permutation& _p) :
            v(_v), p(_p)
        { TMVAssert(v.size()==p.colsize()); }
        inline ptrdiff_t size() const { return p.rowsize(); }
        inline const GenVector<T>& getV() const { return v; }
        inline const Permutation& getP() const { return p; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            p.inverse().apply(v0=v);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            p.inverse().apply(v0=v);
        }
    private:
        const GenVector<T>& v;
        const Permutation& p;
    };

    template <typename T>
    class RQuotVP : public VectorComposite<T>
    {
    public:
        inline RQuotVP(const GenVector<T>& _v, const Permutation& _p) :
            v(_v), p(_p)
        { TMVAssert(v.size()==p.colsize()); }
        inline ptrdiff_t size() const { return p.rowsize(); }
        inline const GenVector<T>& getV() const { return v; }
        inline const Permutation& getP() const { return p; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            p.apply(v0=v);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            p.apply(v0=v);
        }
    private:
        const GenVector<T>& v;
        const Permutation& p;
    };

    template <typename T>
    inline QuotVP<T> operator/(const GenVector<T>& v, const Permutation& p)
    { return QuotVP<T>(v,p); }

    template <typename T>
    inline VectorView<T> operator/=(VectorView<T> v, const Permutation& p)
    {
        TMVAssert(p.colsize() == p.rowsize());
        TMVAssert(p.rowsize() == v.size());
        p.inverse().apply(v);
        return v;
    }

    template <typename T>
    inline RQuotVP<T> operator%(const GenVector<T>& v, const Permutation& p)
    { return RQuotVP<T>(v,p); }

    template <typename T>
    inline VectorView<T> operator%=(VectorView<T> v, const Permutation& p)
    {
        TMVAssert(p.colsize() == p.rowsize());
        TMVAssert(p.rowsize() == v.size());
        p.apply(v);
        return v;
    }

    template <typename T>
    class QuotMP : public MatrixComposite<T>
    {
    public:
        inline QuotMP(const BaseMatrix<T>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert( m.colsize() == p.colsize() ); }
        inline ptrdiff_t colsize() const { return p.rowsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const BaseMatrix<T>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.inverse().applyOnLeft(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            p.inverse().applyOnLeft(m0=m);
        }
    protected:
        const BaseMatrix<T>& m;
        const Permutation& p;
    };


    template <typename T>
    class RQuotMP : public MatrixComposite<T>
    {
    public:
        inline RQuotMP(const BaseMatrix<T>& _m, const Permutation& _p) :
            m(_m), p(_p)
        { TMVAssert( m.rowsize() == p.rowsize() ); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return p.colsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const BaseMatrix<T>& getM() const { return m; }
        inline const Permutation& getP() const { return p; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            p.inverse().applyOnRight(m0=m);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            p.inverse().applyOnRight(m0=m);
        }
    protected:
        const BaseMatrix<T>& m;
        const Permutation& p;
    };

    template <typename T>
    inline QuotMP<T> operator/(const BaseMatrix<T>& m, const Permutation& p)
    { return QuotMP<T>(m,p); }

    template <typename T>
    inline MatrixView<T> operator/=(MatrixView<T> m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.inverse().applyOnLeft(m);
        return m;
    }

    template <typename T>
    inline RQuotMP<T> operator%(const BaseMatrix<T>& m, const Permutation& p)
    { return RQuotMP<T>(m,p); }

    template <typename T>
    inline MatrixView<T> operator%=(MatrixView<T> m, const Permutation& p)
    {
        TMVAssert(p.colsize()==p.rowsize());
        TMVAssert(m.rowsize()==p.rowsize());
        p.inverse().applyOnRight(m);
        return m;
    }

#undef CT
#undef CCT
#undef VCT

} // namespace mv


#endif
