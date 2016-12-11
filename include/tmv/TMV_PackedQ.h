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


#ifndef TMV_PackedQ_H
#define TMV_PackedQ_H

#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

namespace tmv {

    template <typename T>
    class PackedQ : public MatrixComposite<T>
    {
    public :

        inline PackedQ(
            const tmv::GenMatrix<T>& _Q, const tmv::GenVector<T>& _beta) :
            Q(_Q), beta(_beta)
        { TMVAssert(beta.size() == Q.rowsize()); }
        inline ptrdiff_t colsize() const { return Q.colsize(); }
        inline ptrdiff_t rowsize() const { return Q.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline const GenMatrix<T>& getQ() const { return Q; }
        inline const GenVector<T>& getBeta() const { return beta; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            doAssignToM(m0);
        }

        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            doAssignToM(m0);
        }

        // v = Qt v
        template <typename T1>
        void LDivEq(VectorView<T1> v) const;

        // v = v Qt
        template <typename T1>
        void RDivEq(VectorView<T1> v) const;

        // x = Qt v
        template <typename T1, typename T2>
        void LDiv(const GenVector<T1>& v, VectorView<T2> x) const;
        // x = v Qt
        template <typename T1, typename T2>
        void RDiv(const GenVector<T1>& v, VectorView<T2> x) const;

        // m = Qt m
        template <typename T1>
        void LDivEq(MatrixView<T1> m) const;

        // m = m Qt
        template <typename T1>
        void RDivEq(MatrixView<T1> m) const;

        // x = Qt m
        template <typename T1, typename T2>
        void LDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

        // x = m Qt
        template <typename T1, typename T2>
        void RDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

        // v = Q v
        template <typename T1>
        inline void LMultEq(VectorView<T1> v) const
        { RDivEq(v.conjugate()); }

        // v = v Q
        template <typename T1>
        inline void RMultEq(VectorView<T1> v) const
        { LDivEq(v.conjugate()); }

        // x = Q v
        template <typename T1, typename T2>
        inline void LMult(const GenVector<T1>& v, VectorView<T2> x) const
        { RDiv(v.conjugate(),x.conjugate()); }

        // x = v Q
        template <typename T1, typename T2>
        inline void RMult(const GenVector<T1>& v, VectorView<T2> x) const
        { LDiv(v.conjugate(),x.conjugate()); }

        // m = Q m
        template <typename T1>
        inline void LMultEq(MatrixView<T1> m) const
        { RDivEq(m.adjoint()); }

        // m = m Q
        template <typename T1>
        inline void RMultEq(MatrixView<T1> m) const
        { LDivEq(m.adjoint()); }

        // x = Q m
        template <typename T1, typename T2>
        inline void LMult(const GenMatrix<T1>& m, MatrixView<T2> x) const
        { RDiv(m.adjoint(),x.adjoint()); }

        // x = m Q
        template <typename T1, typename T2>
        inline void RMult(const GenMatrix<T1>& m, MatrixView<T2> x) const
        { LDiv(m.adjoint(),x.adjoint()); }

    private :

        const tmv::GenMatrix<T>& Q;
        const tmv::GenVector<T>& beta;

        template <typename T1>
        void doAssignToM(MatrixView<T1> m0) const;

    };

#define CT std::complex<T>

    template <typename T, typename Tm>
    class ProdXpQ : public MatrixComposite<T>
    {
    public:
        inline ProdXpQ(const T _x, const PackedQ<Tm>& _q) : x(_x), q(_q) {}
        inline ptrdiff_t colsize() const { return q.colsize(); }
        inline ptrdiff_t rowsize() const { return q.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline T getX() const { return x; }
        inline const PackedQ<Tm>& getQ() const { return q; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            MultXM(x,m0=q);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            MultXM(x,m0=q);
        }
    private:
        const T x;
        const PackedQ<Tm>& q;
    };

#define GENMATRIX PackedQ
#define PRODXM ProdXpQ
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

    template <typename T, typename T1, typename T2>
    class ProdpQV : public VectorComposite<T>
    {
    public:
        inline ProdpQV(
            const T _x, const PackedQ<T1>& _q, const GenVector<T2>& _v) :
            x(_x), q(_q), v(_v)
        { TMVAssert(v.size()==q.rowsize()); }
        inline ptrdiff_t size() const { return q.colsize(); }
        inline T getX() const { return x; }
        inline const PackedQ<T1>& getQ() const { return q; }
        inline const GenVector<T2>& getV() const { return v; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            q.LMult(v,v0);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            q.LMult(v,v0);
        }
    private:
        const T x;
        const PackedQ<T1>& q;
        const GenVector<T2>& v;
    };

    template <typename T, typename T1, typename T2>
    class ProdVpQ : public VectorComposite<T>
    {
    public:
        inline ProdVpQ(
            const T _x, const GenVector<T1>& _v, const PackedQ<T2>& _q) :
            x(_x), v(_v), q(_q)
        { TMVAssert(v.size()==q.colsize()); }
        inline ptrdiff_t size() const { return q.rowsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            q.RMult(v,v0);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            q.RMult(v,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const PackedQ<T2>& q;
    };

    template <typename T>
    inline VectorView<T> operator*=(
        VectorView<T> v, const PackedQ<T>& m)
    {
        TMVAssert(v.size() == m.colsize());
        MultMV<false>(T(1),m.transpose(),v,v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const PackedQ<T>& m)
    {
        TMVAssert(v.size() == m.colsize());
        MultMV<false>(T(1),m.transpose(),v,v);
        return v;
    }

    template <typename T, typename T2>
    inline VectorView<T> operator*=(
        VectorView<T> v, const ProdXpQ<T,T2>& pxq)
    {
        TMVAssert(pxq.colsize()==pxq.rowsize());
        TMVAssert(v.size()==pxq.rowsize());
        pxq.getQ().RMult(v,v);
        MultXV(pxq.getX(),v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator*=(
        VectorView<CT> v, const ProdXpQ<T,T>& pxq)
    {
        TMVAssert(pxq.colsize()==pxq.rowsize());
        TMVAssert(v.size()==pxq.rowsize());
        pxq.getQ().RMult(v,v);
        MultXV(pxq.getX(),v);
        return v;
    }

#define GENMATRIX1 PackedQ
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXpQ
#define PRODXM2 ProdXV
#define PRODMM ProdpQV
#define GETM1 .getQ()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getQ()
#define GETM2 .getV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM
#define GENMATRIX1 GenVector
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXV
#define PRODXM2 ProdXpQ
#define PRODMM ProdVpQ
#define GETM1 .getV()
#define GETM2 .getQ()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getV()
#define GETM2 .getQ()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM


    template <typename T, typename T1, typename T2>
    class ProdpQM : public MatrixComposite<T>
    {
    public:
        inline ProdpQM(
            const T _x, const PackedQ<T1>& _q, const GenMatrix<T2>& _m) :
            x(_x), q(_q), m(_m)
        { TMVAssert(m.colsize()==q.rowsize()); }
        inline ptrdiff_t colsize() const { return q.colsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline T getX() const { return x; }
        inline const PackedQ<T1>& getQ() const { return q; }
        inline const GenMatrix<T2>& getM() const { return m; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            q.LMult(m,m0);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            q.LMult(m,m0);
        }
    private:
        const T x;
        const PackedQ<T1>& q;
        const GenMatrix<T2>& m;
    };

    template <typename T, typename T1, typename T2>
    class ProdMpQ : public MatrixComposite<T>
    {
    public:
        inline ProdMpQ(const T _x, const GenMatrix<T1>& _m,
                       const PackedQ<T2>& _q) : x(_x), m(_m), q(_q)
        { TMVAssert(m.rowsize()==q.colsize()); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return q.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM() const { return m; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            q.RMult(m,m0);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize());
            TMVAssert(m0.rowsize() == rowsize());
            q.RMult(m,m0);
        }
    private:
        const T x;
        const GenMatrix<T1>& m;
        const PackedQ<T2>& q;
    };

    template <typename T>
    inline MatrixView<T> operator*=(
        MatrixView<T> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize()==q.rowsize());
        TMVAssert(m.rowsize()==q.rowsize());
        q.RMultEq(m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize()==q.rowsize());
        TMVAssert(m.rowsize()==q.rowsize());
        q.RMultEq(m);
        return m;
    }

    template <typename T, typename T2>
    inline MatrixView<T> operator*=(
        MatrixView<T> m, const ProdXpQ<T,T2>& pxq)
    {
        TMVAssert(pxq.colsize()==pxq.rowsize());
        TMVAssert(m.rowsize()==pxq.rowsize());
        pxq.getQ().RMult(m,m);
        MultXM(pxq.getX(),m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator*=(
        MatrixView<CT> m, const ProdXpQ<T,T>& pxq)
    {
        TMVAssert(pxq.colsize()==pxq.rowsize());
        TMVAssert(m.rowsize()==pxq.rowsize());
        pxq.getQ().RMult(m,m);
        MultXM(pxq.getX(),m);
        return m;
    }

#define GENMATRIX1 PackedQ
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXpQ
#define PRODXM2 ProdXM
#define PRODMM ProdpQM
#define GETM1 .getQ()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getQ()
#define GETM2 .getM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXM
#define PRODXM2 ProdXpQ
#define PRODMM ProdMpQ
#define GETM1 .getM()
#define GETM2 .getQ()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .getM()
#define GETM2 .getQ()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM

    template <typename T, typename T1, typename T2>
    class QuotVpQ : public VectorComposite<T>
    {
    public:
        inline QuotVpQ(
            const T _x, const GenVector<T1>& _v, const PackedQ<T2>& _q) :
            x(_x), v(_v), q(_q)
        { TMVAssert(v.size()==q.colsize()); }
        inline ptrdiff_t size() const { return q.rowsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            q.LDiv(v,v0);
            MultXV(x,v0);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            q.LDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const PackedQ<T2>& q;
    };

    template <typename T, typename T1, typename T2>
    class RQuotVpQ : public VectorComposite<T>
    {
    public:
        inline RQuotVpQ(
            const T _x, const GenVector<T1>& _v, const PackedQ<T2>& _q) :
            x(_x), v(_v), q(_q)
        { TMVAssert(v.size()==q.colsize()); }
        inline ptrdiff_t size() const { return q.rowsize(); }
        inline T getX() const { return x; }
        inline const GenVector<T1>& getV() const { return v; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToV(VectorView<TMV_RealType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            TMVAssert(isReal(T()));
            q.RDiv(v,v0);
            MultXV(x,v0);
        }
        inline void assignToV(VectorView<TMV_ComplexType(T)> v0) const
        {
            TMVAssert(v0.size() == size());
            q.RDiv(v,v0);
            MultXV(x,v0);
        }
    private:
        const T x;
        const GenVector<T1>& v;
        const PackedQ<T2>& q;
    };

    template <typename T>
    inline VectorView<T> operator/=(
        VectorView<T> v, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(q.rowsize() == v.size());
        q.LDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator/=(
        VectorView<CT> v, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(q.rowsize() == v.size());
        q.LDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<T> operator%=(
        VectorView<T> v, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(q.rowsize() == v.size());
        q.RDivEq(v);
        return v;
    }

    template <typename T>
    inline VectorView<CT> operator%=(
        VectorView<CT> v, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(q.rowsize() == v.size());
        q.RDivEq(v);
        return v;
    }

#define GENMATRIX1 GenVector
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXV
#define PRODXM2 ProdXpQ
#define QUOTMM QuotVpQ
#define RQUOTMM RQuotVpQ
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTMM
#undef RQUOTMM

    template <typename T, typename T1, typename T2>
    class QuotMpQ : public MatrixComposite<T>
    {
    public:
        inline QuotMpQ(
            const T _x, const GenMatrix<T1>& _m, const PackedQ<T2>& _q) :
            x(_x), m(_m), q(_q)
        { TMVAssert( m.colsize() == q.colsize() ); }
        inline ptrdiff_t colsize() const { return q.rowsize(); }
        inline ptrdiff_t rowsize() const { return m.rowsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM() const { return m; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            q.LDiv(m,m0);
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            q.LDiv(m,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenMatrix<T1>& m;
        const PackedQ<T2>& q;
    };


    template <typename T, typename T1, typename T2>
    class RQuotMpQ : public MatrixComposite<T>
    {
    public:
        inline RQuotMpQ(
            const T _x, const GenMatrix<T1>& _m, const PackedQ<T2>& _q) :
            x(_x), m(_m), q(_q)
        { TMVAssert( m.rowsize() == q.rowsize() ); }
        inline ptrdiff_t colsize() const { return m.colsize(); }
        inline ptrdiff_t rowsize() const { return q.colsize(); }
        inline StorageType stor() const { return ColMajor; }
        inline T getX() const { return x; }
        inline const GenMatrix<T1>& getM() const { return m; }
        inline const PackedQ<T2>& getQ() const { return q; }
        inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            q.RDiv(m,m0);
            MultXM(x,m0);
        }
        inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
        {
            TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
            q.RDiv(m,m0);
            MultXM(x,m0);
        }
    protected:
        const T x;
        const GenMatrix<T1>& m;
        const PackedQ<T2>& q;
    };

    template <typename T>
    inline MatrixView<T> operator/=(
        MatrixView<T> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(m.colsize() == q.rowsize());
        q.LDivEq(m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator/=(
        MatrixView<CT> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(m.colsize() == q.rowsize());
        q.LDivEq(m);
        return m;
    }

    template <typename T>
    inline MatrixView<T> operator%=(
        MatrixView<T> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(m.rowsize() == q.rowsize());
        q.RDivEq(m);
        return m;
    }

    template <typename T>
    inline MatrixView<CT> operator%=(
        MatrixView<CT> m, const PackedQ<T>& q)
    {
        TMVAssert(q.colsize() == q.rowsize());
        TMVAssert(m.rowsize() == q.rowsize());
        q.RDivEq(m);
        return m;
    }


#define GENMATRIX1 GenMatrix
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXM
#define PRODXM2 ProdXpQ
#define QUOTMM QuotMpQ
#define RQUOTMM RQuotMpQ
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#undef CT

} // namespace mv


#endif
