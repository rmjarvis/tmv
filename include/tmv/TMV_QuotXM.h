
#ifndef TMV_QuotXM_H
#define TMV_QuotXM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXV.h"
#include "TMV_ProdXM.h"
#include "TMV_InvertM_Funcs.h"

namespace tmv {
    
    //
    // Scalar * Matrix^-1
    //

    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x,
        const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MakeInverse(x,m1.calc(),m2.mat()); }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <class M1, class M2>
    inline void MakeInverse(
        const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MakeInverse(Scaling<1,typename M2::real_type>(),m1.calc(),m2.mat()); }
    template <class T, class M1, class M2>
    inline void MakeInverse(
        T x, const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MakeInverse(Scaling<0,T>(x),m1.calc(),m2.mat()); }


    template <int ix, class T, class M>
    class QuotXM;

    template <int ix, class T, class M>
    struct Traits<QuotXM<ix,T,M> >
    {
        typedef typename Traits2<T,typename M::value_type>::type value_type;

        enum { _colsize = M::_rowsize };
        enum { _rowsize = M::_colsize };
        enum { _nlo = (
                M::_nlo == 0 ? 0 : 
                IntTraits2<IntTraits<_colsize>::Sm1,0>::max ) };
        enum { _nhi = (
                M::_nhi == 0 ? 0 : 
                IntTraits2<IntTraits<_rowsize>::Sm1,0>::max ) };
        enum { shape1 = ShapeTraits<M::_shape>::inverse_shape };
        enum { _shape = 
            (ix == 1) ?
                int(shape1) :
                int(ShapeTraits<shape1>::nonunit_shape) };
        enum { _fort = M::_fort };
        enum { _calc = false };
        enum { _rowmajor = Traits<typename M::calc_type>::_rowmajor };

        typedef QuotXM<ix,T,M> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) | 
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>:: type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M>
    class QuotXM : 
        public BaseMatrix<QuotXM<ix,T,M> >
    {
    public:

        typedef QuotXM<ix,T,M> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        QuotXM(const T _x, const BaseMatrix<M>& _m) : 
            x(_x), m(_m.mat()) {}
        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const M& getM() const { return m; }

        TMV_INLINE size_t colsize() const { return m.rowsize(); }
        TMV_INLINE size_t rowsize() const { return m.colsize(); }
        TMV_INLINE int nlo() const 
        { return m.nlo() == 0 ? 0 : TMV_MAX(colsize()-1,0); }
        TMV_INLINE int nhi() const 
        { return m.nhi() == 0 ? 0 : TMV_MAX(rowsize()-1,0); }

        template <class M2>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M2::_shape>::assignable)); 
            TMVStaticAssert((Sizes<type::_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<type::_rowsize,M2::_rowsize>::same));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVStaticAssert(type::isreal || M2::iscomplex);
            MakeInverse(x,m.mat(),m2.mat());
        }

    private:
        const Scaling<ix,T> x;
        const M& m;
    };

#define RT typename M::real_type
#define CT typename M::complex_type
#define CCT ConjRef<CT>

    // x / m
    template <class M>
    TMV_INLINE QuotXM<0,RT,M> operator/(const int x, const BaseMatrix<M>& m)
    { return QuotXM<0,RT,M>(RT(x),m); }

    template <class M>
    TMV_INLINE QuotXM<0,RT,M> operator/(const RT x, const BaseMatrix<M>& m)
    { return QuotXM<0,RT,M>(x,m); }

    template <class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const CT x, const BaseMatrix<M>& m)
    { return QuotXM<0,CT,M>(x,m); }

    template <class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const CCT x, const BaseMatrix<M>& m)
    { return CT(x)/m; }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<ix,T,M> operator/(
        const Scaling<ix,T>& x, const BaseMatrix<M>& m)
    { return QuotXM<ix,T,M>(T(x),m); }

    // x / xm
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator/(const int x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,T,M>(RT(x)/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator/(const RT x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,T,M>(x/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const CT x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,CT,M>(x/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const CCT x, const ProdXM<ix,T,M>& pxm)
    { return CT(x)/pxm; }

    template <int ix1, class T1, int ix, class T, class M>
    TMV_INLINE QuotXM<ix1*ix,T,M> operator/(
        const Scaling<ix,T>& x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<ix1*ix,T,M>(T(x)/pxm.getX(),pxm.getM()); }

    // x % m
    // In this context, there is no difference between / and %
    template <class M>
    TMV_INLINE QuotXM<0,RT,M> operator%(const int x, const BaseMatrix<M>& m)
    { return x/m; }

    template <class M>
    TMV_INLINE QuotXM<0,RT,M> operator%(const RT x, const BaseMatrix<M>& m)
    { return x/m; }

    template <class M>
    TMV_INLINE QuotXM<0,CT,M> operator%(const CT x, const BaseMatrix<M>& m)
    { return x/m; }

    template <class M>
    TMV_INLINE QuotXM<0,CT,M> operator%(const CCT x, const BaseMatrix<M>& m)
    { return x/m; }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<ix,T,M> operator%(
        const Scaling<ix,T>& x, const BaseMatrix<M>& m)
    { return x/m; }

    // x % xm
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator%(const int x, const ProdXM<ix,T,M>& pxm)
    { return x/pxm; }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator%(const RT x, const ProdXM<ix,T,M>& pxm)
    { return x/pxm; }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator%(const CT x, const ProdXM<ix,T,M>& pxm)
    { return x/pxm; }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator%(const CCT x, const ProdXM<ix,T,M>& pxm)
    { return x/pxm; }

    template <int ix1, class T1, int ix, class T, class M>
    TMV_INLINE QuotXM<ix1*ix,T,M> operator%(
        const Scaling<ix,T>& x, const ProdXM<ix,T,M>& pxm)
    { return x/pxm; }

#undef RT
#undef CT
#undef CCT

    // Consolidate x*x/m type constructs:

#define RT typename QuotXM<ix,T,M>::real_type
#define CT typename QuotXM<ix,T,M>::complex_type
#define CCT ConjRef<CT>

    // -(x*m)
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<-ix,T,M> operator-(const QuotXM<ix,T,M>& qxm)
    { return QuotXM<-ix,T,M>(-qxm.getX(),qxm.getM()); }

    // x*(x*m)
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator*(const int x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,T,M>(RT(x)*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator*(const RT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,T,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator*(const CT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator*(const CCT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    TMV_INLINE QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const Scaling<ix1,T1>& x, const QuotXM<ix,T,M>& qxm)
    {
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*qxm.getX(),qxm.getM()); 
    }

    // (x*m)*x
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator*(const QuotXM<ix,T,M>& qxm, const int x)
    { return QuotXM<0,T,M>(RT(x)*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,T,M> operator*(const QuotXM<ix,T,M>& qxm, const RT x)
    { return QuotXM<0,T,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator*(const QuotXM<ix,T,M>& qxm, const CT x)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator*(const QuotXM<ix,T,M>& qxm, const CCT x)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    TMV_INLINE QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const QuotXM<ix,T,M>& qxm, const Scaling<ix1,T1>& x)
    {
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*qxm.getX(),qxm.getM()); 
    }

    // (x*m)/x
    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,RT,M> operator/(const QuotXM<ix,T,M>& qxm, const int x)
    { return QuotXM<0,RT,M>(qxm.getX()/RT(x),qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,RT,M> operator/(const QuotXM<ix,T,M>& qxm, const RT x)
    { return QuotXM<0,RT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const QuotXM<ix,T,M>& qxm, const CT x)
    { return QuotXM<0,CT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE QuotXM<0,CT,M> operator/(const QuotXM<ix,T,M>& qxm, const CCT x)
    { return QuotXM<0,CT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    TMV_INLINE QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator/(
        const QuotXM<ix,T,M>& qxm, const Scaling<ix1,T1>& x)
    {
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            qxm.getX()/T1(x),qxm.getM()); 
    }

#undef RT
#undef CT
#undef CCT

#ifdef TMV_TEXT
    template <int ix, class T, class M>
    inline std::string TMV_Text(const QuotXM<ix,T,M>& qxm)
    {
        std::ostringstream s;
        s << "QuotXM< "<<ix<<","<<TMV_Text(T(qxm.getX()));
        s << " , "<<TMV_Text(qxm.getM())<<" >";
        return s.str();
    }
#endif

} // namespace tmv

#endif
