
#ifndef TMV_ProdXM_H
#define TMV_ProdXM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXV.h"
#include "TMV_MultXM_Funcs.h"

namespace tmv {
    
    //
    // Scalar * Matrix
    //

    // These first few are for when an argument is a composite matrix
    // and needs to be calculated before running MultXM.
    template <bool add, int ix, class T, class M1, class M2>
    inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2)
    { MultXM<add>(x,m1.calc(),m2.mat()); }

    // If everything is _Calc, then there should be an overload
    // that says how to do the calculation.  This will give a 
    // compiler error on purpose.
    template <bool add, int ix, class T, class M1, class M2>
    inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2)
    { TMVStaticAssert(ix == 999); }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <bool add, class M1, class M2>
    inline void MultXM(
        const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MultXM<add>(Scaling<1,typename M2::real_type>(),m1.calc(),m2.mat()); }
    template <bool add, class T, class M1, class M2>
    inline void MultXM(
        T x, const BaseMatrix<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MultXM<add>(Scaling<0,T>(x),m1.calc(),m2.mat()); }
    template <class T, class M1>
    inline void Scale(T x, BaseMatrix_Mutable<M1>& m1)
    { Scale(Scaling<0,T>(x),m1.mat()); }



    template <int ix, class T, class M>
    class ProdXM;

    template <int ix, class T, class M>
    struct Traits<ProdXM<ix,T,M> >
    {
        typedef typename Traits2<T,typename M::value_type>::type value_type;

        enum { _colsize = M::_colsize };
        enum { _rowsize = M::_rowsize };
        enum { _nlo = M::_nlo };
        enum { _nhi = M::_nhi };
        enum { _shape = 
            (ix == 1) ?
                int(M::_shape) :
                int(ShapeTraits<M::_shape>::nonunit_shape) };
        enum { _fort = M::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M::calc_type>::_rowmajor };
        enum { _rowmajor = rm1 };

        typedef ProdXM<ix,T,M> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) | 
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>:: type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<M::_calc,const type,calc_type>::type 
            eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M>
    class ProdXM : 
        public BaseMatrix<ProdXM<ix,T,M> >
    {
    public:

        typedef ProdXM<ix,T,M> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ProdXM(const Scaling<ix,T>& _x, const BaseMatrix<M>& _m) : 
            x(_x), m(_m.mat()) {}

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const M& getM() const { return m; }

        TMV_INLINE int colsize() const { return m.colsize(); }
        TMV_INLINE int rowsize() const { return m.rowsize(); }
        TMV_INLINE int nlo() const { return m.nlo(); }
        TMV_INLINE int nhi() const { return m.nhi(); }

        value_type cref(int i, int j) const
        { return x * m.cref(i,j); }

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
            MultXM<false>(x,m.mat(),m2.mat());
        }

    private:
        const Scaling<ix,T> x;
        const M& m;
    };

#define RT typename M::real_type
#define CT typename M::complex_type
#define CCT ConjRef<CT>

    // m *= x
    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const int x)
    { Scale(RT(x),m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const RT x)
    { Scale(x,m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const CT x)
    { Scale(x,m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const CCT x)
    { Scale(CT(x),m.mat()); }

    // m /= x
    template <class M>
    inline void LDivEq(BaseMatrix_Mutable<M>& m, const int x)
    { Scale(RT(1)/RT(x),m.mat()); }

    template <class M>
    inline void LDivEq(BaseMatrix_Mutable<M>& m, const RT x)
    { Scale(RT(1)/x,m.mat()); }

    template <class M>
    inline void LDivEq(BaseMatrix_Mutable<M>& m, const CT x)
    { Scale(RT(1)/x,m.mat()); }

    template <class M>
    inline void LDivEq(BaseMatrix_Mutable<M>& m, const CCT x)
    { Scale(RT(1)/CT(x),m.mat()); }

    // -m
    template <class M>
    TMV_INLINE ProdXM<-1,RT,M> operator-(const BaseMatrix<M>& m)
    { return ProdXM<-1,RT,M>(RT(-1),m); }

    // x * m
    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator*(const int x, const BaseMatrix<M>& m)
    { return ProdXM<0,RT,M>(RT(x),m); }

    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator*(const RT x, const BaseMatrix<M>& m)
    { return ProdXM<0,RT,M>(x,m); }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const CT x, const BaseMatrix<M>& m)
    { return ProdXM<0,CT,M>(x,m); }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const CCT x, const BaseMatrix<M>& m)
    { return CT(x)*m; }

    // m * x
    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const int x)
    { return RT(x)*m; }

    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const RT x)
    { return x*m; }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CT x)
    { return x*m; }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CCT x)
    { return CT(x)*m; }

    // m / x
    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const int x)
    { return (RT(1)/RT(x))*m; }

    template <class M>
    TMV_INLINE ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const RT x)
    { return (RT(1)/x)*m; }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CT x)
    { return (RT(1)/x)*m; }

    template <class M>
    TMV_INLINE ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CCT x)
    { return (RT(1)/CT(x))*m; }

#undef RT
#undef CT
#undef CCT

    // Consolidate x*x*m type constructs:

#define RT typename ProdXM<ix,T,M>::real_type
#define CT typename ProdXM<ix,T,M>::complex_type
#define CCT ConjRef<CT>

    // -(x*m)
    template <int ix, class T, class M>
    TMV_INLINE ProdXM<-ix,T,M> operator-(const ProdXM<ix,T,M>& pxm)
    { return ProdXM<-ix,T,M>(-pxm.getX(),pxm.getM()); }

    // x*(x*m)
    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,T,M> operator*(const int x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,T,M>(RT(x)*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,T,M> operator*(const RT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,T,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const CT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const CCT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    // (x*m)*x
    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const int x)
    { return ProdXM<0,T,M>(RT(x)*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const RT x)
    { return ProdXM<0,T,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CT x)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CCT x)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    // (x*m)/x
    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const int x)
    { return ProdXM<0,RT,M>(pxm.getX()/RT(x),pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const RT x)
    { return ProdXM<0,RT,M>(pxm.getX()/x,pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CT x)
    { return ProdXM<0,CT,M>(pxm.getX()/x,pxm.getM()); }

    template <int ix, class T, class M>
    TMV_INLINE ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CCT x)
    { return ProdXM<0,CT,M>(pxm.getX()/x,pxm.getM()); }

#undef RT
#undef CT
#undef CCT

    // Have SameStorage look into a ProdXM object:
    template <int ix1, class T1, class M1, class M2>
    TMV_INLINE bool SameStorage(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
    { return SameStorage(m1.getM().mat(),m2.mat()); }
    template <class M1, int ix2, class T2, class M2>
    TMV_INLINE bool SameStorage(
        const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SameStorage(m1.mat(),m2.getM().mat()); }
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE bool SameStorage(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SameStorage(m1.getM().mat(),m2.getM().mat()); }

    template <int ix, class T, class M>
    inline std::string TMV_Text(const ProdXM<ix,T,M>& pxm)
    {
        std::ostringstream s;
        s << "ProdXM< "<<ix<<","<<TMV_Text(T(pxm.getX()));
        s << " , "<<TMV_Text(pxm.getM())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
