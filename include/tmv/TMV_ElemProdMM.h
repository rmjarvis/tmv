

#ifndef TMV_ElemProdMM_H
#define TMV_ElemProdMM_H

#include "TMV_BaseMatrix.h"

namespace tmv {

    //
    // ElemProd functions:
    // (Element-wise Matrix * Matrix)
    //

    template <int ix, class T, class M1, class M2>
    class ElemProdMM;

    // I don't currently have any non-generic versions of this.
    // It's probably rare enough that it's not really worth writing
    // a really optimized version, so this is probably fine.

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper;

    // algo 11: ColMajor:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper<11,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            //std::cout<<"ColMajor GenericElemMultMM\n";
            //std::cout<<"add = "<<add<<", x = "<<T(x)<<std::endl;
            //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            //std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<
                Traits<T1>::isreal , typename M3::real_type, 
                typename M3::complex_type>::type T1x;
            typedef typename TypeSelect<
                Traits<T2>::isreal , typename M3::real_type, 
                typename M3::complex_type>::type T2x;
            const int N = m3.rowsize();
            for(int j=0;j<N;++j) {
                const int i1 = m3.colstart(j);
                const int i2 = m3.colend(j);
                //std::cout<<"j = "<<j<<" i = "<<i1<<".."<<i2<<std::endl;
                for(int i=i1;i<i2;++i) {
                    //std::cout<<m1.cref(i,j)<<" * "<<m2.cref(i,j);
                    Maybe<add>::add(
                        m3.ref(i,j) , x * m1.cref(i,j) * m2.cref(i,j));
                    //std::cout<<" -> "<<m3.cref(i,j)<<std::endl;
                }
            }
            //std::cout<<"m3 => "<<m3<<std::endl;
        }
    };

    // algo 12: RowMajor:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper<12,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            //std::cout<<"RowMajor GenericElemMultMM\n";
            //std::cout<<"add = "<<add<<", x = "<<T(x)<<std::endl;
            //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            //std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<
                Traits<T1>::isreal , typename M3::real_type, 
                typename M3::complex_type>::type T1x;
            typedef typename TypeSelect<
                Traits<T2>::isreal , typename M3::real_type, 
                typename M3::complex_type>::type T2x;
            const int M = m3.colsize();
            for(int i=0;i<M;++i) {
                const int j1 = m3.rowstart(i);
                const int j2 = m3.rowend(i);
                //std::cout<<"i = "<<i<<" j = "<<j1<<".."<<j2<<std::endl;
                for(int j=j1;j<j2;++j) {
                    //std::cout<<m1.cref(i,j)<<" * "<<m2.cref(i,j);
                    Maybe<add>::add(
                        m3.ref(i,j) , x * m1.cref(i,j) * m2.cref(i,j));
                    //std::cout<<" -> "<<m3.cref(i,j)<<std::endl;
                }
            }
            //std::cout<<"m3 => "<<m3<<std::endl;
        }
    };

    // algo 99: Check for aliases
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper<99,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3) && !ExactSameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3) && !ExactSameStorage(m2,m3);
            if (!s1 && !s2) {
                // No aliasing
                GenericElemMultMM_Helper<-3,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            } else {
                // Use temporary
                Maybe<add>::add(m3,ElemProdMM<ix,T,M1,M2>(x,m1,m2).calc());
            }
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper<-3,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = M3::_rowmajor ? 12 : 11;
            GenericElemMultMM_Helper<algo,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct GenericElemMultMM_Helper<-1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo =
                M3::_checkalias ? 99 :
                -3;
            GenericElemMultMM_Helper<algo,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void ElemMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M2::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        TMVAssert(m2.colsize() == m3.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        GenericElemMultMM_Helper<-1,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void ElemMultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { ElemMultMM<add>(x,m1.calc(),m2.calc(),m3.mat()); }

    template <bool add, class T, class M1, class M2, class M3>
    inline void ElemMultMM(
        const T& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { ElemMultMM<add>(Scaling<0,T>(x),m1.mat(),m2.mat(),m3.mat()); }

    template <bool add, class M1, class M2, class M3>
    inline void ElemMultMM(
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        ElemMultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.mat(),m2.mat(),m3.mat()); 
    }

    template <int ix, class T, class M1, class M2>
    struct Traits<ElemProdMM<ix,T,M1,M2> >
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        typedef typename Traits2<T,T12>::type value_type;

        enum { _colsize = M1::_colsize };
        enum { _rowsize = M1::_rowsize };
        enum { _nlo = IntTraits2<M1::_nlo,M2::_nlo>::min };
        enum { _nhi = IntTraits2<M1::_nhi,M2::_nlo>::min };
        enum { shape1 = ShapeTraits2<M1::_shape,M2::_shape>::eprod };
        enum { _shape = ix==1 ? shape1 : ShapeTraits<shape1>::nonunit_shape };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { cm2 = Traits<typename M2::calc_type>::_colmajor };
        enum { _rowmajor =
            (rm1 && rm2) || (!cm1 && !cm2 && _rowsize > int(_colsize)) };

        typedef ElemProdMM<ix,T,M1,M2> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>:: type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<M1::_calc&&M2::_calc,
                const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M1, class M2>
    class ElemProdMM : 
        public BaseMatrix<ElemProdMM<ix,T,M1,M2> >
    {
    public:

        typedef ElemProdMM<ix,T,M1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ElemProdMM(
            const T _x, const BaseMatrix<M1>& _m1, const BaseMatrix<M2>& _m2) :
            x(_x), m1(_m1.mat()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
        }

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const M1& getM1() const { return m1; }
        TMV_INLINE const M2& getM2() const { return m2; }

        TMV_INLINE int colsize() const { return m1.colsize(); }
        TMV_INLINE int rowsize() const { return m1.rowsize(); }
        TMV_INLINE int nlo() const { return TMV_MIN(m1.nlo(),m2.nlo()); }
        TMV_INLINE int nhi() const { return TMV_MIN(m1.nhi(),m2.nhi()); }

        value_type cref(int i, int j) const
        { return x * (m1.cref(i,j) * m2.cref(i,j)); }

        template <class M3>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable));
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same));
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same));
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            ElemMultMM<false>(x,m1.calc(),m2.calc(),m3.mat());
        }

    private:
        const Scaling<ix,T> x;
        const M1& m1;
        const M2& m2;
    };

    // [m * m]
#define RT typename M2::real_type
    template <class M1, class M2>
    TMV_INLINE_ND ElemProdMM<1,RT,M1,M2> ElemProd(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return ElemProdMM<1,RT,M1,M2>(RT(1),m1,m2); }
#undef RT

    // [m * xm]
    template <class M1, int ix, class T, class M2>
    TMV_INLINE ElemProdMM<ix,T,M1,M2> ElemProd(
        const BaseMatrix<M1>& m1, const ProdXM<ix,T,M2>& m2)
    { return ElemProdMM<ix,T,M1,M2>(m2.getX(),m1,m2.getM()); }

    // [xm * m]
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<ix,T,M1,M2> ElemProd(
        const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
    { return ElemProdMM<ix,T,M1,M2>(m1.getX(),m1.getM(),m2); }

    // [xm * xm]
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE ElemProdMM<ix1*ix2,PT,M1,M2> ElemProd(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return ElemProdMM<ix1*ix2,PT,M1,M2>(
            m1.getX()*m2.getX(),m1.getM(),m2.getM());
    }
#undef PT

    // m += [mm]
    template <class M3, int ix, class T, class M1, class M2>
    TMV_INLINE void AddEq(
        BaseMatrix_Mutable<M3>& m3, const ElemProdMM<ix,T,M1,M2>& mm)
    {
        ElemMultMM<true>(
            mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
    }

    // m -= [mm]
    template <class M3, int ix, class T, class M1, class M2>
    TMV_INLINE void SubtractEq(
        BaseMatrix_Mutable<M3>& m3, const ElemProdMM<ix,T,M1,M2>& mm)
    {
        ElemMultMM<true>(
            -mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
    }

    // Consolidate x*[xmm] type constructs:

#define RT typename ElemProdMM<ix,T,M1,M2>::real_type
#define CT typename ElemProdMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

    // -[mm]
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<-ix,T,M1,M2> operator-(
        const ElemProdMM<ix,T,M1,M2>& mm)
    { return ElemProdMM<-ix,T,M1,M2>(-mm.getX(),mm.getM1(),mm.getM2()); }

    // x * [mm]
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator*(
        const int x, const ElemProdMM<ix,T,M1,M2>& mm)
    { return ElemProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator*(
        const RT x, const ElemProdMM<ix,T,M1,M2>& mm)
    { return ElemProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator*(
        const CT x, const ElemProdMM<ix,T,M1,M2>& mm)
    { return ElemProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator*(
        const CCT x, const ElemProdMM<ix,T,M1,M2>& mm)
    { return ElemProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    // [mm]*x
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator*(
        const ElemProdMM<ix,T,M1,M2>& mm, const int x)
    { return ElemProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator*(
        const ElemProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ElemProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator*(
        const ElemProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ElemProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator*(
        const ElemProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ElemProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    // [mm]/x
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator/(
        const ElemProdMM<ix,T,M1,M2>& mm, const int x)
    { return ElemProdMM<0,T,M1,M2>(mm.getX()/RT(x),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,T,M1,M2> operator/(
        const ElemProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ElemProdMM<0,T,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator/(
        const ElemProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ElemProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ElemProdMM<0,CT,M1,M2> operator/(
        const ElemProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ElemProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

#undef RT
#undef CT
#undef CCT

    template <int ix, class T, class M1, class M2>
    inline std::string TMV_Text(const ElemProdMM<ix,T,M1,M2>& mm)
    {
        std::ostringstream s;
        s << "ElemProdMM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(mm.getM1())<<" , "<<TMV_Text(mm.getM2())<<" >";
        return s.str();
    }


} // namespace tmv

#endif 
