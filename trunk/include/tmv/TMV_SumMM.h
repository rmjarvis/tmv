

#ifndef TMV_SumMM_H
#define TMV_SumMM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_AddMM_Funcs.h"

namespace tmv {

    template <int algo, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct GenericAddMM_Helper;

    // algo 99: Check for aliases
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct GenericAddMM_Helper<99,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1,
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);

            if (!s1 && !s2) {
                // No aliasing 
                m3.noAlias() = ProdXM<ix1,T1,M1>(x1,m1);
                m3.noAlias() += ProdXM<ix2,T2,M2>(x2,m2);
            } else if (!s2) {
                // Alias with m1 only, do m1 first
                m3.alias() = ProdXM<ix1,T1,M1>(x1,m1);
                m3.noAlias() += ProdXM<ix2,T2,M2>(x2,m2);
            } else if (!s1) {
                // Alias with m2 only, do m2 first
                m3.alias() = ProdXM<ix2,T2,M2>(x2,m2);
                m3.noAlias() += ProdXM<ix1,T1,M1>(x1,m1);
            } else {
                // Need a temporary
                typedef typename M1::copy_type M1c;
                M1c m1c = m1;
                m3.alias() = ProdXM<ix2,T2,M2>(x2,m2);
                m3.noAlias() += ProdXM<ix1,T1,M1c>(x1,m1c);
            }
        }
    };

    // algo -3: No alias. 
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct GenericAddMM_Helper<-3,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1,
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            m3.noAlias() = ProdXM<ix1,T1,M1>(x1,m1);
            m3.noAlias() += ProdXM<ix2,T2,M2>(x2,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct GenericAddMM_Helper<-1,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1,
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = M3::_checkalias ? 99 : -3;
            GenericAddMM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
        BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        GenericAddMM_Helper<-1,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix<M2>& m2,
        BaseMatrix_Mutable<M3>& m3)
    { AddMM(x1,m1.calc(),x2,m2.calc(),m3.mat()); }

    //
    // Matrix + Matrix
    //

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    class SumMM;

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    struct Traits<SumMM<ix1,T1,M1,ix2,T2,M2> >
    {
        typedef typename ProdXM<ix1,T1,M1>::value_type mtype1;
        typedef typename ProdXM<ix2,T2,M2>::value_type mtype2;
        typedef typename Traits2<mtype1,mtype2>::type value_type;

        enum { _colsize = Sizes<M1::_colsize,M2::_colsize>::size };
        enum { _rowsize = Sizes<M1::_rowsize,M2::_rowsize>::size };
        enum { _nlo = IntTraits2<M1::_nlo,M2::_nlo>::max };
        enum { _nhi = IntTraits2<M1::_nhi,M2::_nhi>::max };
        enum { _shape = ShapeTraits2<M1::_shape,M2::_shape>::sum };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { _rowmajor = ( rm1 || (rm2 && !cm1) ) };

        typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (M1::_calc && M2::_calc),const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    class SumMM : 
        public BaseMatrix<SumMM<ix1,T1,M1,ix2,T2,M2> >
    {
    public:

        typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        SumMM(const Scaling<ix1,T1>& _x1, const BaseMatrix<M1>& _m1, 
              const Scaling<ix2,T2>& _x2, const BaseMatrix<M2>& _m2) :
            x1(_x1), m1(_m1.mat()), x2(_x2), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
        }

        TMV_INLINE const Scaling<ix1,T1>& getX1() const { return x1; }
        TMV_INLINE const M1& getM1() const { return m1; }
        TMV_INLINE const Scaling<ix2,T2>& getX2() const { return x2; }
        TMV_INLINE const M2& getM2() const { return m2; }

        TMV_INLINE int colsize() const { return m1.colsize(); }
        TMV_INLINE int rowsize() const { return m1.rowsize(); }
        TMV_INLINE int nlo() const { return TMV_MAX(m1.nlo(),m2.nlo()); }
        TMV_INLINE int nhi() const { return TMV_MAX(m1.nhi(),m2.nhi()); }

        value_type cref(int i, int j) const
        { return x1 * m1.cref(i,j) + x2 * m2.cref(i,j); }

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
            AddMM(x1,m1.mat(),x2,m2.mat(),m3.mat());
        }

    private:
        const Scaling<ix1,T1> x1;
        const M1& m1;
        const Scaling<ix2,T2> x2;
        const M2& m2;
    };

#define RT typename M2::real_type
    // m += m
    template <class M1, class M2>
    inline void AddEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2) 
    { MultXM<true>(m2.mat(),m1.mat()); }

    // m += xm
    template <class M1, int ix2, class T2, class M2>
    inline void AddEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
    { MultXM<true>(m2.getX(),m2.getM().mat(),m1.mat()); }

    // m -= m
    template <class M1, class M2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { MultXM<true>(Scaling<-1,RT>(),m2.mat(),m1.mat()); }

    // m -= xm
    template <class M1, int ix2, class T2, class M2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
    { MultXM<true>(-m2.getX(),m2.getM().mat(),m1.mat()); }

    // m + m
    template <class M1, class M2>
    TMV_INLINE SumMM<1,RT,M1,1,RT,M2> operator+(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return SumMM<1,RT,M1,1,RT,M2>(RT(1),m1,RT(1),m2); }

    // xm + m
    template <int ix1, class T1, class M1, class M2>
    TMV_INLINE SumMM<ix1,T1,M1,1,RT,M2> operator+(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
    { return SumMM<ix1,T1,M1,1,RT,M2>(T1(m1.getX()),m1.getM(),RT(1),m2); }

    // m + xm
    template <class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<1,RT,M1,ix2,T2,M2> operator+(
        const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SumMM<1,RT,M1,ix2,T2,M2>(RT(1),m1,T2(m2.getX()),m2.getM()); }

    // xm + xm
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<ix1,T1,M1,ix2,T2,M2> operator+(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return SumMM<ix1,T1,M1,ix2,T2,M2>(
            T1(m1.getX()),m1.getM(),T2(m2.getX()),m2.getM()); 
    }

    // m - m
    template <class M1, class M2>
    TMV_INLINE SumMM<1,RT,M1,-1,RT,M2> operator-(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return SumMM<1,RT,M1,-1,RT,M2>(RT(1),m1,RT(-1),m2); }

    // xm - m
    template <int ix1, class T1, class M1, class M2>
    TMV_INLINE SumMM<ix1,T1,M1,-1,RT,M2> operator-(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
    { return SumMM<ix1,T1,M1,-1,RT,M2>(T1(m1.getX()),m1.getM(),RT(-1),m2); }

    // m - xm
    template <class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<1,RT,M1,-ix2,T2,M2> operator-(
        const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SumMM<1,RT,M1,-ix2,T2,M2>(RT(1),m1,-T2(m2.getX()),m2.getM()); }

    // xm - xm
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<ix1,T1,M1,-ix2,T2,M2> operator-(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return SumMM<ix1,T1,M1,-ix2,T2,M2>(
            T1(m1.getX()),m1.getM(),-T2(m2.getX()),m2.getM()); 
    }
#undef RT


    // Consolidate x*(xm+xm) type constructs:

#define RT typename SumMM<ix1,T1,M1,ix2,T2,M2>::real_type
#define CT typename SumMM<ix1,T1,M1,ix2,T2,M2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

    // -(xm+xm)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<-ix1,T1,M1,-ix2,T2,M2> operator-(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<-ix1,T1,M1,-ix2,T2,M2>(
            -T1(smm.getX1()),smm.getM1(),-T2(smm.getX2()),smm.getM2());
    }

    // x * (xm+xm)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator*(
        const int x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            RT(x)*smm.getX1(),smm.getM1(), RT(x)*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator*(
        const RT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator*(
        const CT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator*(
        const CCT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            CT(x)*smm.getX1(),smm.getM1(), CT(x)*smm.getX2(),smm.getM2());
    }

    // (xm+xm)*x
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            RT(x)*smm.getX1(),smm.getM1(), RT(x)*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            CT(x)*smm.getX1(),smm.getM1(), CT(x)*smm.getX2(),smm.getM2());
    }

    // (xm+xm)/x
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            ZProd<false,false>::quot(smm.getX1(),RT(x)),smm.getM1(),
            ZProd<false,false>::quot(smm.getX2(),RT(x)),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,T1,M1,0,T2,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            ZProd<false,false>::quot(smm.getX1(),x),smm.getM1(),
            ZProd<false,false>::quot(smm.getX2(),x),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            ZProd<false,false>::quot(smm.getX1(),x),smm.getM1(),
            ZProd<false,false>::quot(smm.getX2(),x),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE SumMM<0,CT,M1,0,CT,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            ZProd<false,false>::quot(smm.getX1(),CT(x)),smm.getM1(),
            ZProd<false,false>::quot(smm.getX2(),CT(x)),smm.getM2());
    }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



    // TMV_Text

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    inline std::string TMV_Text(const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        std::ostringstream s;
        s << "SumMM< "<< ix1<<","<<TMV_Text(T1())
            <<" , "<<TMV_Text(smm.getM1())
            <<" , "<<ix2<<","<<TMV_Text(T2())
            <<" , "<<TMV_Text(smm.getM2())<<" >";
        return s.str();
    }



} // namespace tmv

#endif 
