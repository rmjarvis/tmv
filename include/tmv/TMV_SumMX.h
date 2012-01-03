

#ifndef TMV_SumMX_H
#define TMV_SumMX_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_MultXM_Funcs.h"

namespace tmv {

    //
    // Matrix + Scalar
    //

    template <class T, class M>
    inline void AddMX(const T& x, BaseMatrix_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        m.mat().diag().addToAll(x); 
    }

    template <int ix1, class T1, class M1, class T2, class M3>
    inline void AddMX(
        const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
        const T2& x2, BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        MultXM<false>(x1,m1.mat(),m3.mat());
        AddMX(x2,m3.mat());
    }


    template <int ix1, class T1, class M1, class T2>
    class SumMX;

    template <int ix1, class T1, class M1, class T2>
    struct Traits<SumMX<ix1,T1,M1,T2> >
    {
        typedef typename ProdXM<ix1,T1,M1>::value_type mtype1;
        typedef typename Traits2<mtype1,T2>::type value_type;

        enum { _colsize = M1::_colsize };
        enum { _rowsize = M1::_rowsize };
        enum { _nlo = M1::_nlo };
        enum { _nhi = M1::_nhi };
        enum { _shape = ShapeTraits<M1::_shape>::nonunit_shape };
        enum { _fort = M1::_fort };
        enum { _calc = false };
        enum { _rowmajor = Traits<typename M1::calc_type>::_rowmajor };
        enum { _colmajor = !_rowmajor };

        typedef SumMX<ix1,T1,M1,T2> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            M1::_calc,const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix1, class T1, class M1, class T2>
    class SumMX : 
        public BaseMatrix<SumMX<ix1,T1,M1,T2> >
    {
    public:

        typedef SumMX<ix1,T1,M1,T2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        SumMX(const Scaling<ix1,T1>& _x1, const BaseMatrix<M1>& _m1,
              const T2& _x2) :
            x1(_x1), m1(_m1.mat()), x2(_x2)
        {
            TMVStaticAssert((Sizes<M1::_colsize,M1::_rowsize>::same)); 
            TMVAssert(m1.colsize() == m1.rowsize());
        }

        TMV_INLINE const Scaling<ix1,T1>& getX1() const { return x1; }
        TMV_INLINE const M1& getM1() const { return m1; }
        TMV_INLINE const T2& getX2() const { return x2; }

        TMV_INLINE int colsize() const { return m1.colsize(); }
        TMV_INLINE int rowsize() const { return m1.rowsize(); }
        TMV_INLINE int nlo() const { return m1.nlo(); }
        TMV_INLINE int nhi() const { return m1.nhi(); }
        TMV_INLINE int ls() const { return m1.ls(); }

        value_type cref(int i, int j) const
        { return (i == j) ? (x1 * m1.cref(i,j) + x2) : (x1 * m1.cref(i,j)); }

        template <class M3>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            AddMX(x1,m1.mat(),x2,m3.mat());
        }

    private:
        const Scaling<ix1,T1> x1;
        const M1& m1;
        const T2 x2;
    };

#define RT typename M1::real_type
#define CT typename M1::complex_type
#define CCT ConjRef<CT>
    // m += x
    template <class M1>
    TMV_INLINE void AddEq(BaseMatrix_Mutable<M1>& m1, const int x2)
    { AddMX(RT(x2),m1.mat()); }

    template <class M1>
    TMV_INLINE void AddEq(BaseMatrix_Mutable<M1>& m1, const RT x2)
    { AddMX(x2,m1.mat()); }

    template <class M1>
    TMV_INLINE void AddEq(BaseMatrix_Mutable<M1>& m1, const CT x2)
    { AddMX(x2,m1.mat()); }

    template <class M1>
    TMV_INLINE void AddEq(BaseMatrix_Mutable<M1>& m1, const CCT x2)
    { AddMX(CT(x2),m1.mat()); }

    // m -= x
    template <class M1>
    TMV_INLINE void SubtractEq(BaseMatrix_Mutable<M1>& m1, const int x2)
    { AddMX(-RT(x2),m1.mat()); }

    template <class M1>
    TMV_INLINE void SubtractEq(BaseMatrix_Mutable<M1>& m1, const RT x2)
    { AddMX(-x2,m1.mat()); }

    template <class M1>
    TMV_INLINE void SubtractEq(BaseMatrix_Mutable<M1>& m1, const CT x2)
    { AddMX(-x2,m1.mat()); }

    template <class M1>
    TMV_INLINE void SubtractEq(BaseMatrix_Mutable<M1>& m1, const CCT x2)
    { AddMX(-CT(x2),m1.mat()); }

    // m + x
    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator+(
        const BaseMatrix<M1>& m1, const int x2)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,RT(x2)); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator+(
        const BaseMatrix<M1>& m1, const RT x2)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator+(
        const BaseMatrix<M1>& m1, const CT x2)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator+(
        const BaseMatrix<M1>& m1, const CCT x2)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,CT(x2)); }

    // x + m
    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator+(
        const int x2, const BaseMatrix<M1>& m1)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,RT(x2)); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator+(
        const RT x2, const BaseMatrix<M1>& m1)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator+(
        const CT x2, const BaseMatrix<M1>& m1)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator+(
        const CCT x2, const BaseMatrix<M1>& m1)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,CT(x2)); }

    // m - x
    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator-(
        const BaseMatrix<M1>& m1, const int x2)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,-RT(x2)); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,RT> operator-(
        const BaseMatrix<M1>& m1, const RT x2)
    { return SumMX<1,RT,M1,RT>(RT(1),m1,-x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator-(
        const BaseMatrix<M1>& m1, const CT x2)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,-x2); }

    template <class M1>
    TMV_INLINE SumMX<1,RT,M1,CT> operator-(
        const BaseMatrix<M1>& m1, const CCT x2)
    { return SumMX<1,RT,M1,CT>(RT(1),m1,-CT(x2)); }

    // x - m
    template <class M1>
    TMV_INLINE SumMX<-1,RT,M1,RT> operator-(
        const int x2, const BaseMatrix<M1>& m1)
    { return SumMX<-1,RT,M1,RT>(RT(-1),m1,RT(x2)); }

    template <class M1>
    TMV_INLINE SumMX<-1,RT,M1,RT> operator-(
        const RT x2, const BaseMatrix<M1>& m1)
    { return SumMX<-1,RT,M1,RT>(RT(-1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<-1,RT,M1,CT> operator-(
        const CT x2, const BaseMatrix<M1>& m1)
    { return SumMX<-1,RT,M1,CT>(RT(-1),m1,x2); }

    template <class M1>
    TMV_INLINE SumMX<-1,RT,M1,CT> operator-(
        const CCT x2, const BaseMatrix<M1>& m1)
    { return SumMX<-1,RT,M1,CT>(RT(-1),m1,CT(x2)); }

    // xm + x
    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const int x2)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),RT(x2)); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const RT x2)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const CT x2)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const CCT x2)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),CT(x2)); }

    // x + xm 
    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator+(
        const int x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),int(x2)); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator+(
        const RT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator+(
        const CT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator+(
        const CCT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),CT(x2)); }

    // xm - x
    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const int x2)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),-RT(x2)); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,RT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const RT x2)
    { return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),-x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const CT x2)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),-x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<ix1,T1,M1,CT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const CCT x2)
    { return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),-CT(x2)); }

    // x - xm 
    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<-ix1,T1,M1,RT> operator-(
        const int x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<-ix1,T1,M1,RT>(-T1(m1.getX()),m1.getM(),RT(x2)); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<-ix1,T1,M1,RT> operator-(
        const RT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<-ix1,T1,M1,RT>(-T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<-ix1,T1,M1,CT> operator-(
        const CT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<-ix1,T1,M1,CT>(-T1(m1.getX()),m1.getM(),x2); }

    template <int ix1, class T1, class M1>
    TMV_INLINE SumMX<-ix1,T1,M1,CT> operator-(
        const CCT x2, const ProdXM<ix1,T1,M1>& m1)
    { return SumMX<-ix1,T1,M1,CT>(-T1(m1.getX()),m1.getM(),x2); }

#undef RT
#undef CT
#undef CCT


    // Consolidate x*(xm+x) type constructs:

#define RT typename SumMX<ix1,T1,M1,T2>::real_type
#define CT typename SumMX<ix1,T1,M1,T2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

    // -(xm+x)
    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<-ix1,T1,M1,T2> operator-(
        const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<-ix1,T1,M1,T2>(-T1(smx.getX1()),smx.getM1(),-smx.getX2()); }

    // x * (xm+x)
    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator*(
        const int x, const SumMX<ix1,T1,M1,T2>& smx)
    {
        return SumMX<0,T1,M1,T2>(
            RT(x)*smx.getX1(),smx.getM1(),RT(x)*smx.getX2());
    }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator*(
        const RT x, const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<0,T1,M1,T2>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator*(
        const CT x, const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<0,CT,M1,CT>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator*(
        const CCT x, const SumMX<ix1,T1,M1,T2>& smx)
    {
        return SumMX<0,CT,M1,CT>(
            CT(x)*smx.getX1(),smx.getM1(),CT(x)*smx.getX2());
    }

    // (xm+x)*x
    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const int x)
    {
        return SumMX<0,T1,M1,T2>(
            RT(x)*smx.getX1(),smx.getM1(),RT(x)*smx.getX2());
    }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const RT x)
    { return SumMX<0,T1,M1,T2>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const CT x)
    { return SumMX<0,CT,M1,CT>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const CCT x)
    {
        return SumMX<0,CT,M1,CT>(
            CT(x)*smx.getX1(),smx.getM1(),CT(x)*smx.getX2());
    }

    // (xm+x)/x
    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const int x)
    {
        return SumMX<0,T1,M1,T2>(
            smx.getX1()/RT(x),smx.getM1(),smx.getX2()/RT(x));
    }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,T1,M1,T2> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const RT x)
    {
        return SumMX<0,T1,M1,T2>(smx.getX1()/x,smx.getM1(),smx.getX2()/x);
    }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const CT x)
    { return SumMX<0,CT,M1,CT>(smx.getX1()/x,smx.getM1(), smx.getX2()/x); }

    template <int ix1, class T1, class M1, class T2>
    TMV_INLINE SumMX<0,CT,M1,CT> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const CCT x)
    {
        return SumMX<0,CT,M1,CT>(
            smx.getX1()/CT(x),smx.getM1(),smx.getX2()/CT(x));
    }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



    // TMV_Text

#ifdef TMV_TEXT
    template <int ix1, class T1, class M1, class T2>
    inline std::string TMV_Text(const SumMX<ix1,T1,M1,T2>& smx)
    {
        std::ostringstream s;
        s << "SumMX< "<< ix1<<","<<TMV_Text(T1())
            <<" , "<<TMV_Text(smx.getM1())
            <<" , "<<TMV_Text(T2())<<" >";
        return s.str();
    }
#endif



} // namespace tmv

#endif 
