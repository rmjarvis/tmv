

#ifndef TMV_ProdMM_H
#define TMV_ProdMM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_MultMM_Funcs.h"

//#define XDEBUG_PRODMM

#ifdef XDEBUG_PRODMM
#include <iostream>
#include "TMV_MatrixIO.h"
#endif

namespace tmv {

    //
    // Matrix * Matrix
    //

    // This first one is for matrices that aren't calculated yet.
    // Some BaseMatrix objects are overloaded (e.g. Permutation), but
    // others aren't, in which case they calculate the matrix first
    // and then call MultMM.  If a BaseMatrix_Calc type doesn't overload
    // this function, then an infinite loop will result, so make sure
    // to remember to overload these for any new matrix types
    // Also, this effects the intermediate calculation for things like
    // D = A * B * C, where (A*B) needs a temporary.
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(x,m1.calc(),m2.calc(),m3.mat()); }

    // If everything is a BaseMatrix_Calc, there should be an overload
    // that says how to do the calculation.  If we hit this version,
    // then I'm missing an overload to say how to deal with a particular
    // combination.  This next line gives a compiler error if we ever try 
    // to use this function.
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& , const BaseMatrix_Calc<M1>& , 
        const BaseMatrix_Calc<M2>& , BaseMatrix_Mutable<M3>& )
    { TMVStaticAssert(ix == 999); }

    // On the other hand, if we don't define a special MultEqMM,
    // that just implies that there is no way to do it without the
    // copy.  So just do that here.
    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M2>& m2)
    { 
        MultMM<false>(
            Scaling<1,typename M1::real_type>(),
            (x*m1.mat()).calc(),m2.mat(),m1.mat()); 
    }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <bool add, class M1, class M2, class M3>
    inline void MultMM(
        const BaseMatrix<M1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        MultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.calc(),m2.calc(),m3.mat()); 
    }
    template <bool add, class T, class M1, class M2, class M3>
    inline void MultMM(
        T x, const BaseMatrix<M1>& m1, 
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(Scaling<0,T>(x),m1.calc(),m2.calc(),m3.mat()); }

#ifdef XDEBUG_PRODMM
    template <bool add, int ix, class T, class M1, class M2, class M3>
    static void MultMM_Debug(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();
        Matrix<typename M3::value_type> m3i = m3.mat();
        Matrix<typename M3::value_type> m3c = m3.mat();
        if (!add) m3c.setZero();
        for(int i=0;i<m3.colsize();++i) {
            for(int j=0;j<m3.rowsize();++j) {
                for(int k=0;k<m1.rowsize();++k) {
                    m3c.ref(i,j) += T(x) * m1i.cref(i,k) * m2i.cref(k,j);
                }
            }
        }
        //std::cout<<"Start MultMM XDEBUG"<<std::endl;
        //std::cout<<"add = "<<add<<std::endl;
        //std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3.mat()<<std::endl;
        //std::cout<<"m3i = "<<TMV_Text(m3i)<<"  "<<m3i<<std::endl;
        //std::cout<<"m3c = "<<TMV_Text(m3c)<<"  "<<m3c<<std::endl;
        //std::cout<<"m3c => "<<m3c<<std::endl;

        MultMM<add>(x,m1.mat(),m2.mat(),m3.mat());

        //std::cout<<"m3 => "<<m3.mat()<<std::endl;
        Matrix<typename M3::value_type> m3f = m3.mat();
        //std::cout<<"diff = "<<m3.mat()-m3c<<std::endl;
        //std::cout<<"Norm(diff) = "<<Norm(m3.mat()-m3c)<<std::endl;
        if (Norm(m3f-m3c) > 1.e-6 * std::abs(T(x))*Norm(m1i)*Norm(m2i)) {
            std::cout<<"MultMM:  add = "<<add<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            if (m3.colsize()<100 && m3.rowsize()<100 && m1.colsize()<100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 -> "<<m3f<<std::endl;
                std::cout<<"correct = "<<m3c<<std::endl;
                std::cout<<"diff = "<<(m3c-m3f)<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m3c-m3f)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_PRODMM
    template <class M1, int ix, class T, class M2>
    static void MultEqMM_Debug(
        BaseMatrix_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix<M2>& m2)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M1::value_type> m2i = m2.mat();
        Matrix<typename M1::value_type> m3 = m1.mat();
        m3.setZero();
        for(int i=0;i<m3.colsize();++i) {
            for(int j=0;j<m3.rowsize();++j) {
                for(int k=0;k<m1.rowsize();++k) {
                    m3.ref(i,j) += T(x) * m1i.cref(i,k) * m2i.cref(k,j);
                }
            }
        }

        MultEqMM(m1.mat(),x,m2.mat());

        Matrix<typename M1::value_type> m1f = m1.mat();
        if (Norm(m1f-m3) > 1.e-6 * std::abs(T(x))*Norm(m1i)*Norm(m2i)) {
            std::cout<<"MultEqMM:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            if (m1.colsize() < 100 && m1.rowsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m1 -> "<<m1f<<std::endl;
                std::cout<<"correct = "<<m3<<std::endl;
                std::cout<<"diff = "<<(m3-m1f)<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m3-m1f)<<std::endl;
            exit(1);
        }
    }
#endif

    template <int ix, class T, class M1, class M2>
    class ProdMM;

    template <int ix, class T, class M1, class M2>
    struct Traits<ProdMM<ix,T,M1,M2> >
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        typedef typename Traits2<T,T12>::type value_type;

        enum { _colsize = M1::_colsize };
        enum { _rowsize = M2::_rowsize };
        enum { nlo1 = IntTraits2<M1::_nlo,M2::_nlo>::sum };
        enum { nhi1 = IntTraits2<M1::_nhi,M2::_nlo>::sum };
        enum { csm1 = IntTraits<_colsize>::Sm1 };
        enum { rsm1 = IntTraits<_rowsize>::Sm1 };
        enum { _nlo = IntTraits2<nlo1,csm1>::min };
        enum { _nhi = IntTraits2<nhi1,rsm1>::min };
        enum { shape1 = Traits<ProdXM<ix,T,M2> >::_shape };
        enum { _shape = ShapeTraits2<shape1,M1::_shape>::prod };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { cm2 = Traits<typename M2::calc_type>::_colmajor };
        enum { _rowmajor = 
            (rm1 && rm2) || (!cm1 && !cm2 && _rowsize > int(_colsize)) };

        typedef ProdMM<ix,T,M1,M2> type;
        enum { s = _shape };
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = (
                (_rowmajor ? RowMajor : ColMajor) | 
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename MCopyHelper<value_type,s,cs,rs,A>:: type copy_type;
        typedef const copy_type calc_type;
        enum { caneval = (
                M1::_calc && M2::_calc && 
                M1::_shape == int(Rec) && M2::_shape == int(Rec) ) };
        // Can't use this for non-Rec matrices until I have 
        // get_col(i) and get_row(i) methods for them.
        // This requires having a PartialVector class.
        typedef typename TypeSelect<
            caneval,const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M1, class M2>
    class ProdMM : 
        public BaseMatrix<ProdMM<ix,T,M1,M2> >
    {
    public:

        typedef ProdMM<ix,T,M1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ProdMM(const Scaling<ix,T>& _x, const BaseMatrix<M1>& _m1,
               const BaseMatrix<M2>& _m2) :
            x(_x), m1(_m1.mat()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same)); 
            TMVAssert(m1.rowsize() == m2.colsize());
        }

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const M1& getM1() const { return m1; }
        TMV_INLINE const M2& getM2() const { return m2; }

        TMV_INLINE int colsize() const { return m1.colsize(); }
        TMV_INLINE int rowsize() const { return m2.rowsize(); }
        TMV_INLINE int nlo() const 
        { return TMV_MAX(TMV_MIN(colsize()-1,m1.nlo()+m2.nlo()),0); }
        TMV_INLINE int nhi() const 
        { return TMV_MAX(TMV_MIN(rowsize()-1,m1.nhi()+m2.nhi()),0); }

        value_type cref(int i, int j) const
        { return x * (m1.get_row(i) * m2.get_col(j)); }

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
#ifdef XDEBUG_PRODMM
            MultMM_Debug<false>(x,m1.eval(),m2.eval(),m3.mat());
#else
            MultMM<false>(x,m1.mat(),m2.mat(),m3.mat());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const M1& m1;
        const M2& m2;
    };


    // m * m
#define RT typename M1::real_type
    template <class M1, class M2>
    TMV_INLINE ProdMM<1,RT,M1,M2> operator*(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return ProdMM<1,RT,M1,M2>(RT(1),m1,m2); }
#undef RT

    // m * xm
    template <class M1, int ix, class T, class M2>
    TMV_INLINE ProdMM<ix,T,M1,M2> operator*(
        const BaseMatrix<M1>& m1, const ProdXM<ix,T,M2>& m2)
    { return ProdMM<ix,T,M1,M2>(m2.getX(),m1,m2.getM()); }

    // xm * m
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<ix,T,M1,M2> operator*(
        const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
    { return ProdMM<ix,T,M1,M2>(m1.getX(),m1.getM(),m2); }

    // xm * xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    TMV_INLINE ProdMM<ix1*ix2,PT,M1,M2> operator*(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return ProdMM<ix1*ix2,PT,M1,M2>(
            m1.getX()*m2.getX(),m1.getM(),m2.getM()); 
    }
#undef PT


    // m *= m
#define RT typename M1::real_type
    template <class M1, class M2>
    inline void MultEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_PRODMM
        MultEqMM_Debug(m1.mat(),Scaling<1,RT>(),m2.eval()); 
#else
        MultEqMM(m1.mat(),Scaling<1,RT>(),m2.mat()); 
#endif
    }
#undef RT

    // m *= xm
    template <class M1, int ix2, class T2, class M2>
    inline void MultEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
#ifdef XDEBUG_PRODMM
        MultEqMM_Debug(m1.mat(),m2.getX(),m2.getM().eval()); 
#else
        MultEqMM(m1.mat(),m2.getX(),m2.getM().mat()); 
#endif
    }

    // m += mm
    template <class M3, int ix, class T, class M1, class M2>
    inline void AddEq(
        BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
    {
#ifdef XDEBUG_PRODMM
        MultMM_Debug<true>(
            mm.getX(),mm.getM1().eval(),mm.getM2().eval(),m3.mat()); 
#else
        MultMM<true>(mm.getX(),mm.getM1().mat(),mm.getM2().mat(),m3.mat()); 
#endif
    }

    // m -= mm
    template <class M3, int ix, class T, class M1, class M2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
    {
#ifdef XDEBUG_PRODMM
        MultMM_Debug<true>(
            -mm.getX(),mm.getM1().eval(),mm.getM2().eval(),m3.mat()); 
#else
        MultMM<true>(-mm.getX(),mm.getM1().mat(),mm.getM2().mat(),m3.mat()); 
#endif
    }


    // Consolidate x*(xmm) type constructs:

#define RT typename ProdMM<ix,T,M1,M2>::real_type
#define CT typename ProdMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(mm)
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<-ix,T,M1,M2> operator-(
        const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<-ix,T,M1,M2>(-mm.getX(),mm.getM1(),mm.getM2()); }

    // x * (mm)
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator*(
        const int x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator*(
        const RT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator*(
        const CT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator*(
        const CCT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    // (mm)*x
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const int x)
    { return ProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    // (mm)/x
    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const int x)
    { return ProdMM<0,T,M1,M2>(mm.getX()/RT(x),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,T,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ProdMM<0,T,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    TMV_INLINE ProdMM<0,CT,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }


#undef RT
#undef CT
#undef CCT

    // TMV_Text

    template <int ix, class T, class M1, class M2>
    inline std::string TMV_Text(const ProdMM<ix,T,M1,M2>& mm)
    {
        std::ostringstream s;
        s << "ProdMM< "<<ix<<","<<TMV_Text(T())<<",";
        s << TMV_Text(mm.getM1())<<" , "<<TMV_Text(mm.getM2())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
