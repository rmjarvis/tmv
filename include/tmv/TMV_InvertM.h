

#ifndef TMV_InvertM_H
#define TMV_InvertM_H

#include "TMV_BaseMatrix.h"

#ifdef PRINTALGO_INVM
#include <iostream>
#endif

namespace tmv {

    //
    // minv = x * m^-1
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<0,cs,rs,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {}
    };

    // algo 1: size = 1x1
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<1,1,1,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 1: M,N,cs,rs = "<<M<<','<<N<<','<<
                1<<','<<1<<std::endl;
#endif
            typedef typename M1::value_type T1;
            if (m1.cref(0,0) == T1(0)) ThrowSingular("1x1 Matrix");
            m2.ref(0,0) = ZProd<false,false>::quot(x , m1.cref(0,0)); 
        }
    };

    // algo 2: size = 2x2
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<2,2,2,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 2: M,N,cs,rs = "<<M<<','<<N<<','<<
                2<<','<<2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M1::real_type RT;
            // Find scale to help avoid overflow/underflow.
            RT scale = m1.maxAbs2Element();
            typename MCopyHelper<T1,Rec,2,2>::type m1c = m1/scale;
            T1 det = m1c.det();
            if (det == T1(0)) ThrowSingular("2x2 Matrix");
            T1 invdet = ZProd<false,false>::quot(RT(1) , det);
            // Each of the terms like m1c.cref(0,0) has a 1/scale already
            // so only need one more scale in the xinvdet factor.
            T2 xinvdet = ZProd<false,false>::prod(
                x,ZProd<false,false>::quot(invdet,scale));
            m2.ref(0,0) = ZProd<false,false>::prod(m1c.cref(1,1) , xinvdet);
            m2.ref(0,1) = ZProd<false,false>::prod(-m1c.cref(0,1) , xinvdet);
            m2.ref(1,0) = ZProd<false,false>::prod(-m1c.cref(1,0) , xinvdet);
            m2.ref(1,1) = ZProd<false,false>::prod(m1c.cref(0,0) , xinvdet);
        }
    };

    // algo 3: size = 3x3
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<3,3,3,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 3: M,N,cs,rs = "<<M<<','<<N<<','<<
                3<<','<<3<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M1::real_type RT;
            // Find scale to help avoid overflow/underflow.
            RT scale = m1.maxAbs2Element();
            typename MCopyHelper<T1,Rec,3,3>::type m1c = m1/scale;
            T1 det = m1c.det();
            if (det == T1(0)) ThrowSingular("3x3 Matrix");
            T1 invdet = ZProd<false,false>::quot(RT(1) , det);
            const T1 ae = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(1,1));
            const T1 ai = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(2,2));
            const T1 ei = ZProd<false,false>::prod(m1c.cref(1,1),m1c.cref(2,2));
            const T1 af = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(1,2));
            const T1 ah = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(2,1));
            const T1 fh = ZProd<false,false>::prod(m1c.cref(1,2),m1c.cref(2,1));
            const T1 bf = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(1,2));
            const T1 bg = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(2,0));
            const T1 fg = ZProd<false,false>::prod(m1c.cref(1,2),m1c.cref(2,0));
            const T1 bd = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(1,0));
            const T1 bi = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(2,2));
            const T1 di = ZProd<false,false>::prod(m1c.cref(1,0),m1c.cref(2,2));
            const T1 cd = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(1,0));
            const T1 ch = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(2,1));
            const T1 dh = ZProd<false,false>::prod(m1c.cref(1,0),m1c.cref(2,1));
            const T1 ce = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(1,1));
            const T1 cg = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(2,0));
            const T1 eg = ZProd<false,false>::prod(m1c.cref(1,1),m1c.cref(2,0));
            // Each of the terms like ei-fh has two 1/scales already
            // so only need one more scale in the xinvdet factor.
            T2 xinvdet = ZProd<false,false>::prod(
                x,ZProd<false,false>::quot(invdet,scale));
            m2.ref(0,0) = ZProd<false,false>::prod(ei-fh , xinvdet);
            m2.ref(0,1) = ZProd<false,false>::prod(ch-bi , xinvdet);
            m2.ref(0,2) = ZProd<false,false>::prod(bf-ce , xinvdet);
            m2.ref(1,0) = ZProd<false,false>::prod(fg-di , xinvdet);
            m2.ref(1,1) = ZProd<false,false>::prod(ai-cg , xinvdet);
            m2.ref(1,2) = ZProd<false,false>::prod(cd-af , xinvdet);
            m2.ref(2,0) = ZProd<false,false>::prod(dh-eg , xinvdet);
            m2.ref(2,1) = ZProd<false,false>::prod(bg-ah , xinvdet);
            m2.ref(2,2) = ZProd<false,false>::prod(ae-bd , xinvdet);
        }
    };

    // algo 4: size = 4x4
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<4,4,4,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 4: M,N,cs,rs = "<<M<<','<<N<<','<<
                4<<','<<4<<std::endl;
#endif
            // This algorithm does the direct calculation from 
            // the ajugate matrix divided by the determinant.
            // The letters refer to describing the matrix as:
            // ( a b c d )
            // ( e f g h )
            // ( i j k l )
            // ( m n o p )
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M1::real_type RT;
            
            // Find scale to help avoid overflow/underflow.
            RT scale = m1.maxAbs2Element();
            typename MCopyHelper<T1,Rec,4,4>::type m1c = m1/scale;
            const T1 af = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(1,1));
            const T1 ag = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(1,2));
            const T1 ah = ZProd<false,false>::prod(m1c.cref(0,0),m1c.cref(1,3));
            const T1 be = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(1,0));
            const T1 bg = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(1,2));
            const T1 bh = ZProd<false,false>::prod(m1c.cref(0,1),m1c.cref(1,3));
            const T1 ce = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(1,0));
            const T1 cf = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(1,1));
            const T1 ch = ZProd<false,false>::prod(m1c.cref(0,2),m1c.cref(1,3));
            const T1 de = ZProd<false,false>::prod(m1c.cref(0,3),m1c.cref(1,0));
            const T1 df = ZProd<false,false>::prod(m1c.cref(0,3),m1c.cref(1,1));
            const T1 dg = ZProd<false,false>::prod(m1c.cref(0,3),m1c.cref(1,2));
            const T1 in = ZProd<false,false>::prod(m1c.cref(2,0),m1c.cref(3,1));
            const T1 io = ZProd<false,false>::prod(m1c.cref(2,0),m1c.cref(3,2));
            const T1 ip = ZProd<false,false>::prod(m1c.cref(2,0),m1c.cref(3,3));
            const T1 jm = ZProd<false,false>::prod(m1c.cref(2,1),m1c.cref(3,0));
            const T1 jo = ZProd<false,false>::prod(m1c.cref(2,1),m1c.cref(3,2));
            const T1 jp = ZProd<false,false>::prod(m1c.cref(2,1),m1c.cref(3,3));
            const T1 km = ZProd<false,false>::prod(m1c.cref(2,2),m1c.cref(3,0));
            const T1 kn = ZProd<false,false>::prod(m1c.cref(2,2),m1c.cref(3,1));
            const T1 kp = ZProd<false,false>::prod(m1c.cref(2,2),m1c.cref(3,3));
            const T1 lm = ZProd<false,false>::prod(m1c.cref(2,3),m1c.cref(3,0));
            const T1 ln = ZProd<false,false>::prod(m1c.cref(2,3),m1c.cref(3,1));
            const T1 lo = ZProd<false,false>::prod(m1c.cref(2,3),m1c.cref(3,2));

            // Calculate the Inverse * Det
            const T1 r00 = 
                ZProd<false,false>::prod(m1c.cref(1,1),(kp-lo))
                + ZProd<false,false>::prod(m1c.cref(1,2),(ln-jp))
                + ZProd<false,false>::prod(m1c.cref(1,3),(jo-kn));
            const T1 r10 = 
                ZProd<false,false>::prod(m1c.cref(1,0),(lo-kp))
                + ZProd<false,false>::prod(m1c.cref(1,2),(ip-lm))
                + ZProd<false,false>::prod(m1c.cref(1,3),(km-io));
            const T1 r20 = 
                ZProd<false,false>::prod(m1c.cref(1,0),(jp-ln))
                + ZProd<false,false>::prod(m1c.cref(1,1),(lm-ip))
                + ZProd<false,false>::prod(m1c.cref(1,3),(in-jm));
            const T1 r30 = 
                ZProd<false,false>::prod(m1c.cref(1,0),(kn-jo))
                + ZProd<false,false>::prod(m1c.cref(1,1),(io-km))
                + ZProd<false,false>::prod(m1c.cref(1,2),(jm-in));

            // Interrupt to calculate det now, to make sure it's not 0.
            // We calculate the determinant directly here, since the
            // above intermediate values are useful for calculating the 
            // determinant.  So it saves a few multiplies compared to
            // calling det().
            const T1 det = 
                ZProd<false,false>::prod(m1c.cref(0,0),r00)
                + ZProd<false,false>::prod(m1c.cref(0,1),r10)
                + ZProd<false,false>::prod(m1c.cref(0,2),r20)
                + ZProd<false,false>::prod(m1c.cref(0,3),r30);
            if (det == T1(0)) ThrowSingular("4x4 Matrix");
            T1 invdet = ZProd<false,false>::quot(RT(1) , det);

            const T1 r01 = 
                ZProd<false,false>::prod(m1c.cref(0,1),(lo-kp))
                + ZProd<false,false>::prod(m1c.cref(0,2),(jp-ln))
                + ZProd<false,false>::prod(m1c.cref(0,3),(kn-jo));
            const T1 r11 = 
                ZProd<false,false>::prod(m1c.cref(0,0),(kp-lo))
                + ZProd<false,false>::prod(m1c.cref(0,2),(lm-ip))
                + ZProd<false,false>::prod(m1c.cref(0,3),(io-km));
            const T1 r21 = 
                ZProd<false,false>::prod(m1c.cref(0,0),(ln-jp))
                + ZProd<false,false>::prod(m1c.cref(0,1),(ip-lm))
                + ZProd<false,false>::prod(m1c.cref(0,3),(jm-in));
            const T1 r31 = 
                ZProd<false,false>::prod(m1c.cref(0,0),(jo-kn))
                + ZProd<false,false>::prod(m1c.cref(0,1),(km-io))
                + ZProd<false,false>::prod(m1c.cref(0,2),(in-jm));
            const T1 r02 = 
                ZProd<false,false>::prod(m1c.cref(3,1),(ch-dg))
                + ZProd<false,false>::prod(m1c.cref(3,2),(df-bh))
                + ZProd<false,false>::prod(m1c.cref(3,3),(bg-cf));
            const T1 r12 = 
                ZProd<false,false>::prod(m1c.cref(3,0),(dg-ch))
                + ZProd<false,false>::prod(m1c.cref(3,2),(ah-de))
                + ZProd<false,false>::prod(m1c.cref(3,3),(ce-ag));
            const T1 r22 = 
                ZProd<false,false>::prod(m1c.cref(3,0),(bh-df))
                + ZProd<false,false>::prod(m1c.cref(3,1),(de-ah))
                + ZProd<false,false>::prod(m1c.cref(3,3),(af-be));
            const T1 r32 = 
                ZProd<false,false>::prod(m1c.cref(3,0),(cf-bg))
                + ZProd<false,false>::prod(m1c.cref(3,1),(ag-ce))
                + ZProd<false,false>::prod(m1c.cref(3,2),(be-af));
            const T1 r03 = 
                ZProd<false,false>::prod(m1c.cref(2,1),(dg-ch))
                + ZProd<false,false>::prod(m1c.cref(2,2),(bh-df))
                + ZProd<false,false>::prod(m1c.cref(2,3),(cf-bg));
            const T1 r13 = 
                ZProd<false,false>::prod(m1c.cref(2,0),(ch-dg))
                + ZProd<false,false>::prod(m1c.cref(2,2),(de-ah))
                + ZProd<false,false>::prod(m1c.cref(2,3),(ag-ce));
            const T1 r23 = 
                ZProd<false,false>::prod(m1c.cref(2,0),(df-bh))
                + ZProd<false,false>::prod(m1c.cref(2,1),(ah-de))
                + ZProd<false,false>::prod(m1c.cref(2,3),(be-af));
            const T1 r33 = 
                ZProd<false,false>::prod(m1c.cref(2,0),(bg-cf))
                + ZProd<false,false>::prod(m1c.cref(2,1),(ce-ag))
                + ZProd<false,false>::prod(m1c.cref(2,2),(af-be));

            // Finally divide rij values by det and multiply by x:
            // Each of the terms like r00 has two 1/scales already
            // so only need one more scale in the xinvdet factor.
            T2 xinvdet = ZProd<false,false>::prod(
                x,ZProd<false,false>::quot(invdet,scale));
            m2.ref(0,0) = ZProd<false,false>::prod(r00 , xinvdet);
            m2.ref(1,0) = ZProd<false,false>::prod(r10 , xinvdet);
            m2.ref(2,0) = ZProd<false,false>::prod(r20 , xinvdet);
            m2.ref(3,0) = ZProd<false,false>::prod(r30 , xinvdet);
            m2.ref(0,1) = ZProd<false,false>::prod(r01 , xinvdet);
            m2.ref(1,1) = ZProd<false,false>::prod(r11 , xinvdet);
            m2.ref(2,1) = ZProd<false,false>::prod(r21 , xinvdet);
            m2.ref(3,1) = ZProd<false,false>::prod(r31 , xinvdet);
            m2.ref(0,2) = ZProd<false,false>::prod(r02 , xinvdet);
            m2.ref(1,2) = ZProd<false,false>::prod(r12 , xinvdet);
            m2.ref(2,2) = ZProd<false,false>::prod(r22 , xinvdet);
            m2.ref(3,2) = ZProd<false,false>::prod(r32 , xinvdet);
            m2.ref(0,3) = ZProd<false,false>::prod(r03 , xinvdet);
            m2.ref(1,3) = ZProd<false,false>::prod(r13 , xinvdet);
            m2.ref(2,3) = ZProd<false,false>::prod(r23 , xinvdet);
            m2.ref(3,3) = ZProd<false,false>::prod(r33 , xinvdet);
        }
    };

    // algo 5: size = 4x4 using partitioning
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<5,4,4,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 5: M,N,cs,rs = "<<M<<','<<N<<','<<
                4<<','<<4<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            // This uses an alternate strategy for calculating the inverse.
            // ( P Q )^-1 = ( T U )
            // ( R S )      ( V W )
            // W = (S - R P^-1 Q)^-1
            // V = -W R P^-1
            // U = -P^-1 Q W
            // T = P^-1 - P^-1 Q U
            SmallMatrix<T1,2,2,RowMajor> P; 
            P << m1.cref(0,0) , m1.cref(0,1) , m1.cref(1,0) , m1.cref(1,1);
            SmallMatrix<T1,2,2,RowMajor> S; 
            S << m1.cref(2,2) , m1.cref(2,3) , m1.cref(3,2) , m1.cref(3,3);

            // If P and S are both singular, m1 might not.  
            // Use the more general formula instead.
            T1 detP = P.det();
            T1 detS = S.det();
            if (detP == T1(0) && detS == T1(0)) 
                InvertM_Helper<4,4,4,ix,T,M1,M2>::call(x,m1,m2);

            SmallMatrix<T1,2,2,RowMajor> Q; 
            Q << m1.cref(0,2) , m1.cref(0,3) , m1.cref(1,2) , m1.cref(1,3);
            SmallMatrix<T1,2,2,RowMajor> R; 
            R << m1.cref(2,0) , m1.cref(2,1) , m1.cref(3,0) , m1.cref(3,1);

            SmallMatrix<T1,2,2> W;
            SmallMatrix<T1,2,2> X;
            SmallMatrix<T1,2,2> Y;
            SmallMatrix<T1,2,2> Z;

            if (TMV_ABS(detP) > TMV_ABS(detS)) {
                SmallMatrix<T1,2,2> Pinv = P.inverse();
                SmallMatrix<T1,2,2> PinvQ = Pinv * Q;
                SmallMatrix<T1,2,2> V = S-R*PinvQ;

                T1 detV = V.det();
                if (detV == T1(0)) ThrowSingular("4x4 Matrix");

                Z = V.inverse();
                SmallMatrix<T1,2,2> RPinv = R * Pinv;
                Y = -Z*RPinv;
                X = -PinvQ * Z;
                W = Pinv - PinvQ * Y;
            } else {
                SmallMatrix<T1,2,2> Sinv = S.inverse();
                SmallMatrix<T1,2,2> SinvR = Sinv * R;
                SmallMatrix<T1,2,2> V = P-Q*SinvR;

                T1 detV = V.det();
                if (detV == T1(0)) ThrowSingular("4x4 Matrix");

                W = V.inverse();
                SmallMatrix<T1,2,2> QSinv = Q * Sinv;
                X = -W*QSinv;
                Y = -SinvR * W;
                Z = Sinv - SinvR * X;
            }
            m2.ref(0,0) = W.cref(0,0);
            m2.ref(0,1) = W.cref(0,1);
            m2.ref(1,0) = W.cref(1,0);
            m2.ref(1,1) = W.cref(1,1);
            m2.ref(2,0) = Y.cref(0,0);
            m2.ref(2,1) = Y.cref(0,1);
            m2.ref(3,0) = Y.cref(1,0);
            m2.ref(3,1) = Y.cref(1,1);
            m2.ref(0,2) = X.cref(0,0);
            m2.ref(0,3) = X.cref(0,1);
            m2.ref(1,2) = X.cref(1,0);
            m2.ref(1,3) = X.cref(1,1);
            m2.ref(2,2) = Z.cref(0,0);
            m2.ref(2,3) = Z.cref(0,1);
            m2.ref(3,2) = Z.cref(1,0);
            m2.ref(3,3) = Z.cref(1,1);
        }
    };



    // algo 11: Use Divider 
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<11,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 11: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"divIsSet = "<<m1.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m1.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m1.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m1.getDivType())<<std::endl;
#endif
            m1.setDiv();
            m1.getDiv()->makeInverse(m2);
            m1.doneDiv();
            Scale(x,m2);
        }
    };

    // algo 12: Calculate LU decomposition on the spot.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<12,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 12: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            // This way is slightly faster than going through m1.lud()
            // since it skips the temporary LU matrix.
            Copy(m1,m2);
            Permutation P(m1.rowsize());
            LU_Decompose(m2,P);
            LU_Inverse(m2,P);
            Scale(x,m2);
        }
    };

    // algo 13: Calculate QR decomposition on the spot.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<13,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 13: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m1.qrd().makeInverse(m2);
            Scale(x,m2);
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<14,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 14: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            if (m1.isSquare())
                InvertM_Helper<12,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
            else
                InvertM_Helper<13,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 31: m1 is diagonal.  No alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<31,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 31: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typename M2::diag_type m2d = m2.diag();
            m2.setZero();
            m2d.noAlias() = m1.diag();
            DiagMatrixViewOf(m2d).invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 32: m1 is diagonal.  Safe for aliases.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<32,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 32: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            // For alias versions, we could invert the diagonal and then 
            // zero out the non-diagonal portions, but the non-diagonal part
            // is different for each kind of M2, so to make this code
            // completely generic, we use a temporary copy of m1 which we copy
            // into the diagonal after doing m2.setZero().
            // This is also probably more efficient, since zeroing two
            // triangle matrices is slower than zeroing a contiguous 
            // rectangle matrix.  I don't know if the difference is more or
            // less than the copies required.
            typename M1::copy_type m1c = m1;
            m1c.invertSelf();
            typename M1::copy_type::diag_type m1cd = m1c.diag();
            typename M2::diag_type m2d = m2.diag();
            Scale(x,m1cd);
            m2.setZero();
            m2d.noAlias() = m1cd;
        }
    };

    // algo 33: m1,m2 are both diagonal.  No alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<33,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 33: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typename M2::diag_type m2d = m2.diag();
            m2d.noAlias() = m1.diag();
            m2.invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 34: m1,m2 are both diagonal.  With alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<34,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 34: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typename M2::diag_type m2d = m2.diag();
            Copy(m1.diag(),m2d);
            m2.invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 41: m1 is triangular.  No alias.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<41,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 41: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typename TypeSelect<M1::_upper ,
                     typename M2::uppertri_type ,
                     typename M2::lowertri_type>::type m2u = 
                         Maybe<M1::_upper>::uppertri(m2);
            m2.setZero();
            m2u.noAlias() = m1;
            m2u.invertSelf();
            Scale(x,m2u);
        }
    };

    // algo 42: m1 is triangular.  With alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<42,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const ptrdiff_t N = m1.rowsize();
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            std::cout<<"InvM algo 42: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typename TypeSelect<M1::_upper ,
                     typename M2::uppertri_type ,
                     typename M2::lowertri_type>::type m2u = 
                         Maybe<M1::_upper>::uppertri(m2);
            Copy(m1,m2u);
            m2u.invertSelf();
            Scale(x,m2u);
            if (N > 1) Maybe<!M1::_upper>::uppertri(m2).offDiag().setZero();
        }
    };

    // algo 43: m1,m2 are both uppertri.  No alias
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<43,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 43: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m2.noAlias() = m1;
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 44: m1,m2 are both uppertri.  With alias check
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<44,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvM algo 44: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            Copy(m1,m2);
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<99,cs,rs,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 32 : 34 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 42 : 44 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 42 : 44 ) :
                cs == 2 && rs == 2 ? 2 :
                cs == 3 && rs == 3 ? 3 :
                cs == 4 && rs == 4 ? 4 :
                M1::_hasdivider ? 11 :
                cs == Unknown || rs == Unknown ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_INVM
            std::cout<<"AliasCheck InvertM:\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-3,cs,rs,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 31 : 33 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 41 : 43 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 41 : 43 ) :
                cs == 2 && rs == 2 ? 2 :
                cs == 3 && rs == 3 ? 3 :
                cs == 4 && rs == 4 ? 4 :
                M1::_hasdivider ? 11 :
                cs == Unknown || rs == Unknown ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_INVM
            std::cout<<"Inline InvertM\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"cs,rs, = "<<cs<<','<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -2: No alias
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-2,cs,rs,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { InvertM_Helper<-3,cs,rs,ix,T,M1,M2>::call(x,m1,m2); }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-1,cs,rs,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                M2::_checkalias ? 99 :
                -3;
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const ptrdiff_t cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const ptrdiff_t rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InvertM_Helper<-1,cs,rs,ix,T,M1,M2v>::call(x,m1.mat(),m2v);
    }


    //
    // mata = (at * a)^-1
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<0,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& , M2& )
        {}
    };

    // algo 1: size = 1x1
    template <class M1, class M2>
    struct InverseATA_Helper<1,1,1,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 1: M,N,cs,rs = "<<M<<','<<N<<','<<
                1<<','<<1<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            if (m1.cref(0,0) == T1(0)) ThrowSingular("1x1 Matrix");
            T1 inva00 = ZProd<false,false>::quot(RT(1) , m1.cref(0,0));
            m2.ref(0,0) = ZProd<false,false>::prod(inva00 , inva00);
        }
    };

    // algo 2: Direct (mtm)^-1
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<2,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 2: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(rs != Unknown);
            TMVStaticAssert(cs == Unknown || rs <= cs);
            TMVAssert(m1.rowsize() <= m1.colsize());
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            SmallMatrix<T1,rs,rs,NoAlias> ata;
            MultMM<false>(Scaling<1,RT>(),m1.adjoint(),m1,ata);
            typename M2::noalias_type m2na = m2.noAlias();
            MakeInverse(Scaling<1,RT>(),ata,m2na);
        }
    };

    // algo 3: Direct (mmt)^-1
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<3,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 3: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(cs != Unknown);
            TMVStaticAssert(rs == Unknown || cs < rs);
            TMVAssert(m1.colsize() < m1.rowsize());
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            SmallMatrix<T1,cs,cs,NoAlias> ata;
            MultMM<false>(Scaling<1,RT>(),m1,m1.adjoint(),ata);
            typename M2::noalias_type m2na = m2.noAlias();
            MakeInverse(Scaling<1,RT>(),ata,m2na);
        }
    };

    // algo 11: Use Divider
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<11,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 11: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"divIsSet = "<<m1.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m1.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m1.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m1.getDivType())<<std::endl;
#endif
            m1.setDiv();
            m1.getDiv()->makeInverseATA(m2);
            m1.doneDiv();
        }
    };

    // algo 12: Calculate LU decomposition on the spot. 
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<12,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 12: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m1.lud().makeInverseATA(m2);
        }
    };

    // algo 13: Calculate QR decomposition on the spot. 
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<13,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 13: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m1.qrd().makeInverseATA(m2);
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<14,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 14: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            if (m1.isSquare())
                InverseATA_Helper<12,cs,rs,M1,M2>::call(m1,m2);
            else
                InverseATA_Helper<13,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 21: m1 is diagonal.  No alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<21,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 21: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            if (m1.isSingular()) ThrowSingular("DiagMatrix");
            typename M2::diag_type::noalias_type m2d = m2.diag().noAlias();
            m2.setZero();
            Copy(m1.diag(),m2d);
            ElemInvert(m2d);
            ElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2d);
        }
    };

    // algo 22: m1 is diagonal.  Safe for aliases.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<22,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 22: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            if (m1.isSingular()) ThrowSingular("DiagMatrix");
            typename M1::const_diag_type::copy_type m1d(m1.size());
            typename M2::diag_type::noalias_type m2d = m2.diag();
            Copy(m1.diag(),m1d);
            ElemInvert(m1d);
            ElemMultVV<false>(Scaling<1,RT>(),m1d.conjugate(),m1d,m1d);
            m2.setZero();
            Copy(m1d,m2d);
        }
    };

    // algo 23: m1,m2 are both diagonal.  No alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<23,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 23: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("DiagMatrix");
            typename M2::diag_type::noalias_type m2d = m2.diag().noAlias();
            Copy(m1.diag(),m2d);
            ElemInvert(m2d);
            ElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2d);
        }
    };

    // algo 24: m1,m2 are both diagonal.  With alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<24,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 24: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("DiagMatrix");
            typename M2::diag_type m2d = m2.diag();
            Copy(m1.diag(),m2d);
            ElemInvert(m2d);
            typename M2::diag_type::noalias_type m2dna = m2d.noAlias();
            ElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2dna);
        }
    };

    // algo 31: m1 is uppertri.  No alias.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<31,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 31: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("TriMatrix");
            typename M2::uppertri_type m2u = m2.upperTri();
            m2u.noAlias() = m1;
            InvertSelf(m2u);
            typename M2::noalias_type m2na = m2.noAlias();
            MultMM<false>(Scaling<1,RT>(),m2u,m2u.adjoint(),m2na);
        }
    };

    // algo 32: m1 is uppertri.  With alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<32,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 32: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("TriMatrix");
            typename M2::uppertri_type m2u = m2.upperTri();
            Copy(m1,m2u);
            InvertSelf(m2u);
            typename M2::noalias_type m2na = m2.noAlias();
            MultMM<false>(Scaling<1,RT>(),m2u,m2u.adjoint(),m2na);
        }
    };

    // algo 41: m1 is lowertri.  No alias.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<41,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 41: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("TriMatrix");
            typename M2::lowertri_type m2l = m2.lowerTri();
            m2l.noAlias() = m1;
            InvertSelf(m2l);
            typename M2::noalias_type m2na = m2.noAlias();
            MultMM<false>(Scaling<1,RT>(),m2l,m2l.adjoint(),m2na);
        }
    };

    // algo 42: m1 is lowertri.  With alias check.
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<42,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"InvATA algo 42: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) ThrowSingular("TriMatrix");
            typename M2::lowertri_type m2l = m2.lowerTri();
            Copy(m1,m2l);
            InvertSelf(m2l);
            typename M2::noalias_type m2na = m2.noAlias();
            MultMM<false>(Scaling<1,RT>(),m2l,m2l.adjoint(),m2na);
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 22 : 24 ) :
                !lo1 ? 32 : // m1 is uppertri
                !up1 ? 42 : // m1 is lowertri
                (rs != Unknown && cs != Unknown && rs <= 4 && rs <= cs) ? 2 :
                (rs != Unknown && cs != Unknown && cs <= 4) ? 2 :
                M1::_hasdivider ? 11 :
                cs == Unknown || rs == Unknown ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"AliasCheck InverseATA: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 21 : 23 ) :
                !lo1 ? 31 : // m1 is uppertri
                !up1 ? 41 : // m1 is lowertri
                (rs != Unknown && cs != Unknown && rs <= 4 && rs <= cs) ? 2 :
                (rs != Unknown && cs != Unknown && cs <= 4) ? 2 :
                M1::_hasdivider ? 11 :
                cs == Unknown || rs == Unknown ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_INVM
            const ptrdiff_t M = m1.colsize();
            const ptrdiff_t N = m1.rowsize();
            std::cout<<"Inline InverseATA: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: No alias
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        { InverseATA_Helper<-3,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct InverseATA_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                M2::_checkalias ? 99 :
                -3;
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<IntTraits2<M1::_colsize,M1::_rowsize>::min,
                         M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M2::_colsize>::same));
        TMVAssert(TMV_MIN(m1.colsize(),m1.rowsize()) == m2.rowsize());
        TMVAssert(m2.rowsize() == m2.colsize());
        const ptrdiff_t cs = M1::_colsize;
        const ptrdiff_t rs = M1::_rowsize;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InverseATA_Helper<-1,cs,rs,M1,M2v>::call(m1.mat(),m2v);
    }

} // namespace tmv

#endif 
