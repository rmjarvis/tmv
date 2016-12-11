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

#ifndef SmallMatrixDiv_H
#define SmallMatrixDiv_H

#include "tmv/TMV_SimpleMatrix.h"

namespace tmv {

    //
    // Helper Functions
    // These implement the various Triangle Matrix calculations
    // that are needed for the decompositions and solving routines.
    // We only implement ColMajor versions for L, U, since these are
    // all that are needed in the routines here.
    //

#define TMV_Index(M,N,S,i,j) (S==RowMajor ? (i*N+j) : (j*M+i))
#define TMV_Val(m,M,N,S,i,j) (m[TMV_Index(M,N,S,i,j)])
#define TMV_TransOf(S) (S==RowMajor ? ColMajor : RowMajor)

    template <ptrdiff_t M, ptrdiff_t N, StorageType S, typename T>
    inline void PrintSM(std::ostream& os, const T* m)
    {
        for(int i=0;i<M;++i) {
            os << "( ";
            for(int j=0;j<N;++j) {
                os << " " << TMV_Val(m,M,N,S,i,j) << " ";
            }
            os << " )\n";
        }
    }

#ifndef NOTHROW
    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    class SingularSmallMatrix : public Singular
    {
    public:
        SmallMatrix<T,M,N,A> m;

        inline SingularSmallMatrix(const SmallMatrix<T,M,N,A>& _m) :
            Singular("SmallMatrix."), m(_m) {}
        inline ~SingularSmallMatrix() throw() {}
        inline void write(std::ostream& os) const throw()
        { Singular::write(os); os<<m<<std::endl; }
    };
#endif

    template <ptrdiff_t N2, ptrdiff_t K, StorageType S2, typename T, typename T2, ptrdiff_t N, ptrdiff_t M>
    inline void LDivEq_L(const SmallMatrix<T,N,M,ColMajor>& L, T2* m)
    // SmallMatrix<T2,N2,K,S2>& m
    {
        //std::cout<<"LDivEq_L:\n";
        //std::cout<<"L = "<<TMV_Text(L)<<"  "<<L<<std::endl;
        //std::cout<<"m = "<<m<<std::endl;
        //PrintSM<N2,K,S2>(std::cout,m);
        //std::cout<<"S2 = "<<TMV_Text(S2)<<std::endl;
        //std::cout<<"M,N = "<<M<<','<<N<<std::endl;
        //std::cout<<"N2,K = "<<N2<<','<<K<<std::endl;
        TMVAssert(M >= N);
        TMVAssert(N2 >= N);
        //m.view() /= L.lowerTri(UnitDiag);
        for(ptrdiff_t j=0;j<N;++j) {
            for(ptrdiff_t k=0;k<K;++k) {
                for(ptrdiff_t i=j+1;i<N;++i) {
                    //m.ref(i,k) -= m.cref(j,k) * L.cref(i,j);
                    TMV_Val(m,N2,K,S2,i,k) -=
                        TMV_Val(m,N2,K,S2,j,k) * L.cref(i,j);
                }
            }
        }
    }

    template <ptrdiff_t N2, ptrdiff_t K, StorageType S2, typename T, typename T2, ptrdiff_t N, ptrdiff_t M>
    inline void LDivEq_U(const SmallMatrix<T,M,N,ColMajor>& U, T2* m)
    // SmallMatrix<T2,N2,K,S2>& m
    {
        //std::cout<<"LDivEq_U:\n";
        //std::cout<<"U = "<<TMV_Text(U)<<"  "<<U<<std::endl;
        //std::cout<<"m = "<<m<<std::endl;
        //PrintSM<N2,K,S2>(std::cout,m);
        //std::cout<<"S2 = "<<TMV_Text(S2)<<std::endl;
        //std::cout<<"M,N = "<<M<<','<<N<<std::endl;
        //std::cout<<"N2,K = "<<N2<<','<<K<<std::endl;
        TMVAssert(M >= N);
        TMVAssert(N2 >= N);
        //m.view() /= U.upperTri();
        for(ptrdiff_t j=N-1;j>=0;--j) {
            //std::cout<<"j = "<<j<<std::endl;
            //std::cout<<"U(j,j) = "<<U.cref(j,j)<<std::endl;
            if (U.cref(j,j) == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            //std::cout<<"After check for singular"<<std::endl;
            for(ptrdiff_t k=0;k<K;++k) {
                //std::cout<<"k = "<<k<<std::endl;
                //m.ref(j,k) /= U.cref(j,j);
                TMV_Val(m,N2,K,S2,j,k) /= U.cref(j,j);
                //std::cout<<"After m.ref(j,k) /= U.cref(j,j)"<<std::endl;
                for(ptrdiff_t i=0;i<j;++i) {
                    //std::cout<<"i = "<<i<<std::endl;
                    //std::cout<<"m(i,k) = "<<TMV_Val(m,N2,K,S2,i,k)<<std::endl;
                    //std::cout<<"m(j,k) = "<<TMV_Val(m,N2,K,S2,j,k)<<std::endl;
                    //std::cout<<"U(i,j) = "<<U.cref(i,j)<<std::endl;
                    //std::cout<<"index m(i,k) = "<<TMV_Index(N2,K,S2,i,k)<<std::endl;
                    //std::cout<<"index m(j,k) = "<<TMV_Index(N2,K,S2,j,k)<<std::endl;
                    //std::cout<<"&m(i,k) = "<<&(TMV_Val(m,N2,K,S2,i,k))<<std::endl;
                    //std::cout<<"&m(j,k) = "<<&(TMV_Val(m,N2,K,S2,j,k))<<std::endl;
                    //m.ref(i,k) -= m.cref(j,k) * U.cref(i,j);
                    TMV_Val(m,N2,K,S2,i,k) -=
                        TMV_Val(m,N2,K,S2,j,k) * U.cref(i,j);
                    //std::cout<<"After m.ref(i,k) -= m.cref(j,k) * U.cref(i,j)"<<std::endl;
                }
            }
        }
    }

    template <ptrdiff_t K, ptrdiff_t N2, StorageType S2, typename T, typename T2, ptrdiff_t N, ptrdiff_t M>
    inline void RDivEq_U(const SmallMatrix<T,M,N,ColMajor>& U, T2* m)
    // SmallMatrix<T2,K,N2,S2>& m
    {
        TMVAssert(M >= N);
        TMVAssert(N2 >= N);
        //m.view() %= U.upperTri();
        for(ptrdiff_t i=0;i<N;++i) {
            if (U.cref(i,i) == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            for(ptrdiff_t k=0;k<K;++k) {
                for(ptrdiff_t j=0;j<i;++j) {
                    //m.ref(k,i) -= m.cref(k,j) * U.cref(j,i);
                    TMV_Val(m,K,N2,S2,k,i) -=
                        TMV_Val(m,K,N2,S2,k,j) * U.cref(j,i);
                }
                //m.ref(k,i) /= U.cref(i,i);
                TMV_Val(m,K,N2,S2,k,i) /= U.cref(i,i);
            }
        }
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N>
    inline void InvertU(SmallMatrix<T,M,N,ColMajor>& U)
    {
        // Invert the upper triangle portion of U in place
        TMVAssert(M >= N);
        for(ptrdiff_t j=0;j<N;++j) {
            if (U.cref(j,j) == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            U.ref(j,j) = TMV_RealType(T)(1)/U.cref(j,j);
            //U.col(j,0,j) = -U.SubTriMatrix(0,j) * U.col(j,0,j) * U(j,j);
            for(ptrdiff_t k=0;k<j;++k) {
                for(ptrdiff_t i=0;i<k;++i) U.ref(i,j) += U.cref(k,j) * U.cref(i,k);
                U.ref(k,j) *= U.cref(k,k);
            }
            for(ptrdiff_t i=0;i<j;++i) U.ref(i,j) *= -U.cref(j,j);
        }
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N>
    inline void InvertL(SmallMatrix<T,N,M,ColMajor>& L)
    {
        // Invert the unit lower triangle portion of L in place
        TMVAssert(M >= N);
        for(ptrdiff_t j=N-1;j>=0;--j) {
            //L.col(j,j+1,N) = -L.SubTriMatrix(j+1,N) * L.col(j,j+1,N);
            for(ptrdiff_t k=N-1;k>j;--k) {
                for(ptrdiff_t i=N-1;i>k;--i) L.ref(i,j) += L.cref(k,j) * L.cref(i,k);
            }
            for(ptrdiff_t i=N-1;i>j;--i) L.ref(i,j) = -L.cref(i,j);
        }
    }

    template <ptrdiff_t K, StorageType S2, typename T, typename T2, ptrdiff_t M, ptrdiff_t N>
    inline void Q_LDivEq(const SmallMatrix<T,M,N,ColMajor>& Q,
                         const T* beta, T2* m)
    // SmallMatrix<T2,M,K,S2>& m
    {
        TMVAssert(M > N);
        for(ptrdiff_t j=0;j<N;++j) if (beta[j] != T(0)) {
            //HouseholderLMult(Q.col(j,j+1,M),beta(j),m.rowRange(j,M));
            for(ptrdiff_t k=0;k<K;++k) {
                T2 temp(0);
                for(ptrdiff_t i=j+1;i<M;++i)
                    //temp += TMV_CONJ(Q.cref(i,j)) * m.cref(i,k);
                    temp += TMV_CONJ(Q.cref(i,j)) * TMV_Val(m,M,K,S2,i,k);
                //temp += m.cref(j,k);
                temp += TMV_Val(m,M,K,S2,j,k);
                temp *= beta[j];
                //m.ref(j,k) -= temp;
                TMV_Val(m,M,K,S2,j,k) -= temp;
                for(ptrdiff_t i=j+1;i<M;++i)
                    //m.ref(i,k) -= Q.cref(i,j) * temp;
                    TMV_Val(m,M,K,S2,i,k) -= Q.cref(i,j) * temp;
            }
        }
    }

    template <ptrdiff_t K, StorageType S2, typename T, typename T2, ptrdiff_t M, ptrdiff_t N>
    inline void Q_RDivEq(const SmallMatrix<T,M,N,ColMajor>& Q,
                         const T* beta, T2* m)
    //SmallMatrix<T2,K,M,S2>& m
    {
        TMVAssert(M > N);
        for(ptrdiff_t j=N-1;j>=0;--j) if (beta[j] != T(0)) {
            //HouseholderLMult(Q.col(j,j+1,M).conjugate(),beta(j),
            //    m.colRange(j,M).transpose());
            for(ptrdiff_t k=0;k<K;++k) {
                T2 temp(0);
                for(ptrdiff_t i=j+1;i<M;++i)
                    //temp += Q.cref(i,j) * m.cref(k,i);
                    temp += Q.cref(i,j) * TMV_Val(m,K,M,S2,k,i);
                //temp += m.cref(k,j);
                temp += TMV_Val(m,K,M,S2,k,j);
                temp *= beta[j];
                //m.ref(k,j) -= temp;
                TMV_Val(m,K,M,S2,k,j) -= temp;
                for(ptrdiff_t i=j+1;i<M;++i)
                    //m.ref(k,i) -= TMV_CONJ(Q.cref(i,j)) * temp;
                    TMV_Val(m,K,M,S2,k,i) -= TMV_CONJ(Q.cref(i,j)) * temp;
            }
        }
    }


    //
    // Decompositions
    //

    // This basically just takes the NonBlock LU_Decompose function
    // And inlines everything so the compiler can optimize it for small N
    template <typename T, ptrdiff_t N>
    inline void DoLUD(SmallMatrix<T,N,N,ColMajor>& LU, ptrdiff_t* P)
    {
        // Input as m, output as packed LU
        for (ptrdiff_t j=0; j<N; ++j)
        {
            if (j > 1) {
                //LU.col(j,0,j) /= LU.SubMatrix(0,j,0,j).lowerTri(UnitDiag);
                for(ptrdiff_t k=0;k<j;++k) if (LU.cref(k,j) != T(0)) {
                    for(ptrdiff_t i=k+1;i<j;++i) {
                        LU.ref(i,j) -= LU.cref(k,j) * LU.cref(i,k);
                    }
                }
            }

            if (j > 0) {
                //LU.col(j,j,N) -= LU.SubMatrix(j,N,0,j) * LU.col(j,0,j);
                for(ptrdiff_t k=0;k<j;++k) {
                    for(ptrdiff_t i=j;i<N;++i) {
                        LU.ref(i,j) -= LU.cref(i,k) * LU.cref(k,j);
                    }
                }
            }

            //LU.col(j,j,N).MaxAbsElement(&ip);
            ptrdiff_t ip = j;
            if (j < N-1) {
                TMV_RealType(T) max =
                    isReal(T()) ? TMV_ABS(LU.cref(j,j)) :
                    TMV_NORM(LU.cref(j,j));
                for(ptrdiff_t k=j+1;k<N;++k) {
                    TMV_RealType(T) val =
                        isReal(T()) ? TMV_ABS(LU.cref(k,j)) :
                        TMV_NORM(LU.cref(k,j));
                    if (val > max) {
                        max = val;
                        ip = k;
                    }
                }
                if (ip != j) {
                    LU.swapRows(ip,j);
                }
            }
            P[j] = ip;

            //if (LU(j,j) != T(0)) LU.col(j,j+1,N) /= LU(j,j);
            if (LU.cref(j,j) != T(0)) {
                T invlujj = TMV_RealType(T)(1)/LU.cref(j,j);
                for(ptrdiff_t k=j+1;k<N;++k) LU.ref(k,j) *= invlujj;
            }
        }
    }

    // This isn't really worth inlining.  We just call the regular
    // Reflect function.
    template <typename T>
    T HouseholderReflect(T& x0, VectorView<T> x, T& det);

    // Again, this is just the normal QR_Decompose with everything inlined.
    template <typename T, ptrdiff_t M, ptrdiff_t N>
    inline void DoQRD(SmallMatrix<T,M,N,ColMajor>& QR, T* beta)
    {
        TMVAssert(M >= N);
        T det(0);
        for(ptrdiff_t j=0;j<N;++j) {
            beta[j] = HouseholderReflect(QR.ref(j,j),QR.col(j,j+1,M),det);
            if (beta[j] != T(0) && j<N-1) {
                //HouseholderLMult(v,beta,QR.SubMatrix(j,M,j+1,N));
                for(ptrdiff_t k=j+1;k<N;++k) {
                    T temp(0);
                    for(ptrdiff_t i=j+1;i<M;++i) temp +=
                        TMV_CONJ(QR.cref(i,j)) * QR.cref(i,k);
                    temp += QR.cref(j,k);
                    temp *= beta[j];
                    QR.ref(j,k) -= temp;
                    for(ptrdiff_t i=j+1;i<M;++i) QR.ref(i,k) -= QR.cref(i,j) * temp;
                }
            }
        }
    }

    //
    // Determinant
    //

    template <ptrdiff_t algo, typename T, ptrdiff_t M, ptrdiff_t N, StorageType S>
    struct SMDet;

    template <typename T, StorageType S>
    struct SMDet<1,T,1,1,S> // algo for 1x1 matrix
    { static inline T det(const T* m) { return *m; } };

    template <typename T, StorageType S>
    struct SMDet<2,T,2,2,S> // algo for 2x2 matrix
    { static inline T det(const T* m) { return (m[0]*m[3] - m[1]*m[2]); } };

    template <typename T, StorageType S>
    struct SMDet<3,T,3,3,S> // algo for 3x3 matrix
    {
        static T det(const T* m) {
            return (
                m[0] * (m[4]*m[8] - m[5]*m[7])
                -m[1] * (m[3]*m[8] - m[5]*m[6])
                +m[2] * (m[3]*m[7] - m[4]*m[6]));
        }
    };

    template <typename T, ptrdiff_t N, StorageType S>
    struct SMDet<10,T,N,N,S> // algo for normal calculation using LU
    {
        static T det(const T* m)
        {
            SmallMatrix<T,N,N,ColMajor> LU;
            // LU = m;
            SmallMatrixCopy<N,N,S,ColMajor>(m,LU.ptr());
            ptrdiff_t P[N];
            DoLUD(LU,P);
            T d(1);
            for(ptrdiff_t i=0;i<N;++i) {
                d *= LU.cref(i,i);
                if (P[i] != i) d = -d;
            }
            return d;
        }
    };

    template <typename T, ptrdiff_t N, StorageType S>
    struct SMDet<20,T,N,N,S> // algo for integers without fractions
    {
        template <typename TT>
        struct Helper
        {
            typedef long double longdouble_type;
            static TT convert(longdouble_type x)
            { return x < 0. ? -TT(floor(-x+0.5)) : TT(floor(x+0.5)); }
            static bool isUnity(longdouble_type x)
            { return TMV_ABS2(TMV_ABS2(x)-1.) < 0.1; }
            static bool isZero(longdouble_type x)
            { return TMV_ABS2(x) < 0.1; }
        };
        template <typename TT>
        struct Helper<std::complex<TT> >
        {
            typedef std::complex<TT> CT;
            typedef std::complex<long double> longdouble_type;
            static CT convert(longdouble_type x)
            {
                return CT(Helper<TT>::convert(x.real()),
                          Helper<TT>::convert(x.imag()));
            }
            static bool isUnity(longdouble_type x)
            { return TMV_ABS2(TMV_ABS2(x)-1.) < 0.1; }
            static bool isZero(longdouble_type x)
            { return TMV_ABS2(x) < 0.1; }
        };

        static T det(const T* m)
        {
            typedef typename Helper<T>::longdouble_type DT;
            SimpleMatrix<DT> A(N,N);
            // A = m;
            SmallMatrixCopy<N,N,S,ColMajor>(m,A.ptr());
            // This is the 1x1 Bareiss algorithm
            T det = 1;
            for (ptrdiff_t k=0; k<N-1; ++k) {
                ptrdiff_t imin = k, jmin = k;
                while (imin < N && Helper<T>::isZero(A.cref(imin,jmin))) ++imin;
                if (imin == N) return T(0);
                for (ptrdiff_t j=k; j<N; ++j) {
                    if (Helper<T>::isUnity(A.cref(imin,jmin))) break;
                    for (ptrdiff_t i=k; i<N; ++i) {
                        if (TMV_ABS2(A.cref(i,j))<TMV_ABS2(A.cref(imin,jmin)) &&
                            !Helper<T>::isZero(A.cref(i,j))) {
                            imin = i; jmin = j;
                            if (Helper<T>::isUnity(A.cref(i,j))) break;
                        }
                    }
                }

                if (Helper<T>::isZero(A.cref(imin,jmin))) return T(0);

                if (imin != k) { A.swapRows(imin,k); det *= -1; }
                if (jmin != k) { A.swapCols(jmin,k); det *= -1; }

                for (ptrdiff_t j=k+1; j<N; ++j) for(ptrdiff_t i=k+1; i<N; ++i)
                    A.ref(i,j) =
                        A.cref(k,k)*A.cref(i,j) - A.cref(i,k)*A.cref(k,j);
                if (k > 0) {
                    for (ptrdiff_t j=k+1; j<N; ++j) for(ptrdiff_t i=k+1; i<N; ++i)
                        A.ref(i,j) /= A.cref(k-1,k-1);
                }
            }
            det *= Helper<T>::convert(A.cref(N-1,N-1));
            return det;
        }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A>
    inline T DoDet(const SmallMatrix<T,M,N,A>& m)
    {
        const ptrdiff_t algo =
            N == 1 ? 1 :
            N == 2 ? 2 :
            N == 3 ? 3 :
            Traits<T>::isinteger ? 20 :
            10;
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        return SMDet<algo,T,M,N,S>::det(m.cptr());
    }

    //
    // Matrix Inverse
    //

    template <StorageType S2, typename T, typename T2, ptrdiff_t M, ptrdiff_t N>
    inline void SMQRInv(SmallMatrix<T,M,N,ColMajor>& QR, T2* minv)
    //SmallMatrix<T2,N,M,S2> minv
    {
        TMVAssert(M > N);
        T beta[N];
        DoQRD(QR,beta);
        InvertU(QR);

        //minv.setZero();
        for(ptrdiff_t i=0;i<M*N;++i) minv[i] = T2(0);
        //minv.colRange(0,N)).upperTri() = QR.upperTri()
        for(ptrdiff_t j=0;j<N;++j) for(ptrdiff_t i=0;i<=j;i++)
            //minv.ref(i,j) = QR.cref(i,j);
            TMV_Val(minv,N,M,S2,i,j) = QR.cref(i,j);
        Q_RDivEq<N,S2>(QR,beta,minv);
    }

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, bool MltN, StorageType S, StorageType S2>
    struct SMInv
    {
        static void inv(const T* m, T2* minv)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T2,N,M,S2> minv
        {
            TMVAssert(M > N);
            // (M < N and M == N specialized below)
            SmallMatrix<T,M,N,ColMajor> QR;
            // QR = m;
            SmallMatrixCopy<M,N,S,ColMajor>(m,QR.ptr());
            SMQRInv<S2>(QR,minv);
        }
    };

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, StorageType S, StorageType S2>
    struct SMInv<T,T2,M,N,true,S,S2>
    {
        static void inv(const T* m, T2* minv)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T2,N,M,S2> minv
        {
            TMVAssert(M < N);
            SmallMatrix<T,N,M,ColMajor> QR;
            // QR = m.transpose();
            SmallMatrixCopy<M,N,S,RowMajor>(m,QR.ptr());
            SMQRInv<TMV_TransOf(S2)>(QR,minv);
        }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, bool MltN, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,M,N,MltN,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,M,N,true,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T2, ptrdiff_t N, StorageType S, StorageType S2>
    struct SMInv<T,T2,N,N,false,S,S2>
    {
        static void inv(const T* m, T2* minv)
            //SmallMatrix<T,N,N,S> m
            //SmallMatrix<T2,N,N,S2> minv
        {
            SmallMatrix<T,N,N,ColMajor> LU;
            // LU = m;
            SmallMatrixCopy<N,N,S,ColMajor>(m,LU.ptr());
            ptrdiff_t P[N];
            DoLUD(LU,P);

            InvertU(LU);
            InvertL(LU);

            //minv = LU.upperTri() * LU.lowerTri(UnitDiag);
            //minv.setZero();
            for(ptrdiff_t i=0;i<N*N;++i) minv[i] = T2(0);
            for(ptrdiff_t j=0;j<N;++j) {
                for(ptrdiff_t i=0;i<=j;++i)
                    //minv.ref(i,j) = LU.cref(i,j);
                    TMV_Val(minv,N,N,S2,i,j) = LU.cref(i,j);
                for(ptrdiff_t k=j+1;k<N;++k) {
                    for(ptrdiff_t i=0;i<=k;++i) {
                        //minv.ref(i,j) += LU.cref(i,k) * LU.cref(k,j);
                        TMV_Val(minv,N,N,S2,i,j) += LU.cref(i,k) * LU.cref(k,j);
                    }
                }
            }

            //minv.reversePermuteCols(P);
            for(ptrdiff_t j=N-1;j>=0;--j) if (P[j] != j)
                for(ptrdiff_t i=0;i<N;++i)
                    TMV_SWAP(TMV_Val(minv,N,N,S2,i,j),
                             TMV_Val(minv,N,N,S2,i,P[j]));
        }
    };

    template <typename T, ptrdiff_t N, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,N,N,false,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };


    template <typename T, typename T2, StorageType S, StorageType S2>
    struct SMInv<T,T2,1,1,false,S,S2>
    {
        static void inv(const T* m, T2* minv)
        {
            if (*m == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            *minv = TMV_InverseOf(*m);
        }
    };

    template <typename T, typename T2, StorageType S, StorageType S2>
    struct SMInv<T,T2,2,2,false,S,S2>
    {
        static void inv(const T* m, T2* minv)
        {
            T det = SMDet<2,T,2,2,S>::det(m);
            if (det == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            //minv.ref(0,0) = m.cref(1,1)/d;
            //minv.ref(0,1) = -m.cref(0,1)/d;
            //minv.ref(1,0) = -m.cref(1,0)/d;
            //minv.ref(1,1) = m.cref(0,0)/d;
            minv[0] = m[3]/det;
            minv[3] = m[0]/det;
            if (S == S2) {
                minv[1] = -m[1]/det;
                minv[2] = -m[2]/det;
            } else {
                minv[1] = -m[2]/det;
                minv[2] = -m[1]/det;
            }
        }
    };

    template <typename T, typename T2, StorageType S, StorageType S2>
    struct SMInv<T,T2,3,3,false,S,S2>
    {
        static void inv(const T* m, T2* minv)
        {
            T det = SMDet<3,T,3,3,S>::det(m);
            if (det == T(0))
#ifdef NOTHROW
            { std::cerr<<"Singular SmallMatrix found\n"; exit(1); }
#else
            { throw Singular(); }
#endif
            minv[0] = (m[4]*m[8]-m[7]*m[5])/det;
            minv[4] = (m[0]*m[8]-m[6]*m[2])/det;
            minv[8] = (m[0]*m[4]-m[3]*m[1])/det;
            if (S == S2) {
                minv[1] = -(m[1]*m[8]-m[7]*m[2])/det;
                minv[2] = (m[1]*m[5]-m[4]*m[2])/det;
                minv[3] = -(m[3]*m[8]-m[6]*m[5])/det;
                minv[5] = -(m[0]*m[5]-m[3]*m[2])/det;
                minv[6] = (m[3]*m[7]-m[6]*m[4])/det;
                minv[7] = -(m[0]*m[7]-m[6]*m[1])/det;
            } else {
                minv[1] = -(m[3]*m[8]-m[6]*m[5])/det;
                minv[2] = (m[3]*m[7]-m[6]*m[4])/det;
                minv[3] = -(m[1]*m[8]-m[7]*m[2])/det;
                minv[5] = -(m[0]*m[7]-m[6]*m[1])/det;
                minv[6] = (m[1]*m[5]-m[4]*m[2])/det;
                minv[7] = -(m[0]*m[5]-m[3]*m[2])/det;
            }
        }
    };

    template <typename T, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,1,1,false,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,2,2,false,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, StorageType S, StorageType S2>
    struct SMInv<std::complex<T>,T,3,3,false,S,S2>
    {
        static void inv(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A2>
    inline void DoInverse(
        const SmallMatrix<T,M,N,A>& m, SmallMatrix<T2,N,M,A2>& minv)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMInv<T,T2,M,N,M<N,S,S2>::inv(m.cptr(),minv.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,M,N,A>(m);
        }
#endif
    }

    //
    // Essential Matrix Div Algorithms
    //

    template <ptrdiff_t K, StorageType S1, StorageType S2, typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N>
    inline void SMQRLDiv(SmallMatrix<T,M,N,ColMajor>& QR, const T1* m1, T2* m2)
    //SmallMatrix<T1,M,K,S1> m1
    //SmallMatrix<T2,N,K,S2> m2
    {
        TMVAssert(M > N);
        T beta[N];
        DoQRD(QR,beta);

        // m2 = R^-1 Qt m1
        //SmallMatrix<T,M,K,RowMajor> temp = m1;
        T2 temp[M*K];
        SmallMatrixCopy<M,K,S1,RowMajor>(m1,temp);
        Q_LDivEq<K,RowMajor>(QR,beta,temp);
        //m2 = temp.rowRange(0,N);
        SmallMatrixCopy<N,K,RowMajor,S2>(temp,m2);
        //m2.view() /= QR.upperTri();
        LDivEq_U<N,K,S2>(QR,m2);
    }

    template <ptrdiff_t K, StorageType S1, StorageType S2, typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N>
    inline void SMQRRDiv(SmallMatrix<T,M,N,ColMajor>& QR, const T1* m1, T2* m2)
    //SmallMatrix<T1,K,N,S1> m1
    //SmallMatrix<T2,K,M,S2> m2
    {
        TMVAssert(M > N);
        T beta[N];
        DoQRD(QR,beta);

        // m2 = m1 R^-1 Qt
        //m2.setZero();
        for(ptrdiff_t i=0;i<K*M;++i) m2[i] = T2(0);
        //m2.colRange(0,N) = m1;
        if (S2 == ColMajor)
            SmallMatrixCopy<K,N,S1,S2>(m1,m2);
        else
            for(ptrdiff_t i=0;i<K;++i) for(ptrdiff_t j=0;j<N;j++)
                //m2(i,j) = m1(i,j);
                TMV_Val(m2,K,M,S2,i,j) = TMV_Val(m1,K,N,S1,i,j);
        //m2.colRange(0,N) %= QR.upperTri();
        RDivEq_U<K,M,S2>(QR,m2);
        Q_RDivEq<K,S2>(QR,beta,m2);
    }

    template <ptrdiff_t K, StorageType S2, typename T, typename T2, ptrdiff_t N>
    inline void SMLULDivEq(SmallMatrix<T,N,N,ColMajor>& LU, T2* m2)
    //SmallMatrix<T2,N,K,S2> m2
    {
        //std::cout<<"Start SMLULDivEq:\n";
        //std::cout<<"LU = "<<TMV_Text(LU)<<"  "<<LU<<std::endl;
        //std::cout<<"K = "<<K<<std::endl;
        //std::cout<<"S2 = "<<TMV_Text(S2)<<std::endl;
        //std::cout<<"m2 = "<<m2<<std::endl;
        //PrintSM<N,K,S2>(std::cout,m2);
        ptrdiff_t P[N];
        //std::cout<<"Before DoLUD"<<std::endl;
        DoLUD(LU,P);
        //std::cout<<"After DoLUD"<<std::endl;

        //m2.permuteRows(P);
        for(ptrdiff_t i=0;i<N;++i) if (P[i] != i)
            for(ptrdiff_t j=0;j<K;++j)
                TMV_SWAP(TMV_Val(m2,N,K,S2,i,j),TMV_Val(m2,N,K,S2,P[i],j));
        //std::cout<<"After permuteRows"<<std::endl;

        // m2 = L^-1 m2
        //m2.view() /= LU.lowerTri(UnitDiag);
        LDivEq_L<N,K,S2>(LU,m2);
        //std::cout<<"After LDivEq_L"<<std::endl;

        // m2 = U^-1 m2;
        //m2.view() /= LU.upperTri();
        LDivEq_U<N,K,S2>(LU,m2);
        //std::cout<<"After LDivEq_U"<<std::endl;
    }


    //
    // Matrix LDiv
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM
    {
        static void ldiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T1,M,K,S1> m1
            //SmallMatrix<T2,N,K,S2> m2
        {
            TMVAssert(M > N);
            // (M < N and M == N specialized below)
            SmallMatrix<T,M,N,ColMajor> QR;
            // QR = m;
            SmallMatrixCopy<M,N,S,ColMajor>(m,QR.ptr());
            SMQRLDiv<K,S1,S2>(QR,m1,m2);
        }
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<T,T1,T2,M,N,true,K,S,S1,S2>
    {
        static void ldiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T1,M,K,S1> m1
            //SmallMatrix<T2,N,K,S2> m2
        {
            TMVAssert(M < N);
            // m2 = m^-1 m1
            // m2T = m1T mT^-1
            SmallMatrix<T,N,M,ColMajor> QR;
            // QR = m.transpose();
            SmallMatrixCopy<M,N,S,RowMajor>(m,QR.ptr());
            SMQRRDiv<K,TMV_TransOf(S1),TMV_TransOf(S2)>(QR,m1,m2);
        }
    };

    template <typename T, typename T1, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<std::complex<T>,T1,T,M,N,MltN,K,S,S1,S2>
    {
        static void ldiv(const std::complex<T>*, const T1*, T* )
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<T,std::complex<T>,T,M,N,MltN,K,S,S1,S2>
    {
        static void ldiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T1, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<std::complex<T>,T1,T,M,N,true,K,S,S1,S2>
    {
        static void ldiv(const std::complex<T>*, const T1*, T* )
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<T,std::complex<T>,T,M,N,true,K,S,S1,S2>
    {
        static void ldiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T2, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S2>
    struct SMLDivEqM
    {
        static void ldiveq(const T* m, T2* m2)
            //SmallMatrix<T,N,N,S> m
            //SmallMatrix<T2,N,K,S2> m2
        {
            //std::cout<<"Start SMLDivEqM::ldiveq:\n";
            //std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            //std::cout<<"S2 = "<<TMV_Text(S2)<<std::endl;
            //std::cout<<"N,K = "<<N<<','<<K<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            //PrintSM<N,N,S>(std::cout,m);
            //std::cout<<"m2 = "<<m2<<std::endl;
            //PrintSM<N,K,S2>(std::cout,m2);
            SmallMatrix<T,N,N,ColMajor> LU;
            // LU = m;
            SmallMatrixCopy<N,N,S,ColMajor>(m,LU.ptr());
            //std::cout<<"LU = "<<LU<<std::endl;
            SMLULDivEq<K,S2>(LU,m2);
            //std::cout<<"After LULDivEq\n";
        }
    };

    template <typename T, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S2>
    struct SMLDivEqM<std::complex<T>,T,N,K,S,S2>
    {
        static void ldiveq(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<T,T1,T2,N,N,false,K,S,S1,S2>
    {
        static void ldiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,N,N,S> m
            //SmallMatrix<T1,N,K,S1> m1
            //SmallMatrix<T2,N,K,S2> m2
        {
            //std::cout<<"Start SMLDivM::ldiv:\n";
            //std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            //std::cout<<"S1 = "<<TMV_Text(S1)<<std::endl;
            //std::cout<<"S2 = "<<TMV_Text(S2)<<std::endl;
            //std::cout<<"N,K = "<<N<<','<<K<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
            //PrintSM<N,N,S>(std::cout,m);
            //std::cout<<"m1 = "<<m1<<std::endl;
            //PrintSM<N,K,S1>(std::cout,m1);
            //std::cout<<"m2 = "<<m2<<std::endl;
            //PrintSM<N,K,S2>(std::cout,m2);
            //m2 = m1;
            SmallMatrixCopy<N,K,S1,S2>(m1,m2);
            //std::cout<<"After copy\n";
            SMLDivEqM<T,T2,N,K,S,S2>::ldiveq(m,m2);
            //std::cout<<"After ldiveq\n";
        }
    };

    template <typename T, typename T1, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<std::complex<T>,T1,T,N,N,false,K,S,S1,S2>
    {
        static void ldiv(const std::complex<T>*, const T1*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMLDivM<T,std::complex<T>,T,N,N,false,K,S,S1,S2>
    {
        static void ldiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    // m2 = m^-1 m2
    template <typename T, typename T2, ptrdiff_t N, ptrdiff_t K, int A, int A2>
    inline void DoLDivEq(
        const SmallMatrix<T,N,N,A>& m, SmallMatrix<T2,N,K,A2>& m2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMLDivEqM<T,T2,N,K,S,S2>::ldiveq(m.cptr(),m2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,N,N,A>(m);
        }
#endif
    }

    // m2 = m^-1 m1
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A, int A1, int A2>
    inline void DoLDiv(
        const SmallMatrix<T,M,N,A>& m,
        const SmallMatrix<T1,M,K,A1>& m1, SmallMatrix<T2,N,K,A2>& m2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S1 = static_cast<StorageType>(A1 & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMLDivM<T,T1,T2,M,N,M<N,K,S,S1,S2>::ldiv(m.cptr(),m1.cptr(),m2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,M,N,A>(m);
        }
#endif
    }


    //
    // Vector LDiv
    //

    // v2 = m^-1 v2
    template <typename T, typename T2, ptrdiff_t N, int A, int A2>
    inline void DoLDivEq(
        const SmallMatrix<T,N,N,A>& m, SmallVector<T2,N,A2>& v2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMLDivEqM<T,T2,N,1,S,ColMajor>::ldiveq(m.cptr(),v2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,N,N,A>(m);
        }
#endif
    }

    // v2 = m^-1 v1
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline void DoLDiv(
        const SmallMatrix<T,M,N,A>& m,
        const SmallVector<T1,M,A1>& v1, SmallVector<T2,N,A2>& v2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        //std::cout<<"Start DoLDiv:\n";
        //std::cout<<"m = "<<TMV_Text(m)<<"  "<<m<<std::endl;
        //std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
        //std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        //std::cout<<"S = "<<TMV_Text(S)<<std::endl;
        //std::cout<<"M,N = "<<M<<','<<N<<std::endl;
#ifndef NOTHROW
        try {
#endif
            SMLDivM<T,T1,T2,M,N,M<N,1,S,ColMajor,ColMajor>::ldiv(
                m.cptr(),v1.cptr(),v2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,M,N,A>(m);
        }
#endif
    }


    //
    // Matrix RDiv
    //

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM
    {
        static void rdiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T1,K,N,S1> m1
            //SmallMatrix<T2,K,M,S2> m2
        {
            TMVAssert(M > N);
            // (M < N and M == N specialized below)
            SmallMatrix<T,M,N,ColMajor> QR;
            // QR = m
            SmallMatrixCopy<M,N,S,ColMajor>(m,QR.ptr());
            SMQRRDiv<K,S1,S2>(QR,m1,m2);
        }
    };

    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<T,T1,T2,M,N,true,K,S,S1,S2>
    {
        static void rdiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T1,K,N,S1> m1
            //SmallMatrix<T2,K,M,S2> m2
        {
            TMVAssert(M < N);
            // m2 = m1 m^-1
            // m2t = mt^-1 m1t
            SmallMatrix<T,N,M,ColMajor> QR;
            // QR = m.transpose();
            SmallMatrixCopy<M,N,S,RowMajor>(m,QR.ptr());
            SMQRLDiv<K,TMV_TransOf(S1),TMV_TransOf(S2)>(QR,m1,m2);
        }
    };

    template <typename T, typename T1, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<std::complex<T>,T1,T,M,N,MltN,K,S,S1,S2>
    {
        static void rdiv(const std::complex<T>*, const T1*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, bool MltN, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<T,std::complex<T>,T,M,N,MltN,K,S,S1,S2>
    {
        static void rdiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T1, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<std::complex<T>,T1,T,M,N,true,K,S,S1,S2>
    {
        static void rdiv(const std::complex<T>*, const T1*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<T,std::complex<T>,T,M,N,true,K,S,S1,S2>
    {
        static void rdiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T2, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S2>
    struct SMRDivEqM
    {
        static void rdiveq(const T* m, T2* m2)
            //SmallMatrix<T,N,N,S> m
            //SmallMatrix<T2,K,N,S2> m2
        {
            // m2 = m2 m^-1
            // m2t = mt^-1 m2t
            SmallMatrix<T,N,N,ColMajor> LU;
            // LU = m.transpose()
            SmallMatrixCopy<N,N,S,RowMajor>(m,LU.ptr());
            SMLULDivEq<K,TMV_TransOf(S2)>(LU,m2);
        }
    };

    template <typename T, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S2>
    struct SMRDivEqM<std::complex<T>,T,N,K,S,S2>
    {
        static void rdiveq(const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, typename T1, typename T2, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<T,T1,T2,N,N,false,K,S,S1,S2>
    {
        static void rdiv(const T* m, const T1* m1, T2* m2)
            //SmallMatrix<T,N,N,S> m
            //SmallMatrix<T1,K,N,S1> m1
            //SmallMatrix<T2,K,N,S2> m2
        {
            //m2 = m1;
            SmallMatrixCopy<K,N,S1,S2>(m1,m2);
            //DoRDivEq(m,m2);
            SMRDivEqM<T,T2,N,K,S,S2>::rdiveq(m,m2);
        }
    };

    template <typename T, typename T1, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<std::complex<T>,T1,T,N,N,false,K,S,S1,S2>
    {
        static void rdiv(const std::complex<T>*, const T1*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    template <typename T, ptrdiff_t N, ptrdiff_t K, StorageType S, StorageType S1, StorageType S2>
    struct SMRDivM<T,std::complex<T>,T,N,N,false,K,S,S1,S2>
    {
        static void rdiv(const T*, const std::complex<T>*, T*)
        { TMVAssert(TMV_FALSE); }
    };

    // m2 = m2 m^-1
    template <typename T, typename T2, ptrdiff_t N, ptrdiff_t K, int A, int A2>
    inline void DoRDivEq(
        const SmallMatrix<T,N,N,A>& m, SmallMatrix<T2,K,N,A2>& m2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMRDivEqM<T,T2,N,K,S,S2>::rdiveq(m.cptr(),m2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,N,N,A>(m);
        }
#endif
    }

    // m2 = m^-1 m1
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, int A, int A1, int A2>
    inline void DoRDiv(
        const SmallMatrix<T,M,N,A>& m,
        const SmallMatrix<T1,K,N,A1>& m1, SmallMatrix<T2,K,M,A2>& m2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S1 = static_cast<StorageType>(A1 & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMRDivM<T,T1,T2,M,N,M<N,K,S,S1,S2>::rdiv(
                m.cptr(),m1.cptr(),m2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,M,N,A>(m);
        }
#endif
    }

    //
    // Vector RDiv
    //

    // v2 = v2 m^-1
    template <typename T, typename T2, ptrdiff_t N, int A, int A2>
    inline void DoRDivEq(
        const SmallMatrix<T,N,N,A>& m, SmallVector<T2,N,A2>& v2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMRDivEqM<T,T2,N,1,S,RowMajor>::rdiveq(m.cptr(),v2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,N,N,A>(m);
        }
#endif
    }

    // v2 = m^-1 v1
    template <typename T, typename T1, typename T2, ptrdiff_t M, ptrdiff_t N, int A, int A1, int A2>
    inline void DoRDiv(
        const SmallMatrix<T,M,N,A>& m,
        const SmallVector<T1,N,A1>& v1, SmallVector<T2,M,A2>& v2)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
#ifndef NOTHROW
        try {
#endif
            SMRDivM<T,T1,T2,M,N,M<N,1,S,RowMajor,RowMajor>::rdiv(
                m.cptr(),v1.cptr(),v2.ptr());
#ifndef NOTHROW
        } catch (tmv::Singular) {
            throw SingularSmallMatrix<T,M,N,A>(m);
        }
#endif
    }

    //
    // InverseATA
    //

    template <StorageType S2, typename T, ptrdiff_t M, ptrdiff_t N>
    inline void SMQRATA(SmallMatrix<T,M,N,ColMajor> QR, T* ata)
    {
        TMVAssert(M >= N);
        T beta[N];
        DoQRD(QR,beta);

        InvertU(QR);

        //ata = QR.upperTri() * QR.upperTri().adjoptrdiff_t();
        //ata.setZero();
        for(ptrdiff_t i=0;i<N*N;++i) ata[i] = T(0);
        for(ptrdiff_t j=0;j<N;++j) {
            for(ptrdiff_t k=j;k<N;++k) {
                for(ptrdiff_t i=0;i<=k;++i) {
                    //ata.ref(i,j) += QR.cref(i,k) * TMV_CONJ(QR.cref(j,k));
                    TMV_Val(ata,N,N,S2,i,j) +=
                        QR.cref(i,k) * TMV_CONJ(QR.cref(j,k));
                }
            }
        }
    }

    template <typename T, ptrdiff_t M, ptrdiff_t N, bool MltN, StorageType S, StorageType S2>
    struct SMATA
    {
        static void ata(const T* m, T* ata)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T,N,N,S2>& ata
        {
            TMVAssert(M >= N);
            SmallMatrix<T,M,N,ColMajor> QR;
            // QR = m;
            SmallMatrixCopy<M,N,S,ColMajor>(m,QR.ptr());
            SMQRATA<S2>(QR,ata);
        }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, StorageType S, StorageType S2>
    struct SMATA<T,M,N,true,S,S2>
    {
        static void ata(const T* m, T* ata)
            //SmallMatrix<T,M,N,S> m
            //SmallMatrix<T,N,N,S2> ata
        {
            TMVAssert(M < N);
            SmallMatrix<T,N,M,ColMajor> QR;
            // QR = m.transpose();
            SmallMatrixCopy<M,N,S,RowMajor>(m,QR.ptr());
            SMQRATA<S2>(QR,ata);
        }
    };

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A, int A2>
    inline void DoInverseATA(
        const SmallMatrix<T,M,N,A>& m, SmallMatrix<T,N,N,A2>& ata)
    {
        const StorageType S = static_cast<StorageType>(A & AllStorageType);
        const StorageType S2 = static_cast<StorageType>(A2 & AllStorageType);
        SMATA<T,M,N,M<N,S,S2>::ata(m.cptr(),ata.ptr());
    }

#undef TMV_Index
#undef TMV_Val
#undef TMV_TransOf
}

#endif
