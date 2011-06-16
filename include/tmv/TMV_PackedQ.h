///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
//
// This file defines the PackedQ class.
//
// A PackedQ object represents a unitary matrix that is the 
// result of a QR (or similar) decomposition.  
// The QR decomposition is able to store almost all the information about
// the Q matrix in the lower trapezoid part of the matrix, with R being 
// stored in the upper triangle part.  The only extra information needed
// is a vector of beta values.
//
// Basically each column of the matrix and the corresponding beta value
// define a Householder reflection.  The full Q matrix is the product 
// of these reflections.  The main use of a PackedQ object is to multiply
// Q by some other vector or matrix.  This is implemented by multiplying
// each Householder matrix in turn.  This is more efficient than directly
// constructing the full Q matrix and then multiplying by that.
//
// The constructor takes as arguments a reference to the matrix
// where the columns are stored and a reference to the beta vector.
// The PackedQ object stores views of these.  So if they go out of scope,
// then the PackedQ object will be invalid.  
//
// The PackedQ class is a template with two template parameters.
// These parameters are the types of the Q and beta objects which
// described above.  These templates are called M and V below.
//
// Constructors:
//
//    PackedQ(const M& Q, const V& beta)
//        Create a PackedQ matrix using the values in Q and beta.
//
// Access Functions:
//
//     size_t colsize() const
//     size_t rowsize() const
//         Return the size of the matrix.
//     
//     const M& getQ() const
//     const V& getBeta() const
//         Return the component elements
//
//     size_t operator()(size_t i, size_t j) const
//         This is required for BaseMatrix, but it's extremely inefficient
//         for PackedQ.  (You basically need to construct the full Q
//         matrix up to the j column.)
//         So this function is private to make it a compiler error if anyone
//         tries to use it.
//
// Functions:
//
//     int det() const
//         Returns the determinant of Q.  This is always either 1 or -1.
//     int logDet(int* sign=0) const
//         Returns 0.  If requested, the sign is either 1 or -1.
//
//
// Operators:
//
//     q*v
//     v*q
//     v/q
//     v%q
//     v *= q
//     v /= q
//     v %= q
//
//     q*m
//     m*q
//     m/q
//     m%q
//     m *= q
//     m /= q
//     m %= q
//     
//     These are the main reason for having a PackedQ class.
//     The multiplication and division operations can be performed 
//     without having to construct the actual Q matrix.
//     This is done by successively applying the Householder reflections
//     to the vector or matrix.
//

#ifndef TMV_PackedQ_H
#define TMV_PackedQ_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseVector.h"
#include "TMV_UnpackQ.h"

namespace tmv {

    // TODO: Move this to another file and add block, inst, inline, etc.
    template <class M1, class V1, class M2>
    static void PackedQ_MultEq(
        const M1& Q, const V1& beta, BaseMatrix_Rec_Mutable<M2>& m)
    {
        typedef typename V1::value_type T1;

        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        const int M = Q.colsize();
        const int N = Q.rowsize();
        typedef typename M1::const_col_sub_type M1c;
        typedef typename M2::rowrange_type M2r;
        for(int j=N-1;j>=0;--j) if (beta(j) != T1(0)) {
            Householder<M1c> H(Q.col(j,j+1,M),beta(j));
            M2r mr = m.rowRange(j,M);
            H.multEq(mr);
        }
    }

    template <class M1, class V1, class M2>
    static void PackedQ_LDivEq(
        const M1& Q, const V1& beta, BaseMatrix_Rec_Mutable<M2>& m)
    {
        typedef typename V1::value_type T1;

        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());        
        const int M = Q.colsize();
        const int N = Q.rowsize();
        typedef typename M1::const_col_sub_type M1c;
        typedef typename M2::rowrange_type M2r;
        for(int j=0;j<N;++j) if (beta(j) != T1(0)) {
            Householder<M1c> H(Q.col(j,j+1,M),beta(j));
            M2r mr = m.rowRange(j,M);
            H.multEq(mr);
        }
    }

    template <class M1, class V1, class V2>
    static void PackedQ_MultEq(
        const M1& Q, const V1& beta, BaseVector_Mutable<V2>& v)
    {
        typedef typename V1::value_type T1;

        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        const int M = Q.colsize();
        const int N = Q.rowsize();
        typedef typename M1::const_col_sub_type M1c;
        typedef typename V2::subvector_type V2s;
        for(int j=N-1;j>=0;--j) if (beta(j) != T1(0)) {
            Householder<M1c> H(Q.col(j,j+1,M),beta(j));
            V2s vs = v.subVector(j,M);
            H.multEq(vs);
        }
    }

    template <class M1, class V1, class V2>
    static void PackedQ_LDivEq(
        const M1& Q, const V1& beta, BaseVector_Mutable<V2>& v)
    {
        typedef typename V1::value_type T1;

        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());        
        const int M = Q.colsize();
        const int N = Q.rowsize();
        typedef typename M1::const_col_sub_type M1c;
        typedef typename V2::subvector_type V2s;
        for(int j=0;j<N;++j) if (beta(j) != T1(0)) {
            Householder<M1c> H(Q.col(j,j+1,M),beta(j));
            V2s vs = v.subVector(j,M);
            H.multEq(vs);
        }               
    }

    template <class M, class V> 
    class PackedQ;

    template <class M, class V> 
    struct Traits<PackedQ<M,V> >
    {
        typedef typename M::value_type value_type;
        typedef typename M::real_type real_type;
        typedef typename M::complex_type complex_type;

        enum { isreal = Traits<M>::isreal };
        enum { iscomplex = Traits<M>::iscomplex };

        typedef PackedQ<M,V> type;
        typedef Matrix<value_type> copy_type;
        typedef copy_type calc_type;
        typedef copy_type eval_type;
        typedef InvalidType inverse_type;

        enum { _colsize = M::_colsize };
        enum { _rowsize = M::_rowsize };
        enum { _shape = Rec };
        enum { _fort = M::_fort };
        enum { _calc = false };
    };

    template <class M, class V> 
    class PackedQ : 
        public BaseMatrix<PackedQ<M,V> >
    {
    public :

        typedef typename M::value_type T;
        typedef typename M::real_type RT;

        PackedQ(const BaseMatrix_Rec<M>& _Q, const BaseVector<V>& _beta) :
            Q(_Q.mat()), beta(_beta.vec()) 
        {
            TMVStaticAssert((Traits2<typename V::value_type,RT>::sametype));
            TMVAssert(beta.size() == Q.rowsize()); 
        }

        PackedQ(const PackedQ<M,V>& rhs) : Q(rhs.Q), beta(rhs.beta) {}

        ~PackedQ() {}

        //
        // Accesss
        //

        TMV_INLINE const M& getQ() const { return Q; }
        TMV_INLINE const V& getBeta() const { return beta; }

        int det() 
        { 
            int d=1;
            const int n = beta.size();
            for(int i=0; i<n; ++i) if (beta!=RT(0)) d = -d;
            return d;
        }

        int logDet(int* sign=0) const
        { if (sign) *sign = det(); return 0; }


        //
        // Create matrix version
        //

        template <class M2>
        void assignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            Q.assignTo(m2);
            UnpackQ(m2,beta);
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            Q.newAssignTo(m2);
            UnpackQ(m2,beta);
        }


        // 
        // Auxilliary functions
        //

        TMV_INLINE size_t colsize() const { return Q.colsize(); }
        TMV_INLINE size_t rowsize() const { return Q.rowsize(); }
        
    private : 

        const M& Q;
        const V& beta;

        // Don't allow op=
        void operator=(const PackedQ<M,V>& rhs);
    };

    template <class M, class V>
    static TMV_INLINE int Det(const PackedQ<M,V>& Q)
    { return Q.det(); }
    template <class M, class V>
    static TMV_INLINE int LogDet(const PackedQ<M,V>& Q)
    { return Q.logDet(); }


    //
    // v3 = Q * v2
    //

    template <bool add, int ix, class T, class M1, class V1, class V2, class V3>
    static inline void MultMV(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        if (add) {
            v3 += (x*m1*v2).calc();
        } else if (m1.isSquare()) {
            v3 = v2;
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),v3.vec());
            Scale(x,v3);
        } else {
            TMVAssert(v3.size() > v2.size());
            v3.subVector(0,v2.size()) = v2;
            v3.subVector(v2.size(),v3.size()).setZero();
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),v3.vec());
            Scale(x,v3);
        }
    }
    template <bool add, int ix, class T, class M1, class V1, class V2, class V3>
    static TMV_INLINE void NoAliasMultMV(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m1,v2,v3); }
    template <bool add, int ix, class T, class M1, class V1, class V2, class V3>
    static TMV_INLINE void AliasMultMV(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m1,v2,v3); }

    
    //
    // v3 = v1 * Q
    // v3 = QT * v1
    // v3* = Qt * v1*
    //

    template <bool add, int ix, class T, class V1, class M2, class V2, class V3>
    static inline void MultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typename V3::conjugate_type v3c = v3.conjugate();
        if (add) {
            v3 += (x*v1*m2).calc();
        } else if (m2.isSquare()) {
            v3 = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        } else {
            TMVAssert(v1.size() > v3.size());
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            typename VCopyHelper<T12,V1::_size,false>::type v1c =
                v1.conjugate();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
            v3 = x * v1c.subVector(0,v3.size());
        }
    }
    template <bool add, int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void NoAliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { MultVM<add>(x,v1,m2,v3); }
    template <bool add, int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void AliasMultVM(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { MultVM<add>(x,v1,m2,v3); }
 

    //
    // v *= Q
    //

    template <int ix, class T, class V1, class M2, class V2>
    static inline void MultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    {
        typename V1::conjugate_type v1c = v1.conjugate();
        PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
        m2.divEq(v1c);
        Scale(x,v1);
    }
    template <int ix, class T, class V1, class M2, class V2>
    static TMV_INLINE void NoAliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    { MultEqVM(v1,x,m2); }
    template <int ix, class T, class V1, class M2, class V2>
    static TMV_INLINE void AliasMultEqVM(
        BaseVector_Mutable<V1>& v1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    { MultEqVM(v1,x,m2); }


    //
    // v / Q
    //

    template <int ix, class T, class V1, class M2, class V2, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        if (m2.isSquare()) {
            v3 = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v3);
            Scale(x,v3);
        } else {
            TMVAssert(v1.size() > v3.size());
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            typename VCopyHelper<T12,V1::_size,false>::type v1c = v1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1c);
            v3 = x * v1c.subVector(0,v3.size());
        }
    }
    template <int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }


    //
    // v /= Q
    //

    template <class V1, class M2, class V2>
    static inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { PackedQ_LDivEq(m2.getQ(),m2.getBeta(),v1); }
    template <class V1, class M2, class V2>
    static TMV_INLINE void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { LDivEq(v1,m2); }
    template <class V1, class M2, class V2>
    static TMV_INLINE void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { LDivEq(v1,m2); }


    // 
    // v % Q
    // v3 = v1 * Qt
    // v3 = Q* v1
    // v3* = Q v1*
    //

    template <int ix, class T, class V1, class M2, class V2, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typename V3::conjugate_type v3c = v3.conjugate();
        if (m2.isSquare()) {
            v3 = v1;
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        } else {
            TMVAssert(v3.size() > v1.size());
            v3.subVector(0,v1.size()) = v1;
            v3.subVector(v1.size(),v3.size()).setZero();
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),v3c);
            Scale(x,v3);
        }
    }
    template <int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V2, class V3>
    static TMV_INLINE void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const PackedQ<M2,V2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }


    // 
    // v %= Q
    //

    template <class V1, class M2, class V2>
    static inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    {
        typename V1::conjugate_type v1c = v1.conjugate();
        PackedQ_MultEq(m2.getQ(),m2.getBeta(),v1c);
    }
    template <class V1, class M2, class V2>
    static TMV_INLINE void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { RDivEq(v1,m2); }
    template <class V1, class M2, class V2>
    static TMV_INLINE void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const PackedQ<M2,V2>& m2)
    { RDivEq(v1,m2); }



    //
    // m3 = Q * m2
    //

    template <bool add, int ix, class T, class M1, class V1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else if (m1.isSquare()) {
            m3 = m2;
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),m3.mat());
            Scale(x,m3);
        } else {
            TMVAssert(m3.size() > m2.size());
            m3.rowRange(0,m2.size()) = m2;
            m3.rowRange(m2.size(),m3.size()).setZero();
            PackedQ_MultEq(m1.getQ(),m1.getBeta(),m3.mat());
            Scale(x,m3);
        }
    }
    template <bool add, int ix, class T, class M1, class V1, class M2, class M3>
    static TMV_INLINE void NoAliasMultMM(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMV<add>(x,m1,m2,m3); }
    template <bool add, int ix, class T, class M1, class V1, class M2, class M3>
    static TMV_INLINE void AliasMultMM(
        const Scaling<ix,T>& x, const PackedQ<M1,V1>& m1, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMV<add>(x,m1,m2,m3); }

    
    //
    // m3 = m1 * Q
    // m3t = Qt * m1t
    //

    template <bool add, int ix, class T, class M1, class M2, class V2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (add) {
            m3 += (x*m1*m2).calc();
        } else if (m2.isSquare()) {
            m3 = m1;
            typename M3::adjoint_type m3a = m3.adjoint();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        } else {
            TMVAssert(m1.rowsize() > m3.rowsize());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int cs = M1::_colsize;
            const int rs = M1::_rowsize;
            const int rm = M1::_rowmajor;
            typename MCopyHelper<T12,Rec,cs,rs,rm,false>::type m1a =
                m1.adjoint();
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1a);
            m3.adjoint() = x * m1a.rowRange(0,m3.colsize());
        }
    }
    template <bool add, int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }
    template <bool add, int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { MultMM<add>(x,m1,m2,m3); }
 

    //
    // m *= Q
    //

    template <int ix, class T, class M1, class M2, class V2>
    static inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    {
        typename M1::adjoint_type m1a = m1.adjoint();
        PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1a);
        m2.divEq(m1a);
        Scale(x,m1);
    }
    template <int ix, class T, class M1, class M2, class V2>
    static TMV_INLINE void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    { MultEqMM(m1,x,m2); }
    template <int ix, class T, class M1, class M2, class V2>
    static TMV_INLINE void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1, const Scaling<ix,T>& x,
        const PackedQ<M2,V2>& m2)
    { MultEqMM(m1,x,m2); }


    //
    // m / Q
    //

    template <int ix, class T, class M1, class M2, class V2, class M3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        if (m2.isSquare()) {
            m3 = m1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m3);
            Scale(x,m3);
        } else {
            TMVAssert(m1.colsize() > m3.colsize());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int cs = M1::_colsize;
            const int rs = M1::_rowsize;
            const int rm = M1::_rowmajor;
            typename MCopyHelper<T12,Rec,cs,rs,rm,false>::type m1c = m1;
            PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1c);
            m3 = x * m1c.rowRange(0,m3.colsize());
        }
    }
    template <int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }
    template <int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }


    //
    // m /= Q
    //

    template <class M1, class M2, class V2>
    static inline void LDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { PackedQ_LDivEq(m2.getQ(),m2.getBeta(),m1); }
    template <class M1, class M2, class V2>
    static TMV_INLINE void NoAliasLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { LDivEq(m1,m2); }
    template <class M1, class M2, class V2>
    static TMV_INLINE void AliasLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { LDivEq(m1,m2); }


    // 
    // m % Q
    // m3 = m1 Qt
    // m3t = Q m1t
    //

    template <int ix, class T, class M1, class M2, class V2, class M3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::adjoint_type m3a = m3.adjoint();
        if (m2.isSquare()) {
            m3 = m1;
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        } else {
            TMVAssert(m3.size() > m1.size());
            m3.rowRange(0,m1.size()) = m1;
            m3.rowRange(m1.size(),m3.size()).setZero();
            PackedQ_MultEq(m2.getQ(),m2.getBeta(),m3a);
            Scale(x,m3);
        }
    }
    template <int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { RDiv(x,m1,m2,m3); }
    template <int ix, class T, class M1, class M2, class V2, class M3>
    static TMV_INLINE void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const PackedQ<M2,V2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    { RDiv(x,m1,m2,m3); }


    // 
    // m %= Q
    //

    template <class M1, class M2, class V2>
    static inline void RDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    {
        typename M1::adjoint_type m1a = m1.adjoint();
        PackedQ_MultEq(m2.getQ(),m2.getBeta(),m1a);
    }
    template <class M1, class M2, class V2>
    static TMV_INLINE void NoAliasRDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { RDivEq(m1,m2); }
    template <class M1, class M2, class V2>
    static TMV_INLINE void AliasRDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const PackedQ<M2,V2>& m2)
    { RDivEq(m1,m2); }


    
    //
    // TMV_Text
    //

#ifdef TMV_TEXT
    template <class M, class V>
    static inline std::string TMV_Text(const PackedQ<M,V>& Q)
    {
        std::ostringstream s;
        s << "PackedQ< "<<TMV_Text(Q.getQ())<< " , ";
        s << TMV_Text(Q.getBeta())<<" >";
        return s.str();
    }
#endif


} // namespace mv


#endif
