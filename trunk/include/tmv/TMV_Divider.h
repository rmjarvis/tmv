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
// This file contains the Divider class, which handles the division
// routines for matrices.  It is a virtual base class for things like
// LUDiv, QRDiv, etc.  
//
// There is a bit of a compromise required to set things up so that 
// division can be done through the virtual functions, since virtual 
// functions cannot be templates.  So we can't have the argument of 
// solve (for example) be template matrices.  They need to be explicitly
// something like MatrixView or VectorView with all template parameters
// specified.  This means that we lose any knowledge (at compile time)
// of things like the size and steps if they were known.
//
// It should be ok so long as the resulting call eventually goes through
// in Inst function, since that also (intentionally) loses this information
// to reduce the number of instantiated versions of each function.
// But if you are running with -DTMV_INLINE, then you might want to always
// call things like m.lud().solve(v1,v2) instead to make sure that the 
// compile-time knowledge isn't lost.
//


#ifndef TMV_Divider_H
#define TMV_Divider_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class Divider
    {
    public :

        typedef typename Traits<T>::real_type RT;
        typedef typename Traits<T>::complex_type CT;
        typedef typename Traits<RT>::float_type FT;
        typedef typename Traits<T>::float_type ZFT;

        //
        // Constructors -- none specified.
        // Defaults are ok, since this class doesn't have anything to
        // set up or destroy.
        //
        // We do need to include the destructor though to make it virtual.
        virtual ~Divider() {}

        // MJ: I don't remember why I need this...
        virtual bool isSV() const { return false; }

        // 
        // Perform the division in place
        //
        // Note that we have to explicitly allow both real and complex
        // types for m2, even though RT is invalid if T is complex.
        // We solve this problem later by not compiling code for
        // invalid types, and using an assert statement (below)
        // to make sure that they aren't ever called.
        //
        virtual void doSolveInPlace(MatrixView<RT> m2) const = 0;
        virtual void doSolveInPlace(MatrixView<CT> m2) const = 0;
        virtual void doSolveInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolveInPlace(VectorView<RT> v2) const = 0;
        virtual void doSolveInPlace(VectorView<CT> v2) const = 0;
        virtual void doSolveInPlace(
            VectorView<CT,UNKNOWN,true> v2) const = 0;

        virtual void doSolveTransposeInPlace(MatrixView<RT> m2) const = 0;
        virtual void doSolveTransposeInPlace(MatrixView<CT> m2) const = 0;
        virtual void doSolveTransposeInPlace(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTransposeInPlace(VectorView<RT> v2) const = 0;
        virtual void doSolveTransposeInPlace(VectorView<CT> v2) const = 0;
        virtual void doSolveTransposeInPlace(
            VectorView<CT,UNKNOWN,true> v2) const = 0;

        // The templated versions will take any (valid) argument for 
        // the rhs, and convert it to a MatrixView or VectorView.
        // We also check that the real/complex combo is valid.
        // And, we require that aside from the real/complex-ness,
        // that the underlying types are the same.  If you want to
        // divide a double vector by a float matrix for example, you
        // can, but you need to use something like m.lud().solveInPlace(v)
        // instead of going through the normal division operator.
        template <class M2> 
        void solveInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
        { 
            TMVStaticAssert(Traits<T>::isreal || M2::iscomplex);
            TMVStaticAssert((Traits2<T,typename M2::value_type>::samebase));
            doSolveInPlace(m2.mat().xView()); 
        }
        template <class V2> 
        void solveInPlace(BaseVector_Mutable<V2>& v2) const
        {
            TMVStaticAssert(Traits<T>::isreal || V2::iscomplex);
            TMVStaticAssert((Traits2<T,typename V2::value_type>::samebase));
            doSolveInPlace(v2.vec().xView()); 
        }

        template <class M2> 
        void solveTransposeInPlace(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            TMVStaticAssert(Traits<T>::isreal || M2::iscomplex);
            TMVStaticAssert((Traits2<T,typename M2::value_type>::samebase));
            doSolveTransposeInPlace(m2.mat().xView()); 
        }
        template <class V2>
        void solveTransposeInPlace(BaseVector_Mutable<V2>& v2) const
        {
            TMVStaticAssert(Traits<T>::isreal || V2::iscomplex);
            TMVStaticAssert((Traits2<T,typename V2::value_type>::samebase));
            doSolveTransposeInPlace(v2.vec().xView()); 
        }


        //
        // Next the not-in-place division 
        // 
        virtual void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<RT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT>& m1,
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1,
            MatrixView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1,
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1, VectorView<RT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& m1,
            VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT,UNKNOWN,true>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;

        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1,
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT>& m1, 
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1,
            MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT,UNKNOWN,UNKNOWN,true>& m1,
            MatrixView<CT,UNKNOWN,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1, VectorView<RT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& m1,
            VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT,UNKNOWN,true>& m1,
            VectorView<CT,UNKNOWN,true> m2) const = 0;

        // These two are for matrices that don't fit one of the above.  
        // e.g. DiagMatrix / Matrix
        // If the divider intrinsically likes in-place solves, use that.
        // Otherwise, make a temporary copy of m1.
        template <class M1, class T2, bool C2>
        void doSolve(
            const M1& m1, MatrixView<T2,UNKNOWN,UNKNOWN,C2> m2) const
        {
            if (preferInPlace()) doSolveInPlace(m2=m1);
            else {
                Matrix<typename M1::value_type> m1c(m1);
                doSolve(m1c.xView(),m2);
            }
        }
        template <class M1, class T2, bool C2>
        void doSolveTranspose(
            const M1& m1, MatrixView<T2,UNKNOWN,UNKNOWN,C2> m2) const
        {
            if (preferInPlace()) doSolveTransposeInPlace(m2=m1);
            else {
                Matrix<typename M1::value_type> m1c(m1);
                doSolveTranspose(m1c.xView(),m2);
            }
        }

        template <class M1, class M2> 
        void solve(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolve(m1.mat().xView(),m2.mat().xView()); }
        template <class V1, class V2> 
        void solve(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { doSolve(v1.vec().xView(),v2.vec().xView()); }

        template <class M1, class M2> 
        void solveTranspose(
            const BaseMatrix<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2) const
        { doSolveTranspose(m1.mat().xView(),m2.mat().xView()); }
        template <class V1, class V2> 
        void solveTranspose(
            const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2) const
        { doSolveTranspose(v1.vec().xView(),v2.vec().xView()); }


        //
        // Determinant
        //
        
        virtual T det() const = 0;
        virtual FT logDet(ZFT* sign) const = 0;
        virtual bool isSingular() const = 0;

        //
        // Inverse
        //
        

        virtual void doMakeInverse(MatrixView<RT> minv) const = 0;
        virtual void doMakeInverse(MatrixView<CT> minv) const = 0;
        virtual void doMakeInverse(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> minv) const = 0;

        template <class M2> 
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
        { doMakeInverse(minv.mat().xView()); }


        //
        // InverseATA
        //

        virtual void doMakeInverseATA(MatrixView<RT> ata) const = 0;
        virtual void doMakeInverseATA(MatrixView<CT> ata) const = 0;
        virtual void doMakeInverseATA(
            MatrixView<CT,UNKNOWN,UNKNOWN,true> ata) const = 0;

        template <class M2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
        { doMakeInverseATA(ata.mat().xView()); }


        // 
        // Condition
        // Note: this may be either kappa_inf or kappa_2 depending on the
        // kind of decomposition we are doing.
        // It is kappa_2 if isSV() == true, and kapp_inf if isSV() == false.
        // The input value if m.normInf() of the original matrix,
        // which is ignored if we are doing a kappa_2 calculation.
        //

        virtual RT condition(RT normInf) const = 0;


        //
        // Check the validity of the decomposition 
        // (Used for testing)
        //

        virtual bool doCheckDecomp(
            const ConstMatrixView<T>& m, std::ostream* fout) const = 0;
        virtual bool doCheckDecomp(
            const ConstMatrixView<T,UNKNOWN,UNKNOWN,true>& m,
            std::ostream* fout) const = 0;

        template <class M2>
        bool checkDecomp(
            const BaseMatrix_Rec<M2>& m, std::ostream* fout) const
        { return doCheckDecomp(m.mat().xView(),fout); }

        //
        // Does the solver intrinsically prefer to do the solution in place?
        //
        virtual bool preferInPlace() const = 0;

    private :
        // op= not allowed.
        Divider<T>& operator=(const Divider<T>&);
    };
    
} // namespace tmv

#endif
