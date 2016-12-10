

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
// Also, since I couldn't figure out how to make this work without
// the compiler instantiating a whole bunch of divider functions every
// time I declare a Matrix, I decided to only use the Divider object
// when the type is in the library.  Then I define all the divider functions
// in the compiled library, which keeps the sizes of object files for the 
// user code sane.  Thus, it is ok to drop down to MatrixView and VectorView
// with no Unit, ColMajor, etc attributes since we'll go through an Inst
// function anyway which recovers them.
// 


#ifndef TMV_Divider_H
#define TMV_Divider_H

#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_BaseVector.h"

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
        virtual void doSolveInPlace(MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolveInPlace(VectorView<RT> v2) const = 0;
        virtual void doSolveInPlace(VectorView<CT> v2) const = 0;
        virtual void doSolveInPlace(VectorView<CT,Conj> v2) const = 0;

        virtual void doSolveTransposeInPlace(MatrixView<RT> m2) const = 0;
        virtual void doSolveTransposeInPlace(MatrixView<CT> m2) const = 0;
        virtual void doSolveTransposeInPlace(MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolveTransposeInPlace(VectorView<RT> v2) const = 0;
        virtual void doSolveTransposeInPlace(VectorView<CT> v2) const = 0;
        virtual void doSolveTransposeInPlace(VectorView<CT,Conj> v2) const = 0;

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
            const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1, VectorView<RT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<RT>& m1, VectorView<CT,Conj> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT>& m1, VectorView<CT,Conj> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT,Conj>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolve(
            const ConstVectorView<CT,Conj>& m1, VectorView<CT,Conj> m2) const = 0;

        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<RT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<RT>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstMatrixView<CT,Conj>& m1, MatrixView<CT,Conj> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1, VectorView<RT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<RT>& m1, VectorView<CT,Conj> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT>& m1, VectorView<CT,Conj> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT,Conj>& m1, VectorView<CT> m2) const = 0;
        virtual void doSolveTranspose(
            const ConstVectorView<CT,Conj>& m1, VectorView<CT,Conj> m2) const = 0;

        // These two are for matrices that don't fit one of the above.  
        // e.g. DiagMatrix / Matrix
        // If the divider intrinsically likes in-place solves, use that.
        // Otherwise, make a temporary copy of m1.
        template <class M1, class T2, int A2>
        void doSolve(
            const M1& m1, MatrixView<T2,A2> m2) const
        {
            if (preferInPlace()) doSolveInPlace(m2=m1);
            else {
                Matrix<typename M1::value_type> m1c(m1);
                doSolve(m1c.xView(),m2);
            }
        }
        template <class M1, class T2, int A2>
        void doSolveTranspose(
            const M1& m1, MatrixView<T2,A2> m2) const
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
        virtual void doMakeInverse(MatrixView<CT,Conj> minv) const = 0;

        template <class M2>
        void makeInverse(BaseMatrix_Rec_Mutable<M2>& minv) const
        { doMakeInverse(minv.mat().xView()); }


        //
        // InverseATA
        //

        virtual void doMakeInverseATA(MatrixView<RT> ata) const = 0;
        virtual void doMakeInverseATA(MatrixView<CT> ata) const = 0;
        virtual void doMakeInverseATA(MatrixView<CT,Conj> ata) const = 0;

        template <class M2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<M2>& ata) const
        { doMakeInverseATA(ata.mat().xView()); }


        // 
        // Condition
        // Note: this may be either kappa_inf or kappa_2 depending on the
        // kind of decomposition we are doing.
        // It is kappa_2 if an SV decomposition, and kapp_inf if not.
        // The input value if m.normInf() of the original matrix,
        // which is ignored if we are doing a kappa_2 calculation.
        //

        virtual RT condition(RT normInf) const = 0;


        //
        // Does the solver intrinsically prefer to do the solution in place?
        //
        virtual bool preferInPlace() const = 0;

    private :
        // op= not allowed.
        Divider<T>& operator=(const Divider<T>&);
    };

    // Don't have CheckDecomp be a method of Divider, since it needs
    // to be different for each kind of Divider class, so it should
    // be virtual.  But you can't have virtual functions with template
    // arguments.  So do it here with explicit dynamic_cast checks
    // to determine which kind of divider it is.
    template <class T, class M2>
    static bool CheckDecomp(
        const Divider<T>* div, const BaseMatrix_Calc<M2>& m,
        std::ostream* fout=0)
    {
        TMVAssert(div);
        // Traits<>::type removes any reference that is part of lud_type1, etc.
        typedef typename Traits<typename M2::lud_type>::type lud_type;
        typedef typename Traits<typename M2::qrd_type>::type qrd_type;
        typedef typename Traits<typename M2::qrpd_type>::type qrpd_type;
        typedef typename Traits<typename M2::svd_type>::type svd_type;
        if (dynamic_cast<const lud_type*>(div)) {
            return CheckDecomp(static_cast<const lud_type&>(*div),m,fout);
        } else if (dynamic_cast<const qrd_type*>(div)) {
            return CheckDecomp(static_cast<const qrd_type&>(*div),m,fout);
        } else if (dynamic_cast<const qrpd_type*>(div)) {
            return CheckDecomp(static_cast<const qrpd_type&>(*div),m,fout);
        } else if (dynamic_cast<const svd_type*>(div)) {
            return CheckDecomp(static_cast<const svd_type&>(*div),m,fout);
        } else {
            *fout << "Couldn't cast divider to a known type.\n";
            return false;
        }
    }

    template <class T>
    class DivHelper
    {
    public:

        typedef const Divider<T> div_type;
        typedef const div_type* getdiv_type;

        // Start with XX (= 0) and if we request divtype using getDivType,
        // and it hasn't been manually set yet, then it gets reset to either 
        // LU or QR depending on whether the matrix is square.
        DivHelper() : divider(0), divtype(tmv::XX) {}

        ~DivHelper() { }

        // Set, get the DivType:
        void divideUsing(DivType dt) const
        {
            if (!(divtype & dt)) {
                unsetDiv();
                divtype &= ~tmv::DivTypeFlags;
                divtype |= dt;
            }
        }

        DivType getDivType() const 
        {
            if ((divtype & tmv::DivTypeFlags) == tmv::XX) resetDivType();
            return divtype & tmv::DivTypeFlags; 
        }

        // Set, unset in-place division
        void divideInPlace() const
        { divtype |= tmv::DivInPlaceFlag; saveDiv(); }

        void dontDivideInPlace() const
        { divtype &= ~tmv::DivInPlaceFlag; }

        bool divIsInPlace() const 
        { return divtype & tmv::DivInPlaceFlag; }

        // Set, unset whether to save the divider object
        void saveDiv() const
        { divtype |= tmv::SaveDivFlag; }

        void dontSaveDiv() const
        { divtype &= ~tmv::SaveDivFlag; }

        bool divIsSaved() const 
        { return divtype & tmv::SaveDivFlag; }

        // Set, unset the divider object
        void unsetDiv() const
        { divider.reset(0); }

        void resetDiv() const
        { unsetDiv(); setDiv(); }

        bool divIsSet() const
        { return getDiv(); }

        // This needs to be defined in the derived class:
        virtual void setDiv() const = 0;

        getdiv_type getDiv() const
        { return divider.get(); }

        void doneDiv() const
        { if (!divIsSaved()) unsetDiv(); }

    protected:

        void resetDivType() const
        { divideUsing(mIsSquare() ? tmv::LU : tmv::QR); }

        // Two more that need to be defined in the derived class:
        virtual Matrix<T> getMatrix() const = 0;
        virtual bool mIsSquare() const = 0;

        mutable std::auto_ptr<div_type> divider;
        mutable DivType divtype;

    }; // DivHelper

} // namespace tmv

#endif
