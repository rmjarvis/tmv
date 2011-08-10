
//
// This file contains the code to implement Givens rotations
//
// A Givens rotation is a 2x2 matrix of the form:
//
// G = [  c  s  ]
//     [ -s* c* ]
//
// where |c|^2 + |s|^2 = 1.
// Normally, people consider real Givens matrices where the *'s
// (conjugations) are not necessary.  But they generalize easily
// to the complex case as above.
//
// Givens matrices are used to zero out a single element in a matrix by
// rotating a 2 element vector so that its second element is 0:
//
// [  c  s  ] [ x ] = [ r ]
// [ -s* c* ] [ y ]   [ 0 ]
// 
// [ Note: often (eg. in Golub and Van Loan) Givens matrices are defined
//   as the transpose (or adjoint) of this.  I prefer this definition. ]
//
// (1)  c x + s y = r
// (2)  c* y - s* x = 0
//
// From (2), we get  s* = yc*/x, so 
// |x|^2 |s|^2 = |y|^2 |c|^2 = |y|^2 (1 - |s|^2)
//
// |s| = |y|/sqrt(|x|^2+|y|^2)  and |c| = |x|/sqrt(|x|^2+|y|^2)
//
// From (1), r = cx+sy = cx+yy*c/x* = c/x* (|x|^2+|y|^2)
// |r| = sqrt(|x|^2+|y|^2)
//
// All that remain is to determine the phases of c,s
//
// Let x = |x| e^ip
//     y = |y| e^iq
//     r = |r| e^it
//     c = |x|/|r| e^ia
//     s = |y|/|r| e^ib
//
// Plugging these into (1), we find that  t = a+p = b+q
//
// Since a priori, we should consider all three of (a,b,t)
// to be unknown, we are left with one free parameter.
// There are (at least) 3 reasonable options:
//
// Choice 1: a = -p, b = -q, t = 0
//
// This has the advantage that r is real and positive.
// Also, the equations for c,s look much like the usual
// forms for real Givens matrices:
//
// c = x*/|r|  s = y*/|r|   r = |r|
//
// Choice 2: a = 0, b = p-q, t = p
//
// This has the advantage that c is real, so the Givens
// matrix only needs to store 3 real numbers, not 4.
// Multiplication later is also faster by a factor
// of 4/3.
//
// c = |x|/|r|  s = (x/|x|)y*/|r|  r = (x/|x|)|r|
//
// Choice 3: a = q-p, b = 0, t = q
//
// Similar to choice 2, but s is real:
// 
// c = (y/|y|)x*/|r|  s = |y|/|r|  r = (y/|y|)|r|
//
// 
// The speed of calculation is usually more important than the 
// convenience of a rotating to a real value, so choice 2 or 3
// is preferred.  
// Choice 2 is a bit better when |y|<|x|, since we prefer to 
// calculate sqrt(1+|y/x|^2), and the calculations for 2 work 
// better in this case.
// Since this is more common than |x|<|y|, we always use choice 2.
//

#ifndef TMV_Givens_H
#define TMV_Givens_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_BaseMatrix_Rec.h"

namespace tmv {

    // TODO: Put into the normal Helper structure with Inst, etc.
    
    template <class T> 
    class Givens;

    // Use a Givens matrix G to rotate the vector so y = 0:
    // G [ x ] = [ r ]
    //   [ y ]   [ 0 ]
    // Also, return the Givens matrix used to do so.
    template <class T> 
    static inline Givens<T> GivensRotate(T& x, T& y)
    {
        // Use a Givens matrix G to rotate the vector so y = 0:
        // G [ x ] = [ r ]
        //   [ y ]   [ 0 ]
        // Also, return the Givens matrix used to do so.
        
        //std::cout<<"Start GivensRotate: x,y = "<<x<<','<<y<<std::endl;
        typedef typename Traits<T>::real_type RT;
        const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

        RT absx = TMV_ABS2(x);
        RT absy = TMV_ABS2(y);
        if (absy == RT(0)) {
            y = RT(0);
            //std::cout<<"Simple1: return 1,0\n";
            return Givens<T>(RT(1),T(0));
        } else if (absx == RT(0)) {
            x = RT(0);
            RT absy = TMV_ABS(y);
            T s = TMV_SIGN(TMV_CONJ(y),absy);
            x = absy; y = RT(0);
            //std::cout<<"Simple2: return 0,"<<s<<std::endl;
            return Givens<T>(RT(0),s);
        } else if (absx > absy) {
            if (absy <= sqrteps*absx) {
                // Then everything simplifies:
                // c = 1
                // s = (y/x)* 
                // r = f
                T s = TMV_CONJ(y)/TMV_CONJ(x);
                y = RT(0);
                //std::cout<<"Simple3: return 0,"<<s<<std::endl;
                return Givens<T>(RT(1),s);
            } else {
                // c = 1/sqrt(1+|y/x|^2)
                // s = (y/x)*/sqrt(1+|y/x|^2)
                // r = x sqrt(1+|y/x|^2)
                // We get a slightly more accurate calculation of r
                // if we calculate r-x and add this to x:
                // r-x = x (sqrt(1+|y/x|^2)-1) = x |y/x|^2 / (1 + sqrt(1+|y/x|^2))
                T yoverx = y/x;
                RT n = TMV_NORM(yoverx);
                RT sqrtfactor = TMV_SQRT(RT(1)+n);
                RT c = RT(1)/sqrtfactor;
                T s = TMV_CONJ(yoverx)*c;
                x += x*(n/(RT(1)+sqrtfactor));
                y = RT(0);
                //std::cout<<"x>y: return "<<c<<','<<s<<std::endl;
                return Givens<T>(c,s);
            }
        } else {
            T xovery = x/y;
            RT n = TMV_NORM(xovery);
            RT absxovery = TMV_SQRT(n);
            if (n <= TMV_Epsilon<T>()) {
                // c = |x/y|
                // s = (x/y)/|x/y|
                // r = x/|x/y|
                T s = TMV_SIGN(xovery,absxovery);
                x = s*y;
                y = RT(0);
                //std::cout<<"Simple4: return "<<absxovery<<','<<s<<std::endl;
                return Givens<T>(absxovery,s);
            } else {
                // c = |x/y|/sqrt(1+|x/y|^2)
                // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
                // r = x/|x/y| sqrt(|x|^2+|y|^2)
                RT sqrtfactor = TMV_SQRT(RT(1)+n);
                RT invsqrtfactor = RT(1)/sqrtfactor;
                T signxovery = TMV_SIGN(xovery,absxovery); // (x/y)/|x/y|
                T s = signxovery*invsqrtfactor;
                x = y * signxovery * sqrtfactor;
                y = T(0);
                RT c = absxovery*invsqrtfactor;
                //std::cout<<"x<y: return "<<c<<','<<s<<std::endl;
                return Givens<T>(c,s);
            }
        }
    }

    template <class RTg, class Tg, class T> 
    static inline void GivensMult(RTg c, Tg s, T& x, T& y)
    {
        // [ x' ] = [  c  s ] [ x ] = [ cx+sy  ]
        // [ y' ]   [ -s* c ] [ y ]   [ cy-s*x ]
        T xx = c*x+s*y;
        y = c*y-TMV_CONJ(s)*x;
        x = xx;
    }
         
    template <class RTg, class Tg, class T> 
    static TMV_INLINE void GivensMult(RTg c, Tg s, ConjRef<T> x, ConjRef<T> y)
    { GivensMult(c,TMV_CONJ(s),x.getRef(),y.getRef()); }
         
    template <class RTg, class Tg, class V1, class V2>
    static inline void GivensMultV(
        RTg c, Tg s, BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVAssert(v1.size() == v2.size());
        TMVAssert(!SameStorage(v1,v2));

        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;

        IT1 it1 = v1.begin();
        IT2 it2 = v2.begin();

        for(int i=v1.size();i>0;--i) {
            GivensMult(c,s,*it1++,*it2++);
        }
    }

    //
    // Symmetric Givens Mult
    //

    template <class RTg, class Tg, class T> 
    static inline void GivensHermMult(RTg c, Tg s, T& d0, T& d1, T& e0)
    {
        // [ d0 e0* ] = [  c  s ] [ d0 e0* ] [ c  -s ]
        // [ e0 d1  ]   [ -s* c ] [ e0 d1  ] [ s*  c ]
        // =[c^2 d0 + 2c Re(s e0) + |s|^2 d1   cs(d1-d0) + c^2 e0* - s^2 e0   ]
        //  [cs*(d1-d0) + c^2 e0 - s*^2 e0*    c^2 d1 - 2c Re(s e0) + |s|^2 d0]
        // (using c^2 = 1-|s|^2:)
        // d0' = d0 + 2c Re(s e0) + |s|^2 (d1-d0)
        // e0' = e0 - 2s* Re(s e0) + c s* (d1-d0)
        // d1' = d1 - 2c Re(s e0) - |s|^2 (d1-d0)

        typedef typename Traits<T>::real_type RT;
        RTg s2 = TMV_NORM(s);
        RT Rese0 = TMV_REAL(s*e0);
        T d1md0 = d1-d0;
        T dd = RT(2)*c*Rese0 + s2*d1md0;
        d0 += dd;
        d1 -= dd;
        e0 += TMV_CONJ(s)*(c*d1md0 - RT(2)*Rese0);
    }
          
    template <class RTg, class Tg, class T> 
    static inline void GivensSymMult(RTg c, Tg s, T& d0, T& d1, T& e0)
    {
        // [ d0 e0 ] = [  c  s ] [ d0 e0 ] [ c -s* ]
        // [ e0 d1 ]   [ -s* c ] [ e0 d1 ] [ s  c  ]
        // d0' = d0 + 2c (s e0) + s^2 (d1-d0)
        // e0' = e0 - 2|s|^2 e0 + c (s d1 - s* d0)
        // d1' = d1 - 2c (s e0) - s^2 (d1-d0)

        typedef typename Traits<T>::real_type RT;
        Tg s2 = s*s;
        T se0 = s*e0;
        T d1md0 = d1-d0;
        T dd = RT(2)*c*se0 + s2*d1md0;
        d0 += dd;
        d1 -= dd;
        if (Traits<T>::isreal)
            e0 += s*(c*d1md0 - RT(2)*se0);
        else
            e0 += c*(s*d1-TMV_CONJ(s)*d0) - RT(2)*TMV_CONJ(s)*se0;
    }

    template <class T> 
    class Givens 
    {
        typedef typename Traits<T>::real_type RT;

    public:

        // Constructors
        Givens(const RT& _c, const T& _s) : c(_c), s(_s) {}
        Givens() {} // Uninitialized
        ~Givens() {}
        // Use default copy, op=

        Givens<T> transpose() const 
        { return Givens<T>(c,-TMV_CONJ(s)); }

        Givens<T> conjugate() const 
        { return Givens<T>(c,TMV_CONJ(s)); }

        Givens<T> adjoint() const 
        { return Givens<T>(c,-s); }

        template <class T2> 
        void mult(T2& x, T2& y) const 
        { GivensMult(c,s,x,y); }
        template <class T2> 
        void conjMult(T2& x, T2& y) const 
        { GivensMult(c,TMV_CONJ(s),x,y); }

        template <class T2> 
        void mult(ConjRef<T2> x, ConjRef<T2> y) const 
        { GivensMult(c,TMV_CONJ(s),x.getRef(),y.getRef()); }
        template <class T2> 
        void conjMult(ConjRef<T2> x, ConjRef<T2> y) const 
        { GivensMult(c,s,x.getRef(),y.getRef()); }

        template <class V2> 
        void mult(BaseVector_Mutable<V2>& v) const
        { TMVAssert(v.size()==2); GivensMult(c,s,v.ref(0),v.ref(1)); }
        template <class V2> 
        void conjMult(BaseVector_Mutable<V2>& v) const
        { TMVAssert(v.size()==2); GivensMult(c,TMV_CONJ(s),v.ref(0),v.ref(1)); }

        template <class M2> 
        void mult(BaseMatrix_Rec_Mutable<M2>& m)
        { 
            TMVAssert(m.colsize()==2);
            typedef typename M2::row_type M2r;
            M2r m0 = m.row(0);
            M2r m1 = m.row(1);
            GivensMultV(c,s,m0,m1);
        }
        template <class M2> 
        void conjMult(BaseMatrix_Rec_Mutable<M2>& m)
        { 
            TMVAssert(m.colsize()==2);
            typedef typename M2::row_type M2r;
            M2r m0 = m.row(0);
            M2r m1 = m.row(1);
            GivensMultV(c,TMV_CONJ(s),m0,m1);
        }

        template <class T2> 
        void hermMult(T2& d0, T2& d1, T2& e0) const 
        { GivensHermMult(c,s,d0,d1,e0); }
        template <class T2> 
        void symMult(T2& d0, T2& d1, T2& e0) const 
        { GivensSymMult(c,s,d0,d1,e0); }

    private:

        RT c;
        T s;
    };

} // namespace tmv

#endif
