
#include "tmv/TMV_ScaleB.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    template <class T, class M>
    static void DoScale(const T x, M& m)
    {
        if (x == T(-1))
            InlineScale(Scaling<-1,T>(x),m); 
        else if (x == T(0))
            m.setZero();
        else if (x != T(1))
            InlineScale(Scaling<0,T>(x),m); 
    }

    template <class T, class M>
    static void DoScale(const std::complex<T> x, M& m)
    { 
        if (imag(x) == T(0)) {
            if (real(x) == T(-1))
                InlineScale(Scaling<-1,T>(real(x)),m); 
            else if (real(x) == T(0))
                m.setZero();
            else if (real(x) != T(1))
                InlineScale(Scaling<0,T>(real(x)),m); 
        } else 
            InlineScale(Scaling<0,std::complex<T> >(x),m); 
    }

    template <class T>
    void InstScale(const T x, BandMatrixView<T> m)
    {
        if (m.iscm()) {
            BandMatrixView<T,ColMajor> mcm = m.cmView();
            DoScale(x,mcm); 
        } else if (m.isrm()) {
            BandMatrixView<T,RowMajor> mrm = m.rmView();
            DoScale(x,mrm); 
        } else if (m.isdm()) {
            BandMatrixView<T,DiagMajor> mdm = m.dmView();
            DoScale(x,mdm); 
        } else {
            DoScale(x,m); 
        }
    }


#define InstFile "TMV_ScaleB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


