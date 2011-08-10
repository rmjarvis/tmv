

#include <iostream>
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixIO.h"

namespace tmv {

    template <class T, int C>
    void InstWrite(std::ostream& os, const ConstDiagMatrixView<T,C>& m)
    {
        if (m.step() == 1)
            InlineWrite(os,m.unitView()); 
        else
            InlineWrite(os,m); 
    }

    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,C>& m,
        typename Traits<T>::real_type thresh)
    {
        if (m.step() == 1)
            InlineWrite(os,m.unitView(),thresh); 
        else
            InlineWrite(os,m,thresh); 
    }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


