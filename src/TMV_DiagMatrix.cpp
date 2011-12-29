

#include <iostream>
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixIO.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyV.h"

namespace tmv {

    template <class T, int C>
    void InstWrite(const TMV_Writer& writer, const ConstDiagMatrixView<T,C>& m)
    {
        if (m.step() == 1)
            InlineWrite(writer,m.unitView()); 
        else
            InlineWrite(writer,m); 
    }

    template <class T>
    void InstRead(const TMV_Reader& reader, DiagMatrixView<T> m)
    {
        if (m.step() == 1) {
            DiagMatrixView<T,Unit> mu = m.unitView();
            InlineRead(reader,mu);
        } else
            InlineRead(reader,m);
    }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


