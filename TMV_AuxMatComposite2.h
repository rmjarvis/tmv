// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a Matrix and is of the form:
// (S) op (M)
// (ie. the second operand is the normal Matrix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenBandMatrix
//
// #define SUMMM SumMB
// #define SUMMMa SumBM
// #define PRODMM ProdBM


//
// Matrix + Matrix
//

template <class T, class T1, class T2> class SUMMM :
  public SUMMMa<T,T2,T1>
{
  public:
    SUMMM(T _x1, const GENMATRIX1<T1>& _m1, T _x2,
	const GenMatrix<T2>& _m2) :
      SUMMMa<T,T2,T1>(_x2,_m2,_x1,_m1) {}
};

#include "TMV_AuxSumMM.h"

//
// Matrix * Matrix
//

template <class T, class T1, class T2> class PRODMM : 
  public MatrixComposite<T>
{
  public:
    PRODMM(const T _x, const GENMATRIX1<T1>& _m1,
	const GenMatrix<T2>& _m2) :
      MatrixComposite<T>(BaseStorOf(_m2)),
      x(_x), m1(&_m1), m2(&_m2)
    { TMVAssert(m1->rowsize() == m2->colsize()) ; }
    inline size_t colsize() const { return m1->colsize(); }
    inline size_t rowsize() const { return m2->rowsize(); }
    inline T GetX() const { return x; }
    inline const GENMATRIX1<T1>& GetM1() const { return *m1; }
    inline const GenMatrix<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      MultMM(x, *m1, *m2, T(0), m0);
    }
  private:
    const T x;
    const GENMATRIX1<T1>*const m1;
    const GenMatrix<T2>*const m2;
};

#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"


