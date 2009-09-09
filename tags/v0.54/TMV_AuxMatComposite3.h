// This file sets up the Composite QuotXM classes for Matrix types which
// returns a Matrix for the operations x/m or x%m.
// (Some do not: DiagMatrix for example, returns another DiagMatix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define QUOTXM QuotXB


//
// Scalar / Matrix
//

template <class T, class Tm> class QUOTXM : 
  public MatrixComposite<T> 
{
  public:
    QUOTXM(const T _x, const GENMATRIX2<Tm>& _m) :
      MatrixComposite<T>(ColMajor), x(_x), m(&_m) {}
    inline size_t colsize() const { return m->rowsize(); }
    inline size_t rowsize() const { return m->colsize(); }
    inline T GetX() const { return x; }
    inline const GENMATRIX2<Tm>& GetM() const { return *m; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (x == T(0)) m0.Zero();
      else {
	// Need temporary, since T, Tm are not the same.
	// T = Tm overridden below.
	if (m0.isrm()) {
	  Matrix<Tm,RowMajor> temp(m0.colsize(),m0.rowsize());
	  m->Inverse(temp.View());
	  if (x != T(1)) m0 = temp*x;
	  else m0 = temp;
	} else {
	  Matrix<Tm,ColMajor> temp(m0.colsize(),m0.rowsize());
	  m->Inverse(temp.View());
	  if (x != T(1)) m0 = temp*x;
	  else m0 = temp;
	} 
      }
    }
  private:
    const T x;
    const GENMATRIX2<Tm>*const m;
};

// When T,Tm are the same, don't need temporary.
template <class T> class QUOTXM<T,T> : 
  public MatrixComposite<T> 
{
  public:
    QUOTXM(const T _x, const GENMATRIX2<T>& _m) :
      MatrixComposite<T>(ColMajor), x(_x), m(&_m) {}
    inline size_t colsize() const { return m->rowsize(); }
    inline size_t rowsize() const { return m->colsize(); }
    inline T GetX() const { return x; }
    inline const GENMATRIX2<T>& GetM() const { return *m; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (x == T(0)) m0.Zero();
      else {
	m->Inverse(m0);
	if (x != T(1)) m0 *= x;
      }
    }
  private:
    const T x;
    const GENMATRIX2<T>*const m;
};

#include "TMV_AuxQuotXM.h"
