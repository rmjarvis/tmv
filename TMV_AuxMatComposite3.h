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

template <class T, class T1, class T2> class QUOTXM : 
  public MatrixComposite<T> 
{
  public:
    QUOTXM(const T1 _x1, const GENMATRIX2<T2>& _m2) :
      MatrixComposite<T>(ColMajor), x1(_x1), m2(&_m2) {}
    inline size_t colsize() const { return m2->rowsize(); }
    inline size_t rowsize() const { return m2->colsize(); }
    inline T1 GetX1() const { return x1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (x1 == T1(0)) m0.Zero();
      else {
	Matrix<T2,ColMajor> temp(m0.colsize(),m0.rowsize());
	m2->Inverse(temp.View());
	if (x1 != T1(1)) m0 = temp*x1;
	else m0 = temp;
      }
    }
  private:
    const T1 x1;
    const GENMATRIX2<T2>*const m2;
};

// When T,T2 are the same, don't need temporary.
template <class T, class T1> class QUOTXM<T,T1,T> : 
  public MatrixComposite<T> 
{
  public:
    QUOTXM(const T1 _x1, const GENMATRIX2<T>& _m2) :
      MatrixComposite<T>(ColMajor), x1(_x1), m2(&_m2) {}
    inline size_t colsize() const { return m2->rowsize(); }
    inline size_t rowsize() const { return m2->colsize(); }
    inline T1 GetX1() const { return x1; }
    inline const GENMATRIX2<T>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (x1 == T1(0)) m0.Zero();
      else {
	m2->Inverse(m0);
	if (x1 != T1(1)) m0 *= x1;
      }
    }
  private:
    const T1 x1;
    const GENMATRIX2<T>*const m2;
};
