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
// #define SUMMM1 SumMB
// #define SUMMM SumBM
// #define PRODMM ProdBM


//
// Matrix + Matrix
//

template <class T, class T2, class T4> class SUMMM :
  public SUMMM1<T,T4,T2>
{
  public:
    SUMMM(T _x1, const GENMATRIX1<T2>& _m2, T _x3,
	const GenMatrix<T4>& _m4) :
      SUMMM1<T,T4,T2>(_x3,_m4,_x1,_m2) {}
};


//
// Matrix * Matrix
//

template <class T, class T2, class T3, bool trans> class PRODMM : 
  public MatrixComposite<T>
{
  public:
    PRODMM(const T _x1, const GENMATRIX1<T2>& _m2,
	const GenMatrix<T3>& _m3) :
      MatrixComposite<T>(BaseStorOf(_m3)),
      x1(_x1), m2(&_m2), m3(&_m3)
    { TMVAssert(m2->rowsize() == m3->colsize()) ; }
    inline size_t colsize() const { return m2->colsize(); }
    inline size_t rowsize() const { return m3->rowsize(); }
    inline T GetX1() const { return x1; }
    inline const GENMATRIX1<T2>& GetM2() const { return *m2; }
    inline const GenMatrix<T3>& GetM3() const { return *m3; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (trans)
	MultMM(x1, *m2, m3->Transpose(), T(0), m0);
      else
	MultMM(x1, *m2, *m3, T(0), m0);
    }
  private:
    const T x1;
    const GENMATRIX1<T2>*const m2;
    const GenMatrix<T3>*const m3;
};

