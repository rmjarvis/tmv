// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a Matrix and is of the form:
// (M) op (S)
// (ie. the first operand is the normal Matrix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define SUMMM SumMB
// #define PRODMM ProdMB
// #define QUOTMM QuotMB
// #define RQUOTMM RQuotMB
// #define TQUOTMM TransientQuotMB
// #define TRQUOTMM TransientRQuotMB


//
// Matrix + Matrix
//

template <class T, class T1, class T2> class SUMMM : 
  public MatrixComposite<T> 
{
  public:
    SUMMM(const T _x1, const GenMatrix<T1>& _m1, 
	const T _x2, const GENMATRIX2<T2>& _m2) :
      MatrixComposite<T>(BaseStorOf(_m1)),
      x1(_x1), m1(&_m1), x2(_x2), m2(&_m2)
    { 
      TMVAssert(m1->colsize() == m2->colsize());
      TMVAssert(m1->rowsize() == m2->rowsize()); 
    }
    inline size_t colsize() const { return m1->colsize(); }
    inline size_t rowsize() const { return m1->rowsize(); }
    inline T GetX1() const { return x1; }
    inline const GenMatrix<T1>& GetM1() const { return *m1; }
    inline T GetX2() const { return x2; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (m0.SameStorageAs(*m1)) {
	m0 = *m1;
	AddMM(x2,*m2,x1,m0);
      } else if (m0.SameStorageAs(*m2) || x1 != T(1)) {
	m0 = *m2;
	AddMM(x1,*m1,x2,m0);
      } else {
	m0 = *m1;
	AddMM(x2,*m2,x1,m0);
      }
    }
  private:
    const T x1;
    const GenMatrix<T1>* m1;
    const T x2;
    const GENMATRIX2<T2>* m2;
};

template <class T> inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(T(1),m2,T(1),m1); return m1; 
}

template <class T> inline const MatrixView<CT>& operator+=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2) 
{
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(CT(1),m2,CT(1),m1); 
  return m1; 
}

template <class T> inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(T(-1),m2,T(1),m1); 
  return m1; 
}

template <class T> inline const MatrixView<CT>& operator-=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(CT(-1),m2,CT(1),m1); 
  return m1; 
}


//
// Matrix * Matrix
//

template <class T, class T1, class T2> class PRODMM : 
  public MatrixComposite<T>
{
  public:
    PRODMM(const T _x, const GenMatrix<T1>& _m1, const GENMATRIX2<T2>& _m2) :
      MatrixComposite<T>(BaseStorOf(_m1)),
      x(_x), m1(&_m1), m2(&_m2)
    { TMVAssert(m1->rowsize() == m2->colsize()) ; }
    inline size_t colsize() const { return m1->colsize(); }
    inline size_t rowsize() const { return m2->rowsize(); }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return *m1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      MultMM(x, *m1, *m2, T(0), m0);
    }
  private:
    const T x;
    const GenMatrix<T1>*const m1;
    const GENMATRIX2<T2>*const m2;
};

template <class T> inline const MatrixView<T>& operator*=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize()==m2.rowsize());
  TMVAssert(m1.rowsize()==m2.rowsize());
  MultMM(T(1),m1,m2,T(0),m1); 
  return m1; 
}

template <class T> inline const MatrixView<CT>& operator*=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize()==m2.rowsize());
  TMVAssert(m1.rowsize()==m2.rowsize());
  MultMM(CT(1),m1,m2,CT(0),m1); 
  return m1; 
}

template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
  return m; 
}

template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),T(1),m); 
  return m; 
}

//
// Matrix / % Matrix
//

template <class T, class T1, class T2> class QUOTMM : 
  public MatrixComposite<T>
{
  public:
    QUOTMM(const T _x, const GenMatrix<T1>& _m1, const GENMATRIX2<T2>& _m2) :
      MatrixComposite<T>(BaseStorOf(_m1)),
      x(_x), m1(&_m1), m2(&_m2)
    { TMVAssert( m1->colsize() == m2->colsize() ); }
    inline size_t colsize() const { return m2->rowsize(); }
    inline size_t rowsize() const { return m1->rowsize(); }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return *m1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (m0.SameAs(*m1)) m2->LDivEq(m0);
      else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	if (m0.isrm()) {
	  Matrix<T,RowMajor> temp(m0.colsize(),m0.rowsize());
	  m2->LDiv(*m1,temp.View());
	  m0 = temp;
	} else {
	  Matrix<T,ColMajor> temp(m0.colsize(),m0.rowsize());
	  m2->LDiv(*m1,temp.View());
	  m0 = temp;
	}
      } else {
	m2->LDiv(*m1,m0);
      }
      if (x != T(1)) m0 *= x;
    }
  protected:
    const T x;
    const GenMatrix<T1>*const m1;
    const GENMATRIX2<T2>*const m2;
};

template <class T, class T1, StorageType S1, class T2> class TQUOTMM : 
  public QUOTMM<T,T1,T2>
{
  public :
    TQUOTMM(const Matrix<T1,S1>* m1, const GENMATRIX2<T2>& m2) :
      QUOTMM<T,T1,T2>(T(1),*m1,m2) {}
    ~TQUOTMM() { delete this->m1; }
};

template <class T, class T1, class T2> class RQUOTMM : 
  public MatrixComposite<T>
{
  public:
    RQUOTMM(const T _x, const GenMatrix<T1>& _m1, const GENMATRIX2<T2>& _m2) :
      MatrixComposite<T>(BaseStorOf(_m1)),
      x(_x), m1(&_m1), m2(&_m2)
    { TMVAssert( m1->rowsize() == m2->rowsize() ); }
    inline size_t colsize() const { return m1->colsize(); }
    inline size_t rowsize() const { return m2->colsize(); }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return *m1; }
    inline const GENMATRIX2<T2>& GetM2() const { return *m2; }
    inline void AssignTo(const MatrixView<T>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      if (m0.SameAs(*m1)) m2->RDivEq(m0);
      else if (m0.SameStorageAs(*m1) || m0.SameStorageAs(*m2)) {
	if (m0.isrm()) {
	  Matrix<T,RowMajor> temp(m0.colsize(),m0.rowsize());
	  m2->RDiv(*m1,temp.View());
	  m0 = temp;
	} else {
	  Matrix<T,ColMajor> temp(m0.colsize(),m0.rowsize());
	  m2->RDiv(*m1,temp.View());
	  m0 = temp;
	}
      } else {
	m2->RDiv(*m1,m0);
      }
      if (x != T(1)) m0 *= x;
    }
  protected:
    const T x;
    const GenMatrix<T1>*const m1;
    const GENMATRIX2<T2>*const m2;
};

template <class T, class T1, StorageType S1, class T2> class TRQUOTMM : 
  public RQUOTMM<T,T1,T2>
{
  public :
    TRQUOTMM(const Matrix<T1,S1>* m1, const GENMATRIX2<T2>& m2) :
      RQUOTMM<T,T1,T2>(T(1),*m1,m2) {}
    ~TRQUOTMM() { delete this->m1; }
};

template <class T> inline const MatrixView<T>& operator/=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.colsize() == m2.rowsize());
  m2.LDivEq(m1); 
  return m1; 
}

template <class T> inline const MatrixView<CT>& operator/=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.colsize() == m2.rowsize());
  m2.LDivEq(m1); 
  return m1; 
}

template <class T> inline const MatrixView<T>& operator%=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  m2.RDivEq(m1); 
  return m1; 
}

template <class T> inline const MatrixView<CT>& operator%=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  m2.RDivEq(m1); 
  return m1; 
}


template <class T, class Tm> inline const MatrixView<T>& operator*=(
    const MatrixView<T>& m1, const QUOTXM<T,Tm>& qxm)
{ 
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
  qxm.GetM().RDivEq(m1); 
  m1 *= qxm.GetX();
  return m1; 
}

template <class T> inline const MatrixView<CT>& operator*=(
    const MatrixView<CT>& m1, const QUOTXM<T,T>& qxm)
{ 
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
  qxm.GetM().RDivEq(m1); 
  m1 *= qxm.GetX();
  return m1; 
}


