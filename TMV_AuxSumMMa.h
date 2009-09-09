// Things that need to be #defined on entry:
// (The values for a normal Matrix+Matrix are given)
//
// SUMMM	SumMM
// GENMATRIX1	GenMatrix
// GENMATRIX2	GenMatrix
// PRODXM1	ProdXM
// PRODXM2	ProdXM

// -(x*m+x*m)
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator-(
    const SUMMM<T,T1,T2>& smm)
{ return SUMMM<T,T1,T2>(-smm.GetX1(),smm.GetM1(),-smm.GetX2,smm.GetM2); }

// x*(x*m+x*m)
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator*(
    const T x, const SUMMM<T,T1,T2>& smm)
{ return SUMMM<T,T1,T2>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const T x, const SUMMM<CT,T1,T2>& smm)
{ return SUMMM<CT,T1,T2>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator*(
    const CT x, const SUMMM<T,T,T>& smm)
{ return SUMMM<CT,T,T>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator*(
    const CCT x, const SUMMM<T,T,T>& smm)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T> inline SUMMM<CT,T,T> operator*(
    const VCT x, const SUMMM<T,T,T>& smm)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const CCT x, const SUMMM<CT,T1,T2>& smm)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const VCT x, const SUMMM<CT,T1,T2>& smm)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}


// (x*m+x*m)*x
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator*(
    const SUMMM<T,T1,T2>& smm, const T x)
{ return SUMMM<T,T1,T2>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const SUMMM<CT,T1,T2>& smm, const T x)
{ return SUMMM<CT,T1,T2>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator*(
    const SUMMM<T,T,T>& smm, const CT x)
{ return SUMMM<CT,T,T>(smm.GetX1()*x,smm.GetM1(),smm.GetX2()*x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator*(
    const SUMMM<T,T,T>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T> inline SUMMM<CT,T,T> operator*(
    const SUMMM<T,T,T>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const SUMMM<CT,T1,T2>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator*(
    const SUMMM<CT,T1,T2>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()*CT(x),smm.GetM1(),
      smm.GetX2()*CT(x),smm.GetM2()); 
}

// (x*m+x*m)/x
template <class T, class T1, class T2> inline SUMMM<T,T1,T2> operator/(
    const SUMMM<T,T1,T2>& smm, const T x)
{ return SUMMM<T,T1,T2>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,smm.GetM2()); }

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator/(
    const SUMMM<CT,T1,T2>& smm, const T x)
{ return SUMMM<CT,T1,T2>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator/(
    const SUMMM<T,T,T>& smm, const CT x)
{ return SUMMM<CT,T,T>(smm.GetX1()/x,smm.GetM1(),smm.GetX2()/x,smm.GetM2()); }

template <class T> inline SUMMM<CT,T,T> operator/(
    const SUMMM<T,T,T>& smm, const CCT x)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T> inline SUMMM<CT,T,T> operator/(
    const SUMMM<T,T,T>& smm, const VCT x)
{ 
  return SUMMM<CT,T,T>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator/(
    const SUMMM<CT,T1,T2>& smm, const CCT x)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}

template <class T, class T1, class T2> inline SUMMM<CT,T1,T2> operator/(
    const SUMMM<CT,T1,T2>& smm, const VCT x)
{ 
  return SUMMM<CT,T1,T2>(smm.GetX1()/CT(x),smm.GetM1(),
      smm.GetX2()/CT(x),smm.GetM2()); 
}


