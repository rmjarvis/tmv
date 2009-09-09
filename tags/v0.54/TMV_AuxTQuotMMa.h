// x*(x*m/m)
template <class T, class T1, class T2> inline TQUOTMM<T,T1,T2> operator*(
    T x, const TQUOTMM<T,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    CT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    CCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    VCT x, const TQUOTMM<CT,T,T>& qmm)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    T x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    CCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    VCT x, const TQUOTMM<CT,T1,T2>& qmm)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m/m)*x
template <class T, class T1, class T2> inline TQUOTMM<T,T1,T2> operator*(
    const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator*(
    const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator*(
    const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m/m)/x
template <class T, class T1, class T2> inline TQUOTMM<T,T1,T2> operator/(
    const TQUOTMM<T,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator/(
    const TQUOTMM<CT,T,T>& qmm, CT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator/(
    const TQUOTMM<CT,T,T>& qmm, CCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TQUOTMM<CT,T,T> operator/(
    const TQUOTMM<CT,T,T>& qmm, VCT x)
{ return TQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator/(
    const TQUOTMM<CT,T1,T2>& qmm, T x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator/(
    const TQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TQUOTMM<CT,T1,T2> operator/(
    const TQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

// x*(x*m%m)
template <class T, class T1, class T2> inline TRQUOTMM<T,T1,T2> operator*(
    T x, const TRQUOTMM<T,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    CT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    CCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    VCT x, const TRQUOTMM<CT,T,T>& qmm)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    T x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    CCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    VCT x, const TRQUOTMM<CT,T1,T2>& qmm)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> inline TRQUOTMM<T,T1,T2> operator*(
    const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator*(
    const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator*(
    const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetP(),qmm.GetM2()); }

// (x*m%m)/x
template <class T, class T1, class T2> inline TRQUOTMM<T,T1,T2> operator/(
    const TRQUOTMM<T,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator/(
    const TRQUOTMM<CT,T,T>& qmm, CT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator/(
    const TRQUOTMM<CT,T,T>& qmm, CCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T> inline TRQUOTMM<CT,T,T> operator/(
    const TRQUOTMM<CT,T,T>& qmm, VCT x)
{ return TRQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator/(
    const TRQUOTMM<CT,T1,T2>& qmm, T x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator/(
    const TRQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline TRQUOTMM<CT,T1,T2> operator/(
    const TRQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return TRQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetP(),qmm.GetM2()); }

