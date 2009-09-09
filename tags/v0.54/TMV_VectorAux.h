#ifndef TMV_VECTORAUX_H
#define TMV_VECTORAUX_H

#include <vector>

template <class VI> inline void ConvertIndexToPermute(
          size_t n, const VI& newindex, size_t* P)
{
  // newindex[i]=j means value at original j location needs to go to i.
  std::vector<size_t> currindex(n);
  std::vector<size_t> origindex(n);
  for(size_t i=0;i<n;++i) {
    currindex[i] = i;
    origindex[i] = i;
  } 
  // currindex[i]=j means value at original i location is currently at j.
  // origindex[j]=i means value at original i location is currently at j.
  for(size_t i=0;i<n;++i) {
    size_t ip = currindex[newindex[i]];
    P[i] = ip;
    if (i != ip) { 
      size_t origi = origindex[i];
      size_t origip = origindex[ip];
      currindex[origi] = ip;
      currindex[origip] = i;
      origindex[i] = origip;
      origindex[ip] = origi;
    } 
  } 
}   

#endif
