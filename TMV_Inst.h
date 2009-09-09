// This file is included at the end of each cpp file to instantiate
// the templates for each of the types desired.
// Add more types as desired using the format below:

#define T double
#include InstFile
#undef T

#ifdef INST_FLOAT
#define T float
#include InstFile
#undef T
#endif

#ifdef INST_INT
#define T int
#include InstFile
#undef T
#endif

#ifdef INST_LONGDOUBLE
#define T long double
#include InstFile
#undef T
#endif

