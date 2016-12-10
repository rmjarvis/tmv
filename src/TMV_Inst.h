

// This file is included at the end of each cpp file to instantiate
// the templates for each of the types desired.
// Add more types as desired using the format below:

#ifdef TMV_INST_DOUBLE
#define T double
#include InstFile
#undef T
#endif

#ifdef TMV_INST_FLOAT
#define T float
#include InstFile
#undef T
#endif

#ifdef TMV_INST_SKIP_BLAS
#undef TMV_INST_SKIP_BLAS
#endif

#ifdef TMV_INST_INT
// Define TISINT for any integer type.  e.g. long, short, etc.
#define TISINT 
#define T int
#include InstFile
#undef T
#endif

#ifdef TMV_INST_LONGDOUBLE
#define T long double
#include InstFile
#undef T
#endif

