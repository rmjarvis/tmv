
#ifndef TMV_SafeWrite_H
#define TMV_SafeWrite_H

// This writer structure is both thread safe, and can be turned off
// so you can leave the debug statements in the code without them
// actually compiling into any action.
// (Of course you still nead to beware of things that execute before
// being passed to the writer.  e.g. sw << (a*b) will still calculate a*b
// before ignoring the result.)
// To enable the writing define TMV_USE_WRITER

#ifdef TMV_USE_WRITER
#include <sstream>
#endif

namespace tmv {

    struct SafeWriter
    {
        SafeWriter() {}
        ~SafeWriter()
        {
#ifdef TMV_USE_WRITER
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                std::cout<<s.str();
            }
#endif
        }
        template <class T>
        SafeWriter& operator<<(const T& x)
        {
#ifdef TMV_USE_WRITER
            s << x;
#endif
            return *this; 
        }

        // The next bit is to get dbgcout<<std::endl; working.
        // See: http://stackoverflow.com/questions/1134388/stdendl-is-of-unknown-type-when-overloading-operator

        typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
        typedef CoutType& (*StandardEndLine)(CoutType&);
        SafeWriter& operator<<(StandardEndLine)
        {
#ifdef TMV_USE_WRITER
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                std::cout<<s.str()<<std::endl;
            }
            s.clear();
            s.str(std::string());
#endif
            return *this;
        }

#ifdef TMV_USE_WRITER
        std::stringstream s;
#endif
    };
}


#endif
