///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_IOStyle_H
#define TMV_IOStyle_H

namespace tmv {

    class IOStyle
    {
    public:
        IOStyle()
        { setToDefault(); }

        // Use default copy, destructor, op=

        // Handlers for setting features:
        IOStyle& noPrefix()
        { usecode = false; writesize = false; return *this; }

        IOStyle& useCode()
        { usecode = true; return *this; }

        IOStyle& noCode()
        { usecode = false; return *this; }

        IOStyle& noSize()
        { writesize = false; return *this; }

        IOStyle& simpleSize()
        { writesize = true; simplesize = true; return *this; }

        IOStyle& fullSize()
        { writesize = true; simplesize = false; return *this; }

        IOStyle& markup(
            const std::string& s, const std::string& lp,
            const std::string& sp, const std::string& rp,
            const std::string& re, const std::string& f)
        {
            start = s; lparen = lp; space = sp;
            rparen = rp; rowend = re; final = f;
            return *this;
        }

        IOStyle& fullMatrix()
        { usecompact = false; return *this; }

        IOStyle& compact()
        { usecompact = true; return *this; }

        IOStyle& setThresh(double t)
        { thresh = t; return *this; }

        IOStyle& setPrecision(int p)
        { prec = p; return *this; }

        IOStyle& useDefaultPrecision()
        { prec = -1; return *this; }

        // Revert to the default values.
        IOStyle& setToDefault()
        { *this = getDefaultSingleton(); return *this; }

        // Declare current state to be default IOStyle from now on.
        void makeDefault()
        { getDefaultSingleton() = *this; }

        // Revert the default IO to the original default style.
        static void revertDefault()
        { getDefaultSingleton() = IOStyle(0); }


    private :

        bool usecode;
        bool writesize;
        bool simplesize;
        bool usecompact;
        std::string start;
        std::string lparen;
        std::string space;
        std::string rparen;
        std::string rowend;
        std::string final;
        double thresh;
        int prec; // -1 = don't change precision.

        void write(std::ostream& os)
        {
            os << usecode << " " << writesize << " "
                << simplesize << " " << usecompact << " '"
                << start << "' '" << lparen << "' '" << space << "' '"
                << rparen << "' '" << rowend << "' '" << final << "' "
                << thresh << " " << prec;
        }

        // Helper for dealing with threshold writing.
        template <typename T>
        T outVal(const T& val) const
        { return (thresh > 0. && TMV_ABS(val) < thresh) ? T(0) : val; }

        template <typename T>
        std::complex<T> outVal(const std::complex<T>& val) const
        {
            return thresh > 0. ?
                std::complex<T>(outVal(real(val)),outVal(imag(val))) : val;
        }

        // Private constructor with initial default values.
        // (The int is just to make it easy to resolve on the signature.)
        IOStyle(int) :
            usecode(false), writesize(true), simplesize(true),
            usecompact(false),
            start("\n"), lparen("( "), space("  "),
            rparen(" )"), rowend("\n"), final("\n"),
            thresh(0.), prec(-1) {}

        // Use a singleton idiom for the default IOStyle:
        static inline IOStyle& getDefaultSingleton()
        {
            static IOStyle def(0);
            return def;
        }

        // All actual usage of this class is mediated through a
        // Writer or Reader.
        friend class TMV_Writer;
        friend class TMV_Reader;
    };

    // Some IOStyles that might be typically useful:
    inline IOStyle NormalIO()
    {
        return IOStyle().noCode().simpleSize().fullMatrix().
            markup("\n","( ","  "," )","\n","\n");
    }

    inline IOStyle CompactIO()
    {
        return IOStyle().useCode().fullSize().compact().
            markup("",""," ",""," ","");
    }

    inline IOStyle ThreshIO(double thresh)
    { return IOStyle().setThresh(thresh); }

    inline IOStyle PrecIO(int prec)
    { return IOStyle().setPrecision(prec); }

    inline IOStyle EigenIO()
    { return IOStyle().noPrefix().fullMatrix().markup("",""," ","","\n",""); }

    class TMV_Writer
    {
    public :
        TMV_Writer(std::ostream& _os, const IOStyle& _s) : os(_os), s(_s) {}
        // Use default copy, op=, destr

        void begin() const
        {
            if (s.prec >= 0) {
                oldprec = os.precision(s.prec);
            }
        }

        void end() const
        {
            if (s.prec >= 0) {
                os.precision(oldprec);
            }
        }

        void writeCode(const std::string& code) const
        { if (s.usecode) os << code << s.space; }

        void writeSize(ptrdiff_t n) const
        { if (s.writesize) os << n << s.space; }
        void writeSimpleSize(ptrdiff_t n) const
        { if (s.simplesize) writeSize(n); }
        void writeFullSize(ptrdiff_t n) const
        { if (!s.simplesize) writeSize(n); }

        void writeStart() const
        { os << s.start; }
        void writeLParen() const
        { os << s.lparen; }
        void writeSpace() const
        { os << s.space; }
        void writeRParen() const
        { os << s.rparen; }
        void writeRowEnd() const
        { os << s.rowend; }
        void writeFinal() const
        { os << s.final; }

        template <typename T>
        void writeValue(const T& x) const
        { os << Value(s.outVal(x)); }

        bool isCompact() const
        { return s.usecompact; }

        std::ostream& getos() const { return os; }
        const IOStyle& getstyle() const { return s; }

    private :
        std::ostream& os;
        IOStyle s;

        mutable std::streamsize oldprec;

        // This bit is to workaround a bug in pgCC that was fixed in version 7.
        // I don't know if versions earlier than 6.1 had the bug, but
        // I apply the workaround to all version before 7.
        template <typename T>
        static inline T Value(const T& x) { return x; }
#if defined(__PGI) && (!defined(__PGIC__) || __PGIC__ < 7)
        static inline double Value(const long double& x)
        { return double(x); }
        static inline std::complex<double> Value(
            const std::complex<long double>& x)
        { return std::complex<double>(x); }
#endif
    };

    inline TMV_Writer operator<<(std::ostream& os, const IOStyle& style)
    { return TMV_Writer(os,style); }

    class TMV_Reader
    {
    public :
        TMV_Reader(std::istream& _is, const IOStyle& _s) : is(_is), s(_s) {}
        // Use default copy, op=, destr

        bool readStr(
            const std::string& str, std::string& exp, std::string& got) const
        {
            if (str.size() == 0) return true;
            else {
                skipWhiteSpace();
                std::string getstr(str.size(),' ');
                for(size_t i=0;i<str.size();++i) is.get(getstr[i]);
                if (getstr != str) {
                    exp = str;
                    got = getstr;
                    return false;
                }
                if (!is) return false;
                else return true;
            }
        }

        bool readCode(
            const std::string& code, std::string& exp, std::string& got) const
        {
            if (s.usecode) {
                if (readStr(trim(code),exp,got)) {
                    return readSpace(exp,got);
                } else {
                    return false;
                }
            } else {
                return true;
            }
        }

        // For real SymMatrix, there are two valid codes.
        bool readCode(
            const std::string& code1, const std::string& code2,
            std::string& exp, std::string& got) const
        {
            if (s.usecode) {
                if (readStr(trim(code1),exp,got)) {
                    return readSpace(exp,got);
                } else if (got == code2) {
                    exp = got = "";
                    return readSpace(exp,got);
                } else {
                    return false;
                }
            } else {
                return true;
            }
        }

        bool readStart(std::string& exp, std::string& got) const
        { return readStr(trim(s.start),exp,got); }

        bool readLParen(std::string& exp, std::string& got) const
        { return readStr(trim(s.lparen),exp,got); }

        bool readSpace(std::string& exp, std::string& got) const
        { return readStr(trim(s.space),exp,got); }

        bool readRParen(std::string& exp, std::string& got) const
        { return readStr(trim(s.rparen),exp,got); }

        bool readRowEnd(std::string& exp, std::string& got) const
        { return readStr(trim(s.rowend),exp,got); }

        bool readFinal(std::string& exp, std::string& got) const
        { return readStr(trim(s.final),exp,got); }

        bool readSize(ptrdiff_t& n, std::string& exp, std::string& got) const
        {
            if (s.writesize) {
                skipWhiteSpace();
                is >> n;
                if (!is) return false;
                else return readSpace(exp,got);
            } else {
                return true;
            }
        }

        bool readSimpleSize(ptrdiff_t& n, std::string& exp, std::string& got) const
        { return s.simplesize ? readSize(n,exp,got) : true; }

        bool readFullSize(ptrdiff_t& n, std::string& exp, std::string& got) const
        { return !s.simplesize ? readSize(n,exp,got) : true; }

        template <typename T>
        bool readValue(T& x) const
        {
            skipWhiteSpace();
            is >> x;
            if (!is) return false;
            else return true;
        }

        bool isCompact() const
        { return s.usecompact; }

        std::istream& getis() const { return is; }
        const IOStyle& getstyle() const { return s; }

    private :
        std::istream& is;
        IOStyle s;

        void skipWhiteSpace() const
        {
            static std::string whitespace = " \n\t\v\r\f";
            char c;
            do { is.get(c); } while (whitespace.find(c) != std::string::npos);
            is.unget();
        }

        static std::string trim(std::string s)
        {
            static std::string whitespace = " \n\t\v\r\f";
            size_t i1 = s.find_first_not_of(whitespace);
            if (i1 == std::string::npos) return "";
            else {
                size_t i2 = s.find_last_not_of(whitespace);
                return std::string(s,i1,i2-i1+1);
            }
        }
    };

    inline TMV_Reader operator>>(std::istream& is, const IOStyle& style)
    { return TMV_Reader(is,style); }

} // namespace tmv

#endif
