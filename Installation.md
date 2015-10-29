﻿#summary Instructions for downloading and installing the TMV library

# Obtaining the source code #

### 1. Download the source code ###

You can get the tarball [here](https://googledrive.com/host/0B6hIz9tCW5iZdEcybFNjRHFmOEE/tmv0.72.tar.gz).

Copy it to the directory where you want to put the TMV library.

### 2. Unpack the tarball ###

On linux machines, the commands are:
```
gunzip tmv0.72.tar.gz
tar xf tmv0.72.tar
```
On Windows, you can use a utility like [IZArc](http://www.izarc.org/), or something similar.

This will make a directory called `tmv0.72` with the subdirectories:
`doc`, `examples`, `include`, `lib`, `src` and `test`
along with the files `README`, `INSTALL` and others
in the top directory.

(If you will be using Microsoft Visual C++, please see
[Visual C++ Instructions](Microsoft_Installation.md).  Otherwise continue on with
these instructions.)

### 3. Install SCons ###

Make sure SCons is installed on your system.
It is available for free from http://www.scons.org/.  (It is a very quick installation if you have Python installed.)

### 4. Install the TMV library ###

Type:
```
scons [Optional flags -- see below]
```

This will make the libraries `libtmv.a` and `libtmv_symband.a`
and put them into the directory `lib`.  Like with `make`, you can add
the flag `-j4` to use 4 (or whatever number of) compilers simultaneously.
Also, the command `scons -h` will print some help information, and `scons -H`
will print information about the options specific to TMV.

### 5. Options ###
There are a number of command-line options that you might need (but try it with no flags
first -- it can often find everything automatically).
The options are listed
with their default value.  You change them simply by specifying a different value
on the command line.  For example:
```
scons CXX=icpc INST_LONGDOUBLE=true
```
If you need to run SCons multiple times (for example to compile the test suite or install the libraries as described below), you only need to specify the new parameter values the first time you run SCons. The program automatically saves your options and continues to use them until you change a value again.

  * `CXX=g++` specifies which C++ compiler to use.
  * `FLAGS=''` specifies the basic flags to pass to the compiler.  The default behavior is to automatically choose good flags to use according to which kind of compiler you are using.
  * `DEBUG=false` specifies whether to keep the debugging assert statements in the compiled library code.
  * `PREFIX=/usr/local` specifies where to install the library when running `scons install` (see below).
  * `INST_FLOAT=true` specifies whether to instantiate the `<float>` templates.
  * `INST_DOUBLE=true` specifies whether to instantiate the `<double>` templates.
  * `INST_LONGDOUBLE=false` specifies whether to instantiate the `<long double>` templates.
  * `INST_INT=true` specifies whether to instantiate the `<int>` templates.
  * `WITH_OPENMP=true` specifies whether to use OpenMP to parallelize some parts of the code.
  * `SHARED=false` specifies whether to make the library files shared as opposed to static libraries.
  * `TEST_FLOAT=true` specifies whether to include the `<float>` tests in the test suite.
  * `TEST_DOUBLE=true` specifies whether to include the `<double>` tests in the test suite.
  * `TEST_LONGDOUBLE=false` specifies whether to include the `<long double>` tests in the test suite.
  * `TEST_INT=true` specifies whether to include the `<int>` tests in the test suite.

The next flags set up the paths that SCons will use to try to find your BLAS and LAPACK libraries.
  * `IMPORT_ENV=true` specifies whether to import the entire environment from the calling shell.  The default is for SCons to use the same environment as the shell from which it is called.  However, sometimes it can be useful to start with a clean environment and manually add paths (see below) for various things, in which case you would want to set this to false.
  * `EXTRA_PATH=''` specifies directories in which to search for executables (notably the compiler, although you can also just give the full path in the `CXX` parameter) in addition to the standard locations such as `/usr/bin` and `/usr/local/bin`. If you are giving multiple directories, they should be separated by colons.
  * `EXTRA_INCLUDE_PATH=''` specifies directories in which to search for header files (such as the BLAS or LAPACK header files) in addition to the standard locations such as `/usr/include` and `/usr/local/include`. These directories are specified as `-I` flags to the compiler. If you are giving multiple directories, they should be separated by colons.
  * `EXTRA_LIB_PATH=''` specifies directories in which to search for libraries (such as the BLAS or LAPACK libraries) in addition to the standard locations such as `/usr/lib` and `/usr/local/lib`.   These directories are specified as `-L` flags to the linker. If you are giving multiple directories, they should be separated by colons.
  * `IMPORT_PATHS=false` specifies whether to import extra path directories from the environment variables:  `PATH`, `C_INCLUDE_PATH`, `LD_LIBRARY_PATH` and `LIBRARY_PATH`. The next options can be used to specify what BLAS and/or LAPACK libraries to use (if any), overriding the default of using whatever libraries SCons chooses from searching through your path and trying to link the libraries that it finds.  The `FORCE` options can be useful if SCons finds a library before trying the one that you want, or if SCons fails in the linking step even though the library should link successfully (I'm still not sure why this happens sometimes), or if you want to compile for a library that requires different linking instructions than the ones that SCons tries. The `FORCE` options will try to test linking with the library requested, but if it fails, then it will just give a warning message.

The next options can be used to specify whether to use BLAS and/or LAPACK libraries.

  * `WITH_BLAS=true` specifies whether to look for and try to use a BLAS library.
  * `WITH_LAPACK=false` specifies whether to look for and try to use a LAPACK library.

See also the
[PDF](http://code.google.com/p/tmv-cpp/downloads/detail?name=TMV_Documentation.pdf) documentation for still more options.)

When SCons starts up, it will look through the standard paths, along with any extra paths you have specified with the above options, to find BLAS and LAPACK libraries.  This can sometimes require a few iterations to get working correctly.   You should look at the initial output from SCons to make sure it finds the correct BLAS and LAPACK libraries that you think it should find.  Here is a sample output:

```
$ scons
scons: Reading SConscript files ...
Using the default scons options
Using compiler: /usr/bin/g++
Detected clang++ masquerading as g++
compiler version: 5.1
Determined that there are 4 cpus, so use this many jobs.
You can override this behavior with scons -jN
Debugging turned off
Checking for MKL... no
Checking for ACML... no
Checking for GotoBLAS... no
Checking for CBLAS... no
Checking for ATLAS... no
Checking for Fortran BLAS... yes
Using Fortran BLAS
No OpenMP support for clang++
TMV Version  0.72

scons: done reading SConscript files.
scons: Building targets ...
```

If a "`Checking for`..." line ends with `no`, even though you think that library is installed on your computer, then it probably means that you need to tell SCons which directories to search, in addition to the standard locations.  The most straightforward way to do this is with the parameters `EXTRA_INCLUDE_PATH` and `EXTRA_LIB_PATH`. These are described in detail above.  See also `IMPORT_ENV` and `IMPORT_PATHS`.

### 5. Install the library and header files ###
Type
```
scons install
```
(or possibly `sudo scons install` if you are installing into `/usr/local` or somewhere similar).

This will install the necessary header files into the directory `/usr/local/include` and the libraries into `/usr/local/lib`.  As mentioned above, you can also specify a different prefix with the command line option `PREFIX=`_install-dir_.  A common choice for users without `sudo` privileges is `PREFIX=~` which will install the library in `~/include` and `~/lib`.

At the end of the installation process, you should see a message similar to:
```
The TMV library was successfully installed.  
To link your code against the TMV library, you should use the 
link flags: 

-ltmv -llapack -lblas -lpthread -fopenmp

Or if you are using Band, Sym or SymBand matrices, use: 

-ltmv_symband -ltmv -llapack -lblas -lpthread -fopenmp

These flags (except for the optional -ltmv_symband) have been
saved in the file:

/usr/local/share/tmv-link

so you can automatically use the correct flags in a makefile
(for example) by using lines such as:

TMVLINK := $(shell cat /usr/local/share/tmv-link)
LIBS = $(TMVLINK) [... other libs ...]


scons: done building targets.
```
These instructions tell you what BLAS and LAPACK libraries were found on your system and  are needed for proper linkage of your code.  The linkage flags are stored in a file on your computer so that you can automatically get the linkage correct according to what options you decide to use when installing TMV.  This is especially useful on systems where a system administrator installs the library, which is then used by many users.

You can have your makefiles read this file as described above.  Or if you are using SCons to build your software, you can similarly add the contents of the file into `env['LIBS']`. The `examples` directory has a makefile that uses the above lines (using the local `share` directory rather than the installed location).

### 6. Build the test suite (Optional) ###
Type:
```
scons test
```

This will make three executables called `tmvtest1`, `tmvtest2` and `tmvtest3` in the `test` directory.

Then you should run the three test suites. They should output a bunch of lines reading _Something_` passed all tests`. If one of them ends in a line that starts with `Error`,  then please post a bug report at http://code.google.com/p/tmv-cpp/issues about the problem including what compiler you are using, some details about your system, and what (if any) BLAS and LAPACK libraries you are linking to.

If you specify `SMALL_TESTS=true`, then the smaller test executables `tmvtest1a-d`, `tmvtest2a-c`, and `tmvtest3a-e` (where `a-d` means  four files with each of `a`, `b`, `c` and `d`) will be made instead. These perform the same tests as the larger test executables, but can be easier for some linkers.