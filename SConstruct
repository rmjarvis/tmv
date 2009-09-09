# vim: set filetype=python et ts=4 sw=4:
# to do:  
#   Add more support for other compilers
#   will support g++ and icpc at least
# 
# Always run scons from the root directory of tmv

import os
import sys


# Subdirectories containing SConscript files.  We always process src but
# there are some other optional ones
src_dir = 'src'
#subdirs=['test','examples','doc']
subdirs=['test','examples']

# Configurations will be saved here so command line options don't
# have to be sent more than once
config_file = 'tmv_scons.conf'

# Default directory for installation.  
# This is the only UNIX specific things I am aware
# of in the script.  On the other hand, these are not required for the
# script to work since prefix can be set on the command line and the
# extra paths are not needed, but I wish I knew how to get the default 
# prefix for the system so I didn't have to set this.

default_prefix = '/usr/local'

# Now set up the environment
initial_env = Environment()

# first check for a saved conf file
opts = Variables(config_file)

# Now set up options for the command line
opts.Add('CXX','Name of c++ compiler')
opts.Add('FLAGS','Compile flags to send to the compiler','')
opts.Add(BoolVariable('DEBUG','Turn on debugging statements',False))
opts.Add(PathVariable('PREFIX',
            'prefix for installation','', PathVariable.PathAccept))

opts.Add(EnumVariable('OPT',
            'Set the optimization level for TMV library', '3',
            allowed_values=('0','1','2','3')))
opts.Add(BoolVariable('INST_FLOAT',
            'Instantiate <float> templates in compiled library', True))
opts.Add(BoolVariable('INST_DOUBLE',
            'Instantiate <double> templates in compiled library', True))
opts.Add(BoolVariable('INST_LONGDOUBLE',
            'Instantiate <long double> templates in compiled library', False))
opts.Add(BoolVariable('INST_INT',
            'Instantiate <int> templates in compiled library', False))
opts.Add(BoolVariable('INST_COMPLEX',
            'Instantiate complex<T> templates in compiled library', True))
opts.Add(BoolVariable('INST_MIX',
            'Instantiate functions that mix real with complex', True))

opts.Add(EnumVariable('TEST_OPT',
            'Set the optimization level for TMV test suite', '0',
            allowed_values=('0','1','2','3')))
opts.Add(BoolVariable('TEST_FLOAT',
            'Instantiate <float> in the test suite', True))
opts.Add(BoolVariable('TEST_DOUBLE',
            'Instantiate <double> in the test suite', True))
opts.Add(BoolVariable('TEST_LONGDOUBLE',
            'Instantiate <long double> in the test suite', False))
opts.Add(BoolVariable('TEST_INT',
            'Instantiate <int> in the test suite', True))

opts.Add(PathVariable('EXTRA_PATH',
            'Extra paths for executables (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(PathVariable('EXTRA_LIB_PATH',
            'Extra paths for linking (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(PathVariable('EXTRA_INCLUDE_PATH',
            'Extra paths for header files (separated by : if more than 1)',
            '',PathVariable.PathAccept))
opts.Add(BoolVariable('IMPORT_PATHS',
            'Import PATH, C_INCLUDE_PATH and LIBRARY_PATH/LD_LIBRARY_PATH environment variables',
            False))
opts.Add(BoolVariable('IMPORT_ENV',
            'Import full environment from calling shell', False))

opts.Add(BoolVariable('WITH_BLAS',
            'Look for blas libraries and link if found.', True))
opts.Add(BoolVariable('WITH_LAPACK',
            'Look for lapack libraries and link if found.', True))
opts.Add(BoolVariable('FORCE_MKL',
            'Force scons to use MKL for BLAS and/or LAPACK', False))
opts.Add(BoolVariable('FORCE_ACML',
            'Force scons to use ACML for BLAS and/or LAPACK', False))
opts.Add(BoolVariable('FORCE_GOTO',
            'Force scons to use GOTO BLAS', False))
opts.Add(BoolVariable('FORCE_ATLAS',
            'Force scons to use ATLAS for BLAS', False))
opts.Add(BoolVariable('FORCE_CBLAS',
            'Force scons to use CBLAS', False))
opts.Add(BoolVariable('FORCE_FBLAS',
            'Force scons to use Fortran BLAS', False))
opts.Add(BoolVariable('FORCE_CLAPACK',
            'Force scons to use CLAPACK', False))
opts.Add(BoolVariable('FORCE_ATLAS_LAPACK',
            'Force scons to use ATLAS subset of LAPACK', False))
opts.Add(BoolVariable('FORCE_FLAPACK',
            'Force scons to use Fortran LAPACK', False))
opts.Add('LIBS','Libraries to send to the linker','')

opts.Add(BoolVariable('SMALL_LIB',
            'Avoid optimizations that cause the library to become very large', 
            True))
opts.Add(BoolVariable('WITH_OPENMP',
            'Look for openmp and use if found.', False))
opts.Add(BoolVariable('STATIC',
            'Use static linkage', False))
opts.Add(BoolVariable('XTEST',
            'Do (a LOT of) extra tests in the test suite ', False))
opts.Add(BoolVariable('MEM_TEST',
            'Test for memory leaks', False))
opts.Add(BoolVariable('SMALL_TESTS',
            'Make the small test programs: tmvtest1a, tmvtest1b, etc.', False))
opts.Add(BoolVariable('WARN',
            'Add warning compiler flags, like -Wall', False))
opts.Add(BoolVariable('NOMIX_SMALL',
            'Do not test the mixed Small and regular arithmetic', False))

opts.Update(initial_env)
opts.Save(config_file,initial_env)
Help(opts.GenerateHelpText(initial_env))

# This helps us determine of openmp is available
openmp_mingcc_vers = 4.1
openmp_minicpc_vers = 9.0
openmp_minpgcc_vers = 6.0


def BasicCCFlags(env):
    """
    """
    if env['FLAGS'] == '':
        compiler = env['CXXTYPE']
        version = env['CXXVERSION_NUMERICAL']
    
        if compiler == 'g++':
            env.Replace(CCFLAGS=['-O3'])
            env['TEST_FLAGS'] = ['-O']
            if version <= 4.2:
                env.Append(CCFLAGS=['-fno-strict-aliasing'])
            if env['WARN']:
                env.Append(CCFLAGS=['-ansi','-pedantic-errors','-Wall','-Werror'])
                env['TEST_FLAGS'] = ['-ansi','-pedantic-errors','-Wall','-Werror']
    
        elif compiler == 'icpc':
            env.Replace(CCFLAGS=['-O2'])
            env['TEST_FLAGS'] = []
            if version >= 10:
                env.Append(CCFLAGS=['-vec-report0'])
                env['TEST_FLAGS'] += ['-vec-report0']
            if env['WARN']:
                env.Append(CCFLAGS=['-Wall','-Werror','-wd279,383,810,981'])
                env['TEST_FLAGS'] += ['-Wall','-Werror','-wd279,383,810,981']
                if version >= 9:
                    env.Append(CCFLAGS=['-wd1572'])
                    env['TEST_FLAGS'] += ['-wd1572']
                if version >= 10 and env['WITH_BLAS']:
                    # These warning only show up in mkl.h:
                    env.Append(CCFLAGS=['-wd424,193'])
                if version >= 11:
                    env.Append(CCFLAGS=['-wd2259'])
                    env['TEST_FLAGS'] += ['-wd2259']
            else :
                env.Append(CCFLAGS=['-w'])
                env['TEST_FLAGS'] += ['-w']

        elif compiler == 'pgCC':
            env.Replace(CCFLAGS=['-O2','-fast','-Mcache_align'])
            env['TEST_FLAGS'] = ['-O0']

        elif compiler == 'cl':
            env.Replace(CCFLAGS=['/EHsc','/nologo','/O2','/Oi'])
            env['TEST_FLAGS'] = ['/EHsc','/nologo']
            if env['WARN']:
                env.Append(CCFLAGS=['/W2','/WX'])
                env['TEST_FLAGS'] += ['/W1','/WX']

        else:
            print 'Warning: Unknown compiler.  You should set FLAGS directly.'
            env.Replace(CCFLAGS=[])
            env['TEST_FLAGS'] = []

    else :
        # If flags are specified as an option use them:
        cxx_flags = env['FLAGS'].split(' ')
        env.Replace(CCFLAGS=cxx_flags)
        env['TEST_FLAGS'] = cxx_flags

    # Also parse the LIBS options if present
    if env['LIBS'] == '':
        env.Replace(LIBS=[])
    else:
        libs = env['LIBS'].split(' ')
        env.Replace(LIBS=libs)


def AddOpenMPFlag(env):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    g++ uses -fopemnp
    icpc uses -openmp
    pgCC uses -mp
    
    Other compilers?
    """
    compiler = env['CXXTYPE']
    version = env['CXXVERSION_NUMERICAL']
    if compiler == 'g++':
        if version < openmp_mingcc_vers: 
            print 'No OpenMP support for g++ versions before ',openmp_mingcc_vers
            return
        flag = ['-fopenmp']
        ldflag = ['-fopenmp']
        xlib = ['pthread']
    elif compiler == 'icpc':
        if version < openmp_minicpc_vers:
            print 'No OpenMP support for icpc versions before ',openmp_minicpc_vers
            return
        flag = ['-openmp']
        ldflag = ['-openmp']
        xlib = ['pthread']
    elif compiler == 'pgCC':
        if version < openmp_minpgcc_vers:
            print 'No OpenMP support for pgCC versions before ',openmp_minpgcc_vers
            return
        flag = ['-mp','--exceptions']
        ldflag = ['-mp']
        xlib = ['pthread']
    elif compiler == 'cl':
        flag = ['/openmp']
        ldflag = ['/openmp']
        xlib = []
    else:
        print 'Warning: No OpenMP support for compiler ',compiler

    #print 'Adding openmp support:',flag
    print 'Using OpenMP'
    env['OMP_FLAGS'] = flag 
    env.Append(LINKFLAGS=ldflag)
    env.Append(LIBS=xlib)

def GetCompilerVersion(env):
    """
    """
    compiler = env['CXX']

    # Get the compiler type without suffix or path.  
    # e.g. /sw/bin/g++-4 -> g++
    if 'icpc' in compiler :
        compilertype = 'icpc'
        versionflag = '--version'
        linenum=0
    elif 'pgCC' in compiler :
        compilertype = 'pgCC'
        versionflag = '--version'
        linenum=1
        # pgCC puts the version number on the second line of output.
    elif 'g++' in compiler :
        compilertype = 'g++'
        versionflag = '--version'
        linenum=0
    elif 'cl' in compiler :
        compilertype = 'cl'
        versionflag = ''
        linenum=0
        # With cl, the version seems to be printed in the first line,
        # but the lines read out with popen below seem to skip the
        # first two lines.  So the code below ends up with version = 0.
        # It doesn't really matter though, since we don't use the cl 
        # version for anything. 
    else :
        compilertype = 'unknown'
        version = 0
        vnum = 0

    if compilertype != 'unknown':
        cmd = compiler + ' ' + versionflag
        lines = os.popen(cmd).readlines()
        line = lines[linenum]
    
        import re
        match = re.search(r'[0-9]+(\.[0-9]+)+', line)
    
        if match:
            version = match.group(0)
            # Get the version up to the first decimal
            # e.g. for 4.3.1 we only keep 4.3
            vnum = version[0:version.find('.')+2]
        else:
            version = 0
            vnum = 0

    print '\nUsing compiler:',compiler
    print 'compiler version:',version

    env['CXXTYPE'] = compilertype
    env['CXXVERSION'] = version
    env['CXXVERSION_NUMERICAL'] = float(vnum)


def AddPath(pathlist, newpath):
    """
    Add path(s) to a list of paths.  Check the path exists and that it is
    not already in the list
    """
    if type(newpath) == list:
        for l in newpath:
            AddPath(pathlist, l)
    else:
        # to deal with expansions and possible end / which 
        # messes up uniqueness test
        p = os.path.abspath(newpath) 
        if os.path.exists(p):
            if pathlist.count(p) == 0:
                pathlist.append(p)

def AddExtraPaths(env):
    """
    Add some include and library paths.
    Also merge in $PATH, $C_INCLUDE_PATH and $LIBRARY_PATH/$LD_LIBRARY_PATH 
    environment variables if requested.
    
    The set itself is created in order of appearance here, but then this 
    whole set is prepended.  The order within this list is:

        local lib and include paths
        paths in EXTRA_*PATH parameters
        paths from the user's environment
        paths in PREFIX directory

    Only paths that actually exists are kept.
    """
    # local includes and lib paths
    # The # symbol means to interpret these from the top-level scons
    # directory even when we are in a sub-directory (src,test,etc.)
    bin_paths = []
    cpp_paths = ['#include']
    lib_paths = ['#lib']

    # Paths specified in EXTRA_*
    bin_paths += env['EXTRA_PATH'].split(':')
    lib_paths += env['EXTRA_LIB_PATH'].split(':')
    cpp_paths += env['EXTRA_INCLUDE_PATH'].split(':')

    # Paths found in environment paths
    if env['IMPORT_PATHS'] and os.environ.has_key('PATH'):
        paths=os.environ['PATH']
        paths=paths.split(os.pathsep)
        AddPath(bin_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('C_INCLUDE_PATH'):
        paths=os.environ['C_INCLUDE_PATH']
        paths=paths.split(os.pathsep)
        AddPath(cpp_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('LIBRARY_PATH'):
        paths=os.environ['LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths, paths)

    if env['IMPORT_PATHS'] and os.environ.has_key('LD_LIBRARY_PATH'):
        paths=os.environ['LD_LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths, paths)

    # PREFIX directory
    # If none given, then don't add them to the -L and -I directories.
    # But still use the default /usr/local for installation
    if env['PREFIX'] == '':
        env['INSTALL_PREFIX'] = default_prefix
    else:
        AddPath(bin_paths, os.path.join(env['PREFIX'], 'bin'))
        AddPath(lib_paths, os.path.join(env['PREFIX'], 'lib'))
        AddPath(cpp_paths, os.path.join(env['PREFIX'], 'include'))
        env['INSTALL_PREFIX'] = env['PREFIX']
    

    #print 'bin paths = ',bin_paths
    #print 'cpp paths = ',cpp_paths
    #print 'lib paths = ',lib_paths

    #env.AppendENVPath('PATH', bin_paths)
    #env.Append(LIBPATH= lib_paths)
    #env.Append(CPPPATH= cpp_paths)
    env.PrependENVPath('PATH', bin_paths)
    env.Prepend(LIBPATH= lib_paths)
    env.Prepend(CPPPATH= cpp_paths)

def ReadFileList(fname):
    """
    This reads a list of whitespace separated values from the input file fname
    and stores it as a list.  We will make this part of the environment so
    other SConscripts can use it
    """
    try:
        files=open(fname).read().split()
    except:
        print 'Could not open file:',fname
        sys.exit(45)
    files = [f.strip() for f in files]
    return files


def CheckLibs(context,try_libs,source_file):
    init_libs = context.env['LIBS']
    context.env.Prepend(LIBS=try_libs)
    result = context.TryLink(source_file,'.cpp')
    if not result :
        context.env.Replace(LIBS=init_libs)
    return result
      

def CheckMKL(context):
    mkl_source_file = """
#include "mkl.h"
int main()
{
  char ta='N', tb='N';
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  dgemm(&ta,&tb,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
  return 0;
}
"""

    context.Message('Checking for MKL... ')
    if context.env['CXXTYPE'] == 'icpc':
        threadlib = 'mkl_intel_thread'
    elif context.env['CXXTYPE'] == 'pgCC':
        threadlib = 'mkl_pgi_thread'
    else:
        threadlib = 'mkl_gnu_thread'

    if context.TryCompile(mkl_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],mkl_source_file) or
            CheckLibs(context,['mkl'],mkl_source_file) or
            CheckLibs(context,['mkl','pthread'],mkl_source_file) or
            CheckLibs(context,['mkl','guide','pthread'],mkl_source_file) or
            CheckLibs(context,['mkl_ia32','guide','pthread'],
                        mkl_source_file) or
            CheckLibs(context,['mkl_ia32','mkl_core','mkl_sequential'],
                        mkl_source_file) or
            CheckLibs(context,['mkl_intel_lp64','mkl_core',threadlib,
                        'guide','pthread'],mkl_source_file) or
            CheckLibs(context,['mkl_intel_lp64','mkl_core','mkl_sequential'],
                        mkl_source_file) or
            CheckLibs(context,['mkl_ipf','guide','pthread'],
                        mkl_source_file) or
            CheckLibs(context,['mkl_em64t','guide','pthread'],
                        mkl_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_MKL']:
            print 'Warning: Forced use of MKL even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_MKL']:
            print 'Error: FORCE_MKL, but failed to find or compile with mkl.h'
            Exit(1)
            
    return result

def CheckACML(context):
    acml_source_file = """
#include "acml.h"
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0;
  dbdsdc(uplo,compq,n,d,e,u,ldu,v,ldv,q,iq,info);
  return 0;
}
"""
    context.Message('Checking for ACML... ')

    if context.TryCompile(acml_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],acml_source_file) or
            CheckLibs(context,['acml'],acml_source_file) or
            CheckLibs(context,['acml','pgftnrtl'],acml_source_file) or
            CheckLibs(context,['acml','gfortran'],acml_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_ACML']:
            print 'Warning: Forced use of ACML even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ACML']:
            print 'Error: FORCE_ACML, but failed to find or compile with acml.h'
            Exit(1)
            
    return result


def CheckGOTO(context):
    goto_source_file = """
extern "C" {
#include "util/fblas.h"
}
int main()
{
  char ta='N', tb='N';
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  dgemm_(ta,tb,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1,1);
  return 0;
}
"""

    context.Message('Checking for GotoBLAS... ')

    if context.TryCompile(goto_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],goto_source_file) or
            CheckLibs(context,['goto'],goto_source_file))

        context.Result(result)

        if not result and context.env['FORCE_GOTO']:
            print 'Warning: Forced use of GOTO even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_GOTO']:
            print 'Error: FORCE_GOTO, but failed compile test'
            Exit(1)

    return result



def CheckATLAS(context):
    atlas_source_file = """
extern "C" {
#include "cblas.h"
}
int main()
{
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
      M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
  return 0;
}
"""

    context.Message('Checking for ATLAS... ')

    if context.TryCompile(atlas_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],atlas_source_file) or
            CheckLibs(context,['cblas','atlas'],atlas_source_file))

        context.Result(result)

        if not result and context.env['FORCE_ATLAS']:
            print 'Warning: Forced use of ATLAS even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ATLAS']:
            print 'Error: FORCE_ATLAS, but failed to find or compile with cblas.h'
            Exit(1)
            
    return result




def CheckCBLAS(context):
    cblas_source_file = """
extern "C" {
#include "cblas.h"
}
int main()
{
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
  return 0;
}
"""

    context.Message('Checking for CBLAS... ')

    if context.TryCompile(cblas_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],cblas_source_file) or
            CheckLibs(context,['cblas'],cblas_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_CBLAS']:
            print 'Warning: Forced use of CBLAS even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_CBLAS']:
            print 'Error: FORCE_CBLAS, but failed to find or compile with cblas.h'
            Exit(1)
            
    return result



def CheckFBLAS(context):
    fblas_source_file = """
extern "C" {
#include "util/fblas.h"
}
int main()
{
  char ta='N', tb='N';
  int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
  double alpha=1.,beta=1., *A=0, *B=0, *C=0;
  dgemm_(ta,tb,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1,1);
  return 0;
}
"""

    context.Message('Checking for Fotran BLAS... ')

    if context.TryCompile(fblas_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],fblas_source_file) or
            CheckLibs(context,['blas'],fblas_source_file) or
            CheckLibs(context,['blas','pgftnrtl'],fblas_source_file) or
            CheckLibs(context,['blas','gfortran'],fblas_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_FBLAS']:
            print 'Warning: Forced use of FBLAS even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_FBLAS']:
            print 'Error: FORCE_FBLAS, but failed compile test'
            Exit(1)
            
    return result




def CheckMKL_LAP(context):
    mkl_lap_source_file = """
#include "mkl.h"
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc(&uplo,&compq,&n,d,e,u,&ldu,v,&ldv,q,iq,work,iwork,&info);
  return 0;
}
"""

    context.Message('Checking for MKL LAPACK... ')

    result = (context.TryCompile(mkl_lap_source_file,'.cpp') and
        (CheckLibs(context,[],mkl_lap_source_file) or
         CheckLibs(context,['mkl_lapack'],mkl_lap_source_file) or
         CheckLibs(context,['mkl_lapack','guide'],mkl_lap_source_file)))

    context.Result(result)
    return result



def CheckATLAS_LAP(context):
    atlas_lapack_source_file = """
extern "C" {
#include "clapack.h"
}
int main()
{
  int N=1,lda=1;
  double* A=0;
  clapack_dlauum(CblasRowMajor,CblasUpper,N,A,lda);
  return 0;
}
"""

    context.Message('Checking for ATLAS version of LAPACK... ')

    if context.TryCompile(atlas_lapack_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],atlas_lapack_source_file) or
            CheckLibs(context,['lapack'],atlas_lapack_source_file) or 
            CheckLibs(context,['clapack'],atlas_lapack_source_file) or
            CheckLibs(context,['lapack_atlas'],atlas_lapack_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_ATLAS_LAPACK']:
            print 'Warning: Forced use of ATLAS LAPACK even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ATLAS_LAPACK']:
            print 'Error: FORCE_ATLAS_LAPACK, but failed to find or compile with clapack.h'
            Exit(1)
            
    return result



def CheckCLAPACK(context):
    clapack_source_file = """
extern "C" {
#include "f2c.h"
#include "clapack.h"
}
int main()
{
  char uplo='U', compq='I';
  integer n=1,ldu=1,ldv=1,*iq=0,*iwork=0,info=0,lwork=0,*ipiv=0;
  doublereal *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc_(&uplo,&compq,&n,d,e,u,&ldu,v,&ldv,q,iq,work,iwork,&info);
  dgetri_(&n,u,&ldu,ipiv,work,&lwork,&info);
  return 0;
}
"""

    context.Message('Checking for CLAPACK... ')

    if context.TryCompile(clapack_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],clapack_source_file) or
            CheckLibs(context,['clapack','cblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','cblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack','fblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','fblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack'],clapack_source_file) or
            CheckLibs(context,['lapack'],clapack_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_CLAPACK']:
            print 'Warning: Forced use of CLAPACK even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_CLAPACK']:
            print 'Error: FORCE_CLAPACK, but failed to find or compile with clapack.h and f2c.h'
            Exit(1)
            
    return result




def CheckFLAPACK(context):
    flapack_source_file = """
extern "C" {
#include "util/flapack.h"
}
int main()
{
  char uplo='U', compq='I';
  int n=1,ldu=1,ldv=1,*iq=0,*iwork=0,info=0;
  double *d=0, *e=0, *u=0, *v=0, *q=0, *work=0;
  dbdsdc_(uplo,compq,n,d,e,u,ldu,v,ldv,q,iq,work,iwork,&info);
  return 0;
}
"""

    context.Message('Checking for Fortran LAPACK... ')

    if context.TryCompile(flapack_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],flapack_source_file) or
            CheckLibs(context,['lapack'],flapack_source_file) or
            CheckLibs(context,['lapack','pgftnrtl'],flapack_source_file) or
            CheckLibs(context,['lapack','gfortran'],flapack_source_file) )

        context.Result(result)

        if not result and context.env['FORCE_FLAPACK']:
            print 'Warning: Forced use of FLAPACK even though link test failed.'
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_FLAPACK']:
            print 'Error: FORCE_FLAPACK, but failed compile test'
            Exit(1)
            
    return result




def DoLibraryAndHeaderChecks(config):
    """
    Check for some headers.  
    Mostly we check for a bunch of different BLAS and LAPACK libraries.
    Start with FORCE options if any.
    Otherwise, go through looking for a BLAS library in order:
    MKL, ACML, GOTO, ATLAS, CBLAS, FBLAS
    Once a BLAS library is found, check for a LAPACK library in order:
    (MKL or ACML if that BLAS) CLAPACK, FLAPACK, ATLAS LAPACK
    """

    foundlap = 0
    foundblas = 0
    foundatlasblas = 0

    compiler = config.env['CXXTYPE']
    version = config.env['CXXVERSION_NUMERICAL']

    if not (config.env.has_key('LIBS')) :
      config.env['LIBS'] = []

    if config.env['WITH_BLAS']:
        foundblas = 1  # Set to 0 at end if not found
        # Do FORCE options first:
        if config.env['FORCE_MKL']:
            if config.CheckMKL():
                if (config.env['WITH_LAPACK'] and config.CheckMKL_LAP()) :
                    print 'Using MKL LAPACK'
                    foundlap = 1
                print 'Using MKL BLAS'
                config.env.Append(CPPDEFINES=['MKL'])
                if compiler == 'icpc' and version <= 9.0:
                    env.Append(CPPDEFINES=['NOWORKQUERY'])

        elif config.env['FORCE_ACML']:
            if config.CheckACML():
                if config.env['WITH_LAPACK']:
                    print 'Using ACML LAPACK'
                    foundlap = 1
                print 'Using ACML BLAS'
                config.env.Append(CPPDEFINES=['ACML'])
            else:
                print 'Forced MKL, but failed to find or compile with acml.h'
                Exit(1)
 
        elif config.env['FORCE_GOTO']:
            config.CheckGOTO()
            print 'Using GOTO BLAS'
            config.env.Append(CPPDEFINES=['FBLAS'])
            config.env['NOMIX_SMALL'] = 1

        elif config.env['FORCE_ATLAS']:
            config.CheckATLAS()
            print 'Using ATLAS BLAS'
            config.env.Append(CPPDEFINES=['ATLAS'])
            foundatlasblas = 1

        elif config.env['FORCE_CBLAS']:
            config.CheckCBLAS()
            print 'Using CBLAS'
            config.env.Append(CPPDEFINES=['CBLAS'])

        elif config.env['FORCE_FBLAS']:
            config.CheckFBLAS()
            print 'Using FBLAS'
            config.env.Append(CPPDEFINES=['FBLAS'])

        # If no BLAS is forced, then look for MKL, ACML before more generic
        # (and probably less optimized) BLAS library.
        elif config.CheckMKL() :
            if (config.env['WITH_LAPACK'] and config.CheckMKL_LAP()) :
                print 'Using MKL LAPACK'
                foundlap = 1
            print 'Using MKL BLAS'
            config.env.Append(CPPDEFINES=['MKL'])
 
        elif config.CheckACML() :
            if config.env['WITH_LAPACK']:
                print 'Using ACML LAPACK'
                foundlap = 1
            print 'Using ACML BLAS'
            config.env.Append(CPPDEFINES=['ACML'])

        elif config.CheckGOTO() :
            print 'Using GotoBLAS'
            config.env.Append(CPPDEFINES=['FBLAS'])
            config.env['NOMIX_SMALL'] = 1

        elif config.CheckCBLAS() :
            print 'Using CBLAS'
            config.env.Append(CPPDEFINES=['CBLAS'])

        elif config.CheckATLAS() :
            print 'Using ATLAS'
            foundatlasblas = 1
            config.env.Append(CPPDEFINES=['ATLAS'])

        elif config.CheckFBLAS() :
            print 'Using Fortran BLAS'
            config.env.Append(CPPDEFINES=['FBLAS'])

        else:
            print 'No BLAS libraries found'
            foundblas = 0

    if foundblas and not foundlap and config.env['WITH_LAPACK']:
        foundlap = 1   # Set back to 0 at end if not found.
        if config.env['FORCE_CLAPACK']:
            if config.CheckCLAPACK():
                print 'Using CLAPACK'
                config.env.Append(CPPDEFINES=['CLAPACK'])
            else:
                print 'Forced CLAPACK, but failed to find or compile with clapack.h'
                Exit(1)

        elif config.env['FORCE_FLAPACK']:
            config.CheckFLAPACK()
            print 'Using FLAPACK'
            foundlap = 1
            config.env.Append(CPPDEFINES=['FLAPACK'])

        elif foundatlasblas and config.env['FORCE_ATLAS_LAPACK']:
            config.CheckATLAS_LAP()
            print 'Using ATLAS LAPACK'

        elif config.CheckCLAPACK() :
            config.env.Append(CPPDEFINES=['CLAPACK'])
            print 'Using CLAPACK'

        elif config.CheckFLAPACK() :
            config.env.Append(CPPDEFINES=['FLAPACK'])
            print 'Using Fortran LAPACK'

        elif foundatlasblas and config.CheckATLAS_LAP():
            print 'Using ATLAS LAPACK'

        else :
            print 'No LAPACK libraries found'
            foundlap = 0

    config.env['LAP'] = 0
    if foundlap:
        config.env['LAP'] = 1  # Need this info for test suite compilation
    elif foundblas:
        config.env.Append(CPPDEFINES=['NOLAP'])
    else:
        config.env.Append(CPPDEFINES=['NOBLAS'])
    

def DoConfig(env):
    """
    Configure the system
    """


    # Add some extra paths 
    AddExtraPaths(env)

    # Figure out what kind of compiler we are dealing with
    GetCompilerVersion(env)
   
    # The basic flags for this compiler if not explicitly specified
    BasicCCFlags(env)

    # Some extra flags depending on the options:
    if env['WITH_OPENMP']:
        AddOpenMPFlag(env)
    env.Append(CPPDEFINES=['TMV_OPT=' + env['OPT']])
    if not env['DEBUG']:
        print 'Debugging turned off'
        env.Append(CPPDEFINES=['NDEBUG'])
    if env['MEM_TEST']:
        env.Append(CPPDEFINES=['TMV_MEM_TEST'])
    if env['STATIC'] :
        if env['CXXTYPE'] == 'pgCC':
            env.Append(LINKFLAGS=['-Bstatic'])
        else:
            env.Append(LINKFLAGS=['-static'])
    if not env['INST_FLOAT']:
        env.Append(CPPDEFINES=['TMV_NO_INST_FLOAT'])
    if not env['INST_DOUBLE']:
        env.Append(CPPDEFINES=['TMV_NO_INST_DOUBLE'])
    if env['INST_LONGDOUBLE']:
        env.Append(CPPDEFINES=['TMV_INST_LONGDOUBLE'])
    if env['INST_INT']:
        env.Append(CPPDEFINES=['TMV_INST_INT'])
    if not env['INST_COMPLEX']:
        env.Append(CPPDEFINES=['TMV_NO_INST_COMPLEX'])
    if not env['INST_MIX']:
        env.Append(CPPDEFINES=['TMV_NO_INST_MIX'])
    if not env['TEST_FLOAT']:
        env.Append(CPPDEFINES=['NO_TEST_FLOAT'])
    if not env['TEST_DOUBLE']:
        env.Append(CPPDEFINES=['NO_TEST_DOUBLE'])
    if env['TEST_LONGDOUBLE']:
        env.Append(CPPDEFINES=['TEST_LONGDOUBLE'])
    if env['TEST_INT']:
        env.Append(CPPDEFINES=['TEST_INT'])

    import SCons.SConf

    # Figure out what BLAS and/or LAPACK libraries are on the system
    # MJ: I have had bad luck with scons figuring out when the cache
    #     is invalid.  This just forces a check every time.
    SCons.SConf.SetCacheMode('force')
    config = env.Configure(custom_tests = {
        'CheckMKL' : CheckMKL ,
        'CheckACML' : CheckACML ,
        'CheckGOTO' : CheckGOTO ,
        'CheckATLAS' : CheckATLAS ,
        'CheckCBLAS' : CheckCBLAS ,
        'CheckFBLAS' : CheckFBLAS ,
        'CheckMKL_LAP' : CheckMKL_LAP ,
        'CheckATLAS_LAP' : CheckATLAS_LAP ,
        'CheckCLAPACK' : CheckCLAPACK ,
        'CheckFLAPACK' : CheckFLAPACK })
    DoLibraryAndHeaderChecks(config)
    env = config.Finish()
    # MJ: Turn the cache back on now, since we want it for the
    #     main compilation steps.
    SCons.SConf.SetCacheMode('auto')


#
# main program
#

if not GetOption('help'):

    env = initial_env

    if env['IMPORT_ENV']:
        # I couldn't figure out how to get this option before the 
        # initial constructor.  So this seems a bit inefficient to me.
        # But I think it works, so good enough for now.
        env = Environment(ENV=os.environ)
        # Now repeat the stuff that has already been done to initial_env
        opts.Update(env)
        opts.Save(config_file,env)
        Help(opts.GenerateHelpText(env))

    # Set up the configuration
    DoConfig(env)
 
    # subdirectory SConscript files can use this function
    env['__readfunc'] = ReadFileList

    # subdirectores to process.  We process src by default
    script_files = [os.path.join(src_dir,'SConscript')]
    for d in subdirs:
        #if d in COMMAND_LINE_TARGETS:
            script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports=['env'])


