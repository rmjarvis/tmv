from __future__ import print_function
# vim: set filetype=python et ts=4 sw=4:
# 
# Always run scons from the root directory of tmv

import os
import sys
import SCons
import platform

print('SCons is version',SCons.__version__,'using python version',platform.python_version())


# Subdirectories containing SConscript files.  We always process these, but 
# there are some other optional ones
subdirs = ['src','bin','share','tests']

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

# first check for a saved conf file
opts = Variables(config_file)

# Now set up options for the command line
opts.Add('CXX','Name of c++ compiler','g++')
opts.Add('FLAGS','Flags to send to the compiler','')
opts.Add('EXTRA_FLAGS','Extra flags to send to the compiler','')
opts.Add(BoolVariable('DEBUG',
        'Turn on debugging statements in compilied library',False))
opts.Add(PathVariable('PREFIX',
        'prefix for installation',default_prefix, PathVariable.PathAccept))
opts.Add(PathVariable('FINAL_PREFIX',
        'final installation prefix if different from PREFIX','', PathVariable.PathAccept))

opts.Add(BoolVariable('WITH_OPENMP',
        'Look for openmp and use if found.', True))
opts.Add(BoolVariable('INST_FLOAT',
        'Instantiate <float> templates in compiled library', True))
opts.Add(BoolVariable('INST_DOUBLE',
        'Instantiate <double> templates in compiled library', True))
opts.Add(BoolVariable('INST_INT',
        'Instantiate <int> templates in compiled library', True))
opts.Add(BoolVariable('INST_LONGDOUBLE',
        'Instantiate <long double> templates in compiled library', False))
opts.Add(BoolVariable('SHARED',
        'Build a shared library',True))

opts.Add(BoolVariable('TEST_FLOAT',
        'Instantiate <float> in the test suite', True))
opts.Add(BoolVariable('TEST_DOUBLE',
        'Instantiate <double> in the test suite', True))
opts.Add(BoolVariable('TEST_INT',
        'Instantiate <int> in the test suite', True))
opts.Add(BoolVariable('TEST_LONGDOUBLE',
        'Instantiate <long double> in the test suite', False))

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
        'Import full environment from calling shell', True))
opts.Add(BoolVariable('IMPORT_PREFIX',
        'Use PREFIX/include and PREFIX/lib in search paths', False))

opts.Add(BoolVariable('WITH_BLAS',
        'Look for blas libraries and link if found.', True))
opts.Add(BoolVariable('WITH_LAPACK',
        'Look for lapack libraries and link if found.', False))
opts.Add(BoolVariable('FORCE_MKL',
        'Force scons to use MKL for BLAS and/or LAPACK', False))
opts.Add(BoolVariable('FORCE_ACML',
        'Force scons to use ACML for BLAS and/or LAPACK', False))
#opts.Add(BoolVariable('FORCE_CLAMD',
        #'Force scons to use clAmdBlas', False))
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
opts.Add(BoolVariable('USE_STEGR',
        'Use the LAPACK function ?stegr for finding eigenvectors', False))
opts.Add(BoolVariable('USE_GEQP3',
        'Use the LAPACK function ?qeqp3 for finding strict QRP decomposition', False))
opts.Add('LIBS','Libraries to send to the linker','')
opts.Add('LINKFLAGS','Flags to use when linking','')

opts.Add(BoolVariable('TEST_DEBUG',
        'Turn on debugging statements in the test suite',False))
opts.Add(BoolVariable('STATIC','Use static linkage', False))
opts.Add(BoolVariable('WITH_SSE',
        'Use SSE commands (only necessary for icpc compilations)', True))
opts.Add('XTEST',
        'Do extra tests in the test suite (1=non-unit step, 2=extra sizes/shapes, 4=mix real/complex,  8=degenerate,  16=extra arithmetic, 32=FortranStyle, 64=extreme matrices) ', 0)
opts.Add(BoolVariable('MEM_DEBUG','Do extra tests for memory leaks other memroy problems', False))
opts.Add(BoolVariable('NAN_TEST','Test memory initialized with nans', False))
opts.Add(BoolVariable('SMALL_TESTS',
        'Make the small test programs: tmvtest1a, tmvtest1b, etc.', False))
opts.Add(BoolVariable('WARN',
        'Add warning compiler flags, like -Wall', False))
opts.Add(BoolVariable('PROFILE',
        'Add profiling compiler flags like -pg', False))
opts.Add(BoolVariable('CACHE_LIB',
        'Cache the results of the library checks', True))
opts.Add(BoolVariable('WITH_UPS',
        'Install the ups directory under PREFIX/ups', False))
opts.Add('N_BUILD_THREADS',
        'Number of build threads to use (0 means use ncpus)', 0)
opts.Add(BoolVariable('USE_UNKNOWN_VARS',
        'Allow other parameters besides the ones listed here.',False))


# This helps us determine of openmp is available
openmp_mingcc_vers = 4.1
openmp_minclang_vers = 3.7
openmp_minicpc_vers = 9.1  # 9.0 is supposed to work but has bugs
openmp_minpgcc_vers = 6.0
openmp_mincc_vers = 5.0    # I don't actually know what this should be.

# used only in /bin right now
def RunInstall(env, targets, subdir):
    install_dir = os.path.join(env['INSTALL_PREFIX'], subdir)
    env.Alias(target='install',
          source=env.Install(dir=install_dir, source=targets))

def RunUninstall(env, targets, subdir):
    # There is no env.Uninstall method, we must build our own
    install_dir = os.path.join(env['INSTALL_PREFIX'], subdir)
    deltarget = Delete("$TARGET")

    # delete from $prefix/bin/
    files = []
    for t in targets:
        ifile = os.path.join(install_dir, os.path.basename(str(t))) 
        files.append(ifile)

    for f in files:
        env.Alias('uninstall', env.Command(f, None, deltarget))


def BasicCCFlags(env):
    """
    """

    compiler = env['CXXTYPE']
    version = env['CXXVERSION_NUMERICAL']

    # First parse the LIBS options if present
    if env['LIBS'] == '':
        env.Replace(LIBS=[])
    else:
        libs = env['LIBS'].split(' ')
        env.Replace(LIBS=libs)
    if compiler == 'g++':
        if version >= 4.4:
            # Workaround for a bug in the g++ v4.4 exception handling
            # I don't think 4.5 actually needs it, but keep >= for now
            # just to be safe.
            env.AppendUnique(LIBS='pthread')

    if env['FLAGS'] == '':
        if compiler == 'g++':
            env.Replace(CCFLAGS=['-O2'])
            env.Append(CCFLAGS=['-fno-strict-aliasing'])
            env['TEST_FLAGS'] = ['-O0']
            if env['PROFILE']:
                env.Append(CCFLAGS=['-pg'])
                env['TEST_FLAGS'] += ['-pg']
            if env['WARN']:
                env.Append(CCFLAGS=['-g3','-ansi','-pedantic-errors','-Wall','-Werror','-Wno-long-long'])
                env['TEST_FLAGS'] += ['-g3','-ansi','-pedantic-errors','-Wall','-Werror','-Wno-long-long']

        elif compiler == 'clang++':
            env.Replace(CCFLAGS=['-O2'])
            env['TEST_FLAGS'] = ['-O0']
            if env['PROFILE']:
                env.Append(CCFLAGS=['-pg'])
                env['TEST_FLAGS'] += ['-pg']
            if env['WARN']:
                env.Append(CCFLAGS=['-g3','-ansi','-pedantic-errors','-Wall','-Werror','-Wno-long-long'])
                env['TEST_FLAGS'] += ['-g3','-ansi','-pedantic-errors','-Wall','-Werror','-Wno-long-long']

        elif compiler == 'icpc':
            env.Replace(CCFLAGS=['-O2'])
            if env['WITH_SSE']:
                env.Append(CCFLAGS=['-msse2'])
            env['TEST_FLAGS'] = ['-O0']
            if env['PROFILE']:
                env.Append(CCFLAGS=['-pg'])
                env['TEST_FLAGS'] += ['-pg']
            if env['WARN']:
                env.Append(CCFLAGS=['-g','-Wall','-Werror','-wd279,383,810,981'])
                env['TEST_FLAGS'] += ['-g','-Wall','-Werror','-wd279,383,810,981']
                if version >= 9:
                    env.Append(CCFLAGS=['-wd1572'])
                    env['TEST_FLAGS'] += ['-wd1572']
                if version >= 10 and env['WITH_BLAS']:
                    # These warning only show up in mkl.h:
                    env.Append(CCFLAGS=['-wd424,193'])
                if version >= 11:
                    env.Append(CCFLAGS=['-wd2259'])
                    env['TEST_FLAGS'] += ['-wd2259']

        elif compiler == 'pgCC':
            env.Replace(CCFLAGS=['-O2','-fast','-Mcache_align'])
            env['TEST_FLAGS'] = ['-O0']
            if env['PROFILE']:
                # Not sure if this is right...
                env.Append(CCFLAGS=['-pg'])
                env['TEST_FLAGS'] += ['-pg']
            if env['WARN']:
                env.Append(CCFLAGS=['-g'])
                env['TEST_FLAGS'] += ['-g']

        elif compiler == 'CC':
            env.Replace(CCFLAGS=['-O2','-fast','-instances=semiexplicit'])
            env['TEST_FLAGS'] = ['-O1','-instances=semiexplicit']
            if env['WARN']:
                env.Append(CCFLAGS=['-g','+w'])
                env['TEST_FLAGS'] += ['-g','+w']

        elif compiler == 'cl':
            env.Replace(CCFLAGS=['/EHsc','/nologo','/O2','/Oi'])
            env['TEST_FLAGS'] = ['/EHsc','/nologo','/O1']
            if env['WARN']:
                env.Append(CCFLAGS=['/W2','/WX'])
                env['TEST_FLAGS'] += ['/W1','/WX']

        else:
            print('WARNING: Unknown compiler.  You should set FLAGS directly.')
            env.Replace(CCFLAGS=[])
            env['TEST_FLAGS'] = []

    else :
        # If flags are specified as an option use them:
        cxx_flags = env['FLAGS'].split(' ')
        env.Replace(CCFLAGS=cxx_flags)
        env['TEST_FLAGS'] = cxx_flags

    cxx_flags = env['EXTRA_FLAGS'].split(' ')
    env.AppendUnique(CCFLAGS=cxx_flags)
    env['TEST_FLAGS'] += cxx_flags


def AddOpenMPFlag(env):
    """
    Make sure you do this after you have determined the version of
    the compiler.

    g++ uses -fopemnp
    clang++ uses -fopenmp starting with version 3.7
    icpc uses -qopenmp for >= 17.0, -openmp before that
    pgCC uses -mp
    CC uses -xopenmp
    
    Other compilers?
    """
    compiler = env['CXXTYPE']
    version = env['CXXVERSION_NUMERICAL']
    if compiler == 'g++':
        if version < openmp_mingcc_vers: 
            print('No OpenMP support for g++ versions before ',openmp_mingcc_vers)
            env['WITH_OPENMP'] = False
            return
        flag = ['-fopenmp']
        ldflag = ['-fopenmp']
        xlib = ['pthread']
    elif compiler == 'clang++':
        if version < openmp_minclang_vers: 
            print('No OpenMP support for clang++ versions before ',openmp_minclang_vers)
            env['WITH_OPENMP'] = False
            return
        flag = ['-fopenmp']
        ldflag = ['-fopenmp']
        xlib = ['pthread']
    elif compiler == 'icpc':
        if version < openmp_minicpc_vers:
            print('No OpenMP support for icpc versions before ',openmp_minicpc_vers)
            env['WITH_OPENMP'] = False
            return
        if version >= 17:
            flag = ['-qopenmp']
            ldflag = ['-qopenmp']
        else:
            flag = ['-openmp']
            ldflag = ['-openmp']
        xlib = ['pthread']
    elif compiler == 'pgCC':
        if version < openmp_minpgcc_vers:
            print('No OpenMP support for pgCC versions before ',openmp_minpgcc_vers)
            env['WITH_OPENMP'] = False
            return
        flag = ['-mp','--exceptions']
        ldflag = ['-mp']
        xlib = ['pthread']
    elif compiler == 'CC':
        if version < openmp_mincc_vers:
            print('No OpenMP support for CC versions before ',openmp_mincc_vers)
            env['WITH_OPENMP'] = False
            return
        flag = ['-xopenmp']
        ldflag = ['-xopenmp']
        xlib = ['pthread']
    elif compiler == 'cl':
        #flag = ['/openmp']
        #ldflag = ['/openmp']
        #xlib = []
        # The Express edition, which is the one I have, doesn't come with
        # the file omp.h, which we need.  So I am unable to test TMV's
        # OpenMP with cl.  
        # I believe the Professional edition has full OpenMP support,
        # so if you have that, the above lines might work for you.
        # Just uncomment those, and commend the below three lines.
        print('No OpenMP support for cl')
        env['WITH_OPENMP'] = False
        return
    else:
        print('No OpenMP support for compiler ',compiler)
        env['WITH_OPENMP'] = False
        return

    #print 'Adding openmp support:',flag
    #print 'Using OpenMP'
    env['OMP_FLAGS'] = flag 
    env.AppendUnique(CCFLAGS=flag)
    env.AppendUnique(LINKFLAGS=ldflag)
    env.AppendUnique(LIBS=xlib)

def which(program):
    """
    Mimic functionality of unix which command
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if sys.platform == "win32" and not program.endswith(".exe"):
        program += ".exe"

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def GetCompilerVersion(env):
    """
    """
    compiler = which(env['CXX'])
    print('Using compiler:',compiler)

    compiler_real = os.path.realpath(compiler)
    compiler_base = os.path.basename(compiler)

    # Get the compiler type without suffix or path.  
    # e.g. /sw/bin/g++-4 -> g++
    if 'icpc' in compiler_base :
        compilertype = 'icpc'
        versionflag = '--version'
        linenum=0
    elif 'pgCC' in compiler_base :
        compilertype = 'pgCC'
        versionflag = '--version'
        linenum=1
        # pgCC puts the version number on the second line of output.
    elif 'clang++' in compiler_base :
        compilertype = 'clang++'
        versionflag = '--version'
        linenum=0
    elif 'g++' in compiler_base :
        compilertype = 'g++'
        versionflag = '--version'
        linenum=0
    elif 'CC' in compiler_base :
        compilertype = 'CC'
        versionflag = '-V'
        linenum=0
    elif 'cl' in compiler_base :
        compilertype = 'cl'
        versionflag = ''
        linenum=0
    elif 'c++' in compiler_base :
        compilertype = 'c++'
        versionflag = '--version'
        linenum=0
    else :
        compilertype = 'unknown'
        version = 0
        vnum = 0

    if compilertype != 'unknown':
        cmd = compiler + ' ' + versionflag + ' 2>&1'
        lines = os.popen(cmd).readlines()

        # Check if g++ is a symlink for something else:
        if compilertype == 'g++':
            if 'clang' in lines[0] or 'clang' in lines[1]:
                print('Detected clang++ masquerading as g++')
                compilertype = 'clang++'
                # When it is masquerading, the line with the version is the second line.
                linenum=1

        # Check if c++ is a symlink for something else:
        if compilertype == 'c++':
            if 'clang' in lines[0] or 'clang' in lines[1]:
                print('Detected that c++ is really clang++')
                compilertype = 'clang++'
            elif 'g++' in lines[0] or 'gcc' in lines[0]:
                print('Detected that c++ is really g++')
                compilertype = 'g++'
            else:
                print('Cannot determine what kind of compiler c++ really is')
                compilertype = 'unknown'
            # Any others I should look for?

    # redo this check in case was c++ -> unknown
    if compilertype != 'unknown':
        line = lines[linenum]
        import re
        # For clange, the version can show up in one of two places depending on whether this is
        # Apple clang or regular clang.
        if 'LLVM' in line:
            match = re.search(r'LLVM [0-9]+(\.[0-9]+)+', line)
            match_num = 1
        else:
            match = re.search(r'[0-9]+(\.[0-9]+)+', line)
            match_num = 0
    
        if match:
            version = match.group(match_num)
            # Get the version up to the first decimal
            # e.g. for 4.3.1 we only keep 4.3
            vnum = version[0:version.find('.')+2]
        else:
            version = 0
            vnum = 0

    print('compiler version:',version)

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
        paths in PREFIX directory
        paths in EXTRA_*PATH parameters
        paths from the user's environment

    Only paths that actually exists are kept.
    """
    # local includes and lib paths
    # The # symbol means to interpret these from the top-level scons
    # directory even when we are in a sub-directory (src,tests,etc.)
    bin_paths = []
    cpp_paths = ['#include']
    lib_paths1 = ['#lib']
    lib_paths2 = []

    # PREFIX directory
    env['INSTALL_PREFIX'] = env['PREFIX']

    # FINAL_PREFIX is designed for installations like that done by fink where it installs
    # everything into a temporary directory, and then once it finished successfully, it
    # copies the resulting files to a final location.  This pretty much just matters for the
    # tmv-link file to have the right -L flag.
    if env['FINAL_PREFIX'] == '':
        env['FINAL_PREFIX'] = env['PREFIX']

    if env['IMPORT_PREFIX']:
        AddPath(bin_paths, os.path.join(env['PREFIX'],'bin'))
        AddPath(cpp_paths, os.path.join(env['PREFIX'],'include'))
        AddPath(lib_paths1, os.path.join(env['PREFIX'],'lib'))

    # Paths specified in EXTRA_*
    bin_paths += env['EXTRA_PATH'].split(':')
    cpp_paths += env['EXTRA_INCLUDE_PATH'].split(':')
    lib_paths2 += env['EXTRA_LIB_PATH'].split(':')

    # Paths found in environment paths
    if env['IMPORT_PATHS'] and 'PATH' in os.environ:
        paths=os.environ['PATH']
        paths=paths.split(os.pathsep)
        AddPath(bin_paths, paths)

    if env['IMPORT_PATHS'] and 'C_INCLUDE_PATH' in os.environ:
        paths=os.environ['C_INCLUDE_PATH']
        paths=paths.split(os.pathsep)
        AddPath(cpp_paths, paths)

    if env['IMPORT_PATHS'] and 'LIBRARY_PATH' in os.environ:
        paths=os.environ['LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths2, paths)

    if env['IMPORT_PATHS'] and 'LD_LIBRARY_PATH' in os.environ:
        paths=os.environ['LD_LIBRARY_PATH']
        paths=paths.split(os.pathsep)
        AddPath(lib_paths2, paths)


    #print 'bin paths = ',bin_paths
    #print 'cpp paths = ',cpp_paths
    #print 'lib paths1 = ',lib_paths1
    #print 'lib paths2 = ',lib_paths2

    env.PrependENVPath('PATH', bin_paths)
    env.Prepend(CPPPATH= cpp_paths)
    env.Prepend(LIBPATH= lib_paths2)
    env.Prepend(LIBPATH= lib_paths1)
    env['LIBPATH2'] = lib_paths2    # used for the tmv-link file

def ReadFileList(fname):
    """
    This reads a list of whitespace separated values from the input file fname
    and stores it as a list.  We will make this part of the environment so
    other SConscripts can use it
    """
    try:
        files=open(fname).read().split()
    except:
        print('Could not open file:',fname)
        sys.exit(45)
    files = [f.strip() for f in files]
    return files


def CheckLibs(context,try_libs,source_file):
    init_libs = context.env['LIBS']
    context.env.PrependUnique(LIBS=try_libs)
    result = context.TryLink(source_file,'.cpp')
    if not result :
        context.env.Replace(LIBS=init_libs)
    return result
      

def CheckOMP(context):
    omp_source_file = """
#include <omp.h>
#include <vector>

int main()
{
    double res=0.;
#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int mythread = omp_get_thread_num();
        std::vector<double> x(100);
#pragma omp parallel for
        for (int i=0;i<100;i++) x[i] = i+ num_threads * mythread;
        for (int i=0;i<100;i++) res *= x[i];
    }
    return int(res);
}
"""
    context.Message('Checking for OpenMP library... ')

    AddOpenMPFlag(context.env)

    if context.TryCompile(omp_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],omp_source_file) or
            CheckLibs(context,['gomp'],omp_source_file) or
            CheckLibs(context,['omp'],omp_source_file) or
            False)
        context.Result(result)
        if not result:
            print('Unable to determine correct OpenMP library.')
            print('To enable OpenMP, you may need to explicitly specify the correct library')
            print('with the LIBS flag.')

    else:
        result = 0
        context.Result(result)
        print('Unable to compile code using OpenMP pragmas.')

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

    simple_source_file = "int main() { return 0; }\n"

    context.Message('Checking for MKL... ')

    if context.TryCompile(mkl_source_file,'.cpp'):
        # If guide or iomp5 are available, link them along with pthread
        pthread = ['pthread']
        if CheckLibs(context,['guide']+pthread,simple_source_file):
            pthread = ['guide']+pthread
        if CheckLibs(context,['iomp5']+pthread,simple_source_file):
            pthread = ['iomp5']+pthread
        if CheckLibs(context,['dl']+pthread,simple_source_file):
            pthread = ['dl']+pthread

        # Try to pick out the correct thread library
        if context.env['CXXTYPE'] == 'icpc':
            threadlib = ['mkl_intel_thread']
        elif context.env['CXXTYPE'] == 'pgCC':
            threadlib = ['mkl_pgi_thread']
        else:
            threadlib = ['mkl_gnu_thread']
            if CheckLibs(context,['gomp']+pthread,simple_source_file):
                pthread = ['gomp']+pthread

        result = (
            CheckLibs(context,[],mkl_source_file) or
            CheckLibs(context,['mkl'],mkl_source_file) or
            CheckLibs(context,['mkl']+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_em64t']+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_ipf']+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_ia32']+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_intel_lp64','mkl_core']+threadlib+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_gf_lp64','mkl_core']+threadlib+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_intel_lp64','mkl_core','mkl_sequential'],mkl_source_file) or
            CheckLibs(context,['mkl_intel_lp64','mkl_core','mkl_sequential']+pthread,mkl_source_file) or
            CheckLibs(context,['mkl_gf_lp64','mkl_core','mkl_sequential'],mkl_source_file) or
            CheckLibs(context,['mkl_gf_lp64','mkl_core','mkl_sequential']+pthread,mkl_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_MKL']:
            print('WARNING: Forced use of MKL even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your MKL version.")
            print("         You should figure out what libraries you need to link with here:")
            print("         http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_MKL']:
            print()
            print('Error: FORCE_MKL, but failed to find or compile with mkl.h')
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
            CheckLibs(context,['acml','gfortran'],acml_source_file) or
            CheckLibs(context,['acml','gfortran','pthread'],acml_source_file) or
            CheckLibs(context,['acml','pgftnrtl'],acml_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_ACML']:
            print('WARNING: Forced use of ACML even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your ACML version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ACML']:
            print()
            print('Error: FORCE_ACML, but failed to find or compile with acml.h')
            Exit(1)
            
    return result

def CheckCLAMD(context):
# NB: This doesn't work yet...
    clamd_source_file = """
#include "clAmdBlas.h"
int main()
{
    char ta='N', tb='N';
    int M=1,N=1,K=1,lda=1,ldb=1,ldc=1;
    double alpha=1.,beta=1., *A=0, *B=0, *C=0;
    dgemm_(ta,tb,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1,1);
    return 0;
}
"""
    context.Message('Checking for clAmdBlas... ')

    if context.TryCompile(clamd_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],clamd_source_file) or
            CheckLibs(context,['clAmdBlas'],clamd_source_file) or
            CheckLibs(context,['clAmdBlas','OpenCL'],clamd_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_CLAMD']:
            print('WARNING: Forced use of clAmdBlas even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your CLAMD version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_CLAMD']:
            print()
            print('Error: FORCE_CLAMD, but failed to find or compile with clAmdBlas.h')
            Exit(1)
    
    return result


def CheckGOTO(context):
    fblas_source_file = """
extern "C" {
#include "../src/fblas.h"
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

    if context.TryCompile(fblas_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],fblas_source_file) or
            CheckLibs(context,['goto2'],fblas_source_file) or
            CheckLibs(context,['goto2','pthread'],fblas_source_file) or
            CheckLibs(context,['goto2','gomp','pthread'],fblas_source_file) or
            CheckLibs(context,['goto2','gfortran'],fblas_source_file) or
            CheckLibs(context,['goto2','gfortran','pthread'],fblas_source_file) or
            CheckLibs(context,['goto2','gfortran','gomp','pthread'],fblas_source_file) or
            CheckLibs(context,['goto2','pgftnrtl'],fblas_source_file) or
            CheckLibs(context,['goto'],fblas_source_file) or
            CheckLibs(context,['goto','pthread'],fblas_source_file) or
            CheckLibs(context,['goto','gomp','pthread'],fblas_source_file) or
            CheckLibs(context,['goto','gfortran'],fblas_source_file) or
            CheckLibs(context,['goto','gfortran','pthread'],fblas_source_file) or
            CheckLibs(context,['goto','gfortran','gomp','pthread'],fblas_source_file) or
            CheckLibs(context,['goto','pgftnrtl'],fblas_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_GOTO']:
            print('WARNING: Forced use of GOTO even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your GOTOBLAS version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_GOTO']:
            print()
            print('Error: FORCE_GOTO, but failed compile test')
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
            CheckLibs(context,['ptcblas','atlas'],atlas_source_file) or
            CheckLibs(context,['ptcblas','atlas','pthread'],atlas_source_file) or
            CheckLibs(context,['ptcblas','atlas','guide','pthread'],atlas_source_file) or
            CheckLibs(context,['cblas','atlas'],atlas_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_ATLAS']:
            print('WARNING: Forced use of ATLAS even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your ATLAS version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ATLAS']:
            print()
            print('Error: FORCE_ATLAS, but failed to find or compile with cblas.h')
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
            CheckLibs(context,['cblas'],cblas_source_file) or
            CheckLibs(context,['cblas','gfortran'],cblas_source_file) or
            CheckLibs(context,['cblas','pgftnrtl'],cblas_source_file) or
            CheckLibs(context,['blas','gfortran'],cblas_source_file) or
            CheckLibs(context,['blas','pgftnrtl'],cblas_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_CBLAS']:
            print('WARNING: Forced use of CBLAS even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your BLAS version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_CBLAS']:
            print()
            print('Error: FORCE_CBLAS, but failed to find or compile with cblas.h')
            Exit(1)
            
    return result



def CheckFBLAS(context):
    fblas_source_file = """
extern "C" {
#include "../src/fblas.h"
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

    context.Message('Checking for Fortran BLAS... ')

    if context.TryCompile(fblas_source_file,'.cpp'):
        result = (
            CheckLibs(context,[],fblas_source_file) or
            CheckLibs(context,['blas'],fblas_source_file) or
            CheckLibs(context,['blas','gfortran'],fblas_source_file) or
            CheckLibs(context,['blas','pgftnrtl'],fblas_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_FBLAS']:
            print('WARNING: Forced use of FBLAS even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your BLAS version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_FBLAS']:
            print()
            print('Error: FORCE_FBLAS, but failed compile test')
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
         CheckLibs(context,['mkl_lapack','guide'],mkl_lap_source_file) or
         False))

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
            CheckLibs(context,['lapack_atlas'],atlas_lapack_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_ATLAS_LAPACK']:
            print('WARNING: Forced use of ATLAS LAPACK even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your ATLAS version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_ATLAS_LAPACK']:
            print()
            print('Error: FORCE_ATLAS_LAPACK, but failed to find or compile with clapack.h')
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
            CheckLibs(context,['lapack'],clapack_source_file) or
            CheckLibs(context,['clapack'],clapack_source_file) or
            CheckLibs(context,['lapack','cblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','fblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack','cblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack','fblaswr','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','f2c'],clapack_source_file) or
            CheckLibs(context,['clapack','f2c'],clapack_source_file) or
            CheckLibs(context,['lapack','F77'],clapack_source_file) or
            CheckLibs(context,['clapack','F77'],clapack_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_CLAPACK']:
            print('WARNING: Forced use of CLAPACK even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your CLAPACK version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_CLAPACK']:
            print()
            print('Error: FORCE_CLAPACK, but failed to find or compile with clapack.h and f2c.h')
            Exit(1)
            
    return result




def CheckFLAPACK(context):
    flapack_source_file = """
extern "C" {
#include "../src/flapack.h"
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
            CheckLibs(context,['lapack','gfortran'],flapack_source_file) or
            CheckLibs(context,['lapack','pgftnrtl'],flapack_source_file) or
            CheckLibs(context,['lapack','cblaswr','f2c'],flapack_source_file) or
            CheckLibs(context,['lapack','fblaswr','f2c'],flapack_source_file) or
            CheckLibs(context,['clapack','cblaswr','f2c'],flapack_source_file) or
            CheckLibs(context,['clapack','fblaswr','f2c'],flapack_source_file) or
            CheckLibs(context,['lapack','f2c'],flapack_source_file) or
            CheckLibs(context,['clapack','f2c'],flapack_source_file) or
            CheckLibs(context,['lapack','F77'],flapack_source_file) or
            CheckLibs(context,['clapack','F77'],flapack_source_file) or
            False)

        context.Result(result)

        if not result and context.env['FORCE_FLAPACK']:
            print('WARNING: Forced use of FLAPACK even though link test failed.')
            print("         The compiled library will not have the correct linking commands")
            print("         set up for your LAPACK version.")
            print("         You should figure out what libraries you need to link with")
            print("         and add them explicitly with the LIBS option in SCons.")
            result = 1

    else:
        result = 0
        context.Result(result)
        if context.env['FORCE_FLAPACK']:
            print()
            print('Error: FORCE_FLAPACK, but failed compile test')
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

    if not ('LIBS' in config.env) :
      config.env['LIBS'] = []

    if config.env['WITH_BLAS']:
        foundblas = 1  # Set to 0 at end if not found
        # Do FORCE options first:
        if config.env['FORCE_MKL']:
            config.CheckMKL()
            config.env.Append(CPPDEFINES=['MKL'])
            print('Using MKL BLAS')
            if (config.env['WITH_LAPACK'] and config.CheckMKL_LAP()) :
                foundlap = 1
                print('Using MKL LAPACK')
            if compiler == 'icpc' and version <= 9.0:
                # TODO: A better way to do this would be to check
                # for whether work queries work correctly.
                # For now, I know that icpc 9.0 (and probably earlier?)
                # doesn't do work queries.
                env.Append(CPPDEFINES=['NOWORKQUERY'])

        elif config.env['FORCE_ACML']:
            config.CheckACML()
            config.env.Append(CPPDEFINES=['ACML'])
            print('Using ACML BLAS')
            if config.env['WITH_LAPACK']:
                foundlap = 1
                print('Using ACML LAPACK')
 
        #elif config.env['FORCE_CLAMD']:
            #config.CheckCLAMD()
            #config.env.Append(CPPDEFINES=['CLAMD'])
            #print 'Using clAmdBlas'
 
        elif config.env['FORCE_GOTO']:
            config.CheckGOTO()
            config.env.Append(CPPDEFINES=['FBLAS'])
            print('Using GotoBLAS')
            if config.env['WITH_OPENMP']:
                print('Disabling OpenMP, since GotBLAS becomes very slow with it')
                config.env['WITH_OPENMP'] = False

        elif config.env['FORCE_ATLAS']:
            config.CheckATLAS()
            config.env.Append(CPPDEFINES=['ATLAS'])
            foundatlasblas = 1
            print('Using ATLAS BLAS')

        elif config.env['FORCE_CBLAS']:
            config.CheckCBLAS()
            config.env.Append(CPPDEFINES=['CBLAS'])
            print('Using CBLAS')

        elif config.env['FORCE_FBLAS']:
            config.CheckFBLAS()
            config.env.Append(CPPDEFINES=['FBLAS'])
            print('Using FBLAS')

        # If no BLAS is forced, then look for MKL, ACML before more generic
        # (and probably less optimized) BLAS library.
        elif config.CheckMKL() :
            config.env.Append(CPPDEFINES=['MKL'])
            print('Using MKL BLAS')
            if (config.env['WITH_LAPACK'] and config.CheckMKL_LAP()) :
                foundlap = 1
                print('Using MKL LAPACK')
 
        elif config.CheckACML() :
            config.env.Append(CPPDEFINES=['ACML'])
            print('Using ACML BLAS')
            if config.env['WITH_LAPACK']:
                foundlap = 1
                print('Using ACML LAPACK')

        #elif config.CheckCLAMD() :
            #config.env.Append(CPPDEFINES=['CLAMD'])
            #print 'Using clAmdBlas'

        elif config.CheckGOTO() :
            config.env.Append(CPPDEFINES=['FBLAS'])
            print('Using GotoBLAS')
            if config.env['WITH_OPENMP']:
                print('Disabling OpenMP, since GotBLAS becomes very slow with it')
                config.env['WITH_OPENMP'] = False

        elif config.CheckCBLAS() :
            config.env.Append(CPPDEFINES=['CBLAS'])
            print('Using CBLAS')

        elif config.CheckATLAS() :
            foundatlasblas = 1
            config.env.Append(CPPDEFINES=['ATLAS'])
            print('Using ATLAS')

        elif config.CheckFBLAS() :
            config.env.Append(CPPDEFINES=['FBLAS'])
            print('Using Fortran BLAS')

        else:
            foundblas = 0
            print('No BLAS libraries found')

    if foundblas and not foundlap and config.env['WITH_LAPACK']:
        foundlap = 1   # Set back to 0 at end if not found.
        if config.env['FORCE_FLAPACK']:
            config.CheckFLAPACK()
            config.env.Append(CPPDEFINES=['FLAPACK'])
            print('Using FLAPACK')

        elif config.env['FORCE_CLAPACK']:
            config.CheckCLAPACK()
            config.env.Append(CPPDEFINES=['CLAPACK'])
            print('Using CLAPACK')

        elif foundatlasblas and config.env['FORCE_ATLAS_LAPACK']:
            config.CheckATLAS_LAP()
            print('Using ATLAS LAPACK')

        elif config.CheckFLAPACK() :
            config.env.Append(CPPDEFINES=['FLAPACK'])
            print('Using Fortran LAPACK')

        elif config.CheckCLAPACK() :
            config.env.Append(CPPDEFINES=['CLAPACK'])
            print('Using CLAPACK')

        elif foundatlasblas and config.CheckATLAS_LAP():
            print('Using ATLAS LAPACK')

        else :
            foundlap = 0
            print('No LAPACK libraries found')

    config.env['LAP'] = 0
    if foundlap:
        config.env['LAP'] = 1  # Need this info for test suite compilation
    elif foundblas:
        config.env.Append(CPPDEFINES=['NOLAP'])
    else:
        config.env.Append(CPPDEFINES=['NOBLAS'])
    
def GetNCPU():
    """
    Detects the number of CPUs on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, 'sysconf'):
        if 'SC_NPROCESSORS_ONLN' in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf('SC_NPROCESSORS_ONLN')
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2('sysctl -n hw.ncpu')[1].read())
    # Windows:
    if 'NUMBER_OF_PROCESSORS' in os.environ:
        ncpus = int(os.environ['NUMBER_OF_PROCESSORS'])
        if ncpus > 0:
            return ncpus
    return 1 # Default


def DoConfig(env):
    """
    Configure the system
    """


    # Add some extra paths 
    AddExtraPaths(env)

    # Figure out what kind of compiler we are dealing with
    GetCompilerVersion(env)
   
    # If not explicit, set number of jobs according to number of CPUs
    # Note: this doesn't override an explicit scons -jN
    if env.GetOption('num_jobs') == 1: # not set with -jN
        if int(env['N_BUILD_THREADS']) > 0:
            env.SetOption('num_jobs', int(env['N_BUILD_THREADS']))
            # Do this, because if using scons -jN, it doesn't get updated.
            if (env.GetOption('num_jobs') == int(env['N_BUILD_THREADS']) and
                    env.GetOption('num_jobs') != 1):
                print('Using specified number of jobs =',env.GetOption('num_jobs'))
        elif int(env['N_BUILD_THREADS']) == 0:
            ncpu = GetNCPU()
            if ncpu > 4:
                # Don't use more than 4 cpus for the automatic version.
                env.SetOption('num_jobs', 4)
                if env.GetOption('num_jobs') == 4:
                    print('Determined that there are',ncpu, end=' ')
                    print('cpus, but only using 4 jobs to avoid overly high memory use.')
            elif ncpu > 1:
                env.SetOption('num_jobs', ncpu)
                # Importantly, the above line doesn't do anything if the user has run with
                # scons -j1.  It only works if num_jobs was 1 by default, with no explicit -j1.
                if env.GetOption('num_jobs') == ncpu:
                    print('Determined that there are',ncpu,'cpus, so use this many jobs.')
            if env.GetOption('num_jobs') > 1:
                print('You can override this behavior with scons -jN')

    # The basic flags for this compiler if not explicitly specified
    BasicCCFlags(env)

    # Some extra flags depending on the options:
    if not env['DEBUG']:
        print('Debugging turned off')
        env.Append(CPPDEFINES=['TMV_NDEBUG'])
    if env['MEM_DEBUG']:
        env.Append(CPPDEFINES=['TMV_MEM_DEBUG'])
    if '-m32' in env['CCFLAGS']:
        env.Append(LINKFLAGS=['-m32'])
    if '-m64' in env['CCFLAGS']:
        env.Append(LINKFLAGS=['-m64'])
    if env['STATIC'] :
        if env['CXXTYPE'] == 'pgCC':
            env.Append(LINKFLAGS=['-Bstatic'])
        else:
            env.Append(LINKFLAGS=['-static'])

    import SCons.SConf

    # Figure out what BLAS and/or LAPACK libraries are on the system
    # MJ: I have had bad luck with scons figuring out when the cache
    #     is invalid.  This just forces a check every time.
    if not env['CACHE_LIB']:
        SCons.SConf.SetCacheMode('force')
    config = env.Configure(custom_tests = {
        'CheckOMP' : CheckOMP ,
        'CheckMKL' : CheckMKL ,
        'CheckACML' : CheckACML ,
        #'CheckCLAMD' : CheckCLAMD ,
        'CheckGOTO' : CheckGOTO ,
        'CheckATLAS' : CheckATLAS ,
        'CheckCBLAS' : CheckCBLAS ,
        'CheckFBLAS' : CheckFBLAS ,
        'CheckMKL_LAP' : CheckMKL_LAP ,
        'CheckATLAS_LAP' : CheckATLAS_LAP ,
        'CheckCLAPACK' : CheckCLAPACK ,
        'CheckFLAPACK' : CheckFLAPACK })
    DoLibraryAndHeaderChecks(config)

    # Do this after BLAS checks, since we disable it for GotoBLAS
    if config.env['WITH_OPENMP']:
        if not config.CheckOMP():
            config.env['WITH_OPENMP'] = False

    env = config.Finish()
    # MJ: Turn the cache back on now, since we want it for the
    #     main compilation steps.
    if not env['CACHE_LIB']:
        SCons.SConf.SetCacheMode('auto')


#
# main program
#

env = Environment()
opts.Update(env)

if env['IMPORT_ENV']:
    env = Environment(ENV=os.environ)
    opts.Update(env)

# Check for unknown variables in case something is misspelled
unknown = opts.UnknownVariables()
if unknown and not env['USE_UNKNOWN_VARS']:
    print()
    print("Error: Unknown variables:", unknown.keys())
    print('If you are sure these are right (e.g. you want to set some SCons parameters')
    print('that are not in the list of TMV parameters given by scons -h)')
    print('then you can override this check with USE_UNKNOWN_VARS=true')
    Exit(1)

if any(opt.default != env[opt.key] for opt in opts.options):
    print('Using the following (non-default) scons options:')
    for opt in opts.options:
        if (opt.default != env[opt.key]):
            print('   %s = %s'%(opt.key,env[opt.key]))
    print('These can be edited directly in the file %s.'%config_file)
    print('Type scons -h for a full list of available options.')
else:
    print('Using the default scons options')

opts.Save(config_file,env)
Help(opts.GenerateHelpText(env))

if not GetOption('help'):

    if Dir(env['PREFIX']).abspath == Dir('#').abspath:
        print()
        print('Error: PREFIX=%r is invalid.'%env['PREFIX'])
        if (env['PREFIX'] != Dir(env['PREFIX']).abspath and
            env['PREFIX'] != Dir(env['PREFIX']).abspath + '/'):
            print('(expands to %r)'%Dir(env['PREFIX']).abspath)
        print('You should install into some other location, not the root TMV directory.')
        print("Typical choices are '/usr/local' or '~'")
        Exit(1)

    # Set up the configuration
    DoConfig(env)
 
    # subdirectory SConscript files can use this function
    env['_ReadFileList'] = ReadFileList
    env['_InstallProgram'] = RunInstall
    env['_UninstallProgram'] = RunUninstall

    if env['WITH_UPS']:
        subdirs += ['ups']
    if 'doc' in COMMAND_LINE_TARGETS: 
        subdirs += ['doc']
    if 'examples' in COMMAND_LINE_TARGETS: 
        subdirs += ['examples']

    # subdirectores to process.  We process src by default
    script_files = []
    for d in subdirs:
        script_files.append(os.path.join(d,'SConscript'))

    SConscript(script_files, exports=['env'])

