set(CMAKE_Fortran_COMPILER "/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "19.1.0.20200306")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/ranlib")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED )
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/include;/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/libfabric/1.10.1/include;/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/ucx/1.8.0/include;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/imkl/2020.1.217/mkl/include/fftw;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/imkl/2020.1.217/mkl/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/tbb/include;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/compiler/include/intel64;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/compiler/include/icc;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/compiler/include;/usr/local/include;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/libfabric/1.10.1/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/ucx/1.8.0/lib;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/imkl/2020.1.217/mkl/lib/intel64;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/imkl/2020.1.217/lib/intel64;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/tbb/lib/intel64/gcc4.8;/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin;/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64;/cvmfs/soft.computecanada.ca/gentoo/2020/lib64;/cvmfs/soft.computecanada.ca/gentoo/2020/usr/lib64;/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib;/cvmfs/soft.computecanada.ca/gentoo/2020/lib;/cvmfs/soft.computecanada.ca/gentoo/2020/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
