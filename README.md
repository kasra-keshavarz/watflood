# WATFLOOD/CHARM
# Compile from source
To retrieve the latest version of WATFLOOD/CHARM, you may use the current
commit of the GitHub repository:
```console
foo@bar:~ $ git clone https://github.com/kasra-keshavarz/watflood
```

The [CMake v3.20 or later](https://cmake.org/download/) build system is 
needed to compile WATFLOOD/CHARM source code. Once available, you may
proceed with creating a `build` directory (preferably in the source
code's root directory) and start the build process:
```console
foo@bar:watflood $ mkdir build && cd build
foo@bar:build $ cmake ../
```

## Specifying `ifort` compiler
Currently, WATFLOOD/CHARM can **only** be compiled with
[Intel *ifort* compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.6xmclp).
To specify the compiler to the build process, you may use the following
command on Unix-like systems:
```console
foo@bar:build $ cmake -DCMAKE_Fortran_COMPILER=ifort ../
```
Or, if the path to `ifort` is not located in the `$PATH` environment
variable, you can do the following
```console
foo@bar:build $ cmake -DCMAKE_Fortran_COMPILER:FILEPATH=/path/to/ifort ..
```
Once `cmake` is done, you may continue by using `make`:
```console
foo@bar:build $ make
```

On a Windows machine, if the `ifort` compiler is not found automatically,
you may choose the compiler as follows (your version of MS Visual Studio
may differ - modify accordingly):
```powershell
> cmake .. -G "Visual Studio 15 2017 Win64"
```

## Specifying build type
Currently, there are two build types available for WATFLOOD/CHARM:
`RELEASE` (default) and `DEBUG`. As an example, you may choose each
build type as following:
```console
foo@bar:build $ cmake -DCMAKE_BUILD_TYPE=RELEASE ../
```

# License
WATFLOOD/CHARM is published under the GNU Lesser General Public License
v3.0 or later.
