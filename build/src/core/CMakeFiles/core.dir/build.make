# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/cmake/3.27.7/bin/cmake

# The command to remove a file.
RM = /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/cmake/3.27.7/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kasra545/github-repos/watflood-new

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kasra545/github-repos/watflood-new/build

# Include any dependencies generated for this target.
include src/core/CMakeFiles/core.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/core/CMakeFiles/core.dir/compiler_depend.make

# Include the progress variables for this target.
include src/core/CMakeFiles/core.dir/progress.make

# Include the compile flags for this target's objects.
include src/core/CMakeFiles/core.dir/flags.make

src/core/CMakeFiles/core.dir/area17.f.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/area17.f.o: /home/kasra545/github-repos/watflood-new/src/core/area17.f
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/kasra545/github-repos/watflood-new/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/core/CMakeFiles/core.dir/area17.f.o"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/kasra545/github-repos/watflood-new/src/core/area17.f -o CMakeFiles/core.dir/area17.f.o

src/core/CMakeFiles/core.dir/area17.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/core.dir/area17.f.i"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/kasra545/github-repos/watflood-new/src/core/area17.f > CMakeFiles/core.dir/area17.f.i

src/core/CMakeFiles/core.dir/area17.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/core.dir/area17.f.s"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/kasra545/github-repos/watflood-new/src/core/area17.f -o CMakeFiles/core.dir/area17.f.s

src/core/CMakeFiles/core.dir/Areacg.f.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/Areacg.f.o: /home/kasra545/github-repos/watflood-new/src/core/Areacg.f
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/kasra545/github-repos/watflood-new/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object src/core/CMakeFiles/core.dir/Areacg.f.o"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/kasra545/github-repos/watflood-new/src/core/Areacg.f -o CMakeFiles/core.dir/Areacg.f.o

src/core/CMakeFiles/core.dir/Areacg.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/core.dir/Areacg.f.i"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/kasra545/github-repos/watflood-new/src/core/Areacg.f > CMakeFiles/core.dir/Areacg.f.i

src/core/CMakeFiles/core.dir/Areacg.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/core.dir/Areacg.f.s"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/kasra545/github-repos/watflood-new/src/core/Areacg.f -o CMakeFiles/core.dir/Areacg.f.s

src/core/CMakeFiles/core.dir/area_debug.f90.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/area_debug.f90.o: /home/kasra545/github-repos/watflood-new/src/core/area_debug.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/kasra545/github-repos/watflood-new/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object src/core/CMakeFiles/core.dir/area_debug.f90.o"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/kasra545/github-repos/watflood-new/src/core/area_debug.f90 -o CMakeFiles/core.dir/area_debug.f90.o

src/core/CMakeFiles/core.dir/area_debug.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/core.dir/area_debug.f90.i"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/kasra545/github-repos/watflood-new/src/core/area_debug.f90 > CMakeFiles/core.dir/area_debug.f90.i

src/core/CMakeFiles/core.dir/area_debug.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/core.dir/area_debug.f90.s"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/kasra545/github-repos/watflood-new/src/core/area_debug.f90 -o CMakeFiles/core.dir/area_debug.f90.s

src/core/CMakeFiles/core.dir/area_watflood.f.o: src/core/CMakeFiles/core.dir/flags.make
src/core/CMakeFiles/core.dir/area_watflood.f.o: /home/kasra545/github-repos/watflood-new/src/core/area_watflood.f
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/kasra545/github-repos/watflood-new/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object src/core/CMakeFiles/core.dir/area_watflood.f.o"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/kasra545/github-repos/watflood-new/src/core/area_watflood.f -o CMakeFiles/core.dir/area_watflood.f.o

src/core/CMakeFiles/core.dir/area_watflood.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/core.dir/area_watflood.f.i"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/kasra545/github-repos/watflood-new/src/core/area_watflood.f > CMakeFiles/core.dir/area_watflood.f.i

src/core/CMakeFiles/core.dir/area_watflood.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/core.dir/area_watflood.f.s"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && /cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/kasra545/github-repos/watflood-new/src/core/area_watflood.f -o CMakeFiles/core.dir/area_watflood.f.s

# Object files for target core
core_OBJECTS = \
"CMakeFiles/core.dir/area17.f.o" \
"CMakeFiles/core.dir/Areacg.f.o" \
"CMakeFiles/core.dir/area_debug.f90.o" \
"CMakeFiles/core.dir/area_watflood.f.o"

# External object files for target core
core_EXTERNAL_OBJECTS =

src/core/libcore.a: src/core/CMakeFiles/core.dir/area17.f.o
src/core/libcore.a: src/core/CMakeFiles/core.dir/Areacg.f.o
src/core/libcore.a: src/core/CMakeFiles/core.dir/area_debug.f90.o
src/core/libcore.a: src/core/CMakeFiles/core.dir/area_watflood.f.o
src/core/libcore.a: src/core/CMakeFiles/core.dir/build.make
src/core/libcore.a: src/core/CMakeFiles/core.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/kasra545/github-repos/watflood-new/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking Fortran static library libcore.a"
	cd /home/kasra545/github-repos/watflood-new/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/core.dir/cmake_clean_target.cmake
	cd /home/kasra545/github-repos/watflood-new/build/src/core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/core.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/core/CMakeFiles/core.dir/build: src/core/libcore.a
.PHONY : src/core/CMakeFiles/core.dir/build

src/core/CMakeFiles/core.dir/clean:
	cd /home/kasra545/github-repos/watflood-new/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/core.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/core.dir/clean

src/core/CMakeFiles/core.dir/depend:
	cd /home/kasra545/github-repos/watflood-new/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kasra545/github-repos/watflood-new /home/kasra545/github-repos/watflood-new/src/core /home/kasra545/github-repos/watflood-new/build /home/kasra545/github-repos/watflood-new/build/src/core /home/kasra545/github-repos/watflood-new/build/src/core/CMakeFiles/core.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/core/CMakeFiles/core.dir/depend
