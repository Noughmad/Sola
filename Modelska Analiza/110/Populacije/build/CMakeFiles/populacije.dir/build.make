# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/miha/Sola/Modelska Analiza/110/Populacije"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/miha/Sola/Modelska Analiza/110/Populacije/build"

# Include any dependencies generated for this target.
include CMakeFiles/populacije.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/populacije.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/populacije.dir/flags.make

CMakeFiles/populacije.dir/main.cpp.o: CMakeFiles/populacije.dir/flags.make
CMakeFiles/populacije.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report "/home/miha/Sola/Modelska Analiza/110/Populacije/build/CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/populacije.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/populacije.dir/main.cpp.o -c "/home/miha/Sola/Modelska Analiza/110/Populacije/main.cpp"

CMakeFiles/populacije.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/populacije.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E "/home/miha/Sola/Modelska Analiza/110/Populacije/main.cpp" > CMakeFiles/populacije.dir/main.cpp.i

CMakeFiles/populacije.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/populacije.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S "/home/miha/Sola/Modelska Analiza/110/Populacije/main.cpp" -o CMakeFiles/populacije.dir/main.cpp.s

CMakeFiles/populacije.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/populacije.dir/main.cpp.o.requires

CMakeFiles/populacije.dir/main.cpp.o.provides: CMakeFiles/populacije.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/populacije.dir/build.make CMakeFiles/populacije.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/populacije.dir/main.cpp.o.provides

CMakeFiles/populacije.dir/main.cpp.o.provides.build: CMakeFiles/populacije.dir/main.cpp.o

# Object files for target populacije
populacije_OBJECTS = \
"CMakeFiles/populacije.dir/main.cpp.o"

# External object files for target populacije
populacije_EXTERNAL_OBJECTS =

populacije: CMakeFiles/populacije.dir/main.cpp.o
populacije: CMakeFiles/populacije.dir/build.make
populacije: CMakeFiles/populacije.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable populacije"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/populacije.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/populacije.dir/build: populacije
.PHONY : CMakeFiles/populacije.dir/build

CMakeFiles/populacije.dir/requires: CMakeFiles/populacije.dir/main.cpp.o.requires
.PHONY : CMakeFiles/populacije.dir/requires

CMakeFiles/populacije.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/populacije.dir/cmake_clean.cmake
.PHONY : CMakeFiles/populacije.dir/clean

CMakeFiles/populacije.dir/depend:
	cd "/home/miha/Sola/Modelska Analiza/110/Populacije/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/miha/Sola/Modelska Analiza/110/Populacije" "/home/miha/Sola/Modelska Analiza/110/Populacije" "/home/miha/Sola/Modelska Analiza/110/Populacije/build" "/home/miha/Sola/Modelska Analiza/110/Populacije/build" "/home/miha/Sola/Modelska Analiza/110/Populacije/build/CMakeFiles/populacije.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/populacije.dir/depend

