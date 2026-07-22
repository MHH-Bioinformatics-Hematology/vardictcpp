cmake_minimum_required(VERSION 3.15)
project(vardictcpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")

# htslib: allow an explicit prefix (e.g. a conda env) via -DHTSLIB_ROOT=...
# Falls back to pkg-config / default search paths.
find_path(HTS_INCLUDE_DIR htslib/sam.h
  HINTS ${HTSLIB_ROOT}/include $ENV{CONDA_PREFIX}/include)
find_library(HTS_LIBRARY hts
  HINTS ${HTSLIB_ROOT}/lib $ENV{CONDA_PREFIX}/lib)

if(NOT HTS_INCLUDE_DIR OR NOT HTS_LIBRARY)
  message(FATAL_ERROR "htslib not found. Pass -DHTSLIB_ROOT=<prefix> (containing include/htslib and lib/libhts).")
endif()
message(STATUS "htslib: ${HTS_LIBRARY}  (include ${HTS_INCLUDE_DIR})")

add_executable(vardictcpp
  src/main.cpp
  src/reference.cpp
  src/cigar_parser.cpp
  src/tovars.cpp
  src/realigner.cpp
  src/printer.cpp
)
find_package(Threads REQUIRED)
target_include_directories(vardictcpp PRIVATE ${HTS_INCLUDE_DIR} src)
target_link_libraries(vardictcpp PRIVATE ${HTS_LIBRARY} Threads::Threads)

# Keep the htslib runtime path so it runs without LD_LIBRARY_PATH gymnastics.
get_filename_component(HTS_LIBDIR ${HTS_LIBRARY} DIRECTORY)
set_target_properties(vardictcpp PROPERTIES
  BUILD_RPATH ${HTS_LIBDIR} INSTALL_RPATH ${HTS_LIBDIR})
