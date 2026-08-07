cmake_minimum_required(VERSION 3.15)
project(vardictcpp VERSION 1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Portable optimisation by default: NO -march=native, so the binary runs on any CPU of the target
# architecture and builds cleanly on x86-64 and on ARM / Apple Silicon (AArch64). The source contains
# no architecture-specific SIMD intrinsics; the compiler auto-vectorises to each architecture's
# baseline (SSE2 on x86-64, NEON on AArch64). This is also what distribution (e.g. Bioconda, which
# forbids -march=native) requires. Any CFLAGS the build system injects (CMAKE_CXX_FLAGS) are preserved.
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# Opt in to host-specific tuning for a locally built, non-distributed binary. Uses -march=native where
# the compiler supports it (x86-64 and GCC on ARM) and falls back to -mcpu=native (Apple/clang on ARM);
# each flag is probed first so an unsupported toolchain simply keeps the portable flags.
option(VARDICTCPP_NATIVE "Tune for the build host's CPU (not portable)" OFF)
if(VARDICTCPP_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" VARDICTCPP_HAS_MARCH_NATIVE)
  if(VARDICTCPP_HAS_MARCH_NATIVE)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
  else()
    check_cxx_compiler_flag("-mcpu=native" VARDICTCPP_HAS_MCPU_NATIVE)
    if(VARDICTCPP_HAS_MCPU_NATIVE)
      string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mcpu=native")
    endif()
  endif()
endif()

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
  src/cigar_modifier.cpp
  src/tovars.cpp
  src/realigner.cpp
  src/printer.cpp
  src/fisher.cpp
  src/somatic.cpp
  src/amplicon.cpp
  src/var2vcf.cpp
)
find_package(Threads REQUIRED)
target_include_directories(vardictcpp PRIVATE ${HTS_INCLUDE_DIR} src)
target_link_libraries(vardictcpp PRIVATE ${HTS_LIBRARY} Threads::Threads)
target_compile_definitions(vardictcpp PRIVATE VARDICTCPP_VERSION="${PROJECT_VERSION}")

# Keep the htslib runtime path so it runs without LD_LIBRARY_PATH gymnastics.
get_filename_component(HTS_LIBDIR ${HTS_LIBRARY} DIRECTORY)
set_target_properties(vardictcpp PROPERTIES
  BUILD_RPATH ${HTS_LIBDIR} INSTALL_RPATH ${HTS_LIBDIR})
