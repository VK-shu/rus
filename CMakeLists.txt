cmake_minimum_required(VERSION 3.16)
project(sqct LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#INCLUDE_DIRECTORIES(/usr/local/include)
#LINK_DIRECTORIES(/opt/local/lib)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(OpenMP REQUIRED)
find_package(Boost 1.70 REQUIRED COMPONENTS program_options timer chrono system)

find_package(PkgConfig QUIET)
include(FindPackageHandleStandardArgs)

# GMP / GMPXX
find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_GMP QUIET gmp)
  if(PC_GMP_FOUND)
    if(NOT GMP_INCLUDE_DIR)
      find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h HINTS ${PC_GMP_INCLUDE_DIRS})
    endif()
    if(NOT GMP_LIBRARY)
      find_library(GMP_LIBRARY NAMES gmp HINTS ${PC_GMP_LIBRARY_DIRS})
    endif()
  endif()
endif()
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

# MPFR
find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)
find_library(MPFR_LIBRARY NAMES mpfr)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_MPFR QUIET mpfr)
  if(PC_MPFR_FOUND)
    if(NOT MPFR_INCLUDE_DIR)
      find_path(MPFR_INCLUDE_DIR NAMES mpfr.h HINTS ${PC_MPFR_INCLUDE_DIRS})
    endif()
    if(NOT MPFR_LIBRARY)
      find_library(MPFR_LIBRARY NAMES mpfr HINTS ${PC_MPFR_LIBRARY_DIRS})
    endif()
  endif()
endif()
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIR MPFR_LIBRARY)

# PARI (usually not available via pkg-config)
find_path(PARI_INCLUDE_DIR NAMES pari/pari.h)
find_library(PARI_LIBRARY NAMES pari)
find_package_handle_standard_args(PARI DEFAULT_MSG PARI_INCLUDE_DIR PARI_LIBRARY)
#set( CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc -static")
#set(Boost_USE_STATIC_LIBS ON)

# add_executable(sqctSK
#  hprhelpers.cpp
#  resring.cpp
#  rint.cpp
#  matrix2x2.cpp
#  vector2.cpp
#  output.cpp
#  gatelibrary.cpp
 
#  es/numbersgen.cpp
#  es/optsequencegenerator.cpp
#  es/seqlookupcliff.cpp
#  es/exactdecomposer.cpp

#  sk/vector3hpr.cpp
#  sk/epsilonnet.cpp
#  sk/netgenerator.cpp
#  sk/unitaryapproximator.cpp
#  sk/gcommdecomposer.cpp
#  sk/sk.cpp

#  theory/numbers-stat.cpp
#  theory/theoremverification.cpp
#  theory/toptimalitytest.cpp
#  theory/hoptimalitytest.cpp

#  main.cpp
# )

# add_executable(sqct
#  requestprocessor.cpp
#  hprhelpers.cpp
#  resring.cpp
#  rint.cpp
#  matrix2x2.cpp
#  vector2.cpp
#  output.cpp
#  gatelibrary.cpp
#  symbolic_angle.cpp
#  serializers.h
#  timemeasurement.h

 
#  es/numbersgen.cpp
#  es/optsequencegenerator.cpp
#  es/seqlookupcliff.cpp
#  es/exactdecomposer.cpp
 
#  appr/zrot_cache.cpp
#  appr/toptzrot2.cpp
#  appr/normsolver.cpp
#  appr/approxlist.cpp
#  appr/findhalves.cpp
#  appr/rcup.cpp
#  appr/topt-bfs.cpp
#  appr/cup.cpp
#  appr/observations.cpp

#  test.cpp
#  factorzs2.cpp
#  tcount.cpp
#  solvenormequation.cpp
#  fixedpoint.cpp

#  mainA.cpp
# )


add_executable(rusSyn
 requestprocessor.cpp
 hprhelpers.cpp
 resring.cpp
 rint.cpp
 matrix2x2.cpp
 vector2.cpp
 output.cpp
 gatelibrary.cpp
 symbolic_angle.cpp
 serializers.h
 timemeasurement.h

 
 es/numbersgen.cpp
 es/optsequencegenerator.cpp
 es/seqlookupcliff.cpp
 es/exactdecomposer.cpp
 
 appr/zrot_cache.cpp
 appr/toptzrot2.cpp
 appr/normsolver.cpp
 appr/approxlist.cpp
 appr/findhalves.cpp
 appr/rcup.cpp
 appr/topt-bfs.cpp
 appr/cup.cpp
 appr/observations.cpp

 rus/matrix.cpp
 rus/normalization.cpp
 rus/pslq.cpp
 rus/rus.cpp

 test.cpp
 factorzs2.cpp
 tcount.cpp
 solvenormequation.cpp
 fixedpoint.cpp

 main_rus.cpp
)

target_include_directories(rusSyn PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${GMP_INCLUDE_DIR}
  ${MPFR_INCLUDE_DIR}
  ${PARI_INCLUDE_DIR}
)

target_compile_features(rusSyn PRIVATE cxx_std_17)

target_link_libraries(rusSyn PRIVATE
  ${Boost_LIBRARIES}
  OpenMP::OpenMP_CXX
  pthread
  ${PARI_LIBRARY}
  ${MPFR_LIBRARY}
  ${GMPXX_LIBRARY}
  ${GMP_LIBRARY}
  rt
)

message(STATUS "Used boost libraries: ${Boost_LIBRARIES}")
message(STATUS "GMP include: ${GMP_INCLUDE_DIR}, library: ${GMP_LIBRARY}, gmpxx: ${GMPXX_LIBRARY}")
message(STATUS "MPFR include: ${MPFR_INCLUDE_DIR}, library: ${MPFR_LIBRARY}")
message(STATUS "PARI include: ${PARI_INCLUDE_DIR}, library: ${PARI_LIBRARY}")

