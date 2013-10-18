cmake_minimum_required(VERSION 2.8)
project(vmath)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(vmath_SRC src/lwcoord.c
              src/vmath-cubemap.c
              src/vmath-matvec.c
              src/vmath-convert.c
              src/vmath-geo.c
              src/vmath-quaternions.c)

add_library(vmath STATIC ${vmath_SRC})


# Locate libcheck, we cannot build tests without it.
find_path(LIBCHECK_INCLUDE_DIR check.h)
find_library(LIBCHECK_LIBRARY check)

if (LIBCHECK_LIBRARY AND LIBCHECK_INCLUDE_DIR)
  enable_testing()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
  add_subdirectory(tests)
endif (LIBCHECK_LIBRARY AND LIBCHECK_INCLUDE_DIR)
