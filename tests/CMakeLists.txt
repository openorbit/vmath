# For now there is only one test here. It should perhaps be broken up.
add_executable(vmath-test math-test.c)
include_directories(${LIBCHECK_INCLUDE_DIR})
target_link_libraries(vmath-test vmath m ${LIBCHECK_LIBRARY})
add_test("vmath-test" vmath-test)
add_dependencies(check vmath-test)
