
cmake_minimum_required( VERSION 3.15 )

project( letsberational )

set(SRCS erf_cody.cpp rationalcubic.cpp normaldistribution.cpp lets_be_rational.cpp)
add_library(lbr SHARED ${SRCS})


add_executable(main main.cpp)
target_link_libraries(main lbr)

