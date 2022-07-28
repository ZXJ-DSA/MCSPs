cmake_minimum_required(VERSION 3.20)
project(MCSPs)

set(CMAKE_CXX_STANDARD 14)

### add stxxl
# include the STXXL library in my-project
add_subdirectory(stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

add_definitions(-w) # omit warning information


add_executable(MCSP_EM mcspem.cpp head.h emgraph.h Timer.h emgraph.hpp)


target_link_libraries(MCSP_EM ${STXXL_LIBRARIES})