cmake_minimum_required(VERSION 3.20)
project(MCSPs_master)

set(CMAKE_CXX_STANDARD 14)

### add stxxl
# include the STXXL library in my-project
add_subdirectory(stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

add_executable(MCSPs_master main.cpp head.h emgraph.h Timer.h emgraph.hpp)

target_link_libraries(MCSPs_master ${STXXL_LIBRARIES})
