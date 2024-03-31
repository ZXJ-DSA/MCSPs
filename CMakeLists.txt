cmake_minimum_required(VERSION 3.14)
project(MCSPs)

set(CMAKE_CXX_STANDARD 14)

### add stxxl
# include the STXXL library in my-project
add_subdirectory(stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

# boost
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

add_definitions(-w) # omit warning information

add_executable(Preprocess preprocess.cpp head.h graph.h preprocess.h preprocess.hpp ioxxl.h memory_info.h)
add_executable(MCSP_IM mcspim.cpp head.h graph.h graph.hpp)
add_executable(MCSP_EM mcspem.cpp head.h emgraph.h emgraph.hpp)
#add_executable(Shortcut mainSC.cpp shortcut.hpp)

target_link_libraries(Preprocess ${STXXL_LIBRARIES})
target_link_libraries(MCSP_EM ${STXXL_LIBRARIES})

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(Preprocess pthread boost_system boost_thread)# For linux
        target_link_libraries(MCSP_EM pthread boost_system boost_thread)# For linux
#        target_link_libraries(Shortcut pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
        target_link_libraries(Preprocess boost_system boost_thread-mt)# For mac
        target_link_libraries(MCSP_EM boost_system boost_thread-mt)# For mac
#        target_link_libraries(Shortcut pthread boost_system boost_thread-mt)# For mac
    endif()
endif ()