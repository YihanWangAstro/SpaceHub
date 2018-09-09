cmake_minimum_required(VERSION 3.12)
project(SpaceHub)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(src/vector)
include_directories(src/tools)
include_directories(src/particle_system)
include_directories(src/ode_iterator)
include_directories(src/integrator)
include_directories(src/interaction)
include_directories(src/init_creator)

add_subdirectory(./googletest)
include_directories(./googletest/googletest/include)
link_directories(./googletest/googlemock/include)


set(SOURCE_FILES
        #test/regression_test/simple_test.cpp
        test/unit_test/main.cpp
        test/unit_test/3d_vector.cpp)

add_executable(SpaceHub ${SOURCE_FILES})
target_link_libraries(SpaceHub gtest gtest_main)
target_link_libraries(SpaceHub gmock gmock_main)