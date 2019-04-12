cmake_minimum_required(VERSION 3.12)
project(SpaceHub)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_STANDARD 17)

include_directories(src/)
include_directories(src/vector)
include_directories(src/tools)
include_directories(src/particle-system)
include_directories(src/ode-iterator)
include_directories(src/integrator)
include_directories(src/interaction)
include_directories(src/orbits)

add_subdirectory(./googletest)
include_directories(./googletest/googletest/include)
link_directories(./googletest/googlemock/include)

set(TMP_HEADER_FILES
        src/args-callback/callbacks.h
        src/core-computation.tpp
        src/dev-tools.h
        src/kahan-number.h
        src/macros.h
        src/own-math.h
        src/particle-system/base-system.tpp
        src/particle.h
        src/solver.h
        src/spaceHub.h
        src/type-class.h
        src/orbits/orbits.h
        src/integrator/symplectic/symplectic-integrator.tpp
        src/integrator/Gauss-Dadau.h
        src/ode-iterator/BS-iterator.tpp
        src/ode-iterator/const-iterator.tpp
        src/ode-iterator/IAS15.h
        src/particle-system/chain.tpp
        src/particle-system/regu-system.tpp
        src/tools/timer.h
        src/vector/vector3.h
        src/vector/vector3d.h
        src/vector/vector3pd.h)

set(SOURCE_FILES
        test/unit_test/main.cpp
        test/unit_test/test_3d_vector.cpp
        test/unit_test/test_simple_math.cpp
        src/lazy-evaluation/lazy_expr.h src/lazy-evaluation/lazy_array.h src/lazy-evaluation/slice.h src/lazy-evaluation/mkl-array.h src/multi-thread/multi-thread.h src/multi-thread/thread_pool.h src/multi-thread/concurrent-queue.h src/multi-thread/ciostream.h src/interaction.h src/tools/config-reader.tpp src/particle-system.h src/interaction/newtonian.h src/particle-system/chain-system.tpp src/particle-system/archain.tpp src/integrator/symplectic/symplectic-integrator.h src/ode-iterator/ode-iterator.h src/particle/point-particle.tpp src/particle/finite-size.tpp src/tools/auto-name.tpp src/rand-generator.tpp src/interaction/tidal.h test/unit_test/test_chain.cpp)
set(REGRESSION_TEST
        test/regression_test/simple_test.cpp
        src/interaction/post-newtonian.h)
#add_executable(SpaceHub_unit_test ${SOURCE_FILES} ${TMP_HEADER_FILES} ${REGRESSION_TEST})
#target_link_libraries(SpaceHub_unit_test gtest gtest_main)
#target_link_libraries(SpaceHub_unit_test gmock gmock_main)

add_executable(SpaceHub_reg_test ${TMP_HEADER_FILES} ${REGRESSION_TEST})
