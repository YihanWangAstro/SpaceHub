cmake_minimum_required(VERSION 3.12)
project(SpaceHub)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(./googletest)
include_directories(./googletest/googletest/include)
link_directories(./googletest/googlemock/include)

set(TMP_HEADER_FILES
        src/args-callback/callbacks.hpp

        src/integrator/symplectic/symplectic-integrator.hpp
        src/integrator/Gauss-Dadau.h

        src/interaction/post-newtonian.hpp
        src/interaction/newtonian.hpp
        src/interaction/tidal.hpp

        src/lazy-evaluation/lazy_expr.h
        src/lazy-evaluation/lazy_array.h
        src/lazy-evaluation/slice.h
        src/lazy-evaluation/mkl-array.h

        src/multi-thread/multi-thread.hpp
        src/multi-thread/thread_pool.h
        src/multi-thread/concurrent-queue.h
        src/multi-thread/ciostream.hpp

        src/ode-iterator/ode-iterator.hpp
        src/ode-iterator/BS-iterator.hpp
        src/ode-iterator/const-iterator.hpp
        src/ode-iterator/IAS15.h

        src/orbits/orbits.hpp

        src/particle/point-particles.hpp
        src/particle/finite-size.hpp

        src/particle-system/base-system.hpp
        src/particle-system/chain-system.hpp
        src/particle-system/archain.hpp
        src/particle-system/chain.hpp
        src/particle-system/regu-system.hpp
        src/particle-system/octree.hpp

        src/tools/timer.hpp
        src/tools/config-reader.hpp
        src/tools/auto-name.hpp

        src/vector/vector3.h
        src/vector/vector3d.h
        src/vector/vector3pd.h

        src/coords.hpp
        src/interaction.hpp
        src/particle-system.hpp
        src/rand-generator.hpp
        src/core-computation.hpp
        src/dev-tools.hpp
        src/kahan-number.hpp
        src/macros.hpp
        src/own-math.hpp
        src/particle.hpp
        src/simulator.hpp
        src/spaceHub.h
        src/type-class.hpp)

set(SOURCE_FILES
        test/unit_test/main.cpp
        test/unit_test/test_3d_vector.cpp
        test/unit_test/test_simple_math.cpp
        test/unit_test/test_chain.cpp
        )

set(REGRESSION_TEST
        test/regression_test/simple_test.cpp

        )

set(KOZAI_TEST
        test/regression_test/kozai.cpp)
add_executable(SpaceHub_unit_test ${SOURCE_FILES} ${TMP_HEADER_FILES})
target_link_libraries(SpaceHub_unit_test gtest gtest_main)
target_link_libraries(SpaceHub_unit_test gmock gmock_main)

add_executable(SpaceHub_kozai_test ${TMP_HEADER_FILES} ${KOZAI_TEST})

add_executable(SpaceHub_reg_test ${TMP_HEADER_FILES} ${REGRESSION_TEST})