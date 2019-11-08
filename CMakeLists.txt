cmake_minimum_required(VERSION 3.12)
project(SpaceHub)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++17 -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++17 -O3")


set(CMAKE_CXX_STANDARD 17)

add_subdirectory(./googletest)
include_directories(./googletest/googletest/include)
link_directories(./googletest/googlemock/include)

set(TMP_HEADER_FILES
        src/args-callback/callbacks.hpp

        src/integrator/symplectic/symplectic-integrator.hpp
        src/integrator/Gauss-Dadau.hpp

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
        src/ode-iterator/Burlish-Stoer.hpp
        src/ode-iterator/const-iterator.hpp
        src/ode-iterator/IAS15.hpp

        src/orbits/orbits.hpp

        src/particles/point-particles.hpp
        src/particles/finite-size.hpp

        src/particle-system/base-system.hpp
        src/particle-system/chain-system.hpp
        src/particle-system/archain.hpp
        src/particle-system/chain.hpp
        src/particle-system/regu-system.hpp
        src/particle-system/octree.hpp

        src/tools/timer.hpp
        src/tools/config-reader.hpp
        src/tools/auto-name.hpp

        src/vector/vector3.hpp
        src/vector/vector3d.hpp
        src/vector/vector3pd.hpp

        src/interaction/accelerations.hpp
        src/coords.hpp
        src/interaction/interaction.hpp
        src/particle-system/particle-system.hpp
        src/rand-generator.hpp
        src/core-computation.hpp
        src/dev-tools.hpp
        src/kahan-number.hpp
        src/macros.hpp
        src/math.hpp
        src/particles/particles.hpp
        src/simulator.hpp
        src/spaceHub.hpp
        src/type-class.hpp
        src/ode-iterator/error-checker/error-checker.hpp
        src/ode-iterator/error-checker/worst-offender.hpp
        src/ode-iterator/error-checker/RMS.hpp
        src/ode-iterator/step-controller/stepsize-controller.hpp
        src/ode-iterator/step-controller/PID-controller.hpp src/scattering/cross-section.hpp src/integrator/integrator.hpp src/ode-iterator/error-checker/IAS15-error.hpp src/orbits/particle-manip.hpp src/IO.hpp src/args-callback/collision.hpp test/catch.hpp)

set(SOURCE_FILES
        test/unit_test/main.cpp
        test/unit_test/test_simple_math.cpp
        test/unit_test/test_chain.cpp
        test/unit_test/test_particle.cpp test/unit_test/test_coords.cpp test/unit_test/utest_vector3.cpp test/unit_test/utest_math.cpp)

set(TWOBODY_TEST
        test/regression_test/two-body-test.cpp
        )


set(SOLAR_TEST
        test/regression_test/solar_system.cpp)

set(KOZAI_TEST
        test/regression_test/kozai.cpp)


add_executable(SpaceHub_unit_test ${SOURCE_FILES} ${TMP_HEADER_FILES})
target_link_libraries(SpaceHub_unit_test gtest gtest_main)
target_link_libraries(SpaceHub_unit_test gmock gmock_main)


add_executable(SpaceHub_kozai_test ${TMP_HEADER_FILES} ${KOZAI_TEST})

add_executable(SpaceHub_twobody_test ${TMP_HEADER_FILES} ${TWOBODY_TEST})

add_executable(SpaceHub_solar_test ${TMP_HEADER_FILES} ${SOLAR_TEST})