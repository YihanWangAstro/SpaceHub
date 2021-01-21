
cmake_minimum_required(VERSION 3.12)
project(SpaceHub)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -mavx -pthread")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++20 -O2 -pthread -mavx")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++20 -O3 -pthread -mavx")

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

        src/rand-generator.hpp
        src/core-computation.hpp
        src/dev-tools.hpp
        src/kahan-number.hpp
        src/macros.hpp
        src/math.hpp
        src/simulator.hpp
        src/spaceHub.hpp
        src/type-class.hpp
        src/ode-iterator/error-checker/worst-offender.hpp
        src/ode-iterator/error-checker/RMS.hpp
        src/ode-iterator/step-controller/PID-controller.hpp
        src/scattering/cross-section.hpp
        src/ode-iterator/error-checker/IAS15-error.hpp
        src/orbits/particle-manip.hpp
        src/IO.hpp
        src/args-callback/collision.hpp
        test/catch.hpp)

set(UNIT_TEST
        test/unit_test/utest_vector3.cpp
        test/unit_test/utest_math.cpp
        test/unit_test/utest_point-particle.cpp
        test/unit_test/utest_finite-particle.cpp
        test/unit_test/utest_chain.cpp
        test/unit_test/utest_base-system.cpp)

set(TWOBODY_TEST
        test/regression_test/rtest_two-body.cpp
        )

set(SOLAR_TEST
        test/regression_test/rtest_solar-system.cpp)

set(KOZAI_TEST
        test/regression_test/rtest_kozai.cpp)

set(IAS15_TEST
        test/regression_test/rtest_IAS15.cpp)

set(PN_P_TEST
        test/regression_test/rtest_PN_precession.cpp)

set(PN_R_TEST
        test/regression_test/rtest_PN_radiation.cpp)

set(IMPROVEMENT_TEST
        test/regression_test/rtest_improvement.cpp)


add_executable(SpaceHub_unit_test ${UNIT_TEST} ${TMP_HEADER_FILES})

add_executable(SpaceHub_kozai_test ${TMP_HEADER_FILES} ${KOZAI_TEST})

add_executable(SpaceHub_twobody_test ${TMP_HEADER_FILES} ${TWOBODY_TEST})

add_executable(SpaceHub_IAS15_test ${TMP_HEADER_FILES} ${IAS15_TEST})

add_executable(SpaceHub_improve_test ${TMP_HEADER_FILES} ${IMPROVEMENT_TEST})

add_executable(SpaceHub_PN_p_test ${TMP_HEADER_FILES} ${PN_P_TEST})

add_executable(SpaceHub_PN_r_test ${TMP_HEADER_FILES} ${PN_R_TEST})

add_executable(SpaceHub_solar_test ${TMP_HEADER_FILES} ${SOLAR_TEST} test/regression_test/rtest_samples.hpp)