
cmake_minimum_required(VERSION 3.12)
project(SpaceHub)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -Wall -pthread ")
if(MPFR_VERSION_MAJOR)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lmpfr -lgmp ")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS}  ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ")

set(TMP_HEADER_FILES
        src/args-callback/callbacks.hpp

        src/integrator/symplectic/symplectic-integrator.hpp
        src/integrator/Gauss-Radau.hpp

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

        src/ode-iterator/Bulirsch-Stoer.hpp
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
        src/ode-iterator/error-checker/max-ratio-error.hpp
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

set(EARTH_TEST
        test/regression_test/rtest_earth-system.cpp)

set(KOZAI_TEST
        test/regression_test/rtest_kozai.cpp)

set(PN_P_TEST
        test/regression_test/rtest_PN_precession.cpp)

set(PN_R_TEST
        test/regression_test/rtest_PN_radiation.cpp)


add_executable(SpaceHub_tide_test ${TMP_HEADER_FILES} test/regression_test/rtest_tides.cpp)

add_executable(SpaceHub_simple_test ${TMP_HEADER_FILES} test/regression_test/rtest_basic.cpp)

add_executable(SpaceHub_ABITS_test ${TMP_HEADER_FILES} test/regression_test/rtest_AR_ABITS.cpp)

add_executable(SpaceHub_unit_test ${UNIT_TEST} ${TMP_HEADER_FILES})

add_executable(SpaceHub_kozai_test ${TMP_HEADER_FILES} ${KOZAI_TEST})

add_executable(SpaceHub_earth_test ${TMP_HEADER_FILES} ${EARTH_TEST})

add_executable(SpaceHub_twobody_test ${TMP_HEADER_FILES} ${TWOBODY_TEST})

add_executable(SpaceHub_PN_p_test ${TMP_HEADER_FILES} ${PN_P_TEST})

add_executable(SpaceHub_PN_r_test ${TMP_HEADER_FILES} ${PN_R_TEST})

add_executable(SpaceHub_solar_test ${TMP_HEADER_FILES} ${SOLAR_TEST} test/regression_test/rtest_samples.hpp)

enable_testing()

add_test(UnitTest SpaceHub_unit_test)

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
#if(CMAKE_COMPILER_IS_GNUCXX)
#    include(CodeCoverage)
#    setup_target_for_coverage_gcovr_xml(NAME Utest_coverage UnitTest coverage)
#endif()