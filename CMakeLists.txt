cmake_minimum_required(VERSION 3.12)
project(SpaceHub)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_STANDARD 17)

include_directories(src/)
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

set(TMP_HEADER_FILES
        src/callbacks.h
        src/core_computation.h
        src/dev_tools.h
        src/kahan_number.h
        src/macros.h
        src/own_math.h
        src/particle_system.h
        src/particles.h
        src/solver.h
        src/spaceHub.h
        src/type_class.h
        src/init_creator/init_file_creator.h
        src/init_creator/orbits.h
        src/integrator/symplectic/symplectic2th.h
        src/integrator/symplectic/symplectic4th.h
        src/integrator/symplectic/symplectic6th.h
        src/integrator/symplectic/symplectic8th.h
        src/integrator/symplectic/symplectic10th.h
        src/integrator/Gauss-Dadau.h
        src/interaction/force.h
        src/interaction/interaction.h
        src/interaction/post_newtonian.h
        src/ode_iterator/BSIterator.h
        src/ode_iterator/BSIterator_impl0.h
        src/ode_iterator/const_iterator.h
        src/ode_iterator/IAS15.h
        src/particle_system/chain.h
        src/particle_system/chain_particles.h
        src/particle_system/GAR_system.h
        src/particle_system/regu_system.h
        src/particle_system/regularization.h
        src/tools/timer.h
        src/vector/vector3.h
        src/vector/vector3d.h
        src/vector/vector3pd.h)

set(SOURCE_FILES
        test/unit_test/main.cpp
        test/unit_test/test_3d_vector.cpp
        test/unit_test/test_simple_math.cpp
        test/unit_test/test_particle.cpp
        test/unit_test/test_chain_particle.cpp src/lazy-evaluation/lazy_expr.h src/lazy-evaluation/lazy_array.h src/lazy-evaluation/slice.h src/lazy-evaluation/mkl-array.h src/multi-thread/multi-thread.h src/multi-thread/thread_pool.h src/multi-thread/concurrent-queue.h src/multi-thread/ciostream.h)
set(REGRESSION_TEST
        test/regression_test/simple_test.cpp)
add_executable(SpaceHub_unit_test ${SOURCE_FILES} ${TMP_HEADER_FILES} ${REGRESSION_TEST})
target_link_libraries(SpaceHub_unit_test gtest gtest_main)
target_link_libraries(SpaceHub_unit_test gmock gmock_main)

add_executable(SpaceHub_reg_test ${TMP_HEADER_FILES} ${REGRESSION_TEST})
