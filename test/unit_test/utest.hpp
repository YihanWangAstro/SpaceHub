//
// Created by yihan on 11/8/19.
//

#ifndef SPACEHUB_UTEST_HPP
#define SPACEHUB_UTEST_HPP

#include "../../src/math.hpp"
#include "../../src/rand-generator.hpp"
using utest_scalar = double;
constexpr auto UTEST_EPSILON = space::math::epsilon_v<utest_scalar>;
constexpr auto UTEST_LOW = -1.0;  //-space::math::big_value_v<utest_scalar>;
constexpr auto UTEST_HIGH = 1.0;  // space::math::big_value_v<utest_scalar>;
constexpr size_t RAND_TEST_NUM = 10000;

#define APPROX(x) Approx(x).epsilon(UTEST_EPSILON)
#define UTEST_RAND space::random::Uniform(UTEST_LOW, UTEST_HIGH)
#endif  // SPACEHUB_UTEST_HPP
