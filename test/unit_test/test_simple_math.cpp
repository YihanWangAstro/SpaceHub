#include "../../src/kahan-number.hpp"
#include "../../src/own-math.hpp"
#include "../../src/rand-generator.hpp"
#include "gtest/gtest.h"

namespace UnitTest {

template <size_t N, typename T>
void test_math_sign() {
  T high = space::big_value<T>::value;

  for (size_t i = 0; i < N; ++i) {
    T x = space::randomGen::Uniform<T>::get(0, high);
    ASSERT_EQ(1, space::sign(x));
    x = space::randomGen::Uniform<T>::get(-high, 0);
    ASSERT_EQ(-1, space::sign(x));
  }
  EXPECT_EQ(-1, space::sign(0));
}

template <size_t N, typename T>
void test_math_step() {
  T high = space::big_value<T>::value;

  for (size_t i = 0; i < N; ++i) {
    T x = space::randomGen::Uniform<T>::get(0, high);
    ASSERT_EQ(1, space::stepfunction(x));
    x = space::randomGen::Uniform<T>::get(-high, 0);
    ASSERT_EQ(0, space::stepfunction(x));
  }
  EXPECT_EQ(-1, space::sign(0));
}

template <size_t N, typename T>
void test_math_epsilon() {
  T high = space::big_value<T>::value;

  for (size_t i = 0; i < N; ++i) {
    T x = space::randomGen::Uniform<T>::get(-high, high);
    ASSERT_EQ(x, x + 0.5 * space::epsilon<T>::value);
    ASSERT_EQ(x, x - 0.5 * space::epsilon<T>::value);
  }
}
}  // namespace UnitTest

TEST(MathTest, Sign) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_sign<sample_num, double>();
  UnitTest::test_math_sign<sample_num, float>();
  // UnitTest::test_math_sign<sample_num, space::precise_d>();
  // UnitTest::test_math_sign<sample_num, space::precise_f>();
}

TEST(MathTest, Step) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_step<sample_num, double>();
  UnitTest::test_math_step<sample_num, float>();
  // UnitTest::test_math_step<sample_num, space::precise_d>();
  // UnitTest::test_math_step<sample_num, space::precise_f>();
}

TEST(MathTest, Epsilon) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_epsilon<sample_num, double>();
  UnitTest::test_math_epsilon<sample_num, float>();
  // UnitTest::test_math_epsilon<sample_num, space::precise_d>();
  // UnitTest::test_math_epsilon<sample_num, space::precise_f>();
}