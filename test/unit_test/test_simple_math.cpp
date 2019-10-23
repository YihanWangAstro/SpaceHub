#include "../../src/kahan-number.hpp"
#include "../../src/math.hpp"
#include "../../src/rand-generator.hpp"
#include "gtest/gtest.h"

using namespace space;
using namespace math;
namespace UnitTest {

  template<size_t N, typename T>
  void test_math_sign() {
    T high = big_value<T>::value;

    for (size_t i = 0; i < N; ++i) {
      T x = random::Uniform(0, high);
      ASSERT_EQ(1, sign(x));
      x = random::Uniform(-high, 0);
      ASSERT_EQ(-1, sign(x));
    }
    EXPECT_EQ(-1, sign(0));
  }

  template<size_t N, typename T>
  void test_math_step() {
    T high = big_value<T>::value;

    for (size_t i = 0; i < N; ++i) {
      T x = random::Uniform(0, high);
      ASSERT_EQ(1, step(x));
      x = random::Uniform(-high, 0);
      ASSERT_EQ(0, step(x));
    }
    EXPECT_EQ(-1, sign(0));
  }

  template<size_t N, typename T>
  void test_math_epsilon() {
    T high = big_value<T>::value;

    for (size_t i = 0; i < N; ++i) {
      T x = random::Uniform(-high, high);
      ASSERT_EQ(x, x + 0.5 * epsilon<T>::value);
      ASSERT_EQ(x, x - 0.5 * epsilon<T>::value);
    }
  }
}  // namespace UnitTest

TEST(MathTest, Sign) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_sign<sample_num, double>();
  UnitTest::test_math_sign<sample_num, float>();
  // UnitTest::test_math_sign<sample_num,  precise_d>();
  // UnitTest::test_math_sign<sample_num,  precise_f>();
}

TEST(MathTest, Step) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_step<sample_num, double>();
  UnitTest::test_math_step<sample_num, float>();
  // UnitTest::test_math_step<sample_num,  precise_d>();
  // UnitTest::test_math_step<sample_num,  precise_f>();
}

TEST(MathTest, Epsilon) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_math_epsilon<sample_num, double>();
  UnitTest::test_math_epsilon<sample_num, float>();
  // UnitTest::test_math_epsilon<sample_num,  precise_d>();
  // UnitTest::test_math_epsilon<sample_num,  precise_f>();
}