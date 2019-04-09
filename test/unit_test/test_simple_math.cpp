#include "gtest/gtest.h"
#include "vector3.h"
#include "kahan-number.h"
#include "own-math.h"
#include "rand-generator.tpp"

namespace UnitTest {

    template<size_t N, typename T>
    void test_math_sign() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random::Uniform<T>::get(0, high);
            ASSERT_EQ(1, SpaceH::sign(x));
            x = SpaceH::Random::Uniform<T>::get(-high, 0);
            ASSERT_EQ(-1, SpaceH::sign(x));
        }
        EXPECT_EQ(-1, SpaceH::sign(0));
    }

    template<size_t N, typename T>
    void test_math_step() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random::Uniform<T>::get(0, high);
            ASSERT_EQ(1, SpaceH::stepfunction(x));
            x = SpaceH::Random::Uniform<T>::get(-high, 0);
            ASSERT_EQ(0, SpaceH::stepfunction(x));
        }
        EXPECT_EQ(-1, SpaceH::sign(0));
    }

    template<size_t N, typename T>
    void test_math_epsilon() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random::Uniform<T>::get(-high, high);
            ASSERT_EQ(x, x + 0.5*SpaceH::epsilon<T>::value);
            ASSERT_EQ(x, x - 0.5*SpaceH::epsilon<T>::value);
        }
    }
}

TEST(MathTest, Sign) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_math_sign<sample_num, double>();
    UnitTest::test_math_sign<sample_num, float>();
    UnitTest::test_math_sign<sample_num, SpaceH::precise_d>();
    UnitTest::test_math_sign<sample_num, SpaceH::precise_f>();
}

TEST(MathTest, Step) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_math_step<sample_num, double>();
    UnitTest::test_math_step<sample_num, float>();
    UnitTest::test_math_step<sample_num, SpaceH::precise_d>();
    UnitTest::test_math_step<sample_num, SpaceH::precise_f>();
}

TEST(MathTest, Epsilon) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_math_epsilon<sample_num, double>();
    UnitTest::test_math_epsilon<sample_num, float>();
    UnitTest::test_math_epsilon<sample_num, SpaceH::precise_d>();
    UnitTest::test_math_epsilon<sample_num, SpaceH::precise_f>();
}