#include "gtest/gtest.h"
#include "vector3.h"
#include "kahan_number.h"
#include "own_math.h"


namespace UnitTest {

    template<size_t N, typename T>
    void math_sign() {
        T high = SpaceH::big_value<T>::value;
        SpaceH::Random<T> uniform;
        for(size_t i = 0 ; i < N; ++i) {
            T x = uniform(0,high);
            ASSERT_EQ(1, SpaceH::sign(x));
            x = uniform(-high, 0);
            ASSERT_EQ(-1, SpaceH::sign(x));
        }
        EXPECT_EQ(-1,SpaceH::sign(0));
    }

    template<size_t N, typename T>
    void math_step() {
        T high = SpaceH::big_value<T>::value;
        SpaceH::Random<T> uniform;
        for(size_t i = 0 ; i < N; ++i) {
            T x = uniform(0,high);
            ASSERT_EQ(1, SpaceH::stepfunction(x));
            x = uniform(-high, 0);
            ASSERT_EQ(0, SpaceH::stepfunction(x));
        }
        EXPECT_EQ(-1, SpaceH::sign(0));
    }

    template<size_t N, typename T>
    void math_epsilon() {
        T high = SpaceH::big_value<T>::value;
        SpaceH::Random<T> uniform;
        for(size_t i = 0 ; i < N; ++i) {
            T x = uniform(-high, high);
            ASSERT_EQ(x, x + SpaceH::epsilon<T>::value);
            ASSERT_EQ(x, x - SpaceH::epsilon<T>::value);
        }
    }
}

TEST(MathTest, Sign) {
    constexpr size_t sample_num = 1000000;
    UnitTest::math_sign<sample_num, double>();
    UnitTest::math_sign<sample_num, float>();
    UnitTest::math_sign<sample_num, SpaceH::precise_d>();
    UnitTest::math_sign<sample_num, SpaceH::precise_f>();
}

TEST(MathTest, Step) {
    constexpr size_t sample_num = 1000000;
    UnitTest::math_step<sample_num, double>();
    UnitTest::math_step<sample_num, float>();
    UnitTest::math_step<sample_num, SpaceH::precise_d>();
    UnitTest::math_step<sample_num, SpaceH::precise_f>();
}

TEST(MathTest, Epsilon) {
    constexpr size_t sample_num = 1000000;
    UnitTest::math_epsilon<sample_num, double>();
    UnitTest::math_epsilon<sample_num, float>();
    UnitTest::math_epsilon<sample_num, SpaceH::precise_d>();
    UnitTest::math_epsilon<sample_num, SpaceH::precise_f>();
}