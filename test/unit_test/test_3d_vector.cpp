#include "gtest/gtest.h"
#include "vector3.h"
#include "kahan-number.h"
#include "own-math.h"

namespace UnitTest {
    template<typename T>
    void test_vector_init() {
        SpaceH::Vec3<T> v;
        EXPECT_EQ(0, v.x);
        EXPECT_EQ(0, v.y);
        EXPECT_EQ(0, v.z);
    }

    template<size_t N, typename T>
    void test_vector_construct() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            EXPECT_EQ(x, v.x);
            EXPECT_EQ(y, v.y);
            EXPECT_EQ(z, v.z);
        }

        for (size_t i = 0; i < N; ++i) {
            T s = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(s);
            EXPECT_EQ(s, v.x);
            EXPECT_EQ(s, v.y);
            EXPECT_EQ(s, v.z);
        }

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            SpaceH::Vec3<T> v2(v);
            EXPECT_EQ(x, v2.x);
            EXPECT_EQ(y, v2.y);
            EXPECT_EQ(z, v2.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_negative() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            SpaceH::Vec3<T> nv = -v;
            EXPECT_EQ(-x, nv.x);
            EXPECT_EQ(-y, nv.y);
            ASSERT_EQ(-z, nv.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_abs() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            SpaceH::Vec3<T> abv = v.abs();
            EXPECT_EQ(fabs(x), abv.x);
            EXPECT_EQ(fabs(y), abv.y);
            ASSERT_EQ(fabs(z), abv.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_norm() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            T norm = v.norm();
            ASSERT_EQ(sqrt(x * x + y * y + z * z), norm);
        }
    }

    template<size_t N, typename T>
    void test_vector_norm2() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            T norm2 = v.norm2();
            ASSERT_EQ(x * x + y * y + z * z, norm2);
        }
    }

    template<size_t N, typename T>
    void test_vector_reNorm() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            T rnorm = v.re_norm();
            ASSERT_EQ(1 / sqrt(x * x + y * y + z * z), rnorm);
        }
    }

    template<size_t N, typename T>
    void test_vector_setzero() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x = SpaceH::Random<T>::uniform(-high, high);
            T y = SpaceH::Random<T>::uniform(-high, high);
            T z = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v(x, y, z);
            v.set_zero();
            EXPECT_EQ(0, v.x);
            EXPECT_EQ(0, v.y);
            ASSERT_EQ(0, v.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_add() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            SpaceH::Vec3<T> v3 = v1 + v2;
            EXPECT_EQ(x1 + x2, v3.x);
            EXPECT_EQ(y1 + y2, v3.y);
            ASSERT_EQ(z1 + z2, v3.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_sub() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            SpaceH::Vec3<T> v3 = v1 - v2;
            EXPECT_EQ(x1 - x2, v3.x);
            EXPECT_EQ(y1 - y2, v3.y);
            ASSERT_EQ(z1 - z2, v3.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_mul() {
        T high = sqrt(SpaceH::big_value<T>::value);

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            SpaceH::Vec3<T> v3 = v1 * v2;
            EXPECT_EQ(x1 * x2, v3.x);
            EXPECT_EQ(y1 * y2, v3.y);
            ASSERT_EQ(z1 * z2, v3.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_div() {
        T high = SpaceH::big_value<T>::value;

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            SpaceH::Vec3<T> v3 = v1 / v2;
            EXPECT_EQ(x1 / x2, v3.x);
            EXPECT_EQ(y1 / y2, v3.y);
            ASSERT_EQ(z1 / z2, v3.z);
        }
    }

    template<size_t N, typename T>
    void test_vector_dot() {
        T high = sqrt(SpaceH::big_value<T>::value);

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            T product = dot(v1, v2);
            ASSERT_EQ(product, x1 * x2 + y1 * y2 + z1 * z2);
        }
    }

    template<size_t N, typename T>
    void test_vector_cross() {
        T high = sqrt(SpaceH::big_value<T>::value);

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            SpaceH::Vec3<T> v3 = cross(v1, v2);
            EXPECT_EQ(v3.x, y1 * z2 - y2 * z1);
            EXPECT_EQ(v3.y, z1 * x2 - z2 * x1);
            ASSERT_EQ(v3.z, x1 * y2 - x2 * y1);
        }
    }

    template<size_t N, typename T>
    void test_vector_dist() {
        T high = sqrt(SpaceH::big_value<T>::value);

        for (size_t i = 0; i < N; ++i) {
            T x1 = SpaceH::Random<T>::uniform(-high, high);
            T y1 = SpaceH::Random<T>::uniform(-high, high);
            T z1 = SpaceH::Random<T>::uniform(-high, high);
            T x2 = SpaceH::Random<T>::uniform(-high, high);
            T y2 = SpaceH::Random<T>::uniform(-high, high);
            T z2 = SpaceH::Random<T>::uniform(-high, high);
            SpaceH::Vec3<T> v1(x1, y1, z1);
            SpaceH::Vec3<T> v2(x2, y2, z2);
            T dist = distance(v1, v2);
            ASSERT_EQ(dist, sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
        }
    }
}

TEST(VectorTest, Initialization) {
    UnitTest::test_vector_init<double>();
    UnitTest::test_vector_init<float>();
    UnitTest::test_vector_init<long double>();
    UnitTest::test_vector_init<SpaceH::precise_d>();
    UnitTest::test_vector_init<SpaceH::precise_f>();
}

TEST(VectorTest, Construction) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_construct<sample_num, float>();
    UnitTest::test_vector_construct<sample_num, double>();
    UnitTest::test_vector_construct<sample_num, SpaceH::precise_d>();
    UnitTest::test_vector_construct<sample_num, SpaceH::precise_f>();
}

TEST(VectorTest, Addition) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_add<sample_num, float>();
    UnitTest::test_vector_add<sample_num, double>();
    UnitTest::test_vector_add<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_add<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Subtraction) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_sub<sample_num, float>();
    UnitTest::test_vector_sub<sample_num, double>();
    UnitTest::test_vector_sub<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_sub<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Multiplication) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_mul<sample_num, float>();
    UnitTest::test_vector_mul<sample_num, double>();
    UnitTest::test_vector_mul<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_mul<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Division) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_div<sample_num, float>();
    UnitTest::test_vector_div<sample_num, double>();
    UnitTest::test_vector_div<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_div<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Norm) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_norm<sample_num, float>();
    UnitTest::test_vector_norm<sample_num, double>();
    UnitTest::test_vector_norm<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_norm<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, NormSquare) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_norm2<sample_num, float>();
    UnitTest::test_vector_norm2<sample_num, double>();
    UnitTest::test_vector_norm2<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_norm2<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, ReciprocalNorm) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_reNorm<sample_num, float>();
    UnitTest::test_vector_reNorm<sample_num, double>();
    UnitTest::test_vector_reNorm<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_reNorm<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, AbsoluteValue) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_abs<sample_num, float>();
    UnitTest::test_vector_abs<sample_num, double>();
    UnitTest::test_vector_abs<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_abs<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Opposite) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_negative<sample_num, float>();
    UnitTest::test_vector_negative<sample_num, double>();
    UnitTest::test_vector_negative<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_negative<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, SetZero) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_setzero<sample_num, float>();
    UnitTest::test_vector_setzero<sample_num, double>();
    UnitTest::test_vector_setzero<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_setzero<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, DotProduct) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_dot<sample_num, float>();
    UnitTest::test_vector_dot<sample_num, double>();
    UnitTest::test_vector_dot<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_dot<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, CrossProduct) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_cross<sample_num, float>();
    UnitTest::test_vector_cross<sample_num, double>();
    UnitTest::test_vector_cross<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_cross<sample_num, SpaceH::precise_d>();
}

TEST(VectorTest, Distance) {
    constexpr size_t sample_num = 1000000;
    UnitTest::test_vector_dist<sample_num, float>();
    UnitTest::test_vector_dist<sample_num, double>();
    UnitTest::test_vector_dist<sample_num, SpaceH::precise_f>();
    UnitTest::test_vector_dist<sample_num, SpaceH::precise_d>();
}