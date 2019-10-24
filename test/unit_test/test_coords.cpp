//
// Created by 王艺涵 on 10/6/19.
//
#include "../../src/rand-generator.hpp"
#include "gtest/gtest.h"
#include "../../src/coords.hpp"
using namespace space;
using namespace math;
namespace UnitTest {

  template<typename Container>
  void test_coords_reserve(size_t sample_num) {
    Coords<Container> coords;
    auto high = 10000;

    for (size_t i = 0; i < sample_num; ++i) {
      size_t new_size = static_cast<size_t >(random::Uniform(0, high));

      coords.resize(new_size);

      ASSERT_GE(coords.size(), new_size);

      ASSERT_GE(coords.x.size(), new_size);

      ASSERT_GE(coords.y.size(), new_size);

      ASSERT_GE(coords.z.size(), new_size);
    }
  }

  template<typename Container>
  void test_coords_resize(size_t sample_num) {
    Coords<Container> coords;
    auto high = 10000;

    for (size_t i = 0; i < sample_num; ++i) {
      size_t new_size = static_cast<size_t>(random::Uniform(0, high));

      coords.resize(new_size);

      ASSERT_EQ(coords.size(), new_size);

      ASSERT_EQ(coords.x.size(), new_size);

      ASSERT_EQ(coords.y.size(), new_size);

      ASSERT_EQ(coords.z.size(), new_size);
    }

    coords.clear();
    ASSERT_EQ(coords.size(), 0);

    ASSERT_EQ(coords.x.size(), 0);

    ASSERT_EQ(coords.y.size(), 0);

    ASSERT_EQ(coords.z.size(), 0);
  }

  template<typename Container>
  void test_coords_emplace(size_t sample_num) {
    using Scalar = typename Coords<Container>::Scalar;
    using Vector = typename Coords<Container>::Vector;

    Coords<Container> coords;
    Scalar high = big_value<Scalar>::value;

    coords.reserve(sample_num);

    for (size_t i = 0; i < sample_num; ++i) {
      Scalar x = random::Uniform(-high, high);
      Scalar y = random::Uniform(-high, high);
      Scalar z = random::Uniform(-high, high);

      coords.emplace_back(x, y, z);

      ASSERT_GE(coords.x[i], x);

      ASSERT_GE(coords.y[i], y);

      ASSERT_GE(coords.z[i], z);
    }

    coords.clear();

    for (size_t i = 0; i < sample_num; ++i) {
      Scalar x = random::Uniform(-high, high);
      Scalar y = random::Uniform(-high, high);
      Scalar z = random::Uniform(-high, high);

      coords.emplace_back(Vector(x, y, z));

      ASSERT_GE(coords.x[i], x);

      ASSERT_GE(coords.y[i], y);

      ASSERT_GE(coords.z[i], z);
    }
  }

  template<typename Container>
  void test_coords_load(size_t sample_num) {
    using Scalar = typename Coords<Container>::Scalar;
    using Vector = typename Coords<Container>::Vector;
    Scalar high = big_value<Scalar>::value;

    Container container;

    container.reserve(sample_num * 3);

    for (size_t i = 0; i < sample_num * 3; ++i) {
      container.emplace_back(random::Uniform(-high, high));
    }
    Coords<Container> coords;

    load_to_coords(container.begin(), container.end(), coords);

    auto iter = container.begin();
    for (size_t i = 0; i < sample_num; i++) {
      ASSERT_EQ(coords.x[i], *iter++);
    }

    for (size_t i = 0; i < sample_num; i++) {
      ASSERT_EQ(coords.y[i], *iter++);
    }

    for (size_t i = 0; i < sample_num; i++) {
      ASSERT_EQ(coords.z[i], *iter++);
    }

    Container container1;

    add_coords_to(container1, coords);

    ASSERT_TRUE(container == container1);
  }

  template<typename Container>
  void test_coords_distance(size_t sample_num) {
    using Scalar = typename Coords<Container>::Scalar;
    using Vector = typename Coords<Container>::Vector;

    Coords<Container> coords;
    Scalar high = big_value<Scalar>::value;

    coords.reserve(sample_num);

    for (size_t i = 0; i < sample_num; ++i) {
      Scalar x = random::Uniform(-high, high);
      Scalar y = random::Uniform(-high, high);
      Scalar z = random::Uniform(-high, high);

      coords.emplace_back(x, y, z);
    }

    for (size_t i = 0; i < sample_num; ++i) {
      for (size_t j = 0; j < sample_num; ++j) {
        Scalar dis = distance(coords, i, j);
        Scalar dx = coords.x[i] - coords.x[j];
        Scalar dy = coords.y[i] - coords.y[j];
        Scalar dz = coords.z[i] - coords.z[j];
        ASSERT_EQ(dis, sqrt(dx * dx + dy * dy + dz * dz));
      }
    }
  }
}

TEST(CoordsTest, Reserve) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_coords_reserve<std::vector<double>>(sample_num);
  UnitTest::test_coords_reserve<std::vector<float>>(sample_num);
}

TEST(CoordsTest, Resize) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_coords_resize<std::vector<double>>(sample_num);
  UnitTest::test_coords_resize<std::vector<float>>(sample_num);
}

TEST(CoordsTest, Emplace) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_coords_emplace<std::vector<double>>(sample_num);
  UnitTest::test_coords_emplace<std::vector<float>>(sample_num);
}

TEST(CoordsTest, Load) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_coords_load<std::vector<double>>(sample_num);
  UnitTest::test_coords_load<std::vector<float>>(sample_num);
}

TEST(CoordsTest, distance) {
  constexpr size_t sample_num = 1000;
  UnitTest::test_coords_distance<std::vector<double>>(sample_num);
  UnitTest::test_coords_distance<std::vector<float>>(sample_num);
}