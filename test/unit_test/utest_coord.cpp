#include <vector>
#include "../../src/coords.hpp"
#include "../catch.hpp"
#include "utest.hpp"

TEST_CASE("coords std::vector") {
  using Coord = space::Coords<std::vector<utest_scalar>>;
  using Vector = typename Coord::Vector;
  using Scalar = typename Coord::Scalar;

  Coord coord;

  REQUIRE(coord.x.size() == 0);
  REQUIRE(coord.y.size() == 0);
  REQUIRE(coord.z.size() == 0);
  REQUIRE(coord.size() == 0);

  SECTION("size initialization") {
    {
      Coord coord{0};
      REQUIRE(coord.x.size() == 0);
      REQUIRE(coord.y.size() == 0);
      REQUIRE(coord.z.size() == 0);
      REQUIRE(coord.size() == 0);
    }
    {
      Coord coord{100};
      REQUIRE(coord.x.size() == 100);
      REQUIRE(coord.y.size() == 100);
      REQUIRE(coord.z.size() == 100);
      REQUIRE(coord.size() == 100);
    }
    {
      Coord coord{10000};
      REQUIRE(coord.x.size() == 10000);
      REQUIRE(coord.y.size() == 10000);
      REQUIRE(coord.z.size() == 10000);
      REQUIRE(coord.size() == 10000);
    }
    {
      Coord coord{1000000};
      REQUIRE(coord.x.size() == 1000000);
      REQUIRE(coord.y.size() == 1000000);
      REQUIRE(coord.z.size() == 1000000);
      REQUIRE(coord.size() == 1000000);
    }
  }

  SECTION("reserve & shrink to fit") {
      Coord coord{0};

      coord.reserve(100);
      REQUIRE(coord.x.capacity() >= 100);
      REQUIRE(coord.y.capacity() >= 100);
      REQUIRE(coord.z.capacity() >= 100);
      REQUIRE(coord.capacity() >= 100);

      coord.reserve(10000);
      REQUIRE(coord.x.capacity() >= 10000);
      REQUIRE(coord.y.capacity() >= 10000);
      REQUIRE(coord.z.capacity() >= 10000);
      REQUIRE(coord.capacity() >= 10000);

      coord.reserve(1000000);
      REQUIRE(coord.x.capacity() >= 1000000);
      REQUIRE(coord.y.capacity() >= 1000000);
      REQUIRE(coord.z.capacity() >= 1000000);
      REQUIRE(coord.capacity() >= 1000000);

      coord.clear();
      REQUIRE(coord.x.capacity() >= 1000000);
      REQUIRE(coord.y.capacity() >= 1000000);
      REQUIRE(coord.z.capacity() >= 1000000);
      REQUIRE(coord.capacity() >= 1000000);

      coord.shrink_to_fit();
      REQUIRE(coord.x.capacity() == 0);
      REQUIRE(coord.y.capacity() == 0);
      REQUIRE(coord.z.capacity() == 0);
      REQUIRE(coord.capacity() == 0);
  }

  SECTION("resize & clear") {
      Coord coord{0};

      coord.resize(100);
      REQUIRE(coord.x.size() == 100);
      REQUIRE(coord.y.size() == 100);
      REQUIRE(coord.z.size() == 100);
      REQUIRE(coord.size() == 100);

      coord.resize(10000);
      REQUIRE(coord.x.size() == 10000);
      REQUIRE(coord.y.size() == 10000);
      REQUIRE(coord.z.size() == 10000);
      REQUIRE(coord.size() == 10000);

      coord.resize(1000000);
      REQUIRE(coord.x.size() == 1000000);
      REQUIRE(coord.y.size() == 1000000);
      REQUIRE(coord.z.size() == 1000000);
      REQUIRE(coord.size() == 1000000);

      coord.resize(10000);
      REQUIRE(coord.x.size() == 10000);
      REQUIRE(coord.y.size() == 10000);
      REQUIRE(coord.z.size() == 10000);
      REQUIRE(coord.size() == 10000);

      coord.resize(100);
      REQUIRE(coord.x.size() == 100);
      REQUIRE(coord.y.size() == 100);
      REQUIRE(coord.z.size() == 100);
      REQUIRE(coord.size() == 100);

      coord.clear();

      REQUIRE(coord.x.size() == 0);
      REQUIRE(coord.y.size() == 0);
      REQUIRE(coord.z.size() == 0);
      REQUIRE(coord.size() == 0);
  }

  SECTION("emplace_back Vector") {
      Coord coord{0};

      for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
          REQUIRE(coord.x.size() == i);
          REQUIRE(coord.y.size() == i);
          REQUIRE(coord.z.size() == i);
          REQUIRE(coord.size() == i);

          Vector v{UTEST_RAND, UTEST_RAND, UTEST_RAND};
          coord.emplace_back(v);
          REQUIRE(coord.x[i] == v.x);
          REQUIRE(coord.y[i] == v.y);
          REQUIRE(coord.z[i] == v.z);
      }
  }

  SECTION("emplace_back Scalar") {
      Coord coord{0};

      for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
          REQUIRE(coord.x.size() == i);
          REQUIRE(coord.y.size() == i);
          REQUIRE(coord.z.size() == i);
          REQUIRE(coord.size() == i);
         
          Scalar x = UTEST_RAND;
          Scalar y = UTEST_RAND;
          Scalar z = UTEST_RAND;
        
          coord.emplace_back(x, y, z);
          REQUIRE(coord.x[i] == x);
          REQUIRE(coord.y[i] == y);
          REQUIRE(coord.z[i] == z);
      }
  }

  SECTION("set zero") {
      Coord coord{0};
    for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
        coord.emplace_back(UTEST_RAND, UTEST_RAND, UTEST_RAND);
    }
      coord.set_zero();
      for(size_t i = 0 ; i < coord.size(); ++i){
          REQUIRE(coord.x[i] == 0);
          REQUIRE(coord.y[i] == 0);
          REQUIRE(coord.z[i] == 0);
      }
  }
}