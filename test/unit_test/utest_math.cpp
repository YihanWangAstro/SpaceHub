
#include "../../src/math.hpp"
#include "../catch.hpp"
#include "utest.hpp"

using namespace space;

TEST_CASE("Math") {
  SECTION("min") { REQUIRE(2 == math::min(2, 3)); }
}