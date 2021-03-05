/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "../../src/vector/vector3.hpp"
#include "../catch.hpp"
#include "utest.hpp"

using namespace hub;
using Vector = Vec3<utest_scalar>;

TEST_CASE("Vector3", "[default construction]") {
  Vector vec;

  REQUIRE(vec.x == 0);
  REQUIRE(vec.y == 0);
  REQUIRE(vec.z == 0);
  SECTION("[construct with single scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar r = UTEST_RAND;
      Vector v{r};
      REQUIRE(v.x == APPROX(r));
      REQUIRE(v.y == APPROX(r));
      REQUIRE(v.z == APPROX(r));
    }
  }
  SECTION("[construct with 3 components]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;
      Vector v{x, y, z};
      REQUIRE(v.x == APPROX(x));
      REQUIRE(v.y == APPROX(y));
      REQUIRE(v.z == APPROX(z));
    }
  }
  SECTION("[add vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto v3 = v1 + v2;
      REQUIRE(v3.x == APPROX(v1.x + v2.x));
      REQUIRE(v3.y == APPROX(v1.y + v2.y));
      REQUIRE(v3.z == APPROX(v1.z + v2.z));
    }
  }
  SECTION("[sub vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto v3 = v1 - v2;
      REQUIRE(v3.x == APPROX(v1.x - v2.x));
      REQUIRE(v3.y == APPROX(v1.y - v2.y));
      REQUIRE(v3.z == APPROX(v1.z - v2.z));
    }
  }
  SECTION("[mul vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto v3 = v1 * v2;
      REQUIRE(v3.x == APPROX(v1.x * v2.x));
      REQUIRE(v3.y == APPROX(v1.y * v2.y));
      REQUIRE(v3.z == APPROX(v1.z * v2.z));
    }
  }
  SECTION("[div vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto v3 = v1 / v2;
      REQUIRE(v3.x == APPROX(v1.x / v2.x));
      REQUIRE(v3.y == APPROX(v1.y / v2.y));
      REQUIRE(v3.z == APPROX(v1.z / v2.z));
    }
  }
  SECTION("[add scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = v1 + c;
      REQUIRE(v3.x == APPROX(v1.x + c));
      REQUIRE(v3.y == APPROX(v1.y + c));
      REQUIRE(v3.z == APPROX(v1.z + c));
    }
  }
  SECTION("[sub scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = v1 - c;
      REQUIRE(v3.x == APPROX(v1.x - c));
      REQUIRE(v3.y == APPROX(v1.y - c));
      REQUIRE(v3.z == APPROX(v1.z - c));
    }
  }
  SECTION("[mul scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = v1 * c;
      REQUIRE(v3.x == APPROX(v1.x * c));
      REQUIRE(v3.y == APPROX(v1.y * c));
      REQUIRE(v3.z == APPROX(v1.z * c));
    }
  }
  SECTION("[div scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = v1 / c;
      REQUIRE(v3.x == APPROX(v1.x / c));
      REQUIRE(v3.y == APPROX(v1.y / c));
      REQUIRE(v3.z == APPROX(v1.z / c));
    }
  }
  SECTION("[negative]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      auto v3 = -v1;
      REQUIRE(v3.x == APPROX(-v1.x));
      REQUIRE(v3.y == APPROX(-v1.y));
      REQUIRE(v3.z == APPROX(-v1.z));
    }
  }
  SECTION("[absolute value]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      auto v3 = v1.abs();
      REQUIRE(v3.x == APPROX(fabs(v1.x)));
      REQUIRE(v3.y == APPROX(fabs(v1.y)));
      REQUIRE(v3.z == APPROX(fabs(v1.z)));
    }
  }
  SECTION("[operator += vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;

      Vector v1{x, y, z}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      v1 += v2;
      REQUIRE(v1.x == APPROX(x + v2.x));
      REQUIRE(v1.y == APPROX(y + v2.y));
      REQUIRE(v1.z == APPROX(z + v2.z));
    }
  }
  SECTION("[operator-= vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;

      Vector v1{x, y, z}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      v1 -= v2;
      REQUIRE(v1.x == APPROX(x - v2.x));
      REQUIRE(v1.y == APPROX(y - v2.y));
      REQUIRE(v1.z == APPROX(z - v2.z));
    }
  }
  SECTION("[operator *= vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;

      Vector v1{x, y, z}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      v1 *= v2;
      REQUIRE(v1.x == APPROX(x * v2.x));
      REQUIRE(v1.y == APPROX(y * v2.y));
      REQUIRE(v1.z == APPROX(z * v2.z));
    }
  }
  SECTION("[operator/= vector]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;

      Vector v1{x, y, z}, v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      v1 /= v2;
      REQUIRE(v1.x == APPROX(x / v2.x));
      REQUIRE(v1.y == APPROX(y / v2.y));
      REQUIRE(v1.z == APPROX(z / v2.z));
    }
  }
  SECTION("[operator+= scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;
      utest_scalar r = UTEST_RAND;

      Vector v1{x, y, z};
      v1 += r;
      REQUIRE(v1.x == APPROX(x + r));
      REQUIRE(v1.y == APPROX(y + r));
      REQUIRE(v1.z == APPROX(z + r));
    }
  }
  SECTION("[operator-= scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;
      utest_scalar r = UTEST_RAND;

      Vector v1{x, y, z};
      v1 -= r;
      REQUIRE(v1.x == APPROX(x - r));
      REQUIRE(v1.y == APPROX(y - r));
      REQUIRE(v1.z == APPROX(z - r));
    }
  }
  SECTION("[operator*= scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;
      utest_scalar r = UTEST_RAND;

      Vector v1{x, y, z};
      v1 *= r;
      REQUIRE(v1.x == APPROX(x * r));
      REQUIRE(v1.y == APPROX(y * r));
      REQUIRE(v1.z == APPROX(z * r));
    }
  }
  SECTION("[operator/= scalar]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      utest_scalar x = UTEST_RAND;
      utest_scalar y = UTEST_RAND;
      utest_scalar z = UTEST_RAND;
      utest_scalar r = UTEST_RAND;

      Vector v1{x, y, z};
      v1 /= r;
      REQUIRE(v1.x == APPROX(x / r));
      REQUIRE(v1.y == APPROX(y / r));
      REQUIRE(v1.z == APPROX(z / r));
    }
  }
  SECTION("[operator=]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      Vector v2;
      v2 = v1;
      REQUIRE(v1.x == APPROX(v2.x));
      REQUIRE(v1.y == APPROX(v2.y));
      REQUIRE(v1.z == APPROX(v2.z));
    }
  }
  SECTION("[norm]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      auto m = norm(v1);
      REQUIRE(m == APPROX(sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)));
    }
  }
  SECTION("[norm2]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      auto m = norm2(v1);
      REQUIRE(m == APPROX(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z));
    }
  }
  SECTION("[scalar add]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = c + v1;
      REQUIRE(v3.x == APPROX(v1.x + c));
      REQUIRE(v3.y == APPROX(v1.y + c));
      REQUIRE(v3.z == APPROX(v1.z + c));
    }
  }
  SECTION("[scalar sub]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = c - v1;
      REQUIRE(v3.x == APPROX(c - v1.x));
      REQUIRE(v3.y == APPROX(c - v1.y));
      REQUIRE(v3.z == APPROX(c - v1.z));
    }
  }
  SECTION("[scalar mul]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = c * v1;
      REQUIRE(v3.x == APPROX(v1.x * c));
      REQUIRE(v3.y == APPROX(v1.y * c));
      REQUIRE(v3.z == APPROX(v1.z * c));
    }
  }
  SECTION("[scalar div]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      utest_scalar c = UTEST_RAND;

      auto v3 = c / v1;
      REQUIRE(v3.x == APPROX(c / v1.x));
      REQUIRE(v3.y == APPROX(c / v1.y));
      REQUIRE(v3.z == APPROX(c / v1.z));
    }
  }
  SECTION("[dot product]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      Vector v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto p = dot(v1, v2);

      REQUIRE(p == APPROX(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z));
    }
  }
  SECTION("[dot product]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      Vector v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto p = dot(v1, v2);

      REQUIRE(p == APPROX(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z));
    }
  }
  SECTION("[cross product]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      Vector v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto p = cross(v1, v2);

      REQUIRE(p.x == APPROX(v1.y * v2.z - v1.z * v2.y));
      REQUIRE(p.y == APPROX(v1.z * v2.x - v1.x * v2.z));
      REQUIRE(p.z == APPROX(v1.x * v2.y - v1.y * v2.x));
    }
  }
  SECTION("[distance]") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Vector v1{UTEST_RAND, UTEST_RAND, UTEST_RAND};
      Vector v2{UTEST_RAND, UTEST_RAND, UTEST_RAND};

      auto d = distance(v1, v2);
      auto dx = v1.x - v2.x;
      auto dy = v1.y - v2.y;
      auto dz = v1.z - v2.z;

      REQUIRE(d == APPROX(sqrt(dx * dx + dy * dy + dz * dz)));
    }
  }
}
