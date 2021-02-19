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
#include <vector>

#include "../../src/particles/finite-size.hpp"
#include "../../src/type-class.hpp"
#include "../catch.hpp"
#include "utest.hpp"

TEST_CASE("finite size particle") {
  using namespace space::particle_set;
  using Type = space::Types<utest_scalar>;
  using Particles = SizeParticles<Type>;
  using Particle = typename Particles::Particle;

  Particles ptcs;

  REQUIRE(ptcs.number() == 0);

  SECTION("reserve") {
    ptcs.reserve(100);
    REQUIRE(ptcs.capacity() >= 100);

    ptcs.reserve(10000);
    REQUIRE(ptcs.capacity() >= 10000);

    ptcs.reserve(1000000);
    REQUIRE(ptcs.capacity() >= 1000000);
  }

  SECTION("resize & clear") {
    ptcs.resize(100);
    REQUIRE(ptcs.number() == 100);

    ptcs.resize(10000);
    REQUIRE(ptcs.number() == 10000);

    ptcs.resize(1000000);
    REQUIRE(ptcs.number() == 1000000);

    ptcs.resize(10000);
    REQUIRE(ptcs.number() == 10000);

    ptcs.resize(100);
    REQUIRE(ptcs.number() == 100);

    ptcs.clear();
    REQUIRE(ptcs.number() == 0);
  }

  SECTION("emplace back") {
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      Particle p{UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND};
      REQUIRE(ptcs.number() == i);
      ptcs.emplace_back(p);
      REQUIRE(ptcs.time() == 0);
      REQUIRE(ptcs.idn(i) == i);
      REQUIRE(ptcs.mass(i) == p.mass);
      REQUIRE(ptcs.radius(i) == p.radius);
      REQUIRE(ptcs.pos(i).x == p.pos.x);
      REQUIRE(ptcs.pos(i).y == p.pos.y);
      REQUIRE(ptcs.pos(i).z == p.pos.z);
      REQUIRE(ptcs.vel(i).x == p.vel.x);
      REQUIRE(ptcs.vel(i).y == p.vel.y);
      REQUIRE(ptcs.vel(i).z == p.vel.z);
    }
  }

  SECTION("Initialize from container") {
    std::vector<Particle> init_con;
    init_con.reserve(RAND_TEST_NUM);
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      init_con.emplace_back(
          Particle(UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND));
    }
    double t0 = UTEST_RAND;
    Particles new_ptc{t0, init_con};
    REQUIRE(new_ptc.time() == t0);
    REQUIRE(new_ptc.number() == RAND_TEST_NUM);
    for (size_t i = 0; i < RAND_TEST_NUM; ++i) {
      REQUIRE(new_ptc.idn(i) == i);
      REQUIRE(new_ptc.mass(i) == init_con[i].mass);
      REQUIRE(new_ptc.radius(i) == init_con[i].radius);
      REQUIRE(new_ptc.pos(i).x == init_con[i].pos.x);
      REQUIRE(new_ptc.pos(i).y == init_con[i].pos.y);
      REQUIRE(new_ptc.pos(i).z == init_con[i].pos.z);
      REQUIRE(new_ptc.vel(i).x == init_con[i].vel.x);
      REQUIRE(new_ptc.vel(i).y == init_con[i].vel.y);
      REQUIRE(new_ptc.vel(i).z == init_con[i].vel.z);
    }
  }
}
