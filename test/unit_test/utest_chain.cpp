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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/

#include <iomanip>

#include "../../src/particle-system/chain.hpp"
#include "../../src/type-class.hpp"
#include "../catch.hpp"
#include "utest.hpp"
// TODO test update_chain()
TEST_CASE("particle system chain xy") {
  using type_sys = space::Types<utest_scalar>;
  using VectorArray = typename type_sys::VectorArray;
  using Vector = typename type_sys::Vector;
  using IdxArray = typename type_sys::IdxArray;
  using ScalarArray = typename type_sys::ScalarArray;

  ScalarArray mass{1, 2, 3, 3, 2, 1};
  IdxArray idx{4, 5, 0, 1, 2, 3};
  VectorArray pos;
  pos.emplace_back(0, 0, 0);
  pos.emplace_back(1, 1, 0);
  pos.emplace_back(2, 3, 0);
  pos.emplace_back(-1, 3, 0);
  pos.emplace_back(5, -1, 0);
  pos.emplace_back(-1, -4, 0);

  space::calc::move_to_com(mass, pos);

  SECTION("create index") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    IdxArray expected{4, 5, 0, 1, 2, 3};
    REQUIRE(expected == idx);
  }

  SECTION("calc chain") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    VectorArray chain_pos{pos.size()};
    space::Chain::calc_chain(pos, chain_pos, idx);
    REQUIRE(chain_pos[0].x == APPROX(-6));
    REQUIRE(chain_pos[0].y == APPROX(-3));
    REQUIRE(chain_pos[0].z == APPROX(0));

    REQUIRE(chain_pos[1].x == APPROX(1));
    REQUIRE(chain_pos[1].y == APPROX(4));
    REQUIRE(chain_pos[1].z == APPROX(0));

    REQUIRE(chain_pos[2].x == APPROX(1));
    REQUIRE(chain_pos[2].y == APPROX(1));
    REQUIRE(chain_pos[2].z == APPROX(0));

    REQUIRE(chain_pos[3].x == APPROX(1));
    REQUIRE(chain_pos[3].y == APPROX(2));
    REQUIRE(chain_pos[3].z == APPROX(0));

    REQUIRE(chain_pos[4].x == APPROX(-3));
    REQUIRE(chain_pos[4].y == APPROX(0));
    REQUIRE(chain_pos[4].z == APPROX(0));

    if constexpr (space::Chain::bijective_transfer) {
      REQUIRE(chain_pos[5].x == APPROX(pos[idx[0]].x));
      REQUIRE(chain_pos[5].y == APPROX(pos[idx[0]].y));
      REQUIRE(chain_pos[5].z == APPROX(pos[idx[0]].z));
    } else {
      REQUIRE(chain_pos[5].x == APPROX(0));
      REQUIRE(chain_pos[5].y == APPROX(0));
      REQUIRE(chain_pos[5].z == APPROX(0));
    }
  }

  SECTION("calc cartesian") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    VectorArray chain_pos{pos.size()};
    space::Chain::calc_chain(pos, chain_pos, idx);

    VectorArray cartesian_pos{pos.size()};
    space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);

    for (size_t i = 0; i < pos.size(); ++i) {
      REQUIRE(pos[i].x == APPROX(cartesian_pos[i].x));
      REQUIRE(pos[i].y == APPROX(cartesian_pos[i].y));
      REQUIRE(pos[i].z == APPROX(cartesian_pos[i].z));
    }
  }
}

TEST_CASE("particle system chain yz") {
  using type_sys = space::Types<utest_scalar>;
  using VectorArray = typename type_sys::VectorArray;
  using Vector = typename type_sys::Vector;
  using IdxArray = typename type_sys::IdxArray;
  using ScalarArray = typename type_sys::ScalarArray;

  ScalarArray mass{1, 2, 3, 3, 2, 1};
  IdxArray idx{4, 5, 0, 1, 2, 3};
  VectorArray pos;
  pos.emplace_back(0, 0, 0);
  pos.emplace_back(0, 1, 1);
  pos.emplace_back(0, 2, 3);
  pos.emplace_back(0, -1, 3);
  pos.emplace_back(0, 5, -1);
  pos.emplace_back(0, -1, -4);

  space::calc::move_to_com(mass, pos);

  SECTION("create index") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    IdxArray expected{4, 5, 0, 1, 2, 3};
    REQUIRE(expected == idx);
  }

  SECTION("calc chain") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    VectorArray chain_pos{pos.size()};
    space::Chain::calc_chain(pos, chain_pos, idx);
    REQUIRE(chain_pos[0].y == APPROX(-6));
    REQUIRE(chain_pos[0].z == APPROX(-3));
    REQUIRE(chain_pos[0].x == APPROX(0));

    REQUIRE(chain_pos[1].y == APPROX(1));
    REQUIRE(chain_pos[1].z == APPROX(4));
    REQUIRE(chain_pos[1].x == APPROX(0));

    REQUIRE(chain_pos[2].y == APPROX(1));
    REQUIRE(chain_pos[2].z == APPROX(1));
    REQUIRE(chain_pos[2].x == APPROX(0));

    REQUIRE(chain_pos[3].y == APPROX(1));
    REQUIRE(chain_pos[3].z == APPROX(2));
    REQUIRE(chain_pos[3].x == APPROX(0));

    REQUIRE(chain_pos[4].y == APPROX(-3));
    REQUIRE(chain_pos[4].z == APPROX(0));
    REQUIRE(chain_pos[4].x == APPROX(0));

    if constexpr (space::Chain::bijective_transfer) {
      REQUIRE(chain_pos[5].x == APPROX(pos[idx[0]].x));
      REQUIRE(chain_pos[5].y == APPROX(pos[idx[0]].y));
      REQUIRE(chain_pos[5].z == APPROX(pos[idx[0]].z));
    } else {
      REQUIRE(chain_pos[5].y == APPROX(0));
      REQUIRE(chain_pos[5].z == APPROX(0));
      REQUIRE(chain_pos[5].x == APPROX(0));
    }
  }

  SECTION("calc cartesian") {
    IdxArray idx;
    space::Chain::calc_chain_index(pos, idx);
    VectorArray chain_pos{pos.size()};
    space::Chain::calc_chain(pos, chain_pos, idx);

    VectorArray cartesian_pos{pos.size()};
    space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);

    for (size_t i = 0; i < pos.size(); ++i) {
      REQUIRE(pos[i].x == APPROX(cartesian_pos[i].x));
      REQUIRE(pos[i].y == APPROX(cartesian_pos[i].y));
      REQUIRE(pos[i].z == APPROX(cartesian_pos[i].z));
    }
  }
}