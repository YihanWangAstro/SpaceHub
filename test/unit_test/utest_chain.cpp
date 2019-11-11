#include "../../src/particle-system/chain.hpp"
#include "../../src/type-class.hpp"
#include "../catch.hpp"
#include "utest.hpp"

TEST_CASE("particle system chain xy") {
    using type_sys = space::Types<utest_scalar>;
    using Coord = typename type_sys::Coord;
    using Vector = typename Coord::Vector;
    using IdxArray = typename type_sys::IdxArray;
    using ScalarArray = typename type_sys::ScalarArray;

    ScalarArray mass{1, 2, 3, 3, 2, 1};
    IdxArray idx{4, 5, 0, 1, 2, 3};
    Coord pos;
    pos.emplace_back(0, 0, 0);
    pos.emplace_back(1, 1, 0);
    pos.emplace_back(2, 3, 0);
    pos.emplace_back(-1, 3, 0);
    pos.emplace_back(5, -1, 0);
    pos.emplace_back(-1, -4, 0);

    SECTION("create index") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        IdxArray expected{4, 5, 0, 1, 2, 3};
        REQUIRE(expected == idx);
    }

    SECTION("calc chain") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);
        REQUIRE(chain_pos.x[0] == -6);
        REQUIRE(chain_pos.y[0] == -3);
        REQUIRE(chain_pos.z[0] == 0);

        REQUIRE(chain_pos.x[1] == 1);
        REQUIRE(chain_pos.y[1] == 4);
        REQUIRE(chain_pos.z[1] == 0);

        REQUIRE(chain_pos.x[2] == 1);
        REQUIRE(chain_pos.y[2] == 1);
        REQUIRE(chain_pos.z[2] == 0);

        REQUIRE(chain_pos.x[3] == 1);
        REQUIRE(chain_pos.y[3] == 2);
        REQUIRE(chain_pos.z[3] == 0);

        REQUIRE(chain_pos.x[4] == -3);
        REQUIRE(chain_pos.y[4] == 0);
        REQUIRE(chain_pos.z[4] == 0);

        if constexpr(space::Chain::bijective_transfer) {
            REQUIRE(chain_pos.x[5] == 5);
            REQUIRE(chain_pos.y[5] == -1);
            REQUIRE(chain_pos.z[5] == 0);
        } else {
            REQUIRE(chain_pos.x[5] == 0);
            REQUIRE(chain_pos.y[5] == 0);
            REQUIRE(chain_pos.z[5] == 0);
        }
    }

    SECTION("calc cartesian") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);

        Coord cartesian_pos{pos.size()};
        space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);
        
        REQUIRE(pos.x == cartesian_pos.x);
        REQUIRE(pos.y == cartesian_pos.y);
        REQUIRE(pos.z == cartesian_pos.z);
    }
}

TEST_CASE("particle system chain yz") {
    using type_sys = space::Types<utest_scalar>;
    using Coord = typename type_sys::Coord;
    using Vector = typename Coord::Vector;
    using IdxArray = typename type_sys::IdxArray;
    using ScalarArray = typename type_sys::ScalarArray;

    ScalarArray mass{1, 2, 3, 3, 2, 1};
    IdxArray idx{4, 5, 0, 1, 2, 3};
    Coord pos;
    pos.emplace_back(0, 0, 0);
    pos.emplace_back(0, 1, 1);
    pos.emplace_back(0, 2, 3);
    pos.emplace_back(0, -1, 3);
    pos.emplace_back(0, 5, -1);
    pos.emplace_back(0, -1, -4);

    SECTION("create index") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        IdxArray expected{4, 5, 0, 1, 2, 3};
        REQUIRE(expected == idx);
    }

    SECTION("calc chain") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);
        REQUIRE(chain_pos.y[0] == -6);
        REQUIRE(chain_pos.z[0] == -3);
        REQUIRE(chain_pos.x[0] == 0);

        REQUIRE(chain_pos.y[1] == 1);
        REQUIRE(chain_pos.z[1] == 4);
        REQUIRE(chain_pos.x[1] == 0);

        REQUIRE(chain_pos.y[2] == 1);
        REQUIRE(chain_pos.z[2] == 1);
        REQUIRE(chain_pos.x[2] == 0);

        REQUIRE(chain_pos.y[3] == 1);
        REQUIRE(chain_pos.z[3] == 2);
        REQUIRE(chain_pos.x[3] == 0);

        REQUIRE(chain_pos.y[4] == -3);
        REQUIRE(chain_pos.z[4] == 0);
        REQUIRE(chain_pos.x[4] == 0);

        if constexpr(space::Chain::bijective_transfer) {
            REQUIRE(chain_pos.y[5] == 5);
            REQUIRE(chain_pos.z[5] == -1);
            REQUIRE(chain_pos.x[5] == 0);
        } else {
            REQUIRE(chain_pos.y[5] == 0);
            REQUIRE(chain_pos.z[5] == 0);
            REQUIRE(chain_pos.x[5] == 0);
        }
    }

    SECTION("calc cartesian") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);

        Coord cartesian_pos{pos.size()};
        space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);
        
        REQUIRE(pos.x == cartesian_pos.x);
        REQUIRE(pos.y == cartesian_pos.y);
        REQUIRE(pos.z == cartesian_pos.z);
    }
}