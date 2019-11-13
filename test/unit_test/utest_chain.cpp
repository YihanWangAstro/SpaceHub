#include "../../src/particle-system/chain.hpp"
#include "../../src/type-class.hpp"
#include "../catch.hpp"
#include "utest.hpp"
#include <iomanip>
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

    space::calc::coord_move_to_com(mass, pos);

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
        REQUIRE(chain_pos.x[0] == APPROX(-6));
        REQUIRE(chain_pos.y[0] == APPROX(-3));
        REQUIRE(chain_pos.z[0] == APPROX(0));

        REQUIRE(chain_pos.x[1] == APPROX(1));
        REQUIRE(chain_pos.y[1] == APPROX(4));
        REQUIRE(chain_pos.z[1] == APPROX(0));

        REQUIRE(chain_pos.x[2] == APPROX(1));
        REQUIRE(chain_pos.y[2] == APPROX(1));
        REQUIRE(chain_pos.z[2] == APPROX(0));

        REQUIRE(chain_pos.x[3] == APPROX(1));
        REQUIRE(chain_pos.y[3] == APPROX(2));
        REQUIRE(chain_pos.z[3] == APPROX(0));

        REQUIRE(chain_pos.x[4] == APPROX(-3));
        REQUIRE(chain_pos.y[4] == APPROX(0));
        REQUIRE(chain_pos.z[4] == APPROX(0));

        if constexpr(space::Chain::bijective_transfer) {
            REQUIRE(chain_pos.x[5] == APPROX(pos.x[idx[0]]));
            REQUIRE(chain_pos.y[5] == APPROX(pos.y[idx[0]]));
            REQUIRE(chain_pos.z[5] == APPROX(pos.z[idx[0]]));
        } else {
            REQUIRE(chain_pos.x[5] == APPROX(0));
            REQUIRE(chain_pos.y[5] == APPROX(0));
            REQUIRE(chain_pos.z[5] == APPROX(0));
        }
    }

    SECTION("calc cartesian") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);

        Coord cartesian_pos{pos.size()};
        space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);

        for(size_t i = 0; i < pos.size(); ++i) {
            REQUIRE(pos.x[i] == APPROX(cartesian_pos.x[i]) );
            REQUIRE(pos.y[i] == APPROX(cartesian_pos.y[i]) );
            REQUIRE(pos.z[i] == APPROX(cartesian_pos.z[i]) );
        }
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

    space::calc::coord_move_to_com(mass, pos);

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
        REQUIRE(chain_pos.y[0] == APPROX(-6));
        REQUIRE(chain_pos.z[0] == APPROX(-3));
        REQUIRE(chain_pos.x[0] == APPROX(0));

        REQUIRE(chain_pos.y[1] == APPROX(1));
        REQUIRE(chain_pos.z[1] == APPROX(4));
        REQUIRE(chain_pos.x[1] == APPROX(0));

        REQUIRE(chain_pos.y[2] == APPROX(1));
        REQUIRE(chain_pos.z[2] == APPROX(1));
        REQUIRE(chain_pos.x[2] == APPROX(0));

        REQUIRE(chain_pos.y[3] == APPROX(1));
        REQUIRE(chain_pos.z[3] == APPROX(2));
        REQUIRE(chain_pos.x[3] == APPROX(0));

        REQUIRE(chain_pos.y[4] == APPROX(-3));
        REQUIRE(chain_pos.z[4] == APPROX(0));
        REQUIRE(chain_pos.x[4] == APPROX(0));

        if constexpr(space::Chain::bijective_transfer) {
            REQUIRE(chain_pos.x[5] == APPROX(pos.x[idx[0]]));
            REQUIRE(chain_pos.y[5] == APPROX(pos.y[idx[0]]));
            REQUIRE(chain_pos.z[5] == APPROX(pos.z[idx[0]]));
        } else {
            REQUIRE(chain_pos.y[5] == APPROX(0));
            REQUIRE(chain_pos.z[5] == APPROX(0));
            REQUIRE(chain_pos.x[5] == APPROX(0));
        }
    }

    SECTION("calc cartesian") {
        IdxArray idx;
        space::Chain::calc_chain_index(pos, idx);
        Coord chain_pos{pos.size()};
        space::Chain::calc_chain(pos, chain_pos, idx);

        Coord cartesian_pos{pos.size()};
        space::Chain::calc_cartesian(mass, chain_pos, cartesian_pos, idx);

        for(size_t i = 0; i < pos.size(); ++i) {
            REQUIRE(pos.x[i] == APPROX(cartesian_pos.x[i]));
            REQUIRE(pos.y[i] == APPROX(cartesian_pos.y[i]));
            REQUIRE(pos.z[i] == APPROX(cartesian_pos.z[i]));
        }

    }
}