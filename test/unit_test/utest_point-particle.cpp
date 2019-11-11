#include <vector>
#include "../catch.hpp"
#include "utest.hpp"
#include "../../src/particles/point-particles.hpp"
#include "../../src/type-class.hpp"
#include <vector>

TEST_CASE("point particle") {
    using namespace space::particle_set;
    using Type = space::Types<utest_scalar>;
    using Particles = PointParticles<Type>;
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
        for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
            Particle p{UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND};
            REQUIRE(ptcs.number() == i);
            ptcs.emplace_back(p);
            REQUIRE(ptcs.time() == 0);
            REQUIRE(ptcs.idn()[i] == i);
            REQUIRE(ptcs.mass()[i] == p.mass);
            REQUIRE(ptcs.pos().x[i] == p.pos.x);
            REQUIRE(ptcs.pos().y[i] == p.pos.y);
            REQUIRE(ptcs.pos().z[i] == p.pos.z);
            REQUIRE(ptcs.vel().x[i] == p.vel.x);
            REQUIRE(ptcs.vel().y[i] == p.vel.y);
            REQUIRE(ptcs.vel().z[i] == p.vel.z);
        }
    }

    SECTION("Initialize from container") {
        std::vector<Particle> init_con;
        init_con.reserve(RAND_TEST_NUM);
        for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
            init_con.emplace_back(Particle(UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND, UTEST_RAND));
        }
        double t0 = UTEST_RAND;
        Particles new_ptc{t0, init_con};
        REQUIRE(new_ptc.time() == t0);
        REQUIRE(new_ptc.number() == RAND_TEST_NUM);
        for(size_t i = 0 ; i < RAND_TEST_NUM; ++i) {
            REQUIRE(new_ptc.idn()[i] == i);
            REQUIRE(new_ptc.mass()[i] == init_con[i].mass);
            REQUIRE(new_ptc.pos().x[i] == init_con[i].pos.x);
            REQUIRE(new_ptc.pos().y[i] == init_con[i].pos.y);
            REQUIRE(new_ptc.pos().z[i] == init_con[i].pos.z);
            REQUIRE(new_ptc.vel().x[i] == init_con[i].vel.x);
            REQUIRE(new_ptc.vel().y[i] == init_con[i].vel.y);
            REQUIRE(new_ptc.vel().z[i] == init_con[i].vel.z);
        }
    }

}


