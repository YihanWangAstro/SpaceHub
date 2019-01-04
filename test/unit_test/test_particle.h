#ifndef TEST_PARTICLES_H
#define TEST_PARTICLES_H
#include "gtest/gtest.h"
#include "type_class.h"
#include "own_math.h"
#include "test_tools.h"
namespace UnitTest {

    template<size_t particle_num, typename T>
    class ParticleTestInput {
    protected:
        using type         = SpaceH::TypeSystem<T, particle_num>;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using ScalarBuffer = typename type::ScalarBuffer;
    public:
        ParticleTestInput() {
            Scalar high = 1e6;
            for (size_t i = 0; i < particle_num; ++i) {
                pos.emplace_back(SpaceH::RandomVector<Vector>::uniform(-high, high));
                vel.emplace_back(SpaceH::RandomVector<Vector>::uniform(-high, high));
                mass.emplace_back(SpaceH::Random<Scalar>::uniform(0, high));
                radius.emplace_back(SpaceH::Random<Scalar>::uniform(0, high));
                idn.emplace_back(SpaceH::Random<Scalar>::uniform(0, high));
            }
            time = SpaceH::Random<Scalar>::uniform(0, high);
            Fill_stdbuff();
            Fill_actbuff();
        }

        Scalar time;
        std::vector<Vector> pos;
        std::vector<Vector> vel;
        std::vector<Scalar> mass;
        std::vector<Scalar> radius;
        std::vector<size_t> idn;
        ScalarBuffer stdbuff;
        ScalarBuffer actbuff;
    private:
        void Fill_stdbuff() {
            stdbuff.reserve(particle_num * 8 + 1);
            for (auto &p : pos) {
                stdbuff.emplace_back(p.x);
                stdbuff.emplace_back(p.y);
                stdbuff.emplace_back(p.z);
            }
            for (auto &v : vel) {
                stdbuff.emplace_back(v.x);
                stdbuff.emplace_back(v.y);
                stdbuff.emplace_back(v.z);
            }
            for (const auto &m : mass)
                stdbuff.emplace_back(m);
            for (const auto &r : radius)
                stdbuff.emplace_back(r);
            for (const auto &i : idn)
                stdbuff.emplace_back(i);

            stdbuff.emplace_back(time);
        }

        void Fill_actbuff() {
            actbuff.reserve(particle_num * 6 + 1);
            for (auto &p : pos) {
                actbuff.emplace_back(p.x);
                actbuff.emplace_back(p.y);
                actbuff.emplace_back(p.z);
            }
            for (auto &v : vel) {
                actbuff.emplace_back(v.x);
                actbuff.emplace_back(v.y);
                actbuff.emplace_back(v.z);
            }
            actbuff.emplace_back(time);
        }
    };



    template<typename Particles, size_t particle_num>
    void test_particle_std_read(const ParticleTestInput<particle_num, typename Particles::Scalar> &input) {
        std::unique_ptr<Particles> partc{new Particles};

        if constexpr (Particles::type::arraySize == SpaceH::DYNAMICAL)
            partc->resize(particle_num);
        partc->read(input.stdbuff, SpaceH::IO_flag::STD);

        ASSERT_EQ(input.time, partc->time());
        for (size_t i = 0; i < particle_num; ++i) {
            ASSERT_EQ(input.pos[i].x, partc->pos(i).x);
            ASSERT_EQ(input.pos[i].y, partc->pos(i).y);
            ASSERT_EQ(input.pos[i].z, partc->pos(i).z);
            ASSERT_EQ(input.vel[i].x, partc->vel(i).x);
            ASSERT_EQ(input.vel[i].y, partc->vel(i).y);
            ASSERT_EQ(input.vel[i].z, partc->vel(i).z);
            ASSERT_EQ(input.mass[i], partc->mass(i));
            ASSERT_EQ(input.radius[i], partc->radius(i));
            ASSERT_EQ(input.idn[i], partc->idn(i));
        }
    }

    template<typename Particles, size_t particle_num>
    void test_particle_act_read(const ParticleTestInput<particle_num, typename Particles::Scalar> &input) {
        std::unique_ptr<Particles> partc{new Particles};

        if constexpr (Particles::type::arraySize == SpaceH::DYNAMICAL)
            partc->resize(particle_num);

        partc->read(input.actbuff, SpaceH::IO_flag::EVOLVED);

        ASSERT_EQ(input.time, partc->time());
        for (size_t i = 0; i < particle_num; ++i) {
            ASSERT_EQ(input.pos[i].x, partc->pos(i).x);
            ASSERT_EQ(input.pos[i].y, partc->pos(i).y);
            ASSERT_EQ(input.pos[i].z, partc->pos(i).z);
            ASSERT_EQ(input.vel[i].x, partc->vel(i).x);
            ASSERT_EQ(input.vel[i].y, partc->vel(i).y);
            ASSERT_EQ(input.vel[i].z, partc->vel(i).z);
        }
    }

    template<typename Particles, size_t particle_num>
    void test_particle_std_write(const ParticleTestInput<particle_num, typename Particles::Scalar> &input) {
        std::unique_ptr<Particles> partc{new Particles};

        if constexpr (Particles::type::arraySize == SpaceH::DYNAMICAL)
            partc->resize(particle_num);

        partc->read(input.stdbuff, SpaceH::IO_flag::STD);

        typename Particles::ScalarBuffer out;

        partc->write(out, SpaceH::IO_flag::STD);

        size_t size = input.stdbuff.size();
        for (size_t i = 0; i < particle_num; ++i) {
            ASSERT_EQ(input.stdbuff[i], out[i]);
        }
    }

    template<typename Particles, size_t particle_num>
    void test_particle_act_write(const ParticleTestInput<particle_num, typename Particles::Scalar> &input) {
        std::unique_ptr<Particles> partc{new Particles};
        if constexpr (Particles::type::arraySize == SpaceH::DYNAMICAL)
            partc->resize(particle_num);

        partc->read(input.actbuff, SpaceH::IO_flag::EVOLVED);
        typename Particles::ScalarBuffer out;
        partc->write(out, SpaceH::IO_flag::EVOLVED);

        size_t size = input.stdbuff.size();
        for (size_t i = 0; i < particle_num; ++i) {
            ASSERT_EQ(input.actbuff[i], out[i]);
        }
    }
}
#endif