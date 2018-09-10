#include "test_particle.h"
#include "kahan_number.h"
#include "particles.h"

TEST(ParticleTest, Resize) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_resize<Particles<TypeClass<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeClass<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeClass<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeClass<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ParticleTest, Reserve) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_reserve<Particles<TypeClass<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeClass<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeClass<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeClass<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ParticleTest, StdRead) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_read<Particles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_read<Particles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_read<Particles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_read<Particles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_read<Particles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_read<Particles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_read<Particles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_read<Particles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, ActRead) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_act_read<Particles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_act_read<Particles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_act_read<Particles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_act_read<Particles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, StdWrite) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_write<Particles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_write<Particles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_write<Particles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_write<Particles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_write<Particles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_write<Particles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_write<Particles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_write<Particles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, ActWrite) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_act_write<Particles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_act_write<Particles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_act_write<Particles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_act_write<Particles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_act_write<Particles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_act_write<Particles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_act_read<Particles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}