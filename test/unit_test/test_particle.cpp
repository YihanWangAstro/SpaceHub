#include "test_particle.h"
#include "kahan_number.h"
#include "particles.h"

TEST(ParticleTest, Resize) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_resize<Particles<TypeSystem<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeSystem<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeSystem<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<Particles<TypeSystem<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ParticleTest, Reserve) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_reserve<Particles<TypeSystem<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeSystem<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeSystem<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<Particles<TypeSystem<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ParticleTest, StdRead) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_read<Particles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_read<Particles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_read<Particles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_read<Particles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_read<Particles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_read<Particles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_read<Particles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_read<Particles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, ActRead) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_act_read<Particles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_act_read<Particles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_act_read<Particles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_act_read<Particles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, StdWrite) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_write<Particles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_write<Particles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_write<Particles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_write<Particles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_write<Particles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_write<Particles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_write<Particles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_write<Particles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ParticleTest, ActWrite) {
    constexpr size_t partc_num = 100000;
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_act_write<Particles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_act_write<Particles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_act_write<Particles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_act_write<Particles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_act_write<Particles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_act_write<Particles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_act_read<Particles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}