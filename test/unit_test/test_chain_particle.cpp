#include "chain_particles.h"
#include "test_particle.h"
#include "kahan_number.h"

TEST(ChainParticleTest, Resize) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_resize<ChainParticles<TypeClass<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeClass<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeClass<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeClass<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ChainParticleTest, Reserve) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_reserve<ChainParticles<TypeClass<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeClass<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeClass<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeClass<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ChainParticleTest, StdRead) {
    constexpr size_t partc_num = 1000;//chain algorithm cannot be used in large N(very time consuming).
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_read<ChainParticles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ChainParticleTest, StdWrite) {
    constexpr size_t partc_num = 1000;//chain algorithm cannot be used in large N(very time consuming).
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_write<ChainParticles<TypeClass<precise_d, partc_num>>, partc_num>(inputpd);
}
