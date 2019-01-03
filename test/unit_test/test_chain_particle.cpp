#include "chain_particles.h"
#include "test_particle.h"
#include "kahan_number.h"

TEST(ChainParticleTest, Resize) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_resize<ChainParticles<TypeSystem<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeSystem<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeSystem<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_resize<ChainParticles<TypeSystem<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ChainParticleTest, Reserve) {
    constexpr size_t particle_num = 100000;
    using namespace SpaceH;
    UnitTest::test_particle_reserve<ChainParticles<TypeSystem<float, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeSystem<double, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeSystem<precise_f, DYNAMICAL>>>(particle_num);
    UnitTest::test_particle_reserve<ChainParticles<TypeSystem<precise_d, DYNAMICAL>>>(particle_num);
}

TEST(ChainParticleTest, StdRead) {
    constexpr size_t partc_num = 1000;//chain algorithm cannot be used in large N(very time consuming).
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_read<ChainParticles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}

TEST(ChainParticleTest, StdWrite) {
    constexpr size_t partc_num = 1000;//chain algorithm cannot be used in large N(very time consuming).
    using namespace SpaceH;
    UnitTest::ParticleTestInput<partc_num, float> inputf;
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<float, DYNAMICAL>>, partc_num>(inputf);
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<float, partc_num>>, partc_num>(inputf);

    UnitTest::ParticleTestInput<partc_num, double> inputd;
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<double, DYNAMICAL>>, partc_num>(inputd);
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<double, partc_num>>, partc_num>(inputd);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_f> inputpf;
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<precise_f, DYNAMICAL>>, partc_num>(inputpf);
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<precise_f, partc_num>>, partc_num>(inputpf);

    UnitTest::ParticleTestInput<partc_num, SpaceH::precise_d> inputpd;
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<precise_d, DYNAMICAL>>, partc_num>(inputpd);
    UnitTest::test_particle_std_write<ChainParticles<TypeSystem<precise_d, partc_num>>, partc_num>(inputpd);
}
