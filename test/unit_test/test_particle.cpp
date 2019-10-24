//
// Created by 王艺涵 on 10/5/19.
//
#include "../../src/type-class.hpp"
#include "../../src/particles/point-particles.hpp"
#include "../../src/particles/finite-size.hpp"
#include "../../src/rand-generator.hpp"
#include "gtest/gtest.h"
using namespace space;
using namespace math;

namespace UnitTest {

  template<template<class> class TParticles, typename types>
  void test_particles_size(size_t sample_num) {
    using Particles = TParticles<types>;
    using Particle = typename Particles::Particle;

    particle_set::PointParticles<types> particles;

    EXPECT_EQ(particles.number(), 0);

    for (size_t i = 0; i < sample_num; ++i) {
      size_t particle_num = static_cast<size_t >(space::random::Uniform(0, 10000));

      particles.resize(particle_num);

      ASSERT_EQ(particles.number(), particle_num);
    }

    particles.clear();

    ASSERT_EQ(particles.number(), 0);
  }

  template<template<class> class TParticles, typename types>
  void test_particles_capacity(size_t sample_num) {

    using Particles = TParticles<types>;
    using Particle = typename Particles::Particle;

    particle_set::PointParticles<types> particles;

    EXPECT_EQ(particles.capacity(), 0);

    for (size_t i = 0; i < sample_num; ++i) {
      size_t particle_num = static_cast<size_t>(space::random::Uniform(0, 10000));

      particles.reserve(particle_num);

      ASSERT_GE(particles.capacity(), particle_num);
    }
  }

  template<typename types>
  void test_point_particles_ctor(size_t sample_num) {
    using Scalar = typename types::Scalar;
    using Particles = particle_set::PointParticles<types>;
    using Particle = typename Particles::Particle;
    Scalar high = big_value<Scalar>::value;

    std::vector<Particle> container;

    container.reserve(sample_num);
    Scalar t = random::Uniform(-high, high);
    for (size_t i = 0; i < sample_num; ++i) {
      Scalar mass = random::Uniform(-high, high);
      Scalar px = random::Uniform(-high, high);
      Scalar py = random::Uniform(-high, high);
      Scalar pz = random::Uniform(-high, high);
      Scalar vx = random::Uniform(-high, high);
      Scalar vy = random::Uniform(-high, high);
      Scalar vz = random::Uniform(-high, high);
      container.emplace_back(Particle{mass, px, py, pz, vx, vy, vz});
    }

    Particles particles{t, container};

    ASSERT_EQ(particles.time(), t);
    for (size_t i = 0; i < sample_num; ++i) {
      ASSERT_EQ(particles.idn()[i], i);

      ASSERT_EQ(particles.mass()[i], container[i].mass);

      ASSERT_EQ(particles.pos().x[i], container[i].pos.x);

      ASSERT_EQ(particles.pos().y[i], container[i].pos.y);

      ASSERT_EQ(particles.pos().z[i], container[i].pos.z);

      ASSERT_EQ(particles.vel().x[i], container[i].vel.x);

      ASSERT_EQ(particles.vel().y[i], container[i].vel.y);

      ASSERT_EQ(particles.vel().z[i], container[i].vel.z);
    }
  }

  template<typename types>
  void test_point_particles_emplace_back(size_t sample_num) {
    using Scalar = typename types::Scalar;
    using Particles = particle_set::PointParticles<types>;
    using Particle = typename Particles::Particle;
    Scalar high = big_value<Scalar>::value;

    std::vector<Particle> container;

    container.reserve(sample_num);

    for (size_t i = 0; i < sample_num; ++i) {
      Scalar mass = random::Uniform(-high, high);
      Scalar px = random::Uniform(-high, high);
      Scalar py = random::Uniform(-high, high);
      Scalar pz = random::Uniform(-high, high);
      Scalar vx = random::Uniform(-high, high);
      Scalar vy = random::Uniform(-high, high);
      Scalar vz = random::Uniform(-high, high);
      container.emplace_back(Particle{mass, px, py, pz, vx, vy, vz});
    }

    Particles particles;

    particles.reserve(sample_num);

    for (auto const &p: container) {
      particles.emplace_back(p);
    }

    for (size_t i = 0; i < sample_num; ++i) {
      ASSERT_EQ(particles.idn()[i], i);

      ASSERT_EQ(particles.mass()[i], container[i].mass);

      ASSERT_EQ(particles.pos().x[i], container[i].pos.x);

      ASSERT_EQ(particles.pos().y[i], container[i].pos.y);

      ASSERT_EQ(particles.pos().z[i], container[i].pos.z);

      ASSERT_EQ(particles.vel().x[i], container[i].vel.x);

      ASSERT_EQ(particles.vel().y[i], container[i].vel.y);

      ASSERT_EQ(particles.vel().z[i], container[i].vel.z);
    }
  }

  template<typename types>
  void test_finite_particles_ctor(size_t sample_num) {
    using Scalar = typename types::Scalar;
    using Particles = particle_set::SizeParticles<types>;
    using Particle = typename Particles::Particle;
    Scalar high = big_value<Scalar>::value;

    std::vector<Particle> container;

    container.reserve(sample_num);
    Scalar t = random::Uniform(-high, high);
    for (size_t i = 0; i < sample_num; ++i) {
      Scalar mass = random::Uniform(-high, high);
      Scalar r = random::Uniform(-high, high);
      Scalar px = random::Uniform(-high, high);
      Scalar py = random::Uniform(-high, high);
      Scalar pz = random::Uniform(-high, high);
      Scalar vx = random::Uniform(-high, high);
      Scalar vy = random::Uniform(-high, high);
      Scalar vz = random::Uniform(-high, high);
      container.emplace_back(Particle{mass, r, px, py, pz, vx, vy, vz});
    }

    Particles particles{t, container};

    ASSERT_EQ(particles.time(), t);
    for (size_t i = 0; i < sample_num; ++i) {
      ASSERT_EQ(particles.idn()[i], i);

      ASSERT_EQ(particles.mass()[i], container[i].mass);

      ASSERT_EQ(particles.radius()[i], container[i].radius);

      ASSERT_EQ(particles.pos().x[i], container[i].pos.x);

      ASSERT_EQ(particles.pos().y[i], container[i].pos.y);

      ASSERT_EQ(particles.pos().z[i], container[i].pos.z);

      ASSERT_EQ(particles.vel().x[i], container[i].vel.x);

      ASSERT_EQ(particles.vel().y[i], container[i].vel.y);

      ASSERT_EQ(particles.vel().z[i], container[i].vel.z);
    }
  }

  template<typename types>
  void test_finite_particles_emplace_back(size_t sample_num) {
    using Scalar = typename types::Scalar;
    using Particles = particle_set::SizeParticles<types>;
    using Particle = typename Particles::Particle;
    Scalar high = big_value<Scalar>::value;

    std::vector<Particle> container;

    container.reserve(sample_num);

    for (size_t i = 0; i < sample_num; ++i) {
      Scalar mass = random::Uniform(-high, high);
      Scalar r = random::Uniform(-high, high);
      Scalar px = random::Uniform(-high, high);
      Scalar py = random::Uniform(-high, high);
      Scalar pz = random::Uniform(-high, high);
      Scalar vx = random::Uniform(-high, high);
      Scalar vy = random::Uniform(-high, high);
      Scalar vz = random::Uniform(-high, high);
      container.emplace_back(Particle{mass, r, px, py, pz, vx, vy, vz});
    }
    Particles particles;

    particles.reserve(sample_num);

    for (auto const &p: container) {
      particles.emplace_back(p);
    }
    for (size_t i = 0; i < sample_num; ++i) {
      ASSERT_EQ(particles.idn()[i], i);

      ASSERT_EQ(particles.mass()[i], container[i].mass);

      ASSERT_EQ(particles.radius()[i], container[i].radius);

      ASSERT_EQ(particles.pos().x[i], container[i].pos.x);

      ASSERT_EQ(particles.pos().y[i], container[i].pos.y);

      ASSERT_EQ(particles.pos().z[i], container[i].pos.z);

      ASSERT_EQ(particles.vel().x[i], container[i].vel.x);

      ASSERT_EQ(particles.vel().y[i], container[i].vel.y);

      ASSERT_EQ(particles.vel().z[i], container[i].vel.z);
    }
  }
}
TEST(ParticleTest, Point_Size) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_particles_size<space::particle_set::PointParticles, Types<double>>(sample_num);
  UnitTest::test_particles_size<space::particle_set::PointParticles, Types<float>>(sample_num);
}

TEST(ParticleTest, Point_Capacity) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_particles_capacity<space::particle_set::PointParticles, Types<double>>(sample_num);
  UnitTest::test_particles_capacity<space::particle_set::PointParticles, Types<float>>(sample_num);
}

TEST(ParticleTest, Point_Ctor) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_point_particles_ctor<space::Types<double>>(sample_num);
  UnitTest::test_point_particles_ctor<space::Types<float>>(sample_num);
}

TEST(ParticleTest, Point_emplace) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_point_particles_emplace_back<space::Types<double>>(sample_num);
  UnitTest::test_point_particles_emplace_back<space::Types<float>>(sample_num);
}

TEST(ParticleTest, Finite_Size) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_particles_size<space::particle_set::SizeParticles, Types<double>>(sample_num);
  UnitTest::test_particles_size<space::particle_set::SizeParticles, Types<float>>(sample_num);
}

TEST(ParticleTest, Finite_Capacity) {
  constexpr size_t sample_num = 100000;
  UnitTest::test_particles_capacity<space::particle_set::SizeParticles, Types<double>>(sample_num);
  UnitTest::test_particles_capacity<space::particle_set::SizeParticles, Types<float>>(sample_num);
}

TEST(ParticleTest, Finite_Ctor) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_finite_particles_ctor<space::Types<double>>(sample_num);
  UnitTest::test_finite_particles_ctor<space::Types<float>>(sample_num);
}

TEST(ParticleTest, Finite_emplace) {
  constexpr size_t sample_num = 1000000;
  UnitTest::test_finite_particles_emplace_back<space::Types<double>>(sample_num);
  UnitTest::test_finite_particles_emplace_back<space::Types<float>>(sample_num);
}