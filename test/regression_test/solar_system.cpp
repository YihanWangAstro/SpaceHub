//
// Created by root on 10/9/19.
//
#include "rtest_samples.hpp"
#include "../../src/spaceHub.hpp"

USING_NAMESPACE_ALL;

template<typename simulation>
void run(std::string const &sim_type) {
  auto earth_sys = earth_system<simulation>();

  basic_error_test<simulation>("earth-system-" + sim_type, 100_year, earth_sys);

  auto[rtol, error] = error_scale<simulation>(1e-15, 1e-9, 100_year, earth_sys);

  std::fstream err_stream{"earth-system-" + sim_type + ".scale", std::ios::out};

  err_stream << rtol << '\n' << error;
}

int main(int argc, char **argv) {
  using type = Types<precise_d >;

  using force = interactions::NewtonianGrav;

  using particles = PointParticles<type>;

  using sim_sys = SimpleSystem<particles, force>;

  using regu_sys = RegularizedSystem<particles, force, ReguType::LogH>;

  using chain_sys = ChainSystem<particles, force>;

  using arch_sys = ARchainSystem<particles, force, ReguType::LogH>;

  //using iter = ConstOdeIterator<Symplectic2nd>;

  using iter = BurlishStoer<double, WorstOffender, PIDController>;

  run<Simulator<sim_sys, iter>>("sim");

  run<Simulator<regu_sys, iter>>("regu");

  run<Simulator<chain_sys, iter>>("chain");

  run<Simulator<arch_sys, iter>>("arch");

  return 0;
}

