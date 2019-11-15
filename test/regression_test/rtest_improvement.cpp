
#include <thread>
#include "../../src/spaceHub.hpp"
#include "rtest_samples.hpp"

USING_NAMESPACE_ALL;

template <typename simulation>
void run(std::string const &sim_name) {
  auto earth_sys = earth_system<simulation>();

  basic_error_test<simulation>("improvement-" + sim_name, 1000_year, 1e-13, earth_sys);

  auto [rtol, error] = error_scale<simulation>(1e-14, 1e-8, 1000_year, earth_sys);

  std::fstream err_stream{"improvement-" + sim_name + ".scale", std::ios::out};

  err_stream << rtol << '\n' << error;
}

int main(int argc, char **argv) {
  using type = Types<double_k>;

  using force = interactions::NewtonianGrav;

  using particles = PointParticles<type>;

  using sim_sys = SimpleSystem<particles, force>;

  using regu_sys = RegularizedSystem<particles, force, ReguType::LogH>;

  using arch_sys = ARchainSystem<particles, force, ReguType::LogH>;

  // using iter = ConstOdeIterator<Symplectic2nd>;

  using iter = BurlishStoer<double, WorstOffender, PIDController>;

  run<Simulator<arch_sys, iter>>("arch");

  return 0;
}
