//
// Created by yihan on 4/17/19.
//
#include <iomanip>
#include "../../src/spaceHub.hpp"

using namespace space;
using namespace space::ode_iterator;
using namespace space::integrator;
using namespace space::orbit;
using namespace unit;
using scalar = double;
using type = Types<scalar, std::vector>;

int main() {
  using force = interactions::NewtonianGrav;

  using particles = PointParticles<type>;

  //using particles = SoAFiniteSizeParticles<type>;

  //using sys = SimpleSystem<particles, force>;

  using sys = RegularizedSystem<particles, force, ReguType::LogH>;

  //using sys = ChainSystem<particles, force>;

  //using sys = ARchainSystem<particles, force, ReguType::LogH>;

  //using iter = ConstOdeIterator<symplectic2nd>;

  using iter = BurlishStoer<double, WorstOffender, PIDController>;

  using Simulation = Simulator<sys, iter>;

  using Particle = Simulation::Particle;

  Particle m1{unit::m_solar}, m2{0.5 * unit::m_solar}, m3{0.5 * unit::m_solar};

  auto a1 = 0.5 * unit::au;
  auto a2 = 5 * unit::au;

  move_particles_to(Kepler(m1.mass, m2.mass, a1, 0, 25.01 * unit::deg, 0, 90 * unit::deg, 0.0), m2);

  move_to_com_coord(m1, m2);

  move_particles_to(Kepler(total_mass(m1, m2), m3.mass, a2, 0, -64.99 * unit::deg, 0, 0, 0.0), m3);

  move_to_com_coord(m1, m2, m3);

  Simulation kozai_test{0, m1, m2, m3};

  auto E0 = calc::calc_total_energy(kozai_test.particles());

  Simulation::RunArgs args;

  auto end_time = 3e3 * unit::year;

  std::ofstream output{"kozai.eng"};

  output << std::setprecision(16);

  auto writer = [&](auto &ptc) { output << ptc.time() << ',' << calc::calc_energy_error(ptc, E0) << '\n'; };

  args.add_pre_step_operation(run_operations::TimeSlice(writer, 0, end_time));

  args.add_stop_condition(end_time);

  kozai_test.run(args);
}