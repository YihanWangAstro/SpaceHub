//
// Created by yihan on 4/17/19.
//
#include <iomanip>
#include "../../src/spaceHub.hpp"

USING_NAMESPACE_ALL

int main() {
    using Particle = typename DefaultSolver::Particle;

    Particle  m1{1_Ms}, m2{0.5_Ms}, m3{0.5_Ms};

  auto a1 = 1_AU;
  auto a2 = 5_AU;

  move_particles(Kepler(m1.mass, m2.mass, a1, 0, 25.01 * unit::deg, 0, 90 * unit::deg, 0.0), m2);

  move_to_COM_frame(m1, m2);

  move_particles(Kepler(M_tot(m1, m2), m3.mass, a2, 0, -64.99 * unit::deg, 0, 0, 0.0), m3);

  move_to_COM_frame(m1, m2, m3);

    DefaultSolver kozai_test{0, m1, m2};

  auto E0 = calc::calc_total_energy(kozai_test.particles());

    DefaultSolver::RunArgs args;

  auto end_time = 3e3 * unit::year;

  std::ofstream output{"kozai.eng"};

  output << std::setprecision(16);

  auto writer = [&](auto &ptc) { output << ptc.time() << ',' << calc::calc_energy_error(ptc, E0) << '\n'; };

  args.add_pre_step_operation(run_operations::TimeSlice(writer, 0, end_time));

  args.add_stop_condition(end_time);

  kozai_test.run(args);
}