//
// Created by yihan on 4/17/19.
//
#include <iomanip>
#include "spaceHub.h"
using namespace space;
int main() {
    using Simulation = DefaultSolver<>;
    using Particle = Simulation::Particle;

    Particle m1{unit::m_solar}, m2{0.5*unit::m_solar}, m3{0.05*unit::m_solar};

    auto a1 = 0.5*unit::au;
    auto a2 = 5*unit::au;

    orbit::move_particles_to(orbit::Kepler(m1.mass + m2.mass, a1, 0, 25.01*unit::deg, 0, 90*unit::deg, orbit::thermal), m2);

    orbit::move_to_com_coord(m1, m2);

    orbit::move_particles_to(orbit::Kepler(m1.mass + m2.mass + m3.mass, a2, 0, -64.99*unit::deg, 0, 0, orbit::thermal), m3);

    orbit::move_to_com_coord(m1, m2, m3);

    Simulation kozai_test{0, m1, m2, m3};

    Simulation::RunArgs args;

    auto end_time = 3e4*unit::year;

    std::ofstream output{"kozai.txt"};

    output << std::setprecision(16);

    auto writer = [&](auto& ptc){ output << calc::calc_total_energy(ptc) << ' ' << ptc << '\n';};

    args.add_pre_step_option(argsOpt::TimeSlice(writer,0, end_time));

    args.add_stop_condition(end_time);

    kozai_test.run(args);
}