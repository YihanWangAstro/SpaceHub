//
// Created by yihan on 4/17/19.
//
#include <iomanip>
#include "../../src/spaceHub.h"
using namespace space;
using namespace orbit;
using namespace argsOpt;
int main() {
    using Simulation = DefaultSolver<>;
    using Particle = Simulation::Particle;

    Particle m1{unit::m_solar}, m2{0.5*unit::m_solar}, m3{0.05*unit::m_solar};

    auto a1 = 0.5*unit::au;
    auto a2 = 5*unit::au;

    move_particles_to(Kepler(total_mass(m1, m2), a1, 0, 25.01*unit::deg, 0, 90*unit::deg, thermal), m2);

    move_to_com_coord(m1, m2);

    move_particles_to(Kepler(total_mass(m1, m2, m3), a2, 0, -64.99*unit::deg, 0, 0, thermal), m3);

    move_to_com_coord(m1, m2, m3);

    Simulation kozai_test{0, m1, m2, m3};

    Simulation::RunArgs args;

    auto end_time = 3e4*unit::year;

    std::ofstream output{"kozai.txt"};

    output << std::setprecision(16);

    auto writer = [&](auto& ptc){ output << calc::calc_total_energy(ptc) << ' ' << ptc << '\n';};

    args.add_pre_step_operation(TimeSlice(writer, 0, end_time));

    args.add_stop_condition(end_time);

    kozai_test.run(args);
}