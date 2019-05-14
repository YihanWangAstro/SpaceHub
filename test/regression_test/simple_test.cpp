#include "../../src/spaceHub.h"
#include <vector>
#include <array>
#include <iomanip>

using namespace space;
using namespace space::odeIterator;
using namespace space::integrator;
using namespace space::orbit;
using namespace unit;
using scalar = double;
using type = Types<scalar, std::vector>;

int main(int argc, char **argv) {
    using force = interactions::NewtonianGrav;

    using particles = PointParticles<type>;

    //using particles = SoAFiniteSizeParticles<type>;

    //using sys = SimpleSystem<particles, force>;

    //using sys = RegularizedSystem <particles, force, ReguType::logH>;

    //using sys = ChainSystem <particles, force>;

    using sys = ARchainSystem <particles, force, ReguType::TTL>;

    //using iter = ConstOdeIterator<symplectic2th>;

    using iter = BSIterator<double>;

    using simulation = Simulator<sys, iter>;

    using particle = typename simulation::Particle;

    particle sun{m_solar}, earth{m_earth}, moon{m_moon};

    auto moon_orbit = Kepler{earth.mass, moon.mass, semi_latus_rectum(384748 * km, 0.0549006), 0.0549006, 1.543 * deg,
                             thermal, thermal, thermal};

    auto earth_orbit = Kepler{sun.mass, earth.mass, semi_latus_rectum(au, 0.0167086), 0.0167086, 7.155 * deg,
                              174.9 * deg, 288.1 * deg, thermal};

    move_particles_to(moon_orbit, moon);

    move_particles_to(earth_orbit, earth, moon);

    move_to_com_coord(sun, earth, moon);

    simulation::RunArgs args;

    std::ofstream eng_file("solar.eng");

    eng_file << std::setprecision(16);

    auto end_time = 10*year;

    args.add_pre_step_operation(argsOpt::TimeSlice(argsOpt::DefaultWriter("solar.dat"), 0, end_time));

    args.add_pre_step_operation(
            argsOpt::TimeSlice([&](auto &ptc) { eng_file << calc::calc_total_energy(ptc) << '\n'; }, 0, end_time));

    args.add_stop_condition(end_time);

    simulation nbody{0, sun, earth, moon};

    nbody.run(args);

    return 0;
}
