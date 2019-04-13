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
    using force = Interact::NewtonianGrav;

    using particles = SoAPointParticles<type>;

    //using particles = SoAFiniteSizeParticles<type>;

    //using sys = SimpleSystem<particles, force>;

    //using sys = RegularizedSystem <particles, force, ReguType::logH>;

    using sys = ChainSystem <particles, force>;

    //using sys = ARchainSystem <particles, force, ReguType::TTL>;

    //using iter = ConstOdeIterator<symplectic2th>;

    using iter = BSIterator<double>;

    using simulation = Solver<sys, iter>;

    using particle = typename simulation::Particle;

    particle sun{m_solar}, earth{m_earth}, moon{m_moon};

    auto moon_orbit = Kepler{moon.mass + earth.mass, semi_latus_rectum(384748*km , 0.0549006), 0.0549006, 1.543 * deg, thermal, thermal, thermal};

    auto earth_orbit = Kepler{sun.mass + earth.mass, semi_latus_rectum(au, 0.0167086), 0.0167086, 7.155 * deg, 174.9 * deg, 288.1 * deg, thermal};

    move_particles_to(moon_orbit, moon);

    move_particles_to(earth_orbit, earth, moon);

    move_to_com_coord(sun, earth, moon);

    //print(std::cout, sun,'\n',earth, '\n', moon,'\n',distance(sun.pos, earth.pos), '\n', distance(earth.pos, moon.pos));
    simulation::RunArgs args;

    std::ofstream eng_file("solar.eng");

    eng_file << std::setprecision(16);

    args.add_pre_step_option(ArgsCallBack::DefaultWriter("solar.dat", 0,  10 * year));

    args.add_pre_step_option([&](auto& ptc){eng_file << calc::calc_total_energy(ptc) << '\n';});

    args.add_stop_condition(10* year);

    simulation nbody{0, sun, earth, moon};


    nbody.run(args);

    return 0;
}
