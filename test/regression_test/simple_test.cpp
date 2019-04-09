#include "../../src/spaceHub.h"
#include <vector>
#include <array>
#include <iomanip>

using namespace SpaceH;
using namespace SpaceH::OdeIterator;
using namespace SpaceH::Integrator;
using namespace SpaceH::Orbit;
using namespace Unit;
using scalar = double;
using type = Types<scalar, std::vector>;

int main(int argc, char **argv) {
    using force = Interact::NewtonianGrav;

    using particles = SoAPointParticles<type>;

    //using particles = SoAFiniteSizeParticles<type>;

    //using sys = SimpleSystem<particles, force>;

    //using sys = RegularizedSystem <particles, force, ReguType::logH>;

    using sys = ChainSystem <particles, force>;

    using iter = ConstOdeIterator<symplectic2th>;

    using simulation = Solver<sys, iter>;

    using particle = typename simulation::Particle;

    particle sun{M_SOLAR}, earth{M_EARTH}, moon{M_MOON};

    auto moon_orbit = Kepler{moon.mass + earth.mass, semi_latus_rectum(384748*KM , 0.0549006), 0.0549006, 1.543 * DEG, thermal, thermal, thermal};

    auto earth_orbit = Kepler{sun.mass + earth.mass, semi_latus_rectum(AU, 0.0167086), 0.0167086, 7.155 * DEG, 174.9 * DEG, 288.1 * DEG, thermal};

    move_particles(moon_orbit, moon);

    move_particles(earth_orbit, earth, moon);

    move_to_com_coord(sun, earth, moon);

    //print(std::cout, sun,'\n',earth, '\n', moon,'\n',distance(sun.pos, earth.pos), '\n', distance(earth.pos, moon.pos));
    simulation::RunArgs args;

    args.add_pre_step_option(ArgsCallBack::DefaultWriter("solar.dat", 0,  100 * YEAR));

    args.add_stop_condition(100 * YEAR);

    simulation nbody{0, sun, earth, moon};

    nbody.run(args);

    return 0;
}
