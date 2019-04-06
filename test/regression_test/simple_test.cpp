#include "../../src/spaceHub.h"
#include <vector>
#include <array>
#include <iomanip>

using namespace SpaceH;
using namespace SpaceH::OdeIterator;
using namespace SpaceH::Integrator;
using namespace Unit;
using scalar = double;
using type = Types<scalar, std::vector>;

int main(int argc, char **argv) {

    using force = Interact::NewtonianGrav;

    using particles = SoAPointParticles<type>;

    //using particles = SoAFiniteSizeParticles<type>;

    using sys = SimpleSystem<particles, force>;

    //using sys = RegularizedSystem <particles, force, ReguType::logH>;

    //using sys = ChainSystem <particles, force>;

    using iter = ConstOdeIterator<symplectic2th>;

    using simulation = Solver<sys, iter>;

    using particle = typename simulation::Particle;

    particle sun{1 * M_SOLAR}, earth{3.003e-6 * M_SOLAR};

    auto ecc = 0.0167086;
    auto p = Orbit::semi_latus_rectum(AU, ecc);
    Orbit::move_particles(Orbit::Kepler{sun.mass + earth.mass, p, ecc, 7.155 * DEG, 174.9 * DEG, 288.1 * DEG, Orbit::thermal}, earth);
    Orbit::move_to_com_coord(sun, earth);

    simulation::RunArgs args;

    args.add_pre_step_option(ArgsCallBack::DefaultWriter(Tools::auto_name(), 0, 1000 * YEAR));
    args.add_stop_condition(1000 * Unit::YEAR);

    simulation nbody{0, sun, earth};
    nbody.run(args);

    return 0;
}
