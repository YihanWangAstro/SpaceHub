#include "../../src/spaceHub.h"
#include <vector>
#include <iomanip>

using namespace SpaceH;
using namespace SpaceH::OdeIterator;
using namespace SpaceH::Integrator;

using scalar = double;
using type = Types<scalar, std::vector>;

int main(int argc, char **argv) {

    using force = NewtonianGrav;
    using particles = SoAPointParticles<type>;

    //using particles = SoAFiniteSizeParticles<type>;


    using sys = SimpleSystem<particles, force>;

    //using sys = RegularizedSystem <particles, force, ReguType::logH>;

    //using sys = ChainSystem <particles, force>;

    using iter = ConstOdeIterator<symplectic2th>;

    using simulation = Solver<sys, iter>;

    using particle = typename simulation::Particle;

    std::vector<particle> init;
    init.emplace_back(1, 0, 0, 0, 0, 0, 0);
    init.emplace_back(1e-3, 1, 0, 0, 0, 0, 0);

    simulation nbody{init, 0};

    simulation::RunArgs args;

    args.step_size = 1 * Unit::YEAR;
    auto wtr = ArgsCallBack::DefaultWriter(Tools::auto_name(), 0, 50000 * Unit::YEAR);
    args.add_pre_step_option(wtr);
    args.add_stop_condition(1000 * Unit::YEAR);
    args.add_stop_point_option([&](auto &ptc) { wtr << "!" << ptc; });

    nbody.run(args);

    return 0;
}
