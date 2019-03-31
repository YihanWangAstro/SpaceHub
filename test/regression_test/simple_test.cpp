#include "../../src/spaceHub.h"
#include <vector>
#include <iomanip>
using namespace SpaceH;
using namespace SpaceH::odeIterator;
using namespace SpaceH::integrator;

using scalar = double;
using type = Types<scalar, std::vector>;

int main(int argc, char **argv) {

    using force = NewtonianGrav;
    using particles = SoAPointParticle<type>;

    using sys = SimpleSystem <particles, force>;
    using iter = ConstOdeIterator<symplectic2th>;

    using simulation = Solver<sys,iter>;

    std::vector<PointParticle<scalar>> init;
    init.emplace_back(0,1,0,0,0,0,0,0);
    init.emplace_back(1,1e-3,1,0,0,0,0,0);

    particles ls(init,0.0);
    simulation nbody{init,0};
/*
    std::cout << std::setprecision(16);

    simulation::RunArgs args;

    args.end_time = 50000*Unit::YEAR;
    args.step_size = 1*YEAR;

    nbody.run(args);*/
    return 0;
}
