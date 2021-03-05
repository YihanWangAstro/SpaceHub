
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;  // to save writing unit::
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    Particle star{1 * Ms};
    Particle p1{1 * Me};  // earth mass planet
    Particle p2{1 * Mj};  // Jupiter mass planet
    Particle p3{0.3 * Mj};
    Particle p4{0.1 * Me};
    
    /*--------------------------------------------------New-----------------------------------------------------------*/
    auto orb1 = orbit::EllipOrbit(star.mass, p1.mass, 1 * AU, 0.01, 1 * deg, 2 * deg, 3 * deg, 4 * deg);

    auto orb2 = orbit::EllipOrbit(star.mass, p2.mass, 5 * AU, 0.1, 6 * deg, 7 * deg, 8 * deg, 9 * deg);

    auto orb3 = orbit::EllipOrbit(star.mass, p3.mass, 20 * AU, 0.04, 8 * deg, 7 * deg, 6 * deg, 5 * deg);

    auto orb4 = orbit::EllipOrbit(star.mass, p4.mass, 60 * AU, 0.2, 4 * deg, 3 * deg, 2 * deg, 1 * deg);

    orbit::move_particles(orb1, p1);

    orbit::move_particles(orb2, p2);

    orbit::move_particles(orb3, p3);

    orbit::move_particles(orb4, p4);

    orbit::move_to_COM_frame(star, p1, p2, p3, p4);
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver solver{0, star, p1, p2, p3, p4};

    Solver::RunArgs args;

    args.add_stop_condition(1000 * unit::year);

    solver.run(args);

    print(std::cout, "simulation of planetary system without any output complete!\n");

    return 0;
}
