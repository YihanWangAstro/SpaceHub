
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    Particle p1{1_Ms};  // same as 1*Ms
    Particle p2{1_Ms};
    Particle p3{1_Ms};

    /*--------------------------------------------------New-----------------------------------------------------------*/
    auto inner_orb = orbit::EllipOrbit(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    // use inner binary as the primary object in outer orbit
    auto outer_orb = orbit::EllipOrbit(p1.mass + p2.mass, p3.mass, 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(inner_orb, p2);  // inner binary formed

    orbit::move_to_COM_frame(p1, p2);  // this step is critical

    orbit::move_particles(outer_orb, p3);  // outer binary formed thus hierarchical triple formed

    orbit::move_to_COM_frame(p1, p2, p3);
    // note:for parabolic/hyperbolic orbit creation, move to scattering section.
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver solver{0, p1, p2, p3};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    solver.run(args);

    print(std::cout, "simulation of hierarchical triple without any output complete!\n");

    return 0;
}
