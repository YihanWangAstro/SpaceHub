
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
// using the Newtonian gravity + first order Post-Newtonian correction
using f = force::Interactions<force::NewtonianGrav, force::PN1>;
// pass the force to solver
using Solver = methods::DefaultMethod<f>;// DefaultMethod<> = DefaultMethod<force::Interactions<force::NewtonianGrav>>
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    Particle p1{1_Ms};
    Particle p2{1_Ms};
    Particle p3{1_Ms};

    auto inner_orb = orbit::Elliptic(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    auto outer_orb = orbit::Elliptic(orbit::M_tot(p1, p2), p3.mass, 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    orbit::move_particles(outer_orb, p3);

    orbit::move_to_COM_frame(p1, p2, p3);

    Solver solver{0, p1, p2, p3};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    args.add_operation(DefaultWriter("Post-Newtonian.txt"));

    solver.run(args);

    print(std::cout, "simulation with first order Post-Newtonian term complete!\n");

    return 0;
}
