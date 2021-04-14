
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;
using f = Interactions<NewtonianGrav, Tidal>;  // add static tidal force

using Solver = methods::DefaultMethod<f, particles::TideParticles>;  // use the corresponding tidal particles
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    //(mass, radius, apsidal motion constant, lag time)
    Particle star{1_Ms, 1_Rs, 0, 0};  // apsidal motion const = 0; lag time = 0; thus, no tide on this particle
    Particle p1{1_Mj, 1_Rj, 0.25,
                10_sec};                   // apsidal motion const = 0.25; lag time = 10 seconds;  Tide on this particle
    Particle p2{1_Mj, 1_Rj, 0, 0};         // apsidal motion const = 0; lag time = 0; thus, no tide on this particle
    Particle p3{1_Mj, 1_Rj, 0.3, 20_sec};  // apsidal motion const = 0.3; lag time = 20 seconds;  Tide on this particle
    /*----------------------------------------------------------------------------------------------------------------*/

    auto orb1 = orbit::Elliptic(star.mass, p1.mass, 1 * AU, 0.01, 1 * deg, 2 * deg, 3 * deg, 4 * deg);

    auto orb2 = orbit::Elliptic(star.mass, p2.mass, 5 * AU, 0.1, 6 * deg, 7 * deg, 8 * deg, 9 * deg);

    auto orb3 = orbit::Elliptic(star.mass, p3.mass, 20 * AU, 0.04, 8 * deg, 7 * deg, 6 * deg, 5 * deg);

    orbit::move_particles(orb1, p1);

    orbit::move_particles(orb2, p2);

    orbit::move_particles(orb3, p3);

    orbit::move_to_COM_frame(star, p1, p2, p3);

    Solver solver{0, star, p1, p2, p3};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    args.add_operation(DefaultWriter("static-tidal.txt"));

    solver.run(args);

    print(std::cout, "simulation with static tidal force complete!\n");

    return 0;
}
