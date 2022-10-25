
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace force;
using f = Interactions<NewtonianGrav, StaticGasDrag>;  // add static drag force

using Solver = methods::DefaultMethod<f, particles::DragParticles>;  // use the corresponding drag particles
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    //(mass, sub_sonic_coef, local_sound_speed)
    Particle star{1_Ms, 0, 2};   // sub_sonic_coef = 0, local sound speed = 2;no drag force on this particle
    Particle p1{1_Mj, 1, 100};   // sub_sonic_coef = 1; local sound speed = 100, subsonic drag;
    Particle p2{1_Mj, 2, 0.01};  // sub_sonic_coef = 2; local sound speed = 0.01, supersonic drag;
    Particle p3{1_Mj, 0, 2};     // sub_sonic_coef = 0, local sound speed = 2;no drag force on this particle
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

    args.add_stop_condition(500_year);

    args.add_operation(DefaultWriter("static-drag.txt"));

    solver.run(args);

    print(std::cout, "simulation with static drag force complete!\n");

    return 0;
}
