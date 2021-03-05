
#include "../src/spaceHub.hpp"
using namespace hub;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    Particle p1{1 * unit::Ms};  // 1 solar mass particle rest at origin(x=y=z=vx=vy=vz=default 0)
    Particle p2{1 * unit::Me};  // 1 earth mass particle rest at origin

    // create an ellipse orbit with (primary mass, secondary mass, semi-major axis, eccentricity, inclination,
    // longitude_of_ascending_node,  argument_of_periapsis, true_anomaly)
    auto orb = orbit::EllipOrbit(p1.mass, p2.mass, 1 * unit::AU, 0.2, 0.4 * unit::deg, 25 * unit::deg, 35 * unit::deg,
                                 134 * unit::deg);

    // move particle p2 to the coresponding position/velocity of orbit orb; (p1 stays remain at origin)
    orbit::move_particles(orb, p2);
    /*----------------------------------------------------------------------------------------------------------------*/

    orbit::move_to_COM_frame(p1, p2);

    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    args.add_stop_condition(1000 * unit::year);

    solver.run(args);

    print(std::cout, "simulation of two body keplerian orbit without any output complete!\n");

    return 0;
}
