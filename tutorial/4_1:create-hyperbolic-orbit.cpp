
#include "../src/spaceHub.hpp"
using namespace space;
using namespace unit;
using namespace callback;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    Particle p1{1_Ms};
    Particle p2{1_Ms};

    // create hyperbolic orbit (primary mass, secondary mass, v_inf, impact parameter b, inclination,
    // longitude_of_ascending_node,  argument_of_periapsis, start drop off distance between primary and secondary
    // object, incident in/fly away out)
    auto incident_orb =
        orbit::HyperOrbit(p1.mass, p2.mass, 5_kms, 50_AU, 15_deg, 0_deg, 2_deg, 1000_AU, orbit::Hyper::in);

    // move particle p2 to the coresponding position/velocity of hyperbolic orbit orb; (p1 stays remain at origin)
    orbit::move_particles(incident_orb, p2);

    orbit::move_to_COM_frame(p1, p2);
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver solver{0, p1, p2};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    args.add_operation(TimeSlice(DefaultWriter("single-single.txt"), 0., 1000_year, 5000));

    solver.run(args);

    print(std::cout, "simulation of single single scatterint  complete!\n");

    return 0;
}
