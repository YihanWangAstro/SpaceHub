
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    Particle p1{1_Ms};  // same as 1*Ms
    Particle p2{1_Ms};
    Particle p3{1_Ms};
    Particle p4{1_Ms};

    /*--------------------------------------------------New-----------------------------------------------------------*/
    auto inner_orb1 = orbit::Elliptic(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    orbit::move_particles(inner_orb1, p2);

    orbit::move_to_COM_frame(p1, p2);

    // now the centre of mass of binary (p1, p2) is at origin

    auto inner_orb2 = orbit::Elliptic(p3.mass, p4.mass, 4_AU, 0.1, 10_deg, 2.4_deg, 8.3_deg, 4.8_deg);

    orbit::move_particles(inner_orb2, p4);

    orbit::move_to_COM_frame(p3, p4);

    // now the centre of mass of binary (p3, p4) is at origin

    auto outer_orb =
        orbit::Elliptic(orbit::M_tot(p1, p2), orbit::M_tot(p3, p4), 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(outer_orb, p3, p4);  // move the centre of mass of binary (p3, p4) to the outer orbit. then a
                                               // quadrupole formed ((p1,p2),(p3,p4))

    orbit::move_to_COM_frame(p1, p2, p3, p4);
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver solver{0, p1, p2, p3, p4};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    solver.run(args);

    print(std::cout, "simulation of quadrupole without any output complete!\n");

    return 0;
}
