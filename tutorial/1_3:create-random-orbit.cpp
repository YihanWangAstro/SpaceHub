
#include "../src/spaceHub.hpp"
using namespace space;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    Particle p1{1 * unit::Ms};  // 1 solar mass particle rest at origin(x=y=z=vx=vy=vz=default 0)
    Particle p2{1 * unit::Me};  // 1 earth mass particle rest at origin

    // use place holder orbit::isotherm to generate the inclination, longitude_of_ascending_node,
    // argument_of_periapsis and true_anomaly isothermaly, i.e. cos(i)~[-1, 1], omega~[-pi,pi], Omega~[-pi,pi] and
    // coresponding mean anomaly of true anomaly ~[-pi,pi]
    auto orb = orbit::EllipOrbit(p1.mass, p2.mass, 5 * unit::AU, 0.0, orbit::isotherm, orbit::isotherm, orbit::isotherm,
                                 orbit::isotherm);

    // you can also fix some angles
    // auto orb = orbit::EllipOrbit(p1.mass, p2.mass, 5 * unit::AU, 0.0, 0.0_deg, orbit::isotherm,
    // 2_deg, orbit::isotherm);
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
