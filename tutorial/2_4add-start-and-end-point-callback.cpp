
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using Solver = methods::DefaultMethod<>;
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

    /*--------------------------------------------------New-----------------------------------------------------------*/
    auto greeting = [](auto& particles, auto step_size) {
        // greeting
        print(std::cout, "start integration!\n");
    };

    auto notice = [](auto& particles, auto step_size) {
        // notify
        print(std::cout, "integration finished!\n");
    };

    args.add_start_point_operation(greeting);

    args.add_stop_point_operation(notice);
    /*----------------------------------------------------------------------------------------------------------------*/

    solver.run(args);

    print(std::cout, "simulation of adding start/stop conditions complete!\n");

    return 0;
}
