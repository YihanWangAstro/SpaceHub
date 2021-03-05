
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
    int dummy_variable = 0;

    // capture the variable 'dummy_variable' by reference
    auto print_dummy = [&dummy_variable](auto& particles, auto step_size) {
        // print the value of dummy_variable
        print(std::cout, "dummy_variable:", dummy_variable, '\n');
        dummy_variable += 1;
    };

    args.add_operation(StepSlice(print_dummy, 1000));

    solver.run(args);
    /*----------------------------------------------------------------------------------------------------------------*/

    print(std::cout, "simulation of variable capture complete!\n");

    return 0;
}
