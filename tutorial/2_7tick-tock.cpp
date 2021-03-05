
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

    auto inner_orb = orbit::EllipOrbit(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    auto outer_orb = orbit::EllipOrbit(orbit::M_tot(p1, p2), p3.mass, 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    orbit::move_particles(outer_orb, p3);

    orbit::move_to_COM_frame(p1, p2, p3);

    Solver solver{0, p1, p2, p3};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    /*--------------------------------------------------New-----------------------------------------------------------*/
    tools::Timer timer;  // use the timer in SpaceHub

    // capture the variable 'timer' by reference
    auto tick_tock = [&timer](auto& particles, auto step_size) -> bool {
        // if wall time get by 'get_time' is larger than 2 seconds
        if (timer.get_time() > 2) {
            return true;
        } else {
            return false;
        }
    };

    // set the wall time to 2 seconds; check the condition every 1000 steps.
    args.add_stop_condition(StepSlice(tick_tock, 1000));

    timer.start();  // start the timer

    solver.run(args);  // tick_tok will be invoked during the run.
    /*----------------------------------------------------------------------------------------------------------------*/

    print(std::cout, "simulation of setting wall time complete!\n");

    return 0;
}
