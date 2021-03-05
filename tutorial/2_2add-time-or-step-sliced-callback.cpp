
#include "../src/spaceHub.hpp"
using namespace space;
using namespace unit;
using namespace callback;  // save writing callback::
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
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
    // create a time sliced output writer, the output writer will be invoked every dt between [t_start, t_end]
    // the syntax is 'TimeSlice(any_callback, t_start, t_end, slice_number)' where dt = (t_end - t_start)/slice_number
    auto t_writer = TimeSlice(DefaultWriter("time-sliced-hierarchical.txt"), 0.0, 1000_year, 1000);

    // create a step sliced output writer, the output writer will be invoked every N steps
    // the syntax is 'StepSlice(any_callback, step_interval)' where N = step_interval
    auto s_writer = StepSlice(DefaultWriter("step-sliced-hierarchical.txt"), 10);

    args.add_operation(t_writer);

    args.add_operation(s_writer);

    // now we have two output files, one is time sliced, the other is step sliced.
    /*----------------------------------------------------------------------------------------------------------------*/

    solver.run(args);

    print(std::cout, "simulation of time sliced and step sliced output complete!\n");

    return 0;
}
