
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;

/*--------------------------------------------------New-----------------------------------------------------------*/
template <typename SolverType>
void run_simulation(std::string const& method_name) {
    using Particle = typename SolverType::Particle;  // get particle type from the solver

    Particle p1{1_Ms};
    Particle p2{1_Ms};
    Particle p3{1_Ms};

    auto inner_orb = orbit::Elliptic(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    auto outer_orb = orbit::Elliptic(orbit::M_tot(p1, p2), p3.mass, 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    orbit::move_particles(outer_orb, p3);

    orbit::move_to_COM_frame(p1, p2, p3);

    SolverType solver{0, p1, p2, p3};

    typename SolverType::RunArgs args;  // get run argument type from the solver

    args.add_stop_condition(1000_year);

    args.add_operation(DefaultWriter(method_name + ".txt"));

    solver.run(args);

    print(std::cout, "simulation with method " + method_name + " complete!\n");
}
/*----------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
    run_simulation<methods::AR_Radau_Plus<>>("AR-Radau+");

    run_simulation<methods::AR_Chain<>>("AR-Chain");

    run_simulation<methods::AR_BS_Plus<>>("AR-BS+");

    return 0;
}
