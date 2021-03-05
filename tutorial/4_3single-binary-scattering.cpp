
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using namespace orbit;  // save writing orbit::
/*----------------------------------------------------------------------------------------------------------------*/
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    size_t n = 5000;
    Scalar v_inf = 10_kms;
    Scalar b_max = 10_AU;
    Scalar r_start = 100_AU;

    std::fstream file("scattering2+1.txt", std::ios::out);

    for (size_t i = 0; i < n; ++i) {
        Particle p1{1_Ms};
        Particle p2{1_Ms};
        Particle p3{1_Ms};

        /*--------------------------------------------------New-----------------------------------------------------------*/

        auto binary_orb = Elliptic(p1.mass, p2.mass, 5_AU, 0.0, isotherm, isotherm, isotherm, isotherm);

        move_particles(binary_orb, p2);

        move_to_COM_frame(p1, p2);

        auto orb = scattering::incident_orbit(M_tot(p1, p2), p3.mass, v_inf, b_max, r_start);

        move_particles(orb, p3);

        move_to_COM_frame(p1, p2, p3);
        /*----------------------------------------------------------------------------------------------------------------*/

        Solver solver{0, p1, p2, p3};

        Solver::RunArgs args;

        /*--------------------------------------------------New-----------------------------------------------------------*/
        // orbit::group(p1, p2) is the scattered object, p3 in the incident object
        Scalar t_end = 2 * time_to_periapsis(orbit::group(p1, p2), p3);

        args.add_stop_condition(t_end);
        /*----------------------------------------------------------------------------------------------------------------*/
        // [&] capture all variables in lambda in reference
        args.add_stop_point_operation([&](auto& ptc, auto h) {
            // print the end state of the system into file
            print(file, binary_orb, ',', orb, ',', ptc, '\n');
        });

        solver.run(args);
    }

    print(std::cout, "simulation with binary binary complete!\n");

    return 0;
}
