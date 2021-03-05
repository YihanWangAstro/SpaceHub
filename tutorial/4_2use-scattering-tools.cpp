
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
/*--------------------------------------------------New-----------------------------------------------------------*/
using Scalar = Solver::Scalar;//using the scalar type of the solver
 /*----------------------------------------------------------------------------------------------------------------*/
int main(int argc, char** argv) {
    size_t n = 5000;  // 5000 times scattering
    Scalar v_inf = 10_kms;
    Scalar b_max = 10_AU;
    Scalar r_start = 100_AU;

    std::fstream file("scattering1+1.txt", std::ios::out);

    for (size_t i = 0; i < n; ++i) {
        Particle p1{1_Ms};
        Particle p2{1_Ms};

        /*--------------------------------------------------New-----------------------------------------------------------*/
        // create the incident orbit (total mass of scattered object, total mass of incident object, velocity at
        // infinity, max impact parameter, relative distance at t = 0) The trajectory of r=inf -> r=r_start is
        // calculated analytically and the incident object will be generated uniformly in the area of \pi*b_max^2 at
        // infinity
        auto orb = scattering::incident_orbit(p1.mass, p2.mass, v_inf, b_max, r_start);

        orbit::move_particles(orb, p2);

        orbit::move_to_COM_frame(p1, p2);
        /*----------------------------------------------------------------------------------------------------------------*/

        Solver solver{0, p1, p2};

        Solver::RunArgs args;

        /*--------------------------------------------------New-----------------------------------------------------------*/
        Scalar t_end = 2 * orbit::time_to_periapsis(p1, p2);  // calculate the time from start position to periapsis

        args.add_stop_condition(t_end);
        /*----------------------------------------------------------------------------------------------------------------*/

        args.add_stop_point_operation([&file, &orb](auto& ptc, auto h) {
            // print the end state of the system into file
            print(file, orb, ',', ptc, '\n');
        });

        solver.run(args);
    }

    print(std::cout, "simulation with scattering tools complete!\n");

    return 0;
}
