
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace orbit;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
using Scalar = Solver::Scalar;

int main(int argc, char** argv) {
    size_t n = 5000;
    Scalar v_inf = 10_kms;
    Scalar b_max = 20_AU;
    Scalar r_start = 500_AU;

    std::fstream file("scattering2+2.txt", std::ios::out);

    for (size_t i = 0; i < n; ++i) {
        Particle p1{1_Ms};
        Particle p2{1_Ms};
        Particle p3{1_Ms};
        Particle p4{1_Ms};

        /*--------------------------------------------------New-----------------------------------------------------------*/
        auto binary_orb1 = EllipOrbit(p1.mass, p2.mass, 5_AU, 0.0, isotherm, isotherm, isotherm, isotherm);

        move_particles(binary_orb1, p2);

        move_to_COM_frame(p1, p2);

        auto binary_orb2 = EllipOrbit(p1.mass, p2.mass, 6_AU, 0.0, isotherm, isotherm, isotherm, isotherm);

        move_particles(binary_orb2, p4);

        move_to_COM_frame(p3, p4);

        auto orb = scattering::incident_orbit(M_tot(p1, p2), M_tot(p3, p4), v_inf, b_max, r_start);

        move_particles(orb, p3, p4);
        /*----------------------------------------------------------------------------------------------------------------*/

        move_to_COM_frame(p1, p2, p3, p4);

        Solver solver{0, p1, p2, p3, p4};

        Solver::RunArgs args;

        Scalar t_end = 2 * time_to_periapsis(orbit::group(p1, p2), orbit::group(p3, p4));

        args.add_stop_condition(t_end);

        args.add_stop_point_operation(
            [&](auto& ptc, auto h) { print(file, binary_orb1, ',', binary_orb2, ',', orb, ',', ptc, '\n'); });

        solver.run(args);
    }

    print(std::cout, "simulation with binary binary complete!\n");

    return 0;
}
