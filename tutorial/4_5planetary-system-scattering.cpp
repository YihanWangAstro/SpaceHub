
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
    Scalar b_max = 50_AU;
    Scalar r_start = 500_AU;

    std::fstream file("scattering-planetary.txt", std::ios::out);

    for (size_t i = 0; i < n; ++i) {
        Particle star1{1_Ms};
        Particle star2{1_Ms};
        Particle p1{1_Me};
        Particle p2{1_Mj};
        Particle p3{1_Me};
        Particle p4{1_Me};

        /*--------------------------------------------------New-----------------------------------------------------------*/
        // create first planetary system
        Scalar inclination1 = acos(random::Uniform(-1.0, 1.0));

        auto orb1 = EllipOrbit(star1.mass, p1.mass, 1_AU, 0.0, inclination1, isotherm, isotherm, isotherm);

        auto orb2 = EllipOrbit(star1.mass, p2.mass, 5_AU, 0.0, inclination1, isotherm, isotherm, isotherm);

        move_particles(orb1, p1);

        move_particles(orb2, p2);

        move_to_COM_frame(star1, p1, p2);
        /*----------------------------------------------------------------------------------------------------------------*/
        // create secondary planetary system

        Scalar inclination2 = acos(random::Uniform(-1.0, 1.0));

        auto orb3 = EllipOrbit(star2.mass, p3.mass, 1_AU, 0.0, inclination2, isotherm, isotherm, isotherm);

        auto orb4 = EllipOrbit(star2.mass, p4.mass, 5_AU, 0.0, inclination2, isotherm, isotherm, isotherm);

        move_particles(orb3, p3);

        move_particles(orb4, p4);

        move_to_COM_frame(star2, p3, p4);
        /*----------------------------------------------------------------------------------------------------------------*/
        // create incident orbit
        auto orb = scattering::incident_orbit(M_tot(star1, p1, p2), M_tot(star2, p3, p4), v_inf, b_max, r_start);

        move_particles(orb, star2, p3, p4);

        move_to_COM_frame(star1, star2, p1, p2, p3, p4);
        /*----------------------------------------------------------------------------------------------------------------*/
        Solver solver{0, star1, star2, p1, p2, p3, p4};

        Solver::RunArgs args;

        Scalar t_end = 2 * time_to_periapsis(orbit::group(star1, p1, p2), orbit::group(star2, p3, p4));

        args.add_stop_condition(t_end);

        solver.run(args);
    }

    print(std::cout, "simulation with planetary scattering complete!\n");

    return 0;
}
