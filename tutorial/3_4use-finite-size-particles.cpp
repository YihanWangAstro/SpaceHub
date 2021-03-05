
#include "../src/spaceHub.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
/*--------------------------------------------------New-----------------------------------------------------------*/
using Solver = methods::DefaultMethod<DefaultForce, particles::SizeParticles>;  // default particle are point particle.
                                                                                // check available particles at xx
/*----------------------------------------------------------------------------------------------------------------*/
using Particle = Solver::Particle;

int main(int argc, char** argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    Particle p1{1_Ms, 1_Rs};  // initialize also with radius
    Particle p2{1_Ms, 1_Rs};
    Particle p3{1_Ms, 1_Rs};
    /*----------------------------------------------------------------------------------------------------------------*/

    auto inner_orb = orbit::EllipOrbit(p1.mass, p2.mass, 5_AU, 0.001, 1_deg, 2_deg, 3_deg, 4_deg);

    auto outer_orb = orbit::EllipOrbit(orbit::M_tot(p1, p2), p3.mass, 40_AU, 0.1, 6_deg, 7_deg, 8_deg, 9_deg);

    orbit::move_particles(inner_orb, p2);

    orbit::move_to_COM_frame(p1, p2);

    orbit::move_particles(outer_orb, p3);

    orbit::move_to_COM_frame(p1, p2, p3);

    Solver solver{0, p1, p2, p3};

    Solver::RunArgs args;

    args.add_stop_condition(1000_year);

    args.add_operation(DefaultWriter("size-particles.txt"));

    /*--------------------------------------------------New-----------------------------------------------------------*/
    auto collision_detect = [](auto& ptc, auto h) {
        size_t particle_num = ptc.number();
        for (size_t i = 0; i < particle_num; ++i) {
            for (size_t j = i + 1; j < particle_num; ++j) {
                if (distance(ptc.pos(i), ptc.pos(j)) < ptc.radius(i) + ptc.radius(i)) {
                    return true;
                }
            }
        }
        return false;
    };

    args.add_stop_condition(collision_detect);
    /*----------------------------------------------------------------------------------------------------------------*/

    solver.run(args);

    print(std::cout, "simulation with finite size particles term complete!\n");

    return 0;
}
