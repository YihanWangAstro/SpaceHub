
#include "../src/spaceHub.hpp"
using namespace space;
using namespace unit;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;
using Vector3d = Particle::Vector;  // using Vector 3D type

int main(int argc, char **argv) {
    /*--------------------------------------------------New-----------------------------------------------------------*/
    int particle_number = 20;

    std::vector<Particle> particles(particle_number);

    for (auto &p : particles) {
        p.mass = 1e-10_Me;  // set mass to small value
        p.pos = Vector3d(random::Uniform(-1_AU, 1_AU), random::Uniform(-1_AU, 1_AU),
                         random::Uniform(-1_AU, 1_AU));  // generate x, y, z randomly between [-1,1]AU
        p.vel = Vector3d(random::Uniform(-10_kms, 10_kms), random::Uniform(-10_kms, 10_kms),
                         random::Uniform(-10_kms, 10_kms));  // generate vx, vy, vz randomly between [-10,10]km/s
    }

    Particle centre_obj{1_Ms};  // centre object at origin

    particles.push_back(centre_obj);  // add centre object to 'particles'

    orbit::move_to_COM_frame(particles);

    Solver solver{0, particles};  // create solver with particle array
    /*----------------------------------------------------------------------------------------------------------------*/

    Solver::RunArgs args;

    args.add_stop_condition(100_year);

    solver.run(args);

    print(std::cout, "simulation of multi-bodies without any output complete!\n");

    return 0;
}
