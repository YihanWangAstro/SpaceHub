#include <array>
#include "../../src/spaceHub.hpp"

using namespace space;
using namespace space::ode_iterator;
using namespace space::integrator;
using namespace space::orbit;
using namespace space::particle_system;
using namespace space::particle_set;
using namespace unit;
using scalar = precise_d;
using type = Types<scalar, std::vector>;

template <typename Simulator>
void run_two_body_test(Simulator &nbody, std::string const &file_name) {
  typename Simulator::RunArgs args;

  auto end_time = 1000 * year;

  std::ofstream eng_file(file_name);

  // eng_file << std::setprecision(17);

  auto E0 = calc::calc_total_energy(nbody.particles());

  // args.step_size = 0.001 * unit::year;

  args.add_pre_step_operation(run_operations::TimeSlice(
      [&](auto &ptc) { eng_file << ptc.time() << ',' << calc::calc_energy_error(ptc, E0) << '\n'; }, 0, end_time));

  args.add_stop_condition(end_time);

  args.rtol = 1e-16;
  args.atol = 1e-12;

  nbody.run(args);
}

int main(int argc, char **argv) {
  using force = interactions::NewtonianGrav;

  using particles = PointParticles<type>;

  using particle = typename particles::Particle;

  particle sun{m_solar}, earth{m_earth};

  auto earth_orbit = EllipOrbit{
      sun.mass, earth.mass, semi_latus_rectum(au, 0.0167086), 0.0167086, 7.155 * deg, 174.9 * deg, 288.1 * deg, 0.0};

  move_particles(earth_orbit, earth);

  move_to_COM_frame(sun, earth);

  // using iter = BurlishStoer<double, WorstOffender, PIDController>;
  // using iter = ConstOdeIterator<Symplectic2nd>;
  // using iter = ConstOdeIterator<GaussDadau<particles::Coord>>;
  using iter = IAS15<typename particles::Coord, IAS15Error, PIDController>;
  {
    using sys = SimpleSystem<particles, force>;

    using simulation = Simulator<sys, iter>;

    simulation nbody{0.0, sun, earth};

    run_two_body_test(nbody, "circular-BS-simple.err");
  }

  /*{
    using sys = ChainSystem<particles, force>;

    using simulation = Simulator<sys, iter>;

    simulation nbody{0.0, sun, earth};

    run_two_body_test(nbody, "circular-BS-chain.err");
  }

   {
     using sys = RegularizedSystem<particles, force, ReguType::LogH>;

     using simulation = Simulator<sys, iter>;

     simulation nbody{0.0, sun, earth};

     run_two_body_test(nbody, "circular-BS-regu.err");
   }

   {
     using sys = ARchainSystem<particles, force, ReguType::LogH>;

     using Simulation = Simulator<sys, iter>;

     Simulation nbody{0.0, sun, earth};

     run_two_body_test(nbody, "circular-BS-archain.err");
   }
   particle ecc_sun{m_solar}, ecc_earth{m_earth};

   auto ecc_orbit = Kepler{ecc_sun.mass, ecc_earth.mass, semi_latus_rectum(au, 0.99), 0.99, 7.155 * deg,
                           174.9 * deg, 288.1 * deg, 0.0};

   move_particles_to(ecc_orbit, ecc_earth);

   move_to_com_coord(ecc_sun, ecc_earth);

   {
     using sys = SimpleSystem<particles, force>;

     using simulation = Simulator<sys, iter>;

     simulation nbody{0.0, ecc_sun, ecc_earth};

     run_two_body_test(nbody, "ecc-BS-simple.err");
   }

  {
     using sys = ChainSystem<particles, force>;

     using simulation = Simulator<sys, iter>;

     simulation nbody{0.0, ecc_sun, ecc_earth};

     run_two_body_test(nbody, "ecc-BS-chain.err");
   }

   {
     using sys = RegularizedSystem<particles, force, ReguType::LogH>;

     using simulation = Simulator<sys, iter>;

     simulation nbody{0.0, ecc_sun, ecc_earth};

     run_two_body_test(nbody, "ecc-BS-regu.err");
   }

   {
     using sys = ARchainSystem<particles, force, ReguType::LogH>;

     using simulation = Simulator<sys, iter>;

     simulation nbody{0.0, ecc_sun, ecc_earth};

     run_two_body_test(nbody, "ecc-BS-archain.err");
   }
 */
  return 0;
}
