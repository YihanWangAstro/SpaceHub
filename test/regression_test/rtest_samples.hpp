
#ifndef SPACEHUB_RTEST_SAMPLES_HPP
#define SPACEHUB_RTEST_SAMPLES_HPP
#include "../../src/orbits/orbits.hpp"
#include "../../src/orbits/particle-manip.hpp"
#include "../../src/macros.hpp"

template<typename Solver>
auto two_body(double e) {
  using Particle = typename Solver::Particle;
  using namespace space;
  using namespace space::unit;

  Particle sun{1_Ms}, earth{1_Me};
  auto orbit = orbit::EllipOrbit(sun.mass, earth.mass, 1_AU, e, 0, 0, 0, 0);

  orbit::move_particles(orbit, earth);

  orbit::move_to_COM_frame(sun, earth);

  return std::make_tuple(sun, earth);
}


template<typename Solver>
auto kozai(double e) {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;

    Particle sun{1_Ms}, earth{1_Me};
    auto orbit = orbit::EllipOrbit(sun.mass, earth.mass, 1_AU, e, 0, 0, 0, 0);

    orbit::move_particles(orbit, earth);

    orbit::move_to_COM_frame(sun, earth);

    return std::make_tuple(sun, earth);
}

template<typename Solver>
auto outer_solar() {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;

    Particle sun{1_Ms}, earth{1_Me};
    auto orbit = orbit::EllipOrbit(sun.mass, earth.mass, 1_AU, e, 0, 0, 0, 0);

    orbit::move_particles(orbit, earth);

    orbit::move_to_COM_frame(sun, earth);

    return std::make_tuple(sun, earth);
}

#endif //SPACEHUB_RTEST_SAMPLES_HPP
