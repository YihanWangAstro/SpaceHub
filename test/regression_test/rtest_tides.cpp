/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
#include "../../src/spaceHub.hpp"

using namespace hub::unit;

template <typename Solver>
auto binary(double e) {
    using Particle = typename Solver::Particle;
    using namespace hub;
    using namespace hub::unit;
    using namespace hub::orbit;
    using namespace hub::consts;

    Particle s1{1_Ms, 1_Rs, 0, 0}, s2{1_Ms, 1_Rs, 0.75, 0.25_year};

    auto orbit = Elliptic(s1.mass, s2.mass, 1_AU, e, 0, 0, 0, hub::consts::pi);

    move_particles(orbit, s2);

    move_to_COM_frame(s1, s2);

    return std::vector{s1, s2};
}

template <typename Solver>
void tidal_test(std::string const &fname, double end_time, double rtol,
                std::vector<typename Solver::Particle> const &p) {
    using namespace hub;
    using namespace callback;
    using namespace tools;

    Solver sim{0, p};

    typename Solver::RunArgs args;

    args.rtol = rtol;

    std::cout << std::setprecision(16);

    args.add_operation(TimeSlice(DefaultWriter(fname + ".txt"), 0.0, end_time));

    args.add_stop_condition(end_time);

    Timer timer;

    timer.start();

    sim.run(args);

    std::cout << "time : " << timer.get_time() << " s\n";
}

template <typename simulation>
void run(std::string const &sim_type, double e) {
    auto twobody_sys = binary<simulation>(e);

    tidal_test<simulation>("tidal-" + sim_type, 10000_year, 1e-14, twobody_sys);
}

int main(int argc, char **argv) {
    using namespace hub;
    using f = force::Interactions<force::NewtonianGrav, force::Tidal>;

    using method = methods::AR_Chain_Plus<f, particles::TideParticles>;
    run<method>("e=06", 0.6);
    run<method>("e=089", 0.89);
    run<method>("e=09", 0.9);
    run<method>("e=091", 0.91);
    return 0;
}
