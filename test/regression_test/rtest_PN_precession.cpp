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
#include "rtest_utility.hpp"
using namespace hub::unit;

template <typename Solver>
auto binary(double e = 0) {
    using Particle = typename Solver::Particle;
    using namespace hub;
    using namespace hub::unit;
    using namespace hub::orbit;
    using namespace hub::consts;

    Particle p1{1_Ms}, p2{1_Ms};
    auto orbit = Elliptic(p1.mass, p2.mass, 0.1_AU, e, 0, 0, 0, 0);

    move_particles(orbit, p2);

    move_to_COM_frame(p1, p2);

    return std::vector{p1, p2};
}

template <typename Solver>
void PN_precession_test(std::string const &fname, double end_time, double rtol,
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
void run(std::string const &sim_type, double ecc) {
    auto twobody_sys = binary<simulation>(ecc);

    PN_precession_test<simulation>("PN-precession-" + sim_type, 100_year, 1e-14, twobody_sys);
}

int main(int argc, char **argv) {
    using namespace hub;
    using f = force::Interactions<force::NewtonianGrav, force::PN1>;

    using method = methods::AR_Chain_Plus<f>;

    run<method>("e=06", 0.6);

    run<method>("e=09", 0.9);

    run<method>("e=095", 0.95);

    run<method>("e=099", 0.99);

    return 0;
}
