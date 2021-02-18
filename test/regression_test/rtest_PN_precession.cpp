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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
#include "../../src/spaceHub.hpp"
#include "rtest_utility.hpp"
using namespace space::unit;

template <typename Solver>
auto mercury(double e = 0) {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;
    using namespace space::consts;

    Particle sun{1_Ms}, mercury{0.055_Me};
    auto orbit = EllipOrbit(sun.mass, mercury.mass, a_mercury, e_mercury, i_mercury, LoAN_mercury, AoP_mercury, 0);

    move_particles(orbit, mercury);

    move_to_COM_frame(sun, mercury);

    return std::vector{sun, mercury};
}

template <typename Solver>
void PN_precession_test(std::string const &fname, double end_time, double rtol,
                        std::vector<typename Solver::Particle> const &p) {
    using namespace space;
    using namespace callback;
    using namespace tools;

    Solver sim{0, p};

    typename Solver::RunArgs args;

    args.rtol = rtol;

    std::cout << std::setprecision(16);

    args.add_operation(TimeSlice(DefaultWriter(fname + ".txt"), 0, end_time, 10000));

    args.add_stop_condition(end_time);

    Timer timer;

    timer.start();

    sim.run(args);

    std::cout << "time : " << timer.get_time() << " s\n";
}

template <typename simulation>
void run(std::string const &sim_type) {
    auto twobody_sys = mercury<simulation>();

    PN_precession_test<simulation>("PN-precession-" + sim_type, 100_year, 1e-15, twobody_sys);
}

int main(int argc, char **argv) {
    using f = force::Interactions<force::NewtonianGrav, force::PN1>;

    using method = methods::AR_Chain_Plus<f>;

    return 0;
}
