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
auto BHB(double e = 0.9) {
    using Particle = typename Solver::Particle;
    using namespace hub;
    using namespace hub::unit;
    using namespace hub::orbit;
    using namespace hub::consts;

    Particle bh1{30_Ms}, bh2{50_Ms};
    auto orbit = Elliptic(bh1.mass, bh2.mass, 0.01_AU, e, 0, 0, 0, 0);

    move_particles(orbit, bh2);

    move_to_COM_frame(bh1, bh2);

    return std::vector{bh1, bh2};
}

template <typename Solver>
void PN_radiation_test(std::string const &fname, double end_time, double rtol,
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
    auto twobody_sys = BHB<simulation>(e);

    PN_radiation_test<simulation>("PN-radiation-" + sim_type, 70_year, 1e-14, twobody_sys);
}

int main(int argc, char **argv) {
    using namespace hub;
    using f = force::Interactions<force::NewtonianGrav, force::PN2p5>;

    // using f = force::Interactions<force::NewtonianGrav>;

    using method = methods::AR_Chain_Plus<f>;
    run<method>("e=0894", 0.894);
    run<method>("e=0896", 0.896);
    run<method>("e=0898", 0.898);
    run<method>("e=09", 0.9);
    return 0;
}
