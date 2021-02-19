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

using namespace space::unit;

template <typename Solver>
auto BHB(double e = 0) {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;
    using namespace space::consts;

    Particle bh1{30_Ms}, bh2{50_Ms};
    auto orbit = EllipOrbit(bh1.mass, bh2.mass, 0.01_AU, 0.9, 0, 0, 0, 0);

    move_particles(orbit, bh2);

    move_to_COM_frame(bh1, bh2);

    return std::vector{bh1, bh2};
}

template <typename Solver>
void PN_radiation_test(std::string const &fname, double end_time, double rtol,
                       std::vector<typename Solver::Particle> const &p) {
    using namespace space;
    using namespace callback;
    using namespace tools;

    Solver sim{0, p};

    typename Solver::RunArgs args;

    args.rtol = rtol;

    std::cout << std::setprecision(16);

    args.add_operation(TimeSlice(DefaultWriter(fname + ".txt"), 0.0, end_time, 10000));

    args.add_stop_condition(end_time);

    Timer timer;

    timer.start();

    sim.run(args);

    std::cout << "time : " << timer.get_time() << " s\n";
}

template <typename simulation>
void run(std::string const &sim_type) {
    auto twobody_sys = BHB<simulation>();

    PN_radiation_test<simulation>("PN-radiation-" + sim_type, 70_year, 1e-15, twobody_sys);
}

int main(int argc, char **argv) {
    // run<space::methods::AR_Chain_Plus>("AR-Chain+");
    return 0;
}
