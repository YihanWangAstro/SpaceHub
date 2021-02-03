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

USING_NAMESPACE_SPACEHUB_ALL;

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
    using namespace run_operations;
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
    auto twobody_sys = BHB<simulation>();

    PN_radiation_test<simulation>("PN-radiation-" + sim_type, 70_year, 1e-15, twobody_sys);
}

int main(int argc, char **argv) {
    using namespace interactions;

    using type = Types<double>;

    using force = interactions::Interactions<interactions::NewtonianGrav>;

    using particles = PointParticles<type>;

    using sim_sys = SimpleSystem<particles, force>;

    using regu_sys = RegularizedSystem<particles, force, ReguType::LogH>;

    using chain_sys = ChainSystem<particles, force>;

    using arch_sys = ARchainSystem<particles, force, ReguType::LogH>;

    using base_integrator = LeapFrogDKD<type>;
    //    using iter = ConstOdeIterator<Symplectic2nd>;

    using err_estimator = WorstOffender<type>;

    using step_controller = PIDController<type>;

    using iter = BulirschStoer<base_integrator, err_estimator, step_controller>;

    using ias15_iter = IAS15<integrator::GaussRadau<type>, IAS15Error<type>, step_controller>;

    /*run<Simulator<sim_sys, iter>>("sim");

    run<Simulator<regu_sys, iter>>("regu");

    run<Simulator<chain_sys, iter>>("chain");*/

    run<Simulator<arch_sys, iter>>("arch");

    return 0;
}
