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
#include "rtest_samples.hpp"

USING_NAMESPACE_SPACEHUB_ALL;

template <typename simulation>
void run(std::string const &sim_name, std::fstream &benck_mark_file) {
    auto system = kozai<simulation>();

    double t_end = 10000_year;

    auto rms_err = basic_error_test<simulation>("kozai-" + sim_name, t_end, 1e-15, system);

    auto [rtol, error] = error_scale<simulation>(3e-16, 1e-11, t_end, system);

    std::fstream err_stream{"kozai-" + sim_name + ".scale", std::ios::out};

    err_stream << rtol << '\n' << error;

    benck_mark_file << sim_name << ":" << bench_mark<simulation>(t_end, 1e-15, system) << ":" << rms_err << '\n';
}

int main(int argc, char **argv) {
    using type = Types<double>;

    using adtype = Types<double_k>;

    using force = interactions::Interactions<interactions::NewtonianGrav>;

    using base_integrator = LeapFrogDKD<type>;

    using err_estimator = WorstOffender<type>;

    using step_controller = PIDController<type>;

    using particles = PointParticles<type>;

    using adparticles = PointParticles<adtype>;

    using sim_sys = SimpleSystem<particles, force>;

    using regu_sys = RegularizedSystem<particles, force, ReguType::LogH>;

    using chain_sys = ChainSystem<particles, force>;

    using arch_sys = ARchainSystem<particles, force, ReguType::LogH>;

    using adsim_sys = SimpleSystem<adparticles, force>;

    using adregu_sys = RegularizedSystem<adparticles, force, ReguType::LogH>;

    using adchain_sys = ChainSystem<adparticles, force>;

    using adarch_sys = ARchainSystem<adparticles, force, ReguType::LogH>;

    using iter = BurlishStoer<base_integrator, err_estimator, step_controller>;

    using ias15_iter = IAS15<integrator::GaussDadau<adtype>, IAS15Error<type>, step_controller>;

    using space_iter = BisecOdeIterator<integrator::Symplectic6th<type>, WorstOffender<type>, step_controller>;

    std::fstream benck_mark_file{"kozai-benchmark.txt", std::ios::out};

    run<Simulator<sim_sys, iter>>("BS", benck_mark_file);
    run<Simulator<regu_sys, iter>>("AR", benck_mark_file);
    run<Simulator<chain_sys, iter>>("Chain", benck_mark_file);
    run<Simulator<arch_sys, iter>>("AR-chain", benck_mark_file);
    run<Simulator<adarch_sys, iter>>("AR-chain+", benck_mark_file);
    run<Simulator<adsim_sys, ias15_iter>>("IAS15(SpaceHub)", benck_mark_file);
    run<Simulator<adchain_sys, ias15_iter>>("C-IAS15", benck_mark_file);
    run<Simulator<adregu_sys, ias15_iter>>("AR-IAS15", benck_mark_file);
    run<Simulator<adarch_sys, ias15_iter>>("ARC-IAS15", benck_mark_file);
    run<Simulator<adarch_sys, space_iter>>("ARC-sym6", benck_mark_file);

    return 0;
}
