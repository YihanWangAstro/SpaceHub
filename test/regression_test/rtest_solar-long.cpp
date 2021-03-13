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

#include "rtest_samples.hpp"
#include "rtest_utility.hpp"
using namespace hub::unit;

template <typename Solver>
void job(std::string const &run_name, double rtol, double end_time) {
    using namespace hub;
    using namespace callback;
    using Scalar = typename Solver::Scalar;

    typename Solver::RunArgs args;
    args.rtol = rtol;
    args.atol = 0;

    Scalar E0 = 1;

    std::ofstream err_file;

    err_file.open(run_name + ".err");
    err_file << std::setprecision(16);

    args.add_stop_condition(end_time);

    args.add_start_point_operation([&](auto &ptc, auto step_size) { E0 = calc::calc_total_energy(ptc); });

    args.add_operation(LogTimeSlice(
        [&](auto &ptc, auto step_size) {
            auto err = calc::calc_energy_error(ptc, E0);
            err_file << ptc.time() << ',' << err << std::endl;
        },
        0.0, end_time, 100000));

    Solver sim{0, outer_solar()};

    hub::tools::Timer tick;
    tick.start();
    sim.run(args);
    print(std::cout, run_name, ':', tick.get_time(), '\n');

    return;
}
using namespace hub::unit;
int main(int argc, char **argv) {
    tf::Executor exe;

    double t_end = 11.862_year * 1e9;
    exe.silent_async(job<hub::methods::BS<>>, "BS", 1e-14, t_end);
    exe.silent_async(job<hub::methods::AR_Chain<>>, "AR-Chain", 1e-14, t_end);
    exe.silent_async(job<hub::methods::AR_Chain_Plus<>>, "AR-Chain+", 1e-14, t_end);
    exe.silent_async(job<hub::methods::AR_Sym6_Plus<>>, "AR-Sym6+", 1e-14, t_end);
    exe.silent_async(job<hub::methods::AR_Radau_Plus<>>, "AR-Radau+", 1e-14, t_end);

    exe.wait_for_all();
    return 0;
}
