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
#pragma once

#include <iomanip>
#include <tuple>

#include "../../src/spaceHub.hpp"

template <typename Solver>
auto two_body(double e = 0) {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;

    Particle sun{1_Ms}, earth{1_Me};
    auto orbit = EllipOrbit(sun.mass, earth.mass, 1_AU, e, 0, 0, 0, 0);

    move_particles(orbit, earth);

    move_to_COM_frame(sun, earth);

    return std::vector{sun, earth};
}

template <typename Solver>
auto outer_solar() {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;

    Particle sun{1.00000597682_Ms,     -4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6,
                 +6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9};
    Particle jup{1. / 1047.355,       +3.40546614227466e+0, +3.62978190075864e+0,   +3.42386261766577e-2,
                 -0.3254242234844626, 0.32078376079804843,  -0.00015504584274327233};
    Particle sat{1. / 3501.6,          +6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1,
                 -0.24261807906125546, 0.23236917361652312,  0.0009720111543216125};
    Particle ura{1. / 22869.,         +1.11636331405597e+1, +1.60373479057256e+1,  +3.61783279369958e-1,
                 -0.1894447922304644, 0.12000768830940597,  -0.0012655376696374542};
    Particle nep{1. / 19314.,          -3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1,
                 -0.01264216568441132, -0.18100181375010865, 0.0020831452402001265};
    Particle plu{7.4074074e-09,        -2.13858977531573e+1, +3.20719104739886e+1, +2.49245689556096e+0,
                 -0.10285755114348066, -0.12017192726456442, 0.038256490292647924};

    move_to_COM_frame(sun, jup, sat, ura, nep);

    return std::vector{sun, jup, sat, ura, nep};
}

template <typename Solver>
auto earth_system() {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;

    Particle sun{1_Ms}, earth{1_Me}, moon{1_Mmoon};

    auto moon_orbit = EllipOrbit(earth.mass, moon.mass, 384748_km, 0.0549006, 5.15_deg, 0, 0, 0);

    move_particles(moon_orbit, moon);

    auto orbit = EllipOrbit(sun.mass, earth.mass + moon.mass, 1_AU, 0, 0, 0, 0, 0);

    move_particles(orbit, earth, moon);

    move_to_COM_frame(sun, earth, moon);

    return std::vector{sun, earth, moon};
}

template <typename Solver>
auto kozai() {
    using Particle = typename Solver::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;

    Particle m1{1.4_Ms}, m2{0.1_Ms}, m3{0.8_Ms};

    auto in_orbit = EllipOrbit(m1.mass, m2.mass, 5_AU, 0.5, 0, 0, 120_deg, 0);

    move_particles(in_orbit, m2);

    move_to_COM_frame(m1, m2);

    auto out_orbit = EllipOrbit(m1.mass + m2.mass, m3.mass, 50_AU, 0.5, 88_deg, 0, 0, 0);

    move_particles(out_orbit, m3);

    move_to_COM_frame(m1, m2, m3);

    return std::vector{m1, m2, m3};
}

template <typename Solver>
double basic_error_test(std::string const& fname, double end_time, double rtol,
                        std::vector<typename Solver::Particle> const& p) {
    using namespace space;
    using namespace run_operations;
    using namespace tools;

    Solver sim{0, p};

    typename Solver::RunArgs args;

    std::ofstream err_file(fname + ".err");

    err_file << std::setprecision(16);

    auto E0 = calc::calc_total_energy(sim.particles());

    double tot_error = 0;

    size_t error_num = 0;

    args.rtol = rtol;

    args.atol = 0;

    std::cout << std::setprecision(16);

    // args.atol = args.rtol;
    args.add_operation(TimeSlice(
        [&](auto& ptc, auto step_size) {
            auto err = calc::calc_energy_error(ptc, E0);
            tot_error += err * err;
            error_num++;
            err_file << ptc.time() << ',' << err << '\n';
        },
        0, end_time));

    // args.add_operation(TimeSlice(DefaultWriter(fname + ".txt"), 0, end_time, 10000));

    args.add_stop_condition(end_time);

    sim.run(args);

    double rms_err = sqrt(tot_error / error_num);

    std::cout << "The rms relative error of test: " + fname + " : " << rms_err << "\n";

    return rms_err;
}

template <typename Solver>
double bench_mark(double end_time, double rtol, std::vector<typename Solver::Particle> const& p) {
    using namespace space;
    using namespace run_operations;
    using namespace tools;

    double cpu = 0;
    size_t repeat = 5;
    for (size_t i = 0; i < repeat; ++i) {
        Solver sim{0, p};

        typename Solver::RunArgs args;

        args.rtol = rtol;

        args.atol = 0;

        args.add_stop_condition(end_time);

        Timer timer;

        timer.start();

        sim.run(args);

        cpu += timer.get_time();
    }
    return cpu / repeat;
}

template <typename Solver>
auto error_scale(double rtol_start, double rtol_end, double end_time, std::vector<typename Solver::Particle> const& p) {
    using namespace space;
    using namespace run_operations;
    using namespace tools;

    size_t n = static_cast<size_t>(log(rtol_end / rtol_start) / log(2)) + 1;

    std::vector<double> rtol(n);
    std::vector<double> err(n);

    multi_thread::indexed_multi_thread(n, [&](size_t thid) {
        typename Solver::RunArgs args;

        double tot_error = 0;

        size_t error_num = 0;

        auto E0 = orbit::E_tot(p);  // calc::calc_total_energy(sim.particles());

        args.add_pre_step_operation([&](auto& ptc, auto step_size) {
            auto err = calc::calc_energy_error(ptc, E0);
            tot_error += err * err;
            // std::cout << ' ' << err << ' ' << error_num << ' ' << thid << '\n';
            error_num++;
        });

        args.add_stop_condition(end_time);

        args.rtol = rtol_start * pow(2, thid);

        args.atol = 0;

        Solver sim{0, p};
        sim.run(args);

        rtol[thid] = args.rtol;
        err[thid] = sqrt(tot_error / error_num);
        // std::cout << tot_error << ' ' << error_num << " " << thid << std::endl;
    });

    return std::make_tuple(rtol, err);
}
