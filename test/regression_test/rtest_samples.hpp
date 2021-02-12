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

auto two_body(double e = 0) {
    using Particle = typename space::DefaultMethod::Particle;
    using namespace space;
    using namespace space::unit;
    using namespace space::orbit;

    Particle sun{1_Ms}, earth{1_Me};
    auto orbit = EllipOrbit(sun.mass, earth.mass, 1_AU, e, 0, 0, 0, 0);

    move_particles(orbit, earth);

    move_to_COM_frame(sun, earth);

    return std::vector{sun, earth};
}

auto outer_solar() {
    using Particle = typename space::DefaultMethod::Particle;
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

auto earth_system() {
    using Particle = typename space::DefaultMethod::Particle;
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

auto kozai() {
    using Particle = typename space::DefaultMethod::Particle;
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

template <typename Solver, typename Pt>
auto basic_error_test(std::string const &fname, double end_time, double rtol, std::vector<Pt> const &p) ->
    typename Solver::Scalar {
    using namespace space;
    using namespace callback;
    using namespace tools;
    using Scalar = typename Solver::Scalar;
    Solver sim{0, p};

    typename Solver::RunArgs args;

    std::ofstream err_file(fname + ".err");

    err_file << std::setprecision(16);

    auto E0 = calc::calc_total_energy(sim.particles());

    Scalar tot_error = 0;

    size_t error_num = 0;

    args.rtol = rtol;

    args.atol = 0;

    std::cout << std::setprecision(16);

    // args.atol = args.rtol;
    args.add_operation(TimeSlice(
        [&](auto &ptc, auto step_size) {
            auto err = calc::calc_energy_error(ptc, E0);
            tot_error += err * err;
            error_num++;
            err_file << ptc.time() << ',' << err << '\n';
        },
        0, end_time));

    // args.add_operation(TimeSlice(DefaultWriter(fname + ".txt"), 0, end_time, 10000));

    args.add_stop_condition(end_time);

    sim.run(args);

    Scalar rms_err = sqrt(tot_error / error_num);

    std::cout << "The rms relative error of test: " + fname + " : " << rms_err << "\n";

    return rms_err;
}

template <typename Solver>
double bench_mark(double end_time, double rtol, std::vector<typename Solver::Particle> const &p) {
    using namespace space;
    using namespace callback;
    using namespace tools;

    double cpu = 1e20;  // in seconds
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

        auto t = timer.get_time();

        if (t < cpu) {
            cpu = t;
        }
    }
    return cpu;
}

template <typename Solver>
void error_scale(std::string const &system_name, const std::string &method_name, double rtol_start, double rtol_end,
                 double end_time, std::vector<typename Solver::Particle> const &p) {
    using namespace space;
    using namespace callback;
    using namespace tools;

    size_t n = static_cast<size_t>(log(rtol_end / rtol_start) / log(2)) + 1;

    std::vector<double> rtol(n);
    std::vector<double> err(n);

    multi_thread::indexed_multi_thread(n, [&](size_t thid) {
        typename Solver::RunArgs args;

        double tot_error = 0;

        size_t error_num = 0;

        auto E0 = orbit::E_tot(p);  // calc::calc_total_energy(sim.particles());

        args.add_pre_step_operation([&](auto &ptc, auto step_size) {
            auto r_err = calc::calc_energy_error(ptc, E0);
            tot_error += r_err * r_err;
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

    std::fstream err_stream{system_name + "-" + method_name + ".scale", std::ios::out};

    err_stream << rtol << '\n' << err;
}

template <typename System>
auto fast_err_methods(std::string const &system_name, System const &system, double t_end) {
    using namespace space;
    double rtol = 1e-14;
    std::vector<double> errs;
    errs.reserve(20);
    errs.push_back(basic_error_test<methods::BS<>>(system_name + "-BS", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_BS<>>(system_name + "-AR", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::Chain_BS<>>(system_name + "-Chain", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_Chain<>>(system_name + "-AR-chain", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_Chain_Plus<>>(system_name + "-AR-chain+", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::Radau_Plus<>>(system_name + "-Radau+", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::Chain_Radau_Plus<>>(system_name + "-Radau-chain+", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_Radau_Plus<>>(system_name + "-AR-Radau+", t_end, rtol, system));
    errs.push_back(
        basic_error_test<methods::AR_Radau_Chain_Plus<>>(system_name + "-AR-Radau-chain+", t_end, rtol, system));
    errs.push_back(
        basic_error_test<methods::AR_Sym6_Chain_Plus<>>(system_name + "-AR-sym6-chain+", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_Sym6_Plus<>>(system_name + "-AR-sym6", t_end, rtol, system));
    errs.push_back(
        basic_error_test<methods::AR_Sym8_Chain_Plus<>>(system_name + "-AR-sym8-chain+", t_end, rtol, system));
    errs.push_back(basic_error_test<methods::AR_Sym8_Plus<>>(system_name + "-AR-sym8", t_end, rtol, system));
    return errs;
}

template <typename System>
void bench_mark_methods(std::string const &system_name, System const &system, double t_end) {
    using namespace space;
    double rtol = 1e-14;
    std::ofstream file{system_name + "-benchmark.txt", std::ios::out};

    std::vector<std::string> names{
        "BS",           "AR",        "Chain",           "AR-chain",       "AR-chain+", "Radau+",
        "Radau-chain+", "AR-Radau+", "AR-Radau-chain+", "AR-sym6-chain+", "AR-sym6",   "AR-sym8-chain+",
        "AR-sym8"};
    std::vector<double> errs = fast_err_methods(system_name, system, t_end);
    std::vector<double> cpu_t;
    cpu_t.reserve(20);
    cpu_t.push_back(bench_mark<methods::BS<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_BS<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Chain_BS<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Chain<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Chain_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Radau_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Chain_Radau_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Radau_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Radau_Chain_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym6_Chain_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym6_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym8_Chain_Plus<>>(t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym8_Plus<>>(t_end, rtol, system));

    for (size_t i = 0; i < names.size(); ++i) {
        file << names[i] << ':' << cpu_t[i] << ':' << errs[i] << '\n';
    }
}

template <typename System>
auto err_scale_methods(std::string const &system_name, System const &system, double t_end) {
    using namespace space;
    error_scale<methods::BS<>>(system_name, "BS", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_BS<>>(system_name, "AR", 3e-16, 1e-11, t_end, system);
    error_scale<methods::Chain_BS<>>(system_name, "Chain", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Chain<>>(system_name, "AR-chain", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Chain_Plus<>>(system_name, "AR-chain+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::Radau_Plus<>>(system_name, "Radau+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::Chain_Radau_Plus<>>(system_name, "Radau-chain+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Radau_Plus<>>(system_name, "AR-Radau+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Radau_Chain_Plus<>>(system_name, "AR-Radau-chain+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Sym6_Chain_Plus<>>(system_name, "AR-sym6-chain+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Sym6_Plus<>>(system_name, "AR-sym6", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Sym8_Chain_Plus<>>(system_name, "AR-sym8-chain+", 3e-16, 1e-11, t_end, system);
    error_scale<methods::AR_Sym8_Plus<>>(system_name, "AR-sym8", 3e-16, 1e-11, t_end, system);
}