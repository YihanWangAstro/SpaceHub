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
#pragma once

#include <algorithm>
#include <iomanip>
#include <tuple>

#include "../../src/spaceHub.hpp"
#include "../../src/taskflow/taskflow.hpp"

tf::Executor test_executor;

template <typename Solver, typename Pt>
auto basic_error_test(std::string const &fname, double end_time, double rtol, std::vector<Pt> const &p,
                      bool calc_err = true, bool IO = true, bool output = true) {
    using namespace hub;
    using namespace callback;
    using Scalar = typename Solver::Scalar;

    typename Solver::RunArgs args;
    args.rtol = rtol;
    args.atol = 0;

    set_mpreal_bits_from_rtol(rtol);

    Scalar E0 = 1;
    Scalar tot_error = 0;
    size_t error_num = 0;

    std::ofstream err_file;
    if (calc_err && IO) {
        err_file.open(fname + ".err");
        err_file << std::setprecision(16);
    }

    args.add_stop_condition(end_time);

    if (calc_err) {
        args.add_start_point_operation([&](auto &ptc, auto step_size) { E0 = calc::calc_total_energy(ptc); });
        args.add_operation(TimeSlice(
            [&](auto &ptc, auto step_size) {
                auto err = calc::calc_energy_error(ptc, E0);
                tot_error += err * err;
                error_num++;
                if (IO) {
                    err_file << ptc.time() << ',' << err << '\n';
                }
            },
            decltype(end_time)(0), end_time));
    }

    Solver sim{0, p};

    tools::Timer tick;
    tick.start();
    sim.run(args);
    double duration = tick.get_time();

    Scalar rms_err = error_num ? sqrt(tot_error / error_num) : 0;

    if (output) {
        std::cout << std::scientific << std::setprecision(4) << "  rtol: " << rtol
                  << " | wall time: " << std::defaultfloat << std::setw(8) << duration << " s";
        if (calc_err) {
            std::cout << " | rms/end err: " << std::scientific << rms_err << " / "
                      << calc::calc_energy_error(sim.particles(), E0);
        }
        std::cout << " | " << fname << '\n';
    }
#ifdef MPFR_VERSION_MAJOR
    if constexpr (std::is_same_v<Scalar, mpfr::mpreal>) {
        return std::make_tuple(rms_err.toDouble(), duration);
    } else
#endif
    {
        return std::make_tuple(rms_err, duration);
    }
}

template <typename System>
void fast_test_methods(std::string const &sys_name, System const &system, double t_end, double rtol = 1e-14) {
    using namespace hub;

    std::cout << "Running fast error test...\n";
    basic_error_test<methods::BS<>>(sys_name + "-BS", t_end, rtol, system);
    basic_error_test<methods::AR_BS<>>(sys_name + "-AR", t_end, rtol, system);
    basic_error_test<methods::Chain_BS<>>(sys_name + "-Chain", t_end, rtol, system);
    basic_error_test<methods::AR_Chain<>>(sys_name + "-AR-chain", t_end, rtol, system);
    basic_error_test<methods::AR_Chain_Plus<>>(sys_name + "-AR-chain+", t_end, rtol, system);
    basic_error_test<methods::Radau_Plus<>>(sys_name + "-Radau+", t_end, rtol, system);
    basic_error_test<methods::Chain_Radau_Plus<>>(sys_name + "-Radau-chain+", t_end, rtol, system);
    basic_error_test<methods::AR_Radau_Plus<>>(sys_name + "-AR-Radau+", t_end, rtol, system);
    basic_error_test<methods::AR_Radau_Chain_Plus<>>(sys_name + "-AR-Radau-chain+", t_end, rtol, system);
    basic_error_test<methods::AR_Sym6_Chain_Plus<>>(sys_name + "-AR-sym6-chain+", t_end, rtol, system);
    basic_error_test<methods::AR_Sym6_Plus<>>(sys_name + "-AR-sym6+", t_end, rtol, system);
    basic_error_test<methods::AR_Sym8_Chain_Plus<>>(sys_name + "-AR-sym8-chain+", t_end, rtol, system);
    basic_error_test<methods::AR_Sym8_Plus<>>(sys_name + "-AR-sym8+", t_end, rtol, system);
#ifdef MPFR_VERSION_MAJOR
    basic_error_test<methods::AR_ABITS<>>(sys_name + "-AR-ABITS", t_end, rtol, system);
#else
    std::cout << "mpfr lib is not installed. Skip test for AR_ABITS\n";
#endif
}

template <typename Solver, typename Pt>
auto bench_mark(std::string const &test_name, double end_time, double rtol, std::vector<Pt> const &p,
                size_t repeat = 5) {
    std::vector<double> errs(repeat), ts(repeat);

    auto [err, t] = basic_error_test<Solver>(test_name, end_time, rtol, p, true, false, false);

    for (size_t i = 0; i < repeat; ++i) {
        std::tie(errs[i], ts[i]) = basic_error_test<Solver, Pt>(test_name, end_time, rtol, p, false, false, false);
    }

    double cpu = *(std::min_element(ts.begin(), ts.end()));
    std::cout << "  rtol: " << std::scientific << std::setprecision(4) << rtol << " | wall time: " << std::defaultfloat
              << std::setw(8) << cpu << " s | rms err " << std::scientific << err << " | " << test_name << "\n";
    return std::make_tuple(err, cpu);
}

template <typename System>
void bench_mark_methods(std::string const &sys_name, System const &system, double t_end, double rtol = 1e-14) {
    using namespace hub;
    std::ofstream file{sys_name + "-benchmark.txt", std::ios::out};

    std::vector<std::string> names{
        "BS",           "AR",        "Chain",           "AR-chain",       "AR-chain+", "Radau+",
        "Radau-chain+", "AR-Radau+", "AR-Radau-chain+", "AR-sym6-chain+", "AR-sym6+",  "AR-sym8-chain+",
        "AR-sym8+",     "AR-ABITS"};

    // std::vector<std::string> names{"BS", "AR-chain", "AR-chain+", "AR-Radau+", "AR-sym6+", "AR-ABITS"};

    std::vector<std::tuple<double, double>> cpu_t;
    cpu_t.reserve(20);
    std::cout << "Running benchmark...\n";
    cpu_t.push_back(bench_mark<methods::BS<>>(sys_name + "-BS", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_BS<>>(sys_name + "-AR", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Chain_BS<>>(sys_name + "-Chain", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Chain<>>(sys_name + "-AR-chain", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Chain_Plus<>>(sys_name + "-AR-chain+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Radau_Plus<>>(sys_name + "-Radau+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::Chain_Radau_Plus<>>(sys_name + "-Radau-chain+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Radau_Plus<>>(sys_name + "-AR-Radau+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Radau_Chain_Plus<>>(sys_name + "-AR-Radau-chain+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym6_Chain_Plus<>>(sys_name + "-AR-sym6-chain+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym6_Plus<>>(sys_name + "-AR-sym6+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym8_Chain_Plus<>>(sys_name + "-AR-sym8-chain+", t_end, rtol, system));
    cpu_t.push_back(bench_mark<methods::AR_Sym8_Plus<>>(sys_name + "-AR-sym8+", t_end, rtol, system));
#ifdef MPFR_VERSION_MAJOR
    cpu_t.push_back(bench_mark<methods::AR_ABITS<>>(sys_name + "-AR-ABITS", t_end, rtol, system));
#else
    std::cout << "mpfr lib is not installed. Skip test for AR_ABITS\n";
#endif
    for (size_t i = 0; i < names.size(); ++i) {
        file << names[i] << ':' << std::get<1>(cpu_t[i]) << ':' << std::get<0>(cpu_t[i]) << '\n';
    }
}

template <typename Solver, typename Pt>
void error_scale(std::string const &system_name, const std::string &method_name, double rtol_start, double rtol_end,
                 double end_time, std::vector<Pt> const &p) {
    size_t n = 24;
    double base = POW(10, LOG10(rtol_end / rtol_start) / n);
    std::string test_name = system_name + "-" + method_name;

    std::vector<double> rtols(n);
    std::vector<double> errs(n);
    std::vector<double> wall_times(n);

    std::cout << "Running error scaling of " << test_name << " in range of rtol=[" << std::setw(7) << rtol_start << ","
              << rtol_end << "]\n";
    for (size_t thid = 0; thid < n; ++thid) {
        test_executor.silent_async(
            [&](size_t idx) {
                double rtol = rtol_start * POW(base, idx);
                auto [err, t] = bench_mark<Solver, Pt>(test_name, end_time, rtol, p, 1);
                rtols[idx] = rtol;
                errs[idx] = err;
                wall_times[idx] = t;
            },
            thid);
    }
    test_executor.wait_for_all();
    std::fstream err_stream{system_name + "-" + method_name + ".scale", std::ios::out};

    for (size_t i = 0; i < rtols.size(); ++i) {
        err_stream << rtols[i] << ',' << errs[i] << ',' << wall_times[i] << '\n';
    }
}

template <typename System>
auto err_scale_methods(std::string const &system_name, System const &system, double t_end, double rtol_start = 5e-16,
                       double rtol_end = 1e-6) {
    using namespace hub;

    error_scale<methods::BS<>>(system_name, "BS", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_BS<>>(system_name, "AR", rtol_start, rtol_end, t_end, system);
    error_scale<methods::Chain_BS<>>(system_name, "Chain", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Chain<>>(system_name, "AR-chain", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Chain_Plus<>>(system_name, "AR-chain+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::Radau_Plus<>>(system_name, "Radau+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::Chain_Radau_Plus<>>(system_name, "Radau-chain+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Radau_Plus<>>(system_name, "AR-Radau+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Radau_Chain_Plus<>>(system_name, "AR-Radau-chain+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Sym6_Chain_Plus<>>(system_name, "AR-sym6-chain+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Sym6_Plus<>>(system_name, "AR-sym6", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Sym8_Chain_Plus<>>(system_name, "AR-sym8-chain+", rtol_start, rtol_end, t_end, system);
    error_scale<methods::AR_Sym8_Plus<>>(system_name, "AR-sym8", rtol_start, rtol_end, t_end, system);
#ifdef MPFR_VERSION_MAJOR
    error_scale<methods::AR_ABITS<>>(system_name, "AR-ABITS", rtol_start, rtol_end, t_end, system);
#else
    std::cout << "mpfr lib is not installed. Skip test for AR_ABITS\n";
#endif
}