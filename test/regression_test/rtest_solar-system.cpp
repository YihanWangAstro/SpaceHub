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
using namespace space::unit;
int main(int argc, char **argv) {
    auto system = outer_solar();

    double t_end = 11862_year;

    std::string system_name = "solar";
    double r_tol = 1e-14;

    fast_test_methods(system_name, system, t_end, r_tol);

    bench_mark_methods(system_name, system, t_end, r_tol);

    /*double r_tol_low = 5e-16;
    double r_tol_hi = 1e-6;
    err_scale_methods(system_name, system, t_end, r_tol_low, r_tol_hi);*/
    return 0;
}
