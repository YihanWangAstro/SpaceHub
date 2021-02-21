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
    using method = space::methods::AR_Radau_Plus<>;

    double rtol = 1e-14;

    // basic_error_test<method>("earth", 100_year, rtol, earth_system(), true, false);

    // basic_error_test<method>("ecc", 1000_year, rtol, two_body(0.9999), true, false);

    // basic_error_test<method>("outer", 11862_year, rtol, outer_solar(), true, false);
    // std::cout << std::setprecision(16) << kozai();
    // return 0;
    // basic_error_test<method>("kozai", 150000_year, rtol, kozai(), true, false);

    basic_error_test<method>("kozai", 100000_year, rtol, kozai(), true, false);

    return 0;
}
