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

#include "rtest_samples.hpp"
#include "rtest_utility.hpp"
using namespace space::unit;
int main(int argc, char **argv) {
    using method = space::methods::AR_Chain_Plus<>;

#pragma STDC FENV_ACCESS ON
    std::fesetround(FE_TONEAREST);

    double rtol = 1e-14;

    basic_error_test<method>("outer", 11862_year * 1, rtol, outer_solar(), true, false);

    basic_error_test<method>("earth", 100_year * 1, rtol, earth_system(), true, false);

    basic_error_test<method>("ecc", 1000_year * 1, rtol, two_body(0.9999), true, false);

    return 0;
}
