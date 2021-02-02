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
#include <thread>

#include "rtest_samples.hpp"

USING_NAMESPACE_SPACEHUB_ALL;

int main(int argc, char **argv) {
    Timer t;

    t.start();
    basic_error_test<MethodList::AR_sym6>("outer", 11862_year, 1e-15, outer_solar());

    std::cout << t.get_time() << '\n';
    // basic_error_test<MethodList::AR_sym6>("earth", 100_year, 1e-15, earth_system());

    // basic_error_test<MethodList::AR_sym6>("ecc", 1000_year, 1e-15, two_body(0.9999));

    return 0;
}
