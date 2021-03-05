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
int main(int argc, char **argv) {
#ifdef MPFR_VERSION_MAJOR
    using method = hub::methods::ABITS<>;

    error_scale<method>("short-ecc", "ABITS", 1e-30, 1e-9, 100_year, two_body());
#else
    std::cout << "mpfr lib is not installed. Skip test for AR_ABITS\n";
#endif
    return 0;
}
