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
#ifndef SPACEHUB_UTEST_HPP
#define SPACEHUB_UTEST_HPP

#include "../../src/math.hpp"
#include "../../src/rand-generator.hpp"
using utest_scalar = double;
constexpr auto UTEST_EPSILON = 1 * space::math::epsilon_v<utest_scalar>;
constexpr auto UTEST_LOW = -1.0;  //-space::math::big_value_v<utest_scalar>;
constexpr auto UTEST_HIGH = 1.0;  // space::math::big_value_v<utest_scalar>;
constexpr size_t RAND_TEST_NUM = 10000;

#define APPROX(x) Approx(x).epsilon(UTEST_EPSILON).margin(UTEST_EPSILON)
#define UTEST_RAND space::random::Uniform(UTEST_LOW, UTEST_HIGH)
#endif  // SPACEHUB_UTEST_HPP
