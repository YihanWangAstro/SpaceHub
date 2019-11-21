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
#include "../../src/spaceHub.hpp"
#include <thread>

USING_NAMESPACE_SPACEHUB_ALL;

template<typename simulation>
void run(std::string const &sim_name) {
  auto earth_sys = earth_system<simulation>();

  basic_error_test<simulation>("earth-system-" + sim_name, 1000_year, 1e-13, earth_sys);

  auto[rtol, error] = error_scale<simulation>(1e-14, 1e-8, 1000_year, earth_sys);

  std::fstream err_stream{"earth-system-" + sim_name + ".scale", std::ios::out};

  err_stream << rtol << '\n' << error;
}

int main(int argc, char **argv) {
  using type = Types<double_k>;

  using force = interactions::NewtonianGrav;

  using particles = PointParticles<type>;

  using sim_sys = SimpleSystem<particles, force>;

  using regu_sys = RegularizedSystem<particles, force, ReguType::LogH>;

  using chain_sys = ChainSystem<particles, force>;

  using arch_sys = ARchainSystem<particles, force, ReguType::LogH>;

  //using iter = ConstOdeIterator<Symplectic2nd>;

  using iter = BurlishStoer<double, WorstOffender, PIDController>;

  run<Simulator<sim_sys, iter>>("sim");

  run<Simulator<regu_sys, iter>>("regu");

  run<Simulator<chain_sys, iter>>("chain");

  run<Simulator<arch_sys, iter>>("arch");

  return 0;
}

