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
/**
 * @file rand-generator.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_RAND_GENERATOR_HPP
#define SPACEHUB_RAND_GENERATOR_HPP

#include <mutex>
#include <random>
#include "math.hpp"
#include "multi-thread/multi-thread.hpp"

/**
 * @namespace space::random
 * Documentation for space
 */
namespace space::random {
  static double Uniform(double low, double high) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    std::uniform_real_distribution<double> dist{low, high};
    return dist(generator);
  }

  static double Logarithm(double low, double high) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    double log_low = log10(low);
    double log_high = log10(high);
    std::uniform_real_distribution<double> dist{log_low, log_high};
    return pow(10, dist(generator));
  }

  static double PowerLaw(double power, double low, double high) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    if (!math::iseq(power, -1.0)) {
      double beta = power + 1;
      double f_low = pow(low, beta);
      double f_high = pow(high, beta);
      std::uniform_real_distribution<double> dist{f_low, f_high};
      return pow(dist(generator), 1.0 / beta);
    } else {
      return Logarithm(low, high);
    }
  }

  static double Normal(double mean = 0, double sigma = 1) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    std::normal_distribution<double> dist{mean, sigma};
    return dist(generator);
  }

  static double Maxwellian(double sigma_1d) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    std::normal_distribution<double> dist{0.0, sigma_1d};
    double x = dist(generator);
    double y = dist(generator);
    double z = dist(generator);
    return sqrt(x * x + y * y + z * z);
  }
}  // namespace space::random

#endif  //SPACEHUB_RAND_GENERATOR_HPP
