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
/**
 * @file rand-generator.hpp
 *
 * Header file.
 */
#pragma once

#include <cmath>
#include <functional>
#include <random>

#include "math.hpp"
#include "multi-thread/multi-thread.hpp"

/**
 * @namespace hub::random
 * Documentation for hub
 */
namespace hub::random {

    /**
     * @brief Uniform distributed random number generator
     *
     * @param[in] low The lower limit of the distribution.
     * @param[in] high The higher limit of the distribution.
     * @return double The generated number.
     */
    inline static double Uniform(double low, double high) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        std::uniform_real_distribution<double> dist{low, high};
        return dist(generator);
    }

    /**
     * @brief Logarithmic distributed random number generator
     *
     * @param[in] low The lower limit of the distribution.
     * @param[in] high The higher limit of the distribution.
     * @return double The generated number.
     */
    inline static double Logarithm(double low, double high) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        double log_low = log10(low);
        double log_high = log10(high);
        std::uniform_real_distribution<double> dist{log_low, log_high};
        return pow(10, dist(generator));
    }

    /**
     * @brief Power law distributed random number generator
     *
     * @param[in] power The power of the power law
     * @param[in] low The lower limit of the distribution.
     * @param[in] high The higher limit of the distribution.
     * @return double The generated number.
     */
    inline static double PowerLaw(double power, double low, double high) {
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

    /**
     * @brief Normal distributed random number generator
     *
     * @param[in] mean The mean value of the normal distribution.
     * @param[in] sigma The standard deviation of the normal distribution.
     * @return double The generated number.
     */
    inline static double Normal(double mean = 0, double sigma = 1) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        std::normal_distribution<double> dist{mean, sigma};
        return dist(generator);
    }

    /**
     * @brief Truncated normal distribution random number generator
     *
     * @param[in] low The lower limit of the distribution.
     * @param[in] high The higher limit of the distribution.
     * @param[in] mean The mean value of the untruncated normal distribution.
     * @param[in] sigma The standard deviation of the untruncated normal distribution.
     * @return double The generated number.
     */
    inline static double TruncatedNormal(double low, double high, double mean = 0, double sigma = 1) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        std::normal_distribution<double> dist{mean, sigma};
        // TODO: can be optimized by using inverse transformed method.
        for (;;) {
            auto r = dist(generator);
            if (low <= r && r <= high) {
                return r;
            }
        }
    }

    /**
     * @brief Maxwellian distributed random number generator
     *
     * @param[in] sigma_1d The 1D dispersion of the Maxwellian distribution.
     * @return double The generated number.
     */
    inline static double Maxwellian(double sigma_1d) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        std::normal_distribution<double> dist{0.0, sigma_1d};
        double x = dist(generator);
        double y = dist(generator);
        double z = dist(generator);
        return sqrt(x * x + y * y + z * z);
    }

    /**
     * @brief Truncated Maxwellian distributed random number generator
     *
     * @param[in] low The lower limit of the distribution.
     * @param[in] high The higher limit of the distribution.
     * @param[in] sigma_1d The 1D dispersion of the Maxwellian distribution.
     * @return double The generated number.
     */
    inline static double TruncatedMaxwellian(double low, double high, double sigma_1d) {
        static thread_local std::mt19937 generator{std::random_device{}()};
        std::normal_distribution<double> dist{0.0, sigma_1d};
        // TODO: can be optimized by using inverse transformed method.
        for (;;) {
            double x = dist(generator);
            double y = dist(generator);
            double z = dist(generator);
            double r = sqrt(x * x + y * y + z * z);
            if (low <= r && r <= high) {
                return r;
            }
        }
    }

    inline double Fixed(double x) { return x; }

    class Parameter {
       public:
        template <typename Dist, typename... Args>
        Parameter(Dist distribution, Args... args) : dist_{std::bind(distribution, args...)} {}

           double draw() const {
            double r = dist_();
            for (; r < min_ || r > max_;) {
                r = dist_();
            }
            return r;
        }

        void set_max(double max) { max_ = max; }

        void set_min(double min) { min_ = min; }

        void set_boundary(double min, double max) { min_ = min, max_ = max; }

       private:
        double min_{std::numeric_limits<double>::min()};
        double max_{std::numeric_limits<double>::max()};
        std::function<double()> dist_;
    };

}  // namespace hub::random
