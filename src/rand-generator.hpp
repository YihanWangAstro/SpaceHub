//
// Created by yihan on 4/4/19.
//

#ifndef SPACEHUB_RAND_GENERATOR_H
#define SPACEHUB_RAND_GENERATOR_H

#include <random>
#include "own-math.hpp"

namespace space::randomGen {

    template<typename Dtype>
    class Uniform {
    private:
        Uniform() : gen(rd()), Dist(0, 1) {}

        Uniform(const Uniform &) = default;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<Dtype> Dist;
    public:
        inline static Dtype get(Dtype low = 0, Dtype high = 1) {
            static Uniform singleton;
            return low + (high - low) * singleton.Dist(singleton.gen);
        }

        template <typename RandGen>
        inline static Dtype get(RandGen& gen, Dtype low = 0, Dtype high = 1) {
            static Uniform singleton;
            return low + (high - low) * singleton.Dist(gen);
        }
    };

    template<typename Dtype>
    class Logarithm {
    private:
        Logarithm() : gen(rd()), Dist(0, 1) {}

        Logarithm(const Logarithm &) = default;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<Dtype> Dist;
    public:
        inline static Dtype get(Dtype low, Dtype high) {
            static Logarithm singleton;
            Dtype log_low = log10(low);
            Dtype log_high = log10(high);
            return pow(10, log_low + (log_high - log_low) * singleton.Dist(singleton.gen));
        }

        template <typename RandGen>
        inline static Dtype get(RandGen& gen, Dtype low, Dtype high) {
            static Logarithm singleton;
            Dtype log_low = log10(low);
            Dtype log_high = log10(high);
            return pow(10, log_low + (log_high - log_low) * singleton.Dist(gen));
        }
    };

    template<typename Dtype>
    class PowerLaw {
    private:
        PowerLaw() : gen(rd()), Dist(0, 1) {}

        PowerLaw(const PowerLaw &) = default;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<Dtype> Dist;
    public:

        inline static Dtype get(Dtype alpha, Dtype low, Dtype high) {
            static PowerLaw singleton;
            if(!space::iseq(alpha, -1.0)){
                auto beta = alpha + 1;
                auto f_low = pow(low, beta);
                auto f_high = pow(high, beta);
                return pow(f_low + (f_high - f_low) * singleton.Dist(singleton.gen), 1.0 / beta);
            } else {
                return Logarithm<Dtype>::get(low, high);
            }
        }

        template <typename RandGen>
        inline static Dtype get(RandGen& gen, Dtype alpha, Dtype low, Dtype high) {
            static PowerLaw singleton;
            if(!space::iseq(alpha, -1.0)){
                auto beta = alpha + 1;
                auto f_low = pow(low, beta);
                auto f_high = pow(high, beta);
                return pow(f_low + (f_high - f_low) * singleton.Dist(gen), 1.0 / beta);
            } else {
                return Logarithm<Dtype>::get(gen, low, high);
            }
        }
    };


    template<typename Dtype>
    class Normal {
    private:
        Normal() : gen(rd()), Dist(0, 1) {}

        Normal(const Normal &) = default;

        std::random_device rd;
        std::mt19937 gen;
        std::normal_distribution<Dtype> Dist;
    public:
        inline static Dtype get(Dtype mean = 0, Dtype sigma = 1) {
            static Normal singleton;
            return mean + sigma * singleton.Dist(singleton.gen);
        }

        template <typename RandGen>
        inline static Dtype get(RandGen& gen, Dtype mean = 0, Dtype sigma = 1) {
            static Normal singleton;
            return mean + sigma * singleton.Dist(gen);
        }
    };

    template<typename Dtype>
    class Maxwell {
    private:
        Maxwell() : gen(rd()), Dist(0, 1) {}

        Maxwell(const Maxwell &) = default;

        std::random_device rd;
        std::mt19937 gen;
        std::normal_distribution<Dtype> Dist;
    public:
        inline static Dtype get(Dtype sigma = 1) {
            static Maxwell singleton;
            auto x = singleton.Dist(singleton.gen);
            auto y = singleton.Dist(singleton.gen);
            auto z = singleton.Dist(singleton.gen);

            return sigma * sqrt(x * x + y * y + z * z);
        }

        template <typename RandGen>
        inline static Dtype get(RandGen& gen, Dtype sigma = 1) {
            static Maxwell singleton;
            auto x = singleton.Dist(gen);
            auto y = singleton.Dist(gen);
            auto z = singleton.Dist(gen);

            return sigma * sqrt(x * x + y * y + z * z);
        }
    };
}
#endif //SPACEHUB_RAND_GENERATOR_H
