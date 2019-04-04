//
// Created by yihan on 4/4/19.
//

#ifndef SPACEHUB_RAND_GENERATOR_H
#define SPACEHUB_RAND_GENERATOR_H

#include <random>

namespace SpaceH::Random {
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
    };
}
#endif //SPACEHUB_RAND_GENERATOR_H
