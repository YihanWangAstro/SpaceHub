
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++17-extensions"
#ifndef POSTNEWTONIAN_H
#define POSTNEWTONIAN_H

#include "../protoType.h"
#include "../macros.h"
#include "../ownMath.h"

namespace SpaceH {

    namespace Const {
        constexpr double INV_C = 1 / C;
        constexpr double INV_C2 = INV_C * INV_C;
        constexpr double INV_C3 = INV_C2 * INV_C;
        constexpr double INV_C4 = INV_C3 * INV_C;
        constexpr double INV_C5 = INV_C4 * INV_C;
    }


    template<typename Particles>
    struct NewtonianForce {
        /* Typedef */
        using type = typename Particles::type;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */

        inline void operator()(const Particles& partc, Vector& acc) {

            for (auto& a : acc)
                a.setZero();

            size_t size = partc.particleNumber();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 1; j < size; ++j) {
                    Vector dr = partc.pos(j) - partc.pos(i);
                    Scalar re_r = dr.reNorm();
                    Scalar re_r3 = re_r * re_r * re_r;
                    acc[i] += dr * (re_r3 * partc.mass(j));
                    acc[j] -= dr * (re_r3 * partc.mass(i));
                }
            }
        }
    };


    template<typename Particles>
    struct KarmackNewtonian {
        /* Typedef */
        using type    = typename Particles::type;
        using Scalar  = typename type::Scalar;
        using Vector  = typename type::Vector;
        using Buildin = typename  SpaceH::get_value_type<Scalar>::type;
        /* Typedef */

        inline void operator()(const Particles& partc, Vector& acc) {

            for (auto& a : acc)
                a.setZero();

            size_t size = partc.particleNumber();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 1; j < size; ++j) {
                    Vector dr = partc.pos(j) - partc.pos(i);
                    Scalar re_r = KarmackFastInverseSquareRoot<Buildin>(dr.norm2());
                    Scalar re_r3 = re_r * re_r * re_r;
                    acc[i] += dr * (re_r3 * partc.mass(j));
                    acc[j] -= dr * (re_r3 * partc.mass(i));
                }
            }
        }
    };


    /**
     *
     * @tparam Particles
     * @tparam First
     * @tparam Second
     * @tparam Radiative
     */
    template<typename Particles, bool First = true, bool Second = false, bool Radiative = false>
    class PostNewtonianForce {
    public:
        /* Typedef */
        using type   = typename Particles::type;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */

        inline void operator()(const Particles& partc, Vector& acc) {

            for (auto &a : acc)
                a.setZero();
            
            if constexpr (Radiative) {
                RadiativeOrder(partc, acc);
            }
            if constexpr (Second) {
                SecondOrder(partc, acc);
            }
            if constexpr (First) {
                FirstOrder(partc, acc);
            }
        }
    private:
        inline void FirstOrder(const Particles& partc, Vector& acc) {
            size_t size = partc.particleNumber();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 1; j < size; ++j) {
                    Vector dr = partc.pos(j) - partc.pos(i);
                    Vector dv = partc.vel(j) - partc.vel(i);

                }
            }
        }

        inline void SecondOrder(const Particles& partc, Vector& acc) {
            size_t size = partc.particleNumber();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 1; j < size; ++j) {
                    Vector dr = partc.pos(j) - partc.pos(i);
                    Vector dv = partc.vel(j) - partc.vel(i);

                }
            }
        }

        inline void RadiativeOrder(const Particles& partc, Vector& acc) {
            size_t size = partc.particleNumber();
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 1; j < size; ++j) {
                    Vector dr = partc.pos(j) - partc.pos(i);
                    Vector dv = partc.vel(j) - partc.vel(i);

                }
            }
        }
    };
}

#endif

#pragma clang diagnostic pop