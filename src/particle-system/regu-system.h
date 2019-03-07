
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.h"
#include "core-computation.h"
#include "dev-tools.h"
#include <memory>

namespace SpaceH {
    enum class ReguType {
        logH, TTL, none
    };

    template<typename Scalar, ReguType Type = ReguType::logH>
    class Regularization {
    public:
        SPACEHUB_STD_ACCESSOR(omega, omega_);

        SPACEHUB_STD_ACCESSOR(bindE, bindE_);

        template<typename Particles>
        explicit Regularization(Particles const &partc) {
            omega_ = capital_omega(partc);
            bindE_ = -calc::calc_total_energy(partc);
        }

        template<typename Particles>
        inline auto eval_pos_phy_time(Particles const &partc, Scalar stepSize) {
            if constexpr (Type == ReguType::logH) {
                return stepSize / (bindE_ + calc::calc_kinetic_energy(partc));
            } else if (Type == ReguType::TTL) {
                return stepSize / omega_;
            } else if (Type == ReguType::none) {
                return stepSize;
            } else {
                SPACEHUB_ABORT("Undefined regularization type!");
            }
        }

        template<typename Particles>
        inline auto eval_vel_phy_time(Particles const &partc, Scalar stepSize) {
            if constexpr (Type == ReguType::logH) {
                return stepSize / -calc::calc_potential_energy(partc);
            } else if (Type == ReguType::TTL) {
                return stepSize / capital_omega();
            } else if (Type == ReguType::none) {
                return stepSize;
            } else {
                SPACEHUB_ABORT("Undefined regularization type!");
            }
        }

    private:
        template<typename Particles>
        inline auto capital_omega(Particles const &partc) {
            return -calc::calc_potential_energy(partc);
        }

        Scalar omega_;
        Scalar bindE_;
    };

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
     * http://iopscience.iop.org/article/10.1086/301102/meta and
     * https://link.springer.com/article/10.1023%2FA%3A1021149313347.
     * @tparam Particles
     * @tparam Interactions
     */
    template<typename Particles, typename Interactions, ReguType ReguType>
    class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interactions, ReguType>> {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        SPACEHUB_STD_ACCESSOR(impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(impl_time, ptc_.time());
        /* Typedef */

        RegularizedSystem() = delete;

        template<typename Container>
        RegularizedSystem(Container const &partc, Scalar t) : ptc_(partc, t), acc_(partc.size()), regu_(partc) {}

        size_t impl_number() {
            return ptc_.number();
        }
        void impl_advance_time(Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            ptc_.time() += phyTime;
        }

        void impl_advance_pos(Coord const &velocity, Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(velocity, phyTime);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            pure_advance_vel(acceleration, phyTime);
        }

        void impl_evaluate_acc(Coord const &acceleration) {
            eom_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar stepSize) {
            auto phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(vel(), phyTime);
            ptc_.time() += phyTime;
        }

        void impl_kick(Scalar stepSize) {
            auto phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            auto halfTime = 0.5 * phyTime;

            if constexpr (Interactions::is_vel_dep) {
                eom_.eval_vel_indep_acc(ptc_, acc_.vel_indep_acc());
                kick_pseu_vel(halfTime);
                kick_real_vel(phyTime);
                advance_omega(ptc_.aux_vel(), acc_.vel_indep_acc(), phyTime);
                advance_bindE(ptc_.aux_vel(), phyTime);
                kick_pseu_vel(halfTime);
            } else {
                eom_.eval_acc(ptc_, acc());
                pure_advance_vel(acc(), halfTime);
                advance_omega(vel(), acc(), phyTime);
                pure_advance_vel(acc(), halfTime);
            }
        }

        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem const &ps) {
            os << ps.ptc_;
        }

    private:
        void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar stepSize) {
            Scalar d_omega = calc::coord_contract_to_scalar(mass(), velocity, d_omega_dr);
            regu_.omega() += d_omega * stepSize;
        }

        void advance_bindE(Coord const &velocity, Scalar stepSize) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(mass(), velocity, acc_.vel_dep_acc());
            regu_.bindE() += d_bindE * stepSize;
        }

        void pure_advance_pos(Coord const &velocity, Scalar stepSize) {
            calc::array_advance(pos(), velocity, stepSize);
        }

        void pure_advance_vel(Coord const &acceleration, Scalar stepSize) {
            calc::array_advance(vel(), acceleration, stepSize);
        }

        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::array_add(acc(), acc_.vel_indep_acc(), acc_.vel_dep_acc());
            calc::array_advance(ptc_.aux_vel(), acc(), stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            eom_.eval_aux_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::array_add(acc(), acc_.vel_indep_acc(), acc_.vel_dep_acc());
            pure_advance_vel(acc(), stepSize);
        }

        Particles ptc_;
        Interactions eom_;
        Regularization<Scalar, ReguType> regu_;
    };
}

#endif
