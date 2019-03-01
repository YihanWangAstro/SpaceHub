
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.h"
#include "../accelerations.h"
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
        SPACEHUB_STD_SCALAR_INTERFACES(omega, omega_);
        SPACEHUB_STD_SCALAR_INTERFACES(bindE, bindE_);

        template<typename Particles>
        explicit Regularization(Particles const &partc) {
            omega_ = capital_omega(partc);
            bindE_ = -calc::calc_total_energy(partc);
        }

        template<typename Particles>
        inline auto eval_pos_phy_time(Particles const &partc, Scalar stepSize) {
            if constexpr (Type == ReguType::logH) {
                return stepSize / (bindE_ + get_kinetic_energy(partc));
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
                return stepSize / -get_potential_energy(partc);
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

        SPACEHUB_STD_ARRAY_INTERFACES(mass, ptc_.mass());
        SPACEHUB_STD_ARRAY_INTERFACES(idn, ptc_.idn());
        SPACEHUB_STD_SCALAR_INTERFACES(pos, ptc_.pos());
        SPACEHUB_STD_SCALAR_INTERFACES(vel, ptc_.vel());
        SPACEHUB_STD_SCALAR_INTERFACES(time, ptc_.time());
        SPACEHUB_STD_SCALAR_INTERFACES(acc, acc_.acc());

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

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(vel(), stepSize);
        }

        void impl_advance_pos(Coord const& velocity, Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(velocity, phyTime);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(acc(), stepSize);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            pure_advance_vel(acceleration, phyTime);
        }

        void impl_evaluate_acc(Coord const &acceleration) {
            eom_.eval_acc(ptc_, acceleration);
        }

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void impl_drift(Scalar stepSize) {
            auto phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(vel(), phyTime);
            ptc_.time() += phyTime;
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void impl_kick(Scalar stepSize) {
            auto phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            auto halfTime = 0.5 * phyTime;

            if constexpr (Interactions::isVelDependent) {
                eom_.eval_vel_indep_acc(ptc_, acc_.v_indep_acc());
                kick_pseu_vel(halfTime);
                kick_real_vel(phyTime);
                advance_omega(ptc_.aux_vel(), acc_.v_indep_acc(), phyTime);
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
        void
        advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar stepSize) {
            auto prodx = calc::array_dot(mass(), v_x, d_x);
            auto prody = calc::array_dot(mass(), v_y, d_y);
            auto prodz = calc::array_dot(mass(), v_z, d_z);
            regu_.omega() += (prodx + prody + prodz) * stepSize;
        }

        void advance_bindE(Coord const &velocity, Scalar stepSize) {
            auto prodx = calc::array_dot(mass(), v_x, acc_.vd_ax());
            auto prody = calc::array_dot(mass(), v_y, acc_.vd_ay());
            auto prodz = calc::array_dot(mass(), v_z, acc_.vd_az());
            regu_.bindE() -= (prodx + prody + prodz) * stepSize;
        }

        void pure_advance_pos(Coord const& velocity, Scalar stepSize) {
            calc::array_advance(pos(), velocity, stepSize);
        }

        void pure_advance_vel(Coord const& acceleration, Scalar stepSize) {
            calc::array_advance(vel(), acceleration, stepSize);
        }

        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_vel_dep_acc(ptc_, acc_.v_dep_acc());
            calc::array_add(acc(), acc_.v_indep_acc(), acc_.v_dep_acc());
            calc::array_advance(ptc_.aux_vel(), acc(), stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            eom_.eval_aux_vel_dep_acc(ptc_, acc_.v_dep_acc());
            calc::array_add(acc(), acc_.v_indep_acc(), acc_.v_dep_acc());
            pure_advance_vel(acc(), stepSize);
        }

        Particles ptc_;
        Interactions eom_;
        Accelerations<ScalarArray, Interactions::isVelDependent> acc_;
        Regularization<Scalar, ReguType> regu_;
    };
}

#endif
