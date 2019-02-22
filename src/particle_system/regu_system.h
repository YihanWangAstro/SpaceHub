
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle_system.h"
#include "../core_computation.h"
#include "../dev_tools.h"

namespace SpaceH {

    enum class ReguType {
        logH, TTL, none
    };

    template<typename Scalar, ReguType Type = ReguType::logH>
    class Regularization {
    public:
        template<typename Particles>
        Regularization(Particles const &partc) {
            omega = capital_omega(partc);
            bindE = -comp::get_total_energy(partc);
        }

        template<typename Particles>
        inline Scalar get_physical_pos_time(Particles const &partc, Scalar stepSize) {
            if constexpr (Type == ReguType::logH) {
                return stepSize / (bindE + get_kinetic_energy(partc));
            } else if (Type == ReguType::TTL) {
                return stepSize / omega;
            } else if (Type == ReguType::none) {
                return stepSize;
            } else {
                SPACEHUB_ABORT("Undefined regularization type!");
            }
        }

        template<typename Particles>
        inline Scalar get_physical_vel_time(Particles const &partc, Scalar stepSize) {
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

        Scalar omega;
        Scalar bindE;
    private:
        template<typename Particles>
        inline Scalar capital_omega(Particles const &partc) {
            return -comp::get_potential_energy(partc);
        }
    };

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
     * http://iopscience.iop.org/article/10.1086/301102/meta and
     * https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
     * @tparam Particles
     * @tparam Interaction
     */
    template<typename Particles, typename Interaction, ReguType ReguType>
    class RegularizedSystem : public ParticleSystem<Particles, Interaction> {
    public:
        /* Typedef */
        using Base = ParticleSystem<Particles, Interaction>;
        SPACEHUB_USING_TYPE_SYSTEM_OF(Base);
        /* Typedef */

        using Base::partc_;
        using Base::eom_;

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void drift(Scalar stepSize) {
            Scalar physicalTime = regular_.get_physical_pos_time(partc_, stepSize);
            Base::drift(physicalTime);
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            Scalar physicalTime = regular_.get_physical_vel_time(partc_, stepSize);

            if constexpr (Interaction::isVelDependent) {
                eom_.eval_vel_indep_acc(partc_);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                comp::array_advance(eom_.auxi_vx, eom_.ax, 0.5 * physicalTime);
                comp::array_advance(eom_.auxi_vy, eom_.ay, 0.5 * physicalTime);
                comp::array_advance(eom_.auxi_vz, eom_.az, 0.5 * physicalTime);

                eom_.eval_auxi_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                advance_vel(physicalTime);
                advance_omega(partc_, eom_, physicalTime);
                advance_bindE(partc_, eom_, physicalTime);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                comp::array_advance(eom_.auxi_vx, eom_.ax, 0.5 * physicalTime);
                comp::array_advance(eom_.auxi_vy, eom_.ay, 0.5 * physicalTime);
                comp::array_advance(eom_.auxi_vz, eom_.az, 0.5 * physicalTime);
            } else {
                eom_.eval_acc(partc_);
                advance_vel(0.5 * physicalTime);
                advance_omega(partc_, eom_, physicalTime);
                advance_vel(0.5 * physicalTime);
            }
        }


    private:
        Scalar advance_omega(Particles const &partc, Interaction const &eom, Scalar stepSize) {

        }

        Scalar advance_bindE(Particles const &partc, Interaction const &eom, Scalar stepSize) {

        }

        Regularization<Scalar, ReguType> regular_;
    };
}

#endif
