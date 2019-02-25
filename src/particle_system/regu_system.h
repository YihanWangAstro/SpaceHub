
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.h"
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
        explicit Regularization(Particles const &partc) {
            omega = capital_omega(partc);
            bindE = -calcu::get_total_energy(partc);
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
            return -calcu::get_potential_energy(partc);
        }
    };

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
     * http://iopscience.iop.org/article/10.1086/301102/meta and
     * https://link.springer.com/article/10.1023%2FA%3A1021149313347.
     * @tparam Particles
     * @tparam Interaction
     */
    template<typename Particles, typename Interaction, ReguType ReguType>
    class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interaction, ReguType>> {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        /* Typedef */

        RegularizedSystem() = delete;

        template<typename Container>
        RegularizedSystem(Container const &partc, Scalar t) : partc_(partc, t) {
            if constexpr (Interaction::isVelDependent) {
                eom_.auxi_vx = partc_.vx;
                eom_.auxi_vy = partc_.vy;
                eom_.auxi_vz = partc_.vz;
            }
        }

        size_t impl_number() {
            return partc_.number();
        }

        void impl_advance_time(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(partc_, stepSize);
            partc_.time += phyTime;
        }

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(partc_.vx, partc_.vy, partc_.vz, stepSize);
        }

        void impl_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(partc_, stepSize);
            pure_advance_pos(vx, vy, vz, phyTime);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(eom_.ax, eom_.ay, eom_.az, stepSize);
        }

        void impl_advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_vel_time(partc_, stepSize);
            pure_advance_vel(ax, ay, az, phyTime);
        }

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void impl_drift(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(partc_, stepSize);
            pure_advance_pos(partc_.vx, partc_.vy, partc_.vz, phyTime);
            partc_.time += phyTime;
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void impl_kick(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_vel_time(partc_, stepSize);

            if constexpr (Interaction::isVelDependent) {
                eom_.eval_vel_indep_acc(partc_);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                pure_advance_auxi_vel(eom_.ax, eom_.ay, eom_.az, 0.5 * phyTime);

                eom_.eval_auxi_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                pure_advance_vel(eom_.ax, eom_.ay, eom_.az, phyTime);
                advance_omega(eom_.vel_indep_ax, eom_.vel_indep_ay, eom_.vel_indep_az, phyTime);
                advance_bindE(eom_.vel_dep_ax, eom_.vel_dep_ay, eom_.vel_dep_az, phyTime);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                pure_advance_auxi_vel(eom_.ax, eom_.ay, eom_.az, 0.5 * phyTime);
            } else {
                eom_.eval_acc(partc_);
                pure_advance_vel(eom_.ax, eom_.ay, eom_.az, 0.5 * phyTime);
                advance_omega(eom_.ax, eom_.ay, eom_.az, phyTime);
                pure_advance_vel(eom_.ax, eom_.ay, eom_.az, 0.5 * phyTime);
            }
        }

    private:
        Scalar
        advance_omega(ScalarArray const &dodrx, ScalarArray const &dodry, ScalarArray const &dodrz, Scalar stepSize) {
            auto prodx = calcu::array_dot(partc_.mass, eom_.auxi_vx, dodrx);
            auto prody = calcu::array_dot(partc_.mass, eom_.auxi_vy, dodry);
            auto prodz = calcu::array_dot(partc_.mass, eom_.auxi_vz, dodrz);
            regular_.omega += (prodx + prody + prodz) * stepSize;
        }

        Scalar
        advance_bindE(ScalarArray const &dbdrx, ScalarArray const &dbdry, ScalarArray const &dbdrz, Scalar stepSize) {
            auto prodx = calcu::array_dot(partc_.mass, eom_.auxi_vx, dbdrx);
            auto prody = calcu::array_dot(partc_.mass, eom_.auxi_vy, dbdry);
            auto prodz = calcu::array_dot(partc_.mass, eom_.auxi_vz, dbdrz);
            regular_.bindE -= (prodx + prody + prodz) * stepSize;
        }

        void pure_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            calcu::array_advance(partc_.px, vx, stepSize);
            calcu::array_advance(partc_.py, vy, stepSize);
            calcu::array_advance(partc_.pz, vz, stepSize);
        }

        void pure_advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            calcu::array_advance(partc_.vx, ax, stepSize);
            calcu::array_advance(partc_.vy, ay, stepSize);
            calcu::array_advance(partc_.vz, az, stepSize);
        }

        void
        pure_advance_auxi_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            calcu::array_advance(eom_.auxi_vx, ax, stepSize);
            calcu::array_advance(eom_.auxi_vy, ay, stepSize);
            calcu::array_advance(eom_.auxi_vz, az, stepSize);
        }

        Particles partc_;
        Interaction eom_;
        Regularization<Scalar, ReguType> regular_;
    };
}

#endif
