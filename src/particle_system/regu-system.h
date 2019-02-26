
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
        SPACEHUB_STD_SCALAR_INTERFACES(omega, omega_);
        SPACEHUB_STD_SCALAR_INTERFACES(bindE, bindE_);

        template<typename Particles>
        explicit Regularization(Particles const &partc) {
            omega_ = capital_omega(partc);
            bindE_ = -calcu::get_total_energy(partc);
        }

        template<typename Particles>
        inline Scalar get_physical_pos_time(Particles const &partc, Scalar stepSize) {
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

    private:
        template<typename Particles>
        inline Scalar capital_omega(Particles const &partc) {
            return -calcu::get_potential_energy(partc);
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
     * @tparam Interaction
     */
    template<typename Particles, typename Interaction, ReguType ReguType>
    class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interaction, ReguType>> {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        SPACEHUB_STD_ARRAY_INTERFACES(px, partc_.px());

        SPACEHUB_STD_ARRAY_INTERFACES(py, partc_.py());

        SPACEHUB_STD_ARRAY_INTERFACES(pz, partc_.pz());

        SPACEHUB_STD_ARRAY_INTERFACES(vx, partc_.vx());

        SPACEHUB_STD_ARRAY_INTERFACES(vy, partc_.vy());

        SPACEHUB_STD_ARRAY_INTERFACES(vz, partc_.vz());

        SPACEHUB_STD_ARRAY_INTERFACES(ax, action_.ax());

        SPACEHUB_STD_ARRAY_INTERFACES(ay, action_.ay());

        SPACEHUB_STD_ARRAY_INTERFACES(az, action_.az());

        SPACEHUB_STD_ARRAY_INTERFACES(mass, partc_.mass());

        SPACEHUB_STD_ARRAY_INTERFACES(idn, partc_.idn());

        SPACEHUB_CONDITIONAL_ARRAY_INTERFACES(Interaction::isVelDependent, auxi_vx, (*auxi_vx_ptr));

        SPACEHUB_CONDITIONAL_ARRAY_INTERFACES(Interaction::isVelDependent, auxi_vy, (*auxi_vy_ptr));

        SPACEHUB_CONDITIONAL_ARRAY_INTERFACES(Interaction::isVelDependent, auxi_vz, (*auxi_vz_ptr));
        /* Typedef */

        RegularizedSystem() = delete;

        template<typename Container>
        RegularizedSystem(Container const &partc, Scalar t) : partc_(partc, t), action_(partc.size()), regular_(partc) {
            if constexpr (Interaction::isVelDependent) {
                auxi_vx_ptr = std::make_unique<ScalarArray>(vx());
                auxi_vy_ptr = std::make_unique<ScalarArray>(vy());
                auxi_vz_ptr = std::make_unique<ScalarArray>(vz());
            } else {
                auxi_vx_ptr = auxi_vy_ptr = auxi_vz_ptr = nullptr;
            }
        }

        size_t impl_number() {
            return partc_.number();
        }

        void impl_advance_time(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(*this, stepSize);
            partc_.time() += phyTime;
        }

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(vx(), vy(), vz(), stepSize);
        }

        void impl_advance_pos(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(*this, stepSize);
            pure_advance_pos(v_x, v_y, v_z, phyTime);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(ax(), ay(), az(), stepSize);
        }

        void impl_advance_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_vel_time(*this, stepSize);
            pure_advance_vel(a_x, a_y, a_z, phyTime);
        }

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void impl_drift(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_pos_time(*this, stepSize);
            pure_advance_pos(vx(), vy(), vz(), phyTime);
            partc_.time() += phyTime;
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void impl_kick(Scalar stepSize) {
            Scalar phyTime = regular_.get_physical_vel_time(*this, stepSize);

            if constexpr (Interaction::isVelDependent) {
                action_.eval_vel_indep_acc(*this);

                action_.eval_vel_dep_acc(*this);
                action_.sum_tot_acc();
                pure_advance_auxi_vel(ax(), ay(), az(), 0.5 * phyTime);

                action_.eval_auxi_vel_dep_acc(*this);
                action_.sum_tot_acc();
                pure_advance_vel(ax(), ay(), az(), phyTime);
                advance_omega(action_.vel_indep_ax(), action_.vel_indep_ay(), action_.vel_indep_az(), phyTime);
                advance_bindE(action_.vel_dep_ax(), action_.vel_dep_ay(), action_.vel_dep_az(), phyTime);

                action_.eval_vel_dep_acc(*this);
                action_.sum_tot_acc();
                pure_advance_auxi_vel(ax(), ay(), az(), 0.5 * phyTime);
            } else {
                action_.eval_acc(*this);
                pure_advance_vel(ax(), ay(), az(), 0.5 * phyTime);
                advance_omega(ax(), ay(), az(), phyTime);
                pure_advance_vel(ax(), ay(), az(), 0.5 * phyTime);
            }
        }

        friend std::ostream&operator<<(std::ostream& os, RegularizedSystem const & ps) {
            os << ps.partc_;
        }
    private:
        void
        advance_omega(ScalarArray const &dodrx, ScalarArray const &dodry, ScalarArray const &dodrz, Scalar stepSize) {
            auto prodx = calcu::array_dot(mass(), auxi_vx(), dodrx);
            auto prody = calcu::array_dot(mass(), auxi_vy(), dodry);
            auto prodz = calcu::array_dot(mass(), auxi_vz(), dodrz);
            regular_.omega() += (prodx + prody + prodz) * stepSize;
        }

        void
        advance_bindE(ScalarArray const &dbdrx, ScalarArray const &dbdry, ScalarArray const &dbdrz, Scalar stepSize) {
            auto prodx = calcu::array_dot(mass(), auxi_vx(), dbdrx);
            auto prody = calcu::array_dot(mass(), auxi_vy(), dbdry);
            auto prodz = calcu::array_dot(mass(), auxi_vz(), dbdrz);
            regular_.bindE() -= (prodx + prody + prodz) * stepSize;
        }

        void pure_advance_pos(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, Scalar stepSize) {
            calcu::array_advance(px(), v_x, stepSize);
            calcu::array_advance(py(), v_y, stepSize);
            calcu::array_advance(pz(), v_z, stepSize);
        }

        void pure_advance_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            calcu::array_advance(vx(), a_x, stepSize);
            calcu::array_advance(vy(), a_y, stepSize);
            calcu::array_advance(vz(), a_z, stepSize);
        }

        void
        pure_advance_auxi_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            calcu::array_advance(auxi_vx(), a_x, stepSize);
            calcu::array_advance(auxi_vy(), a_y, stepSize);
            calcu::array_advance(auxi_vz(), a_z, stepSize);
        }

        Particles partc_;
        Interaction action_;
        Regularization<Scalar, ReguType> regular_;
        std::unique_ptr<ScalarArray> auxi_vx_ptr;
        std::unique_ptr<ScalarArray> auxi_vy_ptr;
        std::unique_ptr<ScalarArray> auxi_vz_ptr;
    };
}

#endif
