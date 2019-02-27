
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

        SPACEHUB_STD_ARRAY_INTERFACES(px, ptc_.px());
        SPACEHUB_STD_ARRAY_INTERFACES(py, ptc_.py());
        SPACEHUB_STD_ARRAY_INTERFACES(pz, ptc_.pz());
        SPACEHUB_STD_ARRAY_INTERFACES(vx, ptc_.vx());
        SPACEHUB_STD_ARRAY_INTERFACES(vy, ptc_.vy());
        SPACEHUB_STD_ARRAY_INTERFACES(vz, ptc_.vz());
        SPACEHUB_STD_ARRAY_INTERFACES(mass, ptc_.mass());
        SPACEHUB_STD_ARRAY_INTERFACES(idn, ptc_.idn());
        SPACEHUB_STD_ARRAY_INTERFACES(ax, acc_.ax());
        SPACEHUB_STD_ARRAY_INTERFACES(ay, acc_.ay());
        SPACEHUB_STD_ARRAY_INTERFACES(az, acc_.az());

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
            impl_advance_pos(vx(), vy(), vz(), stepSize);
        }

        void impl_advance_pos(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(v_x, v_y, v_z, phyTime);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(ax(), ay(), az(), stepSize);
        }

        void impl_advance_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            pure_advance_vel(a_x, a_y, a_z, phyTime);
        }

        void impl_evaluate_acc(ScalarArray &a_x, ScalarArray &a_y, ScalarArray &a_z) {
            eom_.eval_acc(ptc_, a_x, a_y, a_z);
        }

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void impl_drift(Scalar stepSize) {
            auto phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            pure_advance_pos(vx(), vy(), vz(), phyTime);
            ptc_.time() += phyTime;
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void impl_kick(Scalar stepSize) {
            auto phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            auto halfTime = 0.5 * phyTime;

            if constexpr (Interactions::isVelDependent) {
                eom_.eval_vel_indep_acc(ptc_, acc_.vid_ax(), acc_.vid_ay(), acc_.vid_az());
                kick_pseu_vel(halfTime);
                kick_real_vel(phyTime);
                advance_omega(ptc_.aux_vx(), ptc_.aux_vy(), ptc_.aux_vz(), acc_.vid_ax(), acc_.vid_ay(), acc_.vid_az(), phyTime);
                advance_bindE(ptc_.aux_vx(), ptc_.aux_vy(), ptc_.aux_vz(), phyTime);
                kick_pseu_vel(halfTime);
            } else {
                eom_.eval_acc(ptc_, ax(), ay(), az());
                pure_advance_vel(ax(), ay(), az(), halfTime);
                advance_omega(vx(), vy(), vz(), ax(), ay(), az(), phyTime);
                pure_advance_vel(ax(), ay(), az(), halfTime);
            }
        }

        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem const &ps) {
            os << ps.ptc_;
        }

    private:
        void
        advance_omega(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, ScalarArray const &d_x,
                      ScalarArray const &d_y, ScalarArray const &d_z, Scalar stepSize) {
            auto prodx = calc::array_dot(mass(), v_x, d_x);
            auto prody = calc::array_dot(mass(), v_y, d_y);
            auto prodz = calc::array_dot(mass(), v_z, d_z);
            regu_.omega() += (prodx + prody + prodz) * stepSize;
        }

        void advance_bindE(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, Scalar stepSize) {
            auto prodx = calc::array_dot(mass(), v_x, acc_.vd_ax());
            auto prody = calc::array_dot(mass(), v_y, acc_.vd_ay());
            auto prodz = calc::array_dot(mass(), v_z, acc_.vd_az());
            regu_.bindE() -= (prodx + prody + prodz) * stepSize;
        }

        void pure_advance_pos(ScalarArray const &v_x, ScalarArray const &v_y, ScalarArray const &v_z, Scalar stepSize) {
            calc::array_advance(px(), v_x, stepSize);
            calc::array_advance(py(), v_y, stepSize);
            calc::array_advance(pz(), v_z, stepSize);
        }

        void pure_advance_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            calc::array_advance(vx(), a_x, stepSize);
            calc::array_advance(vy(), a_y, stepSize);
            calc::array_advance(vz(), a_z, stepSize);
        }

        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_vel_dep_acc(ptc_, acc_.vd_ax(), acc_.vd_ay(), acc_.vd_az());
            sum_all_acc(acc_);
            calc::array_advance(ptc_.aux_vx(), ax(), stepSize);
            calc::array_advance(ptc_.aux_vy(), ay(), stepSize);
            calc::array_advance(ptc_.aux_vz(), az(), stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            eom_.eval_aux_vel_dep_acc(ptc_, acc_.vd_ax(), acc_.vd_ay(), acc_.vd_az());
            sum_all_acc(acc_);
            pure_advance_vel(ax(), ay(), az(), stepSize);
        }

        Particles ptc_;
        Interactions eom_;
        Accelerations<ScalarArray, Interactions::isVelDependent> acc_;
        Regularization<Scalar, ReguType> regu_;
    };
}

#endif
