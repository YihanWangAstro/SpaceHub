
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.h"
#include "core-computation.tpp"
#include "../dev-tools.h"

namespace SpaceH {
    enum class ReguType {
        logH, TTL, none
    };

    template<typename Scalar, ReguType Type = ReguType::logH>
    class Regularization {
    public:
        SPACEHUB_STD_ACCESSOR(auto, omega, omega_);

        SPACEHUB_STD_ACCESSOR(auto, bindE, bindE_);

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
    template<typename Particles, typename Interactions, ReguType RegType>
    class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interactions, RegType>> {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(auto, impl_time, ptc_.time());

        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());
        /* Typedef */

        RegularizedSystem() = delete;

        static constexpr ReguType regu_type{RegType};

        template<typename STL>
        RegularizedSystem(STL const &ptc, Scalar t) : ptc_(ptc, t), acc_(ptc.size()), newtonian_acc_(ptc.size()), regu_(ptc_) {

            if constexpr (Interactions::has_extra_vel_indep_acc) {
                extra_vel_indep_acc_.resize(ptc.size());
            }

            if constexpr (Interactions::has_extra_vel_dep_acc) {
                extra_vel_dep_acc_.resize(ptc.size());
                aux_vel_ = ptc_.vel();
            }
        }

        size_t impl_number() const {
            return ptc_.number();
        }
        void impl_advance_time(Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            ptc_.time() += phyTime;
        }

        void impl_advance_pos(Coord const &velocity, Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            calc::array_advance(ptc_.pos(), velocity, stepSize);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            pure_advance_vel(acceleration, phyTime);
        }

        void impl_evaluate_acc(Coord const &acceleration) const {
            eom_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            calc::coord_advance(ptc_.pos(), ptc_.vel(), stepSize);
            ptc_.time() += phyTime;
        }

        void impl_kick(Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            Scalar halfTime = 0.5 * phyTime;


            eval_vel_indep_acc();

            if constexpr (Interactions::has_extra_vel_dep_acc) {
                kick_pseu_vel(halfTime);
                kick_real_vel(phyTime);
                kick_pseu_vel(halfTime);
            } else {
                if constexpr (Interactions::has_extra_vel_indep_acc) {
                    calc::coord_add(acc_, newtonian_acc_, extra_vel_indep_acc_);

                    calc::coord_advance(ptc_.vel(), acc_, halfTime);
                    advance_omega(ptc_.vel(), newtonian_acc_, phyTime);
                    advance_bindE(ptc_.vel(), extra_vel_indep_acc_, phyTime);
                    calc::coord_advance(ptc_.vel(), acc_, halfTime);
                } else {
                    calc::coord_advance(ptc_.vel(), newtonian_acc_, halfTime);
                    advance_omega(ptc_.vel(), newtonian_acc_, phyTime);
                    calc::coord_advance(ptc_.vel(), newtonian_acc_, halfTime);
                }
            }
        }

        void impl_pre_iter_process() {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                aux_vel_ = ptc_.vel();
            }
        }

        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem const &ps) {
            os << ps.ptc_;
            return os;
        }
        friend std::istream &operator>>(std::istream &is, RegularizedSystem &ps) {
            is >> ps.ptc_;
            return is;
        }

    private:
        void eval_vel_indep_acc() {
            eom_.eval_newtonian_acc(ptc_, newtonian_acc_);
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                eom_.eval_extra_vel_indep_acc(ptc_, extra_vel_indep_acc_);
            }
        }

        void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar stepSize) {
            Scalar d_omega = calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_omega_dr);
            regu_.omega() += d_omega * stepSize;
        }

        void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar stepSize) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_bindE_dr);
            regu_.bindE() += d_bindE * stepSize;
        }

        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::coord_add(acc_, newtonian_acc_, extra_vel_dep_acc_);
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
            calc::coord_advance(aux_vel_, acc_, stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            std::swap(aux_vel_, ptc_.vel());
            eom_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            std::swap(aux_vel_, ptc_.vel());

            calc::coord_add(acc_, newtonian_acc_, extra_vel_dep_acc_);
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
            calc::coord_advance(ptc_.vel(), acc_, stepSize);

            advance_omega(aux_vel_, newtonian_acc_, stepSize);
            advance_bindE(aux_vel_, extra_vel_dep_acc_, stepSize);
        }

        Particles ptc_;
        Interactions eom_;
        Regularization<Scalar, RegType> regu_;

        Coord aux_vel_{0};
        Coord acc_{0};
        Coord newtonian_acc_{0};
        Coord extra_vel_indep_acc_{0};
        Coord extra_vel_dep_acc_{0};

    };
}

#endif
