
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.h"
#include "../core-computation.tpp"
#include "../dev-tools.h"

namespace space {
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
        inline auto eval_pos_phy_time(Particles const &partc, Scalar step_size) {
            if constexpr (Type == ReguType::logH) {
                return step_size / (bindE_ + calc::calc_kinetic_energy(partc));
            } else if constexpr (Type == ReguType::TTL) {
                return step_size / omega_;
            } else if constexpr (Type == ReguType::none) {
                return step_size;
            } else {
                SPACEHUB_ABORT("Undefined regularization type!");
            }
        }

        template<typename Particles>
        inline auto eval_vel_phy_time(Particles const &partc, Scalar step_size) {
            if constexpr (Type == ReguType::logH) {
                return step_size / -calc::calc_potential_energy(partc);
            } else if constexpr (Type == ReguType::TTL) {
                return step_size / capital_omega(partc);
            } else if constexpr (Type == ReguType::none) {
                return step_size;
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

        RegularizedSystem(RegularizedSystem const &) = default;

        RegularizedSystem(RegularizedSystem &&) = default;

        RegularizedSystem &operator=(RegularizedSystem const &) = default;

        RegularizedSystem &operator=(RegularizedSystem &&) = default;

        static constexpr ReguType regu_type{RegType};

        template<typename STL>
        RegularizedSystem(Scalar t, STL const &ptc)
                : ptc_(t, ptc),
                  acc_(ptc.size()),
                  newtonian_acc_(ptc.size()),
                  regu_(ptc_) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
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

        void impl_advance_time(Scalar step_size) {
            Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
            ptc_.time() += phy_time;
        }

        void impl_advance_pos(Coord const &velocity, Scalar step_size) {
            Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
            calc::coord_advance(ptc_.pos(), velocity, phy_time);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar step_size) {
            Scalar phy_time = regu_.eval_vel_phy_time(ptc_, step_size);
            calc::coord_advance(ptc_.vel(), acceleration, phy_time);
        }

        void impl_evaluate_acc(Coord const &acceleration) const {
            eom_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar step_size) {
            Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
            calc::coord_advance(ptc_.pos(), ptc_.vel(), phy_time);
            ptc_.time() += phy_time;
        }

        void impl_kick(Scalar step_size) {
            Scalar phy_time = regu_.eval_vel_phy_time(ptc_, step_size);
            Scalar half_time = 0.5 * phy_time;

            eval_vel_indep_acc();

            if constexpr (Interactions::has_extra_vel_dep_acc) {
                kick_pseu_vel(half_time);
                kick_real_vel(phy_time);
                kick_pseu_vel(half_time);
            } else {
                if constexpr (Interactions::has_extra_vel_indep_acc) {
                    calc::coord_add(acc_, newtonian_acc_, extra_vel_indep_acc_);

                    calc::coord_advance(ptc_.vel(), acc_, half_time);
                    advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                    advance_bindE(ptc_.vel(), extra_vel_indep_acc_, phy_time);
                    calc::coord_advance(ptc_.vel(), acc_, half_time);
                } else {
                    calc::coord_advance(ptc_.vel(), newtonian_acc_, half_time);
                    advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                    calc::coord_advance(ptc_.vel(), newtonian_acc_, half_time);
                }
            }
        }

        void impl_pre_iter_process() {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                aux_vel_ = ptc_.vel();
            }
        }

        template <typename STL>
        void impl_to_linear_container(STL& stl){
            stl.clear();
            stl.reserve(impl_number()*6 +3);
            stl.emplace_back(impl_time());
            stl.emplace_back(omega());
            stl.emplace_back(bindE());
            add_coords_to(stl, impl_pos());
            add_coords_to(stl, impl_vel());
        }

        template <typename STL>
        void impl_load_from_linear_container(STL const& stl){
            auto i = stl.begin();
            impl_time() = *i, ++i;
            omega() = *i, ++i;
            bindE() = *i, ++i;
            load_to_coords(i, impl_pos());
            load_to_coords(i, impl_vel());
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

        void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar phy_time) {
            Scalar d_omega = calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_omega_dr);
            regu_.omega() += d_omega * phy_time;
        }

        void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar phy_time) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_bindE_dr);
            regu_.bindE() += d_bindE * phy_time;
        }

        void kick_pseu_vel(Scalar phy_time) {
            eom_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::coord_add(acc_, newtonian_acc_, extra_vel_dep_acc_);
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
            calc::coord_advance(aux_vel_, acc_, phy_time);
        }

        void kick_real_vel(Scalar phy_time) {
            std::swap(aux_vel_, ptc_.vel());
            eom_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            std::swap(aux_vel_, ptc_.vel());

            calc::coord_add(acc_, newtonian_acc_, extra_vel_dep_acc_);
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
            calc::coord_advance(ptc_.vel(), acc_, phy_time);

            advance_omega(aux_vel_, newtonian_acc_, phy_time);
            advance_bindE(aux_vel_, extra_vel_dep_acc_, phy_time);
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
