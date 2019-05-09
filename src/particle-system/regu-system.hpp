
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle-system.hpp"
#include "../core-computation.hpp"
#include "../dev-tools.hpp"
#include <type_traits>

namespace space {
    enum class ReguType {
        logH, TTL, none
    };

    /*---------------------------------------------------------------------------*\
        Class Regularization Declaration
    \*---------------------------------------------------------------------------*/
    template<typename Scalar, ReguType Type = ReguType::logH>
    class Regularization {
    public:
        //Constructors
        template<typename Particles>
        explicit Regularization(Particles const &partc);

        //Public methods
        SPACEHUB_STD_ACCESSOR(auto, omega, omega_);

        SPACEHUB_STD_ACCESSOR(auto, bindE, bindE_);

        template<typename Particles>
        auto eval_pos_phy_time(Particles const &partc, Scalar step_size);

        template<typename Particles>
        auto eval_vel_phy_time(Particles const &partc, Scalar step_size);

    private:
        //Private methods
        template<typename Particles>
        inline auto capital_omega(Particles const &partc);

        //Private members
        Scalar omega_;
        Scalar bindE_;
    };


    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Declaration
    \*---------------------------------------------------------------------------*/

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
     * http://iopscience.iop.org/article/10.1086/301102/meta and
     * https://link.springer.com/article/10.1023%2FA%3A1021149313347.
     * @tparam Particles
     * @tparam Forces
     */
    template<typename Particles, typename Forces, ReguType RegType>
    class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Forces, RegType>> {
    public:
        //Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Base = ParticleSystem<RegularizedSystem<Particles, Forces, RegType>>;

        using Particle = typename Particles::Particle;

        //Constructors
        RegularizedSystem() = delete;

        RegularizedSystem(RegularizedSystem const &) = default;

        RegularizedSystem(RegularizedSystem &&) noexcept = default;

        RegularizedSystem &operator=(RegularizedSystem const &) = default;

        RegularizedSystem &operator=(RegularizedSystem &&) noexcept = default;

        template<typename STL>
        RegularizedSystem(Scalar t, STL const &ptc);

        //Static members
        static constexpr ReguType regu_type{RegType};


        //Public methods
        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

        //Friend functions
        template <typename P, typename F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem<P, F, R> const &ps);

        template <typename P, typename F, ReguType R>
        friend std::istream &operator>>(std::istream &is, RegularizedSystem<P, F, R> &ps);

    CRTP_impl:
        //CRTP implementation
        SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(auto, impl_time, ptc_.time());

        size_t impl_number() const;

        void impl_advance_time(Scalar step_size);

        void impl_advance_pos(Coord const &velocity, Scalar step_size);

        void impl_advance_vel(Coord const &acceleration, Scalar step_size);

        void impl_evaluate_acc(Coord const &acceleration) const;

        void impl_drift(Scalar step_size);

        void impl_kick(Scalar step_size);

        void impl_pre_iter_process();

        template<typename STL>
        void impl_to_linear_container(STL &stl);

        template<typename STL>
        void impl_load_from_linear_container(STL const &stl);

    private:
        //Private methods
        void eval_vel_indep_acc();

        void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar phy_time);

        void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        //Private members
        Particles ptc_;
        Forces forces_;
        Regularization<Scalar, RegType> regu_;

        Coord acc_;
        Coord newtonian_acc_;

        std::conditional_t<Forces::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;
        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;
        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> aux_vel_;
    };


    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    RegularizedSystem<Particles, Forces, RegType>::RegularizedSystem(Scalar t, const STL &ptc)
            : ptc_(t, ptc),
              acc_(ptc.size()),
              newtonian_acc_(ptc.size()),
              regu_(ptc_) {
        static_assert(is_container_v<STL>, "Only STL-like container can be used");
        if constexpr (Forces::ext_vel_indep) {
            ext_vel_indep_acc_.resize(ptc.size());
        }

        if constexpr (Forces::ext_vel_dep) {
            ext_vel_dep_acc_.resize(ptc.size());
            aux_vel_ = ptc_.vel();
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    std::istream &operator>>(std::istream &is, RegularizedSystem<Particles, Forces, RegType> &ps) {
        is >> ps.ptc_;
        return is;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    std::ostream &operator<<(std::ostream &os, const RegularizedSystem<Particles, Forces, RegType> &ps) {
        os << ps.ptc_;
        return os;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    size_t RegularizedSystem<Particles, Forces, RegType>::impl_number() const {
        return ptc_.number();
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_advance_time(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        ptc_.time() += phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        calc::coord_advance(ptc_.pos(), velocity, phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
        Scalar phy_time = regu_.eval_vel_phy_time(ptc_, step_size);
        calc::coord_advance(ptc_.vel(), acceleration, phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_evaluate_acc(const Coord &acceleration) const {
        forces_.eval_acc(ptc_, acceleration);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        calc::coord_advance(ptc_.pos(), ptc_.vel(), phy_time);
        ptc_.time() += phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_kick(Scalar step_size) {
        Scalar phy_time = regu_.eval_vel_phy_time(ptc_, step_size);
        Scalar half_time = 0.5 * phy_time;

        eval_vel_indep_acc();

        if constexpr (Forces::ext_vel_dep) {
            kick_pseu_vel(half_time);
            kick_real_vel(phy_time);
            kick_pseu_vel(half_time);
        } else {
            if constexpr (Forces::ext_vel_indep) {
                calc::coord_add(acc_, newtonian_acc_, ext_vel_indep_acc_);

                calc::coord_advance(ptc_.vel(), acc_, half_time);
                advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                advance_bindE(ptc_.vel(), ext_vel_indep_acc_, phy_time);
                calc::coord_advance(ptc_.vel(), acc_, half_time);
            } else {
                calc::coord_advance(ptc_.vel(), newtonian_acc_, half_time);
                advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                calc::coord_advance(ptc_.vel(), newtonian_acc_, half_time);
            }
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::impl_pre_iter_process() {
        if constexpr (Forces::ext_vel_dep) {
            aux_vel_ = ptc_.vel();
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    void RegularizedSystem<Particles, Forces, RegType>::impl_to_linear_container(STL &stl) {
        stl.clear();
        stl.reserve(impl_number() * 6 + 3);
        stl.emplace_back(impl_time());
        stl.emplace_back(omega());
        stl.emplace_back(bindE());
        add_coords_to(stl, impl_pos());
        add_coords_to(stl, impl_vel());
    }

    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    void RegularizedSystem<Particles, Forces, RegType>::impl_load_from_linear_container(const STL &stl) {
        auto i = stl.begin();
        impl_time() = *i, ++i;
        omega() = *i, ++i;
        bindE() = *i, ++i;
        load_to_coords(i, impl_pos());
        load_to_coords(i, impl_vel());
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::eval_vel_indep_acc() {
        forces_.eval_newtonian_acc(ptc_, newtonian_acc_);
        if constexpr (Forces::ext_vel_dep) {
            forces_.eval_extra_vel_indep_acc(ptc_, ext_vel_indep_acc_);
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::advance_omega(const Coord &velocity, const Coord &d_omega_dr,
                                                                      Scalar phy_time) {
        Scalar d_omega = calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_omega_dr);
        regu_.omega() += d_omega * phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::advance_bindE(const Coord &velocity, const Coord &d_bindE_dr,
                                                                      Scalar phy_time) {
        Scalar d_bindE = -calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_bindE_dr);
        regu_.bindE() += d_bindE * phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::kick_pseu_vel(Scalar phy_time) {
        forces_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
        calc::coord_add(acc_, newtonian_acc_, ext_vel_dep_acc_);
        if constexpr (Forces::ext_vel_indep) {
            calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
        }
        calc::coord_advance(aux_vel_, acc_, phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void RegularizedSystem<Particles, Forces, RegType>::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, ptc_.vel());
        forces_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
        std::swap(aux_vel_, ptc_.vel());

        calc::coord_add(acc_, newtonian_acc_, ext_vel_dep_acc_);
        if constexpr (Forces::ext_vel_indep) {
            calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
        }
        calc::coord_advance(ptc_.vel(), acc_, phy_time);

        advance_omega(aux_vel_, newtonian_acc_, phy_time);
        advance_bindE(aux_vel_, ext_vel_dep_acc_, phy_time);
    }

    /*---------------------------------------------------------------------------*\
        Class Regularization Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Scalar, ReguType Type>
    template<typename Particles>
    Regularization<Scalar, Type>::Regularization(const Particles &partc) {
        omega_ = capital_omega(partc);
        bindE_ = -calc::calc_total_energy(partc);
    }

    template<typename Scalar, ReguType Type>
    template<typename Particles>
    auto Regularization<Scalar, Type>::eval_pos_phy_time(const Particles &partc, Scalar step_size) {
        if constexpr (Type == ReguType::logH) {
            return step_size / (bindE_ + calc::calc_kinetic_energy(partc));
        } else if constexpr (Type == ReguType::TTL) {
            return step_size / omega_;
        } else if constexpr (Type == ReguType::none) {
            return step_size;
        } else {
            spacehub_abort("Undefined regularization type!");
        }
    }

    template<typename Scalar, ReguType Type>
    template<typename Particles>
    auto Regularization<Scalar, Type>::eval_vel_phy_time(const Particles &partc, Scalar step_size) {
        if constexpr (Type == ReguType::logH) {
            return step_size / -calc::calc_potential_energy(partc);
        } else if constexpr (Type == ReguType::TTL) {
            return step_size / capital_omega(partc);
        } else if constexpr (Type == ReguType::none) {
            return step_size;
        } else {
            spacehub_abort("Undefined regularization type!");
        }
    }

    template<typename Scalar, ReguType Type>
    template<typename Particles>
    auto Regularization<Scalar, Type>::capital_omega(const Particles &partc) {
        return -calc::calc_potential_energy(partc);
    }
}

#endif
