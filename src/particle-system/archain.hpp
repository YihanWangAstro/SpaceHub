//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_ARCHAIN_H
#define SPACEHUB_ARCHAIN_H

#include "regu-system.hpp"
#include "chain.hpp"

namespace space {

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Declaration
    \*---------------------------------------------------------------------------*/
    template<typename Particles, typename Forces, ReguType RegType>
    class ARchainSystem : public ParticleSystem<ARchainSystem<Particles, Forces, RegType>> {
    public:
        //Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Base = ParticleSystem<ARchainSystem<Particles, Forces, RegType>>;

        using Particle = typename Particles::Particle;

        //Static public members
        static constexpr ReguType regu_type{RegType};

        ARchainSystem() = delete;

        ARchainSystem(ARchainSystem const &) = default;

        ARchainSystem(ARchainSystem &&) noexcept = default;

        ARchainSystem &operator=(ARchainSystem const &) = default;

        ARchainSystem &operator=(ARchainSystem &&) noexcept = default;

        //Constructors
        template<typename STL>
        ARchainSystem(Scalar t, STL const &ptc);

        //Public methods
        SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(auto, index, index_);

        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

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

        void impl_evaluate_acc(Coord &acceleration) const;

        void impl_drift(Scalar step_size);

        void impl_kick(Scalar step_size);

        void impl_pre_iter_process();

        void impl_post_iter_process();

        template<typename STL>
        void impl_to_linear_container(STL &stl);

        template<typename STL>
        void impl_load_from_linear_container(STL const &stl);

        //Friend functions
        template<typename P, typename F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, ARchainSystem<P, F, R> const &ps);

        template<typename P, typename F, ReguType R>
        friend std::istream &operator>>(std::istream &is, ARchainSystem<P, F, R> &ps);

    private:
        //Private methods
        void chain_advance(Coord &var, Coord &ch_var, Coord const &ch_inc, Scalar phy_time);

        void eval_vel_indep_acc();

        void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar phy_time);

        void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        //Private members
        Particles ptc_;

        Forces forces_;

        Regularization <Scalar, RegType> regu_;

        Coord chain_pos_;

        Coord chain_vel_;

        Coord acc_;

        Coord chain_acc_;

        Coord newtonian_acc_;

        IdxArray index_;

        IdxArray new_index_;

        std::conditional_t<Forces::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;

        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;

        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> aux_vel_;

        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> chain_aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    ARchainSystem<Particles, Forces, RegType>::ARchainSystem(Scalar t, const STL &ptc)
            : ptc_(t, ptc),
              regu_(ptc_),
              chain_pos_(ptc.size()),
              chain_vel_(ptc.size()),
              index_(ptc.size()),
              new_index_(ptc.size()),
              acc_(ptc.size()),
              newtonian_acc_(ptc.size()),
              chain_acc_(ptc.size()) {
        static_assert(is_container_v<STL>, "Only STL-like container can be used");
        Chain::calc_chain_index(ptc_.pos(), index_);
        Chain::calc_chain(ptc_.pos(), chain_pos(), index());
        Chain::calc_chain(ptc_.vel(), chain_vel(), index());
        if constexpr (Forces::ext_vel_indep) {
            ext_vel_indep_acc_.resize(ptc.size());
        }

        if constexpr (Forces::ext_vel_dep) {
            ext_vel_dep_acc_.resize(ptc.size());
            aux_vel_ = ptc_.vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    size_t ARchainSystem<Particles, Forces, RegType>::impl_number() const {
        return ptc_.number();
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_advance_time(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        ptc_.time() += phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        Chain::calc_chain(velocity, chain_vel(), index());
        chain_advance(ptc_.pos(), chain_pos(), chain_vel(), phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
        Scalar phy_time = regu_.eval_vel_phy_time(ptc_, step_size);
        Chain::calc_chain(acceleration, chain_acc_, index());
        chain_advance(ptc_.vel(), chain_vel(), chain_acc_, phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_evaluate_acc(Coord &acceleration) const {
        forces_.eval_acc(*this, acceleration);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(ptc_, step_size);
        chain_advance(ptc_.pos(), chain_pos(), chain_vel(), phy_time);
        ptc_.time() += phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_kick(Scalar step_size) {
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
                Chain::calc_chain(acc_, chain_acc_, index());
                chain_advance(ptc_.vel(), chain_vel(), chain_acc_, half_time);
                advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                advance_bindE(ptc_.vel(), ext_vel_indep_acc_, phy_time);
                chain_advance(ptc_.vel(), chain_vel(), chain_acc_, half_time);
            } else {
                Chain::calc_chain(newtonian_acc_, chain_acc_, index());
                chain_advance(ptc_.vel(), chain_vel(), chain_acc_, half_time);
                advance_omega(ptc_.vel(), newtonian_acc_, phy_time);
                chain_advance(ptc_.vel(), chain_vel(), chain_acc_, half_time);
            }
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_pre_iter_process() {
        if constexpr (Forces::ext_vel_dep) {
            aux_vel_ = ptc_.vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::impl_post_iter_process() {
        Chain::calc_chain_index(ptc_.pos(), new_index_);
        if (new_index_ != index_) {
            Chain::update_chain(chain_pos_, index_, new_index_);
            Chain::calc_cartesian(ptc_.mass(), chain_pos_, ptc_.pos(), new_index_);
            Chain::update_chain(chain_vel_, index_, new_index_);
            Chain::calc_cartesian(ptc_.mass(), chain_vel_, ptc_.vel(), new_index_);
            index_ = new_index_;
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    void ARchainSystem<Particles, Forces, RegType>::impl_to_linear_container(STL &stl) {
        stl.clear();
        stl.reserve(impl_number() * 6 + 3);
        stl.emplace_back(impl_time());
        stl.emplace_back(omega());
        stl.emplace_back(bindE());
        add_coords_to(stl, chain_pos_);
        add_coords_to(stl, chain_vel_);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    template<typename STL>
    void ARchainSystem<Particles, Forces, RegType>::impl_load_from_linear_container(const STL &stl) {
        auto i = stl.begin();
        impl_time() = *i, ++i;
        omega() = *i, ++i;
        bindE() = *i, ++i;
        load_to_coords(i, chain_pos_);
        load_to_coords(i, chain_vel_);
        Chain::calc_cartesian(ptc_.mass(), chain_pos_, impl_pos(), index());
        Chain::calc_cartesian(ptc_.mass(), chain_vel_, impl_vel(), index());
    }

    template<typename Particles, typename Forces, ReguType RegType>
    std::ostream &operator<<(std::ostream &os, ARchainSystem<Particles, Forces, RegType> const &ps) {
        os << ps.ptc_;
        return os;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    std::istream &operator>>(std::istream &is, ARchainSystem<Particles, Forces, RegType> &ps) {
        is >> ps.ptc_;
        return is;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::chain_advance(Coord &var, Coord &ch_var, const Coord &ch_inc,
                                                                  Scalar phy_time) {
        calc::coord_advance(ch_var, ch_inc, phy_time);
        Chain::calc_cartesian(ptc_.mass(), ch_var, var, index());
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::eval_vel_indep_acc() {
        forces_.eval_newtonian_acc(ptc_, newtonian_acc_);
        if constexpr (Forces::ext_vel_dep) {
            forces_.eval_extra_vel_indep_acc(ptc_, ext_vel_indep_acc_);
        }
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::advance_omega(const Coord &velocity, const Coord &d_omega_dr,
                                                                  Scalar phy_time) {
        Scalar d_omega = calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_omega_dr);
        regu_.omega() += d_omega * phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::advance_bindE(const Coord &velocity, const Coord &d_bindE_dr,
                                                                  Scalar phy_time) {
        Scalar d_bindE = -calc::coord_contract_to_scalar(ptc_.mass(), velocity, d_bindE_dr);
        regu_.bindE() += d_bindE * phy_time;
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::kick_pseu_vel(Scalar phy_time) {
        forces_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
        calc::coord_add(acc_, newtonian_acc_, ext_vel_dep_acc_);
        if constexpr (Forces::ext_vel_indep) {
            calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
        }
        Chain::calc_chain(acc_, chain_acc_, index());
        chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, phy_time);
    }

    template<typename Particles, typename Forces, ReguType RegType>
    void ARchainSystem<Particles, Forces, RegType>::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, ptc_.vel());
        std::swap(chain_aux_vel_, chain_vel());
        forces_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
        std::swap(aux_vel_, ptc_.vel());
        std::swap(chain_aux_vel_, chain_vel());
        calc::coord_add(acc_, newtonian_acc_, ext_vel_dep_acc_);
        if constexpr (Forces::ext_vel_indep) {
            calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
        }

        Chain::calc_chain(acc_, chain_acc_, index());
        chain_advance(ptc_.vel(), chain_vel(), chain_acc_, phy_time);
        advance_omega(aux_vel_, newtonian_acc_, phy_time);
        advance_bindE(aux_vel_, ext_vel_dep_acc_, phy_time);
    }
}
#endif //SPACEHUB_ARCHAIN_H
