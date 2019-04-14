//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_CHAINSYSTEM_H
#define SPACEHUB_CHAINSYSTEM_H

#include "../core-computation.tpp"
#include "../dev-tools.h"
#include "../particle-system.h"
#include "chain.tpp"
#include <type_traits>

namespace space {

    template<typename Particles, typename Forces>
    class ChainSystem : public ParticleSystem<ChainSystem<Particles, Forces>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(auto, impl_time, ptc_.time());

        SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(auto, index, index_);

        ChainSystem() = delete;

        ChainSystem(ChainSystem const &) = default;

        ChainSystem(ChainSystem &&) = default;

        ChainSystem &operator=(ChainSystem const &) = default;

        ChainSystem &operator=(ChainSystem &&) = default;

        template<typename STL>
        ChainSystem(Scalar t, STL const &ptc)
                : ptc_(t, ptc),
                  chain_pos_(ptc.size()),
                  chain_vel_(ptc.size()),
                  index_(ptc.size()),
                  new_index_(ptc.size()),
                  acc_(ptc.size()),
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

        size_t impl_number() const {
            return ptc_.number();
        }

        void impl_advance_time(Scalar dt) {
            ptc_.time() += dt;
        }

        void impl_advance_pos(Coord const &velocity, Scalar step_size) {
            Chain::calc_chain(velocity, chain_vel(), index());
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), step_size);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar step_size) {
            Chain::calc_chain(acceleration, chain_acc_, index());
            chain_advance(ptc_.vel(), chain_vel(), chain_acc_, step_size);
        }

        void impl_evaluate_acc(Coord &acceleration) const {
            forces_.eval_acc(*this, acceleration);
        }

        void impl_drift(Scalar step_size) {
            ptc_.time() += step_size;
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), step_size);
        }

        void impl_kick(Scalar step_size) {
            if constexpr (Forces::ext_vel_dep) {
                Scalar half_step = 0.5 * step_size;
                eval_vel_indep_acc();
                kick_pseu_vel(half_step);
                kick_real_vel(step_size);
                kick_pseu_vel(half_step);
            } else {
                forces_.eval_acc(*this, acc_);
                impl_advance_vel(acc_, step_size);
            }
        }

        void impl_pre_iter_process() {
            if constexpr (Forces::ext_vel_dep) {
                aux_vel_ = ptc_.vel();
                chain_aux_vel_ = chain_vel_;
            }
        }

        void impl_post_iter_process() {
            Chain::calc_chain_index(ptc_.pos(), new_index_);
            if (new_index_ != index_) {
                Chain::update_chain(chain_pos_, index_, new_index_);
                Chain::calc_cartesian(ptc_.mass(), chain_pos_, ptc_.pos(), new_index_);
                Chain::update_chain(chain_vel_, index_, new_index_);
                Chain::calc_cartesian(ptc_.mass(), chain_vel_, ptc_.vel(), new_index_);
                index_ = new_index_;
            }
        }

        template<typename STL>
        void impl_to_linear_container(STL &stl) {
            stl.clear();
            stl.reserve(impl_number() * 6 + 1);
            stl.emplace_back(impl_time());
            add_coords_to(stl, chain_pos_);
            add_coords_to(stl, chain_vel_);
        }

        template<typename STL>
        void impl_load_from_linear_container(STL const &stl) {
            auto i = stl.begin();
            impl_time() = *i, ++i;
            load_to_coords(i, chain_pos_);
            load_to_coords(i, chain_vel_);
            Chain::calc_cartesian(ptc_.mass(), chain_pos_, ptc_.pos(), index_);
            Chain::calc_cartesian(ptc_.mass(), chain_vel_, ptc_.vel(), index_);
        }

        friend std::ostream &operator<<(std::ostream &os, ChainSystem const &ps) {
            os << ps.ptc_;
            return os;
        }

        friend std::istream &operator>>(std::istream &is, ChainSystem &ps) {
            is >> ps.ptc_;
            return is;
        }

    private:
        void chain_advance(Coord &var, Coord &ch_var, Coord &ch_inc, Scalar step_size) {
            calc::coord_advance(ch_var, ch_inc, step_size);
            Chain::calc_cartesian(ptc_.mass(), ch_var, var, index());
        }

        void eval_vel_indep_acc() {
            forces_.eval_newtonian_acc(*this, acc_);
            if constexpr (Forces::ext_vel_dep) {
                forces_.eval_extra_vel_indep_acc(*this, ext_vel_indep_acc_);
                calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
            }
        }

        void kick_pseu_vel(Scalar step_size) {
            forces_.eval_extra_vel_dep_acc(*this, ext_vel_dep_acc_);
            calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
            Chain::calc_chain(acc_, chain_acc_, index());
            chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, step_size);
        }

        void kick_real_vel(Scalar step_size) {
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());
            forces_.eval_extra_vel_dep_acc(*this, ext_vel_dep_acc_);
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());

            calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
            Chain::calc_chain(acc_, chain_acc_, index());
            chain_advance(ptc_.vel(), chain_vel_(), chain_acc_, step_size);
        }

        Particles ptc_;
        Forces forces_;

        Coord chain_pos_;
        Coord chain_vel_;
        Coord acc_;
        Coord chain_acc_;

        IdxArray index_;
        IdxArray new_index_;

        std::conditional_t<Forces::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;
        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;
        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> aux_vel_;
        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> chain_aux_vel_;
    };

}
#endif //SPACEHUB_ARCHAIN_H