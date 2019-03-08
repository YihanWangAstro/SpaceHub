//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_ARCHAIN_H
#define SPACEHUB_ARCHAIN_H

#include "../core-computation.h"
#include "../dev-tools.h"
#include "../particle-system.h"
#include "chain.tpp"

namespace SpaceH {

    template<typename Particles, typename Interactions>
    class ChainSystem : public ParticleSystem<ChainSystem<Particles, Interactions>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        SPACEHUB_STD_ACCESSOR(impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(impl_time, ptc_.time());

        SPACEHUB_STD_ACCESSOR(chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(index, index_);

        ChainSystem() = delete;

        template<typename Container>
        ChainSystem(Container const &ptc, Scalar t) : ptc_(ptc, t), chain_pos_(ptc.size()), chain_vel_(ptc.size()), index_(ptc.size()), new_index_(ptc.size()), acc_(ptc.size()) , chain_acc_(ptc.size()){
            chain::calc_chain_index(ptc_.pos(), index_);
            chain::coord_calc_chain(ptc_.pos(), chain_pos(), index());
            chain::coord_calc_chain(ptc_.vel(), chain_vel(), index());
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                extra_vel_indep_acc_.resize(ptc.size());
            }

            if constexpr (Interactions::has_extra_vel_dep_acc) {
                extra_vel_dep_acc_.resize(ptc.size());
                aux_vel_ = ptc_.vel();
                chain_aux_vel_ = chain_vel_;
            }
        }

        size_t impl_number() {
            return ptc_.number();
        }

        void impl_advance_time(Scalar dt) {
            ptc_.time() += dt;
        }

        void impl_advance_pos(Coord const &velocity, Scalar stepSize) {
            chain::coord_calc_chain(velocity, chain_vel(), index());
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), stepSize);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            chain::coord_calc_chain(acceleration, chain_acc_, index());
            chain_advance(ptc_.vel(), chain_vel(), chain_acc_, stepSize);
        }

        void impl_evaluate_acc(Coord &acceleration) {
            eom_.eval_acc(*this, acceleration);
        }

        void impl_drift(Scalar stepSize) {
            ptc_.time() += stepSize;
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                Scalar halfStep = 0.5 * stepSize;
                eval_vel_indep_acc();
                kick_pseu_vel(halfStep);
                kick_real_vel(stepSize);
                kick_pseu_vel(halfStep);
            } else {
                eom_.eval_acc(*this, acc_);
                impl_advance_vel(acc_, stepSize);
            }
        }

        void impl_pre_iter_process() {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                aux_vel_ = ptc_.vel();
                chain_aux_vel_ = chain_vel_;
            }
        }

        void impl_post_iter_process() {
            chain::calc_chain_index(ptc_.pos(), new_index_);
            if(new_index_ != index_){
                chain::update_chain(chain_pos_, index_, new_index_);
                chain::coord_calc_cartesian(chain_pos_, ptc_.pos(), new_index_);
                calc::coord_move_to_com(ptc_.mass(), ptc_.pos());
                chain::update_chain(chain_vel_, index_, new_index_);
                chain::coord_calc_cartesian(chain_vel_, ptc_.vel(), new_index_);
                calc::coord_move_to_com(ptc_.mass(), ptc_.vel());
                index_ = new_index_;
            }
        }

        friend std::ostream &operator<<(std::ostream &os, ChainSystem const &ps) {
            os << ps.ptc_;
        }
    private:
        void chain_advance(Coord &var, Coord& ch_var, Coord & ch_inc, Scalar stepSize) {
            calc::coord_advance(ch_var, ch_inc, stepSize);
            chain::coord_calc_cartesian(ch_var, var, index());
            calc::coord_move_to_com(ptc_.mass(), var);
        }
        void eval_vel_indep_acc() {
            eom_.eval_newtonian_acc(*this, acc_);
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                eom_.eval_extra_vel_indep_acc(*this, extra_vel_indep_acc_);
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
        }

        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_extra_vel_dep_acc(*this, extra_vel_dep_acc_);
            calc::coord_add(acc_, acc_, extra_vel_dep_acc_);
            chain::coord_calc_chain(acc_, chain_acc_, index());
            chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());
            eom_.eval_extra_vel_dep_acc(*this, extra_vel_dep_acc_);
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());

            calc::coord_add(acc_, acc_, extra_vel_dep_acc_);
            chain::coord_calc_chain(acc_, chain_acc_, index());
            chain_advance(ptc_.vel(), chain_vel_(), chain_acc_, stepSize);
        }

        Particles ptc_;
        Interactions eom_;

        Coord chain_pos_{0};
        Coord chain_vel_{0};
        Coord acc_{0};
        Coord chain_acc_{0};

        Coord extra_vel_indep_acc_{0};
        Coord extra_vel_dep_acc_{0};
        Coord aux_vel_{0};
        Coord chain_aux_vel_{0};
        IdxArray index_{0};
        IdxArray new_index_{0};
    };

}
#endif //SPACEHUB_ARCHAIN_H
