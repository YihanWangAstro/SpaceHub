//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_ARCHAIN_H
#define SPACEHUB_ARCHAIN_H

#include "regu-system.tpp"
#include "chain.tpp"

namespace SpaceH {

    template<typename Particles, typename Interactions, ReguType ReguType>
    class ARchainSystem : public ParticleSystem<ARchainSystem<Particles, Interactions, ReguType>> {
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

        ARchainSystem() = delete;

        template<typename Container>
        ARchainSystem(Container const &ptc, Scalar t) : ptc_(ptc, t), regu_(ptc), chain_pos_(ptc.size()), chain_vel_(ptc.size()), index_(ptc.size()), new_index_(ptc.size()), acc_(ptc.size()) , newtonian_acc_(ptc.size()),chain_acc_(ptc.size()){
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

        void impl_advance_time(Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            ptc_.time() += phyTime;
        }

        void impl_advance_pos(Coord const &velocity, Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            chain::coord_calc_chain(velocity, chain_vel(), index());
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), phyTime);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            Scalar phyTime = regu_.eval_vel_phy_time(ptc_, stepSize);
            chain::coord_calc_chain(acceleration, chain_acc_, index());
            chain_advance(ptc_.vel(), chain_vel(), chain_acc_, phyTime);
        }

        void impl_evaluate_acc(Coord &acceleration) {
            eom_.eval_acc(*this, acceleration);
        }

        void impl_drift(Scalar stepSize) {
            Scalar phyTime = regu_.eval_pos_phy_time(ptc_, stepSize);
            ptc_.time() += phyTime;
            chain_advance(ptc_.pos(), chain_pos(), chain_vel(), phyTime);
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

                    chain_advance(ptc_.vel(), chain_vel(), acc_, halfTime);
                    advance_omega(ptc_.vel(), newtonian_acc_, phyTime);
                    advance_bindE(ptc_.vel(), extra_vel_indep_acc_, phyTime);
                    chain_advance(ptc_.vel(), chain_vel(), acc_, halfTime);
                } else {
                    chain_advance(ptc_.vel(), chain_vel(), newtonian_acc_, halfTime);
                    advance_omega(ptc_.vel(), newtonian_acc_, phyTime);
                    chain_advance(ptc_.vel(), chain_vel(), newtonian_acc_, halfTime);
                }
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

        friend std::ostream &operator<<(std::ostream &os, ARchainSystem const &ps) {
            os << ps.ptc_;
        }
    private:
        void chain_advance(Coord &var, Coord& ch_var, Coord & ch_inc, Scalar stepSize) {
            calc::coord_advance(ch_var, ch_inc, stepSize);
            chain::coord_calc_cartesian(ch_var, var, index());
            calc::coord_move_to_com(ptc_.mass(), var);
        }

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
            chain_advance(aux_vel_, chain_aux_vel_, acc_, stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());
            eom_.eval_extra_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            std::swap(aux_vel_, ptc_.vel());
            std::swap(chain_aux_vel_, chain_vel());

            calc::coord_add(acc_, newtonian_acc_, extra_vel_dep_acc_);
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }

            chain_advance(ptc_.vel(), chain_vel(), acc_, stepSize);

            advance_omega(aux_vel_, newtonian_acc_, stepSize);
            advance_bindE(aux_vel_, extra_vel_dep_acc_, stepSize);
        }


        Particles ptc_;
        Interactions eom_;
        Regularization<Scalar, ReguType> regu_;

        Coord chain_pos_{0};
        Coord chain_vel_{0};
        Coord acc_{0};
        Coord chain_acc_{0};

        Coord newtonian_acc_{0};
        Coord extra_vel_indep_acc_{0};
        Coord extra_vel_dep_acc_{0};
        Coord aux_vel_{0};
        Coord chain_aux_vel_{0};
        IdxArray index_{0};
        IdxArray new_index_{0};
    };
}
#endif //SPACEHUB_ARCHAIN_H
