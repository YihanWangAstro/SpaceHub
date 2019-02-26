//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_ARCHAIN_H
#define SPACEHUB_ARCHAIN_H

#include "regu-system.h"

namespace SpaceH {

    template<typename Particles, typename Interaction>
    class ChainSystem : public ParticleSystem<ChainSystem<Particles, Interaction>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        ChainSystem() = delete;

        template<typename Container>
        ChainSystem(Container const &partc, Scalar t) : partc_(partc, t), action_(partc.size()) {
            if constexpr (Interaction::isVelDependent) {
                action_.auxi_vx = partc_.vx;
                action_.auxi_vy = partc_.vy;
                action_.auxi_vz = partc_.vz;
            }
        }

        size_t impl_number() {
            return partc_.number();
        }

        void impl_advance_time(Scalar dt) {
            partc_.time_ += dt;
        }

        void impl_advance_pos(Scalar stepSize) {
            advance_pos(partc_.vx, partc_.vy, partc_.vz, stepSize);
        }

        void impl_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            calcu::array_advance(partc_.px, partc_.vx, stepSize);
            calcu::array_advance(partc_.py, partc_.vy, stepSize);
            calcu::array_advance(partc_.pz, partc_.vz, stepSize);
        }

        void impl_advance_vel(Scalar stepSize) {
            advance_vel(action_.ax, action_.ay, action_.az, stepSize);
        }

        void impl_advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            calcu::array_advance(partc_.vx, ax, stepSize);
            calcu::array_advance(partc_.vy, ay, stepSize);
            calcu::array_advance(partc_.vz, az, stepSize);
        }

        void impl_evaluate_acc() {
            action_.eval_acc(partc_);
        }

        void impl_drift(Scalar stepSize) {
            advance_time(stepSize);
            advance_pos(stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interaction::isVelDependent) {
                action_.eval_vel_indep_acc(partc_);

                action_.eval_vel_dep_acc(partc_);
                action_.sum_tot_acc();
                calcu::array_advance(action_.auxi_vx, action_.ax, 0.5 * stepSize);
                calcu::array_advance(action_.auxi_vy, action_.ay, 0.5 * stepSize);
                calcu::array_advance(action_.auxi_vz, action_.az, 0.5 * stepSize);

                action_.eval_auxi_vel_dep_acc(partc_);
                action_.sum_tot_acc();
                advance_vel(stepSize);

                action_.eval_vel_dep_acc(partc_);
                action_.sum_tot_acc();
                calcu::array_advance(action_.auxi_vx, action_.ax, 0.5 * stepSize);
                calcu::array_advance(action_.auxi_vy, action_.ay, 0.5 * stepSize);
                calcu::array_advance(action_.auxi_vz, action_.az, 0.5 * stepSize);
            } else {
                action_.eval_acc(partc_);
                advance_vel(stepSize);
            }
        }

        friend std::ostream&operator<<(std::ostream& os, BaseSystem const & ps) {
            os << ps.partc_;
        }


        Particles partc_;
        ScalarArray chain_px;
        ScalarArray chain_py;
        ScalarArray chain_pz;
        ScalarArray chain_vx;
        ScalarArray chain_vy;
        ScalarArray chain_vz;
        IndexArray index;
        Interaction action_;
        Regularization<Scalar, ReguType> regular_;
    };

}
#endif //SPACEHUB_ARCHAIN_H
