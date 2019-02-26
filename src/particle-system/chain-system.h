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

        ChainSystem() = delete;

        template<typename Container>
        ChainSystem(Container const &partc, Scalar t) : partc_(partc, t), action_(partc.size()) {
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

        void impl_advance_time(Scalar dt) {
            partc_.time() += dt;
        }

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(vx(), vy(), vz(), stepSize);
        }

        void impl_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            calcu::array_advance(px(), vx, stepSize);
            calcu::array_advance(py(), vy, stepSize);
            calcu::array_advance(pz(), vz, stepSize);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(action_.ax(), action_.ay(), action_.az(), stepSize);
        }

        void impl_advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            calcu::array_advance(vx(), ax, stepSize);
            calcu::array_advance(vy(), ay, stepSize);
            calcu::array_advance(vz(), az, stepSize);
        }

        void impl_evaluate_acc() {
            action_.eval_acc(*this);
        }

        void impl_drift(Scalar stepSize) {
            impl_advance_time(stepSize);
            impl_advance_pos(stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interaction::isVelDependent) {
                action_.eval_vel_indep_acc(*this);

                action_.eval_vel_dep_acc(*this);
                action_.sum_tot_acc();
                calcu::array_advance(auxi_vx(), ax(), 0.5 * stepSize);
                calcu::array_advance(auxi_vy(), ay(), 0.5 * stepSize);
                calcu::array_advance(auxi_vz(), az(), 0.5 * stepSize);

                action_.eval_auxi_vel_dep_acc(*this);
                action_.sum_tot_acc();
                impl_advance_vel(stepSize);

                action_.eval_vel_dep_acc(*this);
                action_.sum_tot_acc();
                calcu::array_advance(auxi_vx(), ax(), 0.5 * stepSize);
                calcu::array_advance(auxi_vy(), ay(), 0.5 * stepSize);
                calcu::array_advance(auxi_vz(), az(), 0.5 * stepSize);
            } else {
                action_.eval_acc(*this);
                impl_advance_vel(stepSize);
            }
        }

        friend std::ostream&operator<<(std::ostream& os, ChainSystem const & ps) {
            os << ps.partc_;
        }


        Particles partc_;
        Interaction action_;
        ScalarArray chain_px;
        ScalarArray chain_py;
        ScalarArray chain_pz;
        ScalarArray chain_vx;
        ScalarArray chain_vy;
        ScalarArray chain_vz;
        IndexArray index;
    };

}
#endif //SPACEHUB_ARCHAIN_H
