
#ifndef GARPARTICLESYSTEM_H
#define GARPARTICLESYSTEM_H

#include "regu_system.h"
#include "../core_computation.h"
#include "../dev_tools.h"

namespace SpaceH {
/**  @brief General midpoint method class of particle System. Used in velocity dependent interaction system to create
 *   symmetry leapfrog scheme. See details in http://adsabs.harvard.edu/abs/2006MNRAS.372..219M
 *   @note The call of afterIterProcess() after each iteraction is required to synchronize the auxiliary velocity array.
 */
    template<typename Particles, typename Interaction, typename Regularitor>
    class GAR_system : public RegularizedSystem<Particles, Interaction, Regularitor> {
    public:
        /* Typedef */
        using Base         = RegularizedSystem<Particles, Interaction, Regularitor>;
        using type         = typename Base::type;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarArray  = typename type::ScalarArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        using State        = typename Particles::State;
        using ParticleType = Particles;
        /* Typedef */

        using Base::partc;
        using Base::act;
        using Base::regular;

        /*Template parameter check*/
        CHECK_TYPE(Particles, Interaction);
        /*Template parameter check*/

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            Scalar physicalTime = regular.getPhysicalVelTime(partc, stepSize);

            act.zeroTotalAcc();
            act.calcuVelIndepAcc(partc);//evaluate velocity independent acc

            if constexpr (!Interaction::isVelDep) {
                act.sumTotalAcc();
                partc.advanceVel(act.acc(), 0.5*physicalTime);
                regular.advanceOmega(act.pairVelIndepAcc(), partc.vel(), partc.mass(), physicalTime);
                partc.advanceVel(act.acc(), 0.5*physicalTime);
                //regular.advanceOmega(act.pairVelIndepAcc(), partc.vel(), partc.mass(), 0.5*physicalTime);
            } else {
                act.calcuVelDepAcc(partc);//evaluate velocity dependent acc with velocity
                act.sumTotalAcc();

                //evaluate velocity dependent acc with auxiliary velocity
                /*swap for std::vector is O(1), however, for std::array is O(n). This is kinda overhead. Opt in the future.*/
                partc.swap_vel_state(auxi_vel);//swap vel and auxi_vel. after this line, auxi_vel is vel.
                partc.advanceVel(act.acc(), 0.5*physicalTime);//advance auxiliary velocity
                act.calcuVelDepAcc(partc);
                act.sumTotalAcc();
                partc.swap_vel_state(auxi_vel);//swap vel and auxi_vel. after this line, everything is normal.

                partc.advanceVel(act.acc(), physicalTime);//advance velocity with auxi_vel evaluated acc
                regular.advanceOmega(act.pairVelIndepAcc(), auxi_vel.cartesian(), partc.mass(), physicalTime);
                regular.advanceBindE(act.pairVelDepAcc(),   auxi_vel.cartesian(), partc.mass(), physicalTime);

                act.calcuVelDepAcc(partc);
                act.sumTotalAcc();
                partc.swap_vel_state(auxi_vel);
                partc.advanceVel(act.acc(), 0.5*physicalTime);//advance auxiliary velocity
                partc.swap_vel_state(auxi_vel);
            }
        }

        /** @brief After process after iteration*/
        void afterIterProcess() {
            Base::afterIterProcess();
            if constexpr (Interaction::isVelDep) {
                auxi_vel = partc.vel_state();
            }
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &is, GAR_system &sys) {
            is >> static_cast<Base &>(sys);
            if constexpr (Interaction::isVelDep) {
                sys.auxi_vel = sys.partc.vel_state();
            }
            return is;
        }
    protected:
        /** @brief auxiliary velocity array used in General midpoint method*/
        State auxi_vel;
    };
}

#endif
