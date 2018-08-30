
#ifndef GARPARTICLESYSTEM_H
#define GARPARTICLESYSTEM_H

#include "../../particleSystemC++17.h"
#include "coreComputation.h"
#include "devTools.h"

namespace SpaceH {
/**  @brief General midpoint method class of particle System. Used in velocity dependent interaction system to create
 *   symmetry leapfrog scheme. See details in http://adsabs.harvard.edu/abs/2006MNRAS.372..219M
 */
    template<typename Particles, typename Interaction>
    class GMParticleSystem : public ParticleSystem<Particles, Interaction> {
    public:
        /* Typedef */
        using Base         = typename ParticleSystem<Particles, Interaction>;
        using type         = typename Base::type;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarArray  = typename type::ScalarArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        using ParticleType = Particles;
        /* Typedef */

        /*Template parameter check*/
        CHECK_TYPE(Particles, Interaction);
        static_assert(!Interaction::isVelDep, "General midpoint method is unecessary for velocity independent system.")
        /*Template parameter check*/

        constexpr static bool isVelDep{Interaction::isVelDep};
        constexpr static size_t arraySize{type::arraySize};

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            act.zeroTotalAcc();
            //evaluate velocity independent acc
            act.calcuVelIndepAcc(partc);
            act.calcuExtVelIndepAcc(partc);

            //evaluate velocity dependent acc with velocity
            act.calcuVelDepAcc(partc);
            act.calcuExtVelDepAcc(partc);
            act.sumTotalAcc();

            SpaceH::advanceVector(auxi_vel, act.totalAcc(), 0.5*stepSize);//advance auxiliary velocity

            //evaluate velocity dependent acc with auxiliary velocity
            /*swap for std::vector is O(1), however, for std::array is O(n). This is kinda overhead. Opt in the future.*/
            partc.vel.swap(auxi_vel);//swap vel and auxi_vel. after this line, auxi_vel is vel.
            act.calcuVelDepAcc(partc);
            act.calcuExtVelDepAcc(partc);
            act.sumTotalAcc();
            partc.vel.swap(auxi_vel);//swap vel and auxi_vel. after this line, every is normal.

            SpaceH::advanceVector(partc.vel, act.totalAcc(), stepSize);//advance velocity with auxi_vel evaluated acc

            act.calcuVelDepAcc(partc);
            act.calcuExtVelDepAcc(partc);
            act.sumTotalAcc();
            SpaceH::advanceVector(auxi_vel, act.totalAcc(), 0.5*stepSize);//advance auxiliary velocity
        }

        /** @brief After process after iteration*/
        void afterIterProcess() {
            auxi_vel = partc.vel;
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &is, GMParticleSystem &sys) {
            is >> static_cast<Base &>(sys);
            auxi_vel = partc.vel;
            return is;
        }
    protected:
        /** @brief auxiliary velocity array used in General midpoint method*/
        Vector auxi_vel;
    };
}

#endif
