
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../particle_system.h"
#include "../core_computation.h"
#include "../dev_tools.h"

namespace SpaceH {

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
     * http://iopscience.iop.org/article/10.1086/301102/meta and
     * https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
     * @tparam Particles
     * @tparam Interaction
     */
    template<typename Particles, typename Interaction, typename Regularitor>
    class RegularizedSystem : public ParticleSystem<Particles, Interaction> {
    public:
        /* Typedef */
        using State        = typename Particles::State;
        using ParticleType = Particles;
        using Base         = ParticleSystem<Particles, Interaction>;
        SPACEHUB_USING_TYPE_SYSTEM_OF(Base);
        /* Typedef */

        using Base::partc;
        using Base::act;

        /*Template parameter check*/
        CHECK_TYPE(Particles, Interaction);
        CHECK_TYPE(Particles, Regularitor);
        /*Template parameter check*/

        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_SCALAR(omega, Scalar, regular, omega);
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_SCALAR(bindE, Scalar, regular, bindE);

        void advanceOmega(const VectorArray& velIndepAcc, const VectorArray& vel, const ScalarArray& mass, Scalar physicalTime){
            regular.advanceOmegaO(velIndepAcc, vel, mass, physicalTime);
        };
        void advanceBindE(const VectorArray& velDepAcc, const VectorArray& vel, const ScalarArray& mass, Scalar physicalTime){
            regular.advanceBindE(velDepAcc, vel, mass, physicalTime);
        }
        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void drift(Scalar stepSize) {
            Scalar physicalTime = regular.getPhysicalPosTime(partc, stepSize);
            Base::drift(physicalTime);
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            Scalar physicalTime = regular.getPhysicalVelTime(partc, stepSize);

            act.zeroTotalAcc();
            act.calcuVelIndepAcc(partc);
            if constexpr (!Interaction::isVelDep){
                act.sumTotalAcc();
                partc.advenceVel(act.acc(), 0.5*physicalTime);
                regular.advanceOmega(act.pairVelIndepAcc(), partc.vel(), partc.mass(), physicalTime);
                partc.advanceVel(act.acc(), 0.5*physicalTime);
            } else {
                act.calcuVelDepAcc(partc);
                act.sumTotalAcc();
                State v0 = partc.vel_state();
                partc.advenceVel(act.acc(), 0.5 * physicalTime);
                Base::iterateVeltoConvergent(v0, 0.5 * physicalTime);
                regular.advanceOmega(act.pairVelIndepAcc(), partc.vel(), partc.mass(), physicalTime);
                regular.advanceBindE(act.pairVelDepAcc(),   partc.vel(), partc.mass(), physicalTime);
                partc.advenceVel(act.acc(), 0.5 * physicalTime);
            }
        }

        /** @brief Interface to rescale the time.
             *
             *  Interace used by dynamic system. Transfer integration time to physical time.
             *  @return The phsyical time.
             */
        Scalar timeScale() {
            return Base::timeScale() / regular.getPhysicalPosTime(partc, 1);
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &is, RegularizedSystem &sys) {
            is >> static_cast<Base &>(sys);
            if constexpr (Regularitor::nonTrivial) {
                sys.regular.init(sys.partc);
            }
            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {
            size_t loc = static_cast<Base &>(*this).read(data, flag);
            if constexpr (Regularitor::nonTrivial) {
                if (flag == IO_flag::EVOLVED ) {
                    regular.set_omega(data[loc++]);
                    regular.set_bindE(data[loc++]);
                }
            }
            return loc;
        }

        /** @brief Output variables to plain scalar array.*/
        size_t write(ScalarBuffer &data, const IO_flag flag = IO_flag::STD) const {
            size_t loc = static_cast<const Base &>(*this).write(data, flag);
            if constexpr (Regularitor::nonTrivial) {
                if (flag == IO_flag::EVOLVED ) {
                    data.emplace_back(regular.omega());
                    data.emplace_back(regular.bindE());
                }
            }
            return data.size();
        }

    protected:
        /** @brief Regularization interface.*/
        Regularitor regular;
    };
}

#endif
