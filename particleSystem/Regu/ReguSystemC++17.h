
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include "../../particleSystemC++17.h"
#include "coreComputation.h"
#include "devTools.h"

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
    class ReguParticleSystem : public ParticleSystem<Particles, Interaction> {
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
        CHECK_TYPE(Particles, Regularitor);
        /*Template parameter check*/

        constexpr static bool isVelDep{Interaction::isVelDep};
        constexpr static size_t arraySize{type::arraySize};

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void drift(Scalar stepSize) {
            Scalar physicalTime = regular.getPhysicalPosTime(partc, stepSize);
            Base::drift(physicalTime);
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            Scalar physicalTime = regular.getPhysicalVelTime(partc, stepSize);

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
        friend std::istream &operator>>(std::istream &is, ReguParticleSystem &sys) {
            is >> static_cast<Base &>(sys);
            omega = regular.Omega(partc);
            bindE = regular.bindE(partc);
            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {
            size_t loc = static_cast<Base &>(*this).read(data, flag);
            if (flag == IO_flag::EVOLVED) {
                omega = data[loc++];
                bindE = data[loc++];
            }
            return loc;
        }

        /** @brief Output variables to plain scalar array.*/
        size_t write(ScalarBuffer &data, const IO_flag flag = IO_flag::STD) const {
            size_t loc = static_cast<const Base &>(*this).write(data, flag);
            if (flag == IO_flag::EVOLVED) {
                data.reserve(loc + 2);
                data.emplace_back(omega);
                data.emplace_back(bindE);
            }
            return data.size();
        }

    private:
        /** @brief Regularied variables*/
        Scalar omega;

        /** @brief Regularied variables*/
        Scalar bindE;

        /** @brief Regularization interface.*/
        Regularitor regular;
    }
}

#endif
