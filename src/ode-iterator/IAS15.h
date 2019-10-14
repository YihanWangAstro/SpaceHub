
#ifndef IAS15ITERATOR_H
#define IAS15ITERATOR_H

#include "../dev-tools.hpp"
#include "../own-math.hpp"
namespace space {

    namespace Radau {
        constexpr double maxStepCof = 1.74867862159014;// = 1/pow(0.02,1/7)
        constexpr double minStepCof = 0.5718603679678214;// = pow(0.02,1/7)
        constexpr double maxFloat = 1e100;
    }

/** @brief IAS15 iterator see details in https://arxiv.org/abs/1409.4779 .
 *
 */
    template<typename ParticSys, typename Integrator>
    class IAS15 {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
        using RadauArray   = typename Integrator::RadauArray;
        using RadauTab     = typename Integrator::RadauTab;
        /* Typedef */
        /* Typedef */

        /*Template parameter check*/
        CHECK_TYPE(ParticSys, Integrator);
        /*Template parameter check*/

        /** @brief interface to iterate particle system for one step
         *  @param particles  Particle system needs evolution.
         *  @param integrator Integrator to integrate the particle system.
         *  @param stepLength Macro step length for iteration(Here, the step length of the integrator).
         *  @return step length for next iteration.
         */
        Scalar iterate(ParticSys &particles, Scalar stepLength) {
            integrator_.checkTabVolume(particles.particleNumber());

            Scalar iterH = stepLength;
            RadauTab iterBTab = integrator_.getBTab();//get the b value of the last step_sequence.

            resetLastConvergence();
            for (size_t k = 0; k < max_iter_; ++k) {
                integrator_.calcuBTab(particles, iterH);

                if (isConvergent(iterBTab, integrator_.getBTab(), integrator_.localAcc())) {
                    Scalar error = calcuBError(integrator_.getBTab(), integrator_.localAcc());
                    Scalar stepQ = optimalStepCoef(error);

                    if (error < 1) {
                        DEBUG_MSG(true, "accept: error= ", error, "iter = ", k, "stepSize = ", iterH);
                        integrator_.evaluateSystemAt(particles, iterH, Integrator::finalPoint);
                        integrator_.predictNewB(stepQ);
                        iterH *= stepQ;
                        return iterH;
                    } else//current stepSize is too large, restart the iteration with smaller iterH that has been determined by current error.
                    {
                        DEBUG_MSG(true, "reject: error= ", error, "iter = ", k, "stepSize = ", iterH);
                        iterH *= stepQ;
                        integrator_.predictNewB(stepQ);
                        k = 0;
                    }
                }
                iterBTab = integrator_.getBTab();
            }
            SPACEHUB_ABORT("IAS15: iteration exceed the max iteration depth!");
            return 0;
        }

        /**
         *
         * @param rel_tol
         */
        void setRelativeError(Scalar rel_tol) {
            relativeError_ = rel_tol;
        }

        /**
         *
         * @param limit
         */
        void SetConvergentLimit(Scalar limit) {
            convergent_limit_ = limit;
        }

    private:
        /**
         *
         * @param BTab
         * @param newBTab
         * @param acc
         * @return
         */
        bool isConvergent(const RadauTab &BTab, const RadauTab &newBTab, const VectorArray &acc) {
            size_t size = acc.size();

            Scalar diff = 0;
            Scalar scale = 0;
            for (size_t i = 0; i < size; ++i) {
                diff = space::max(diff, (BTab[i][6] - newBTab[i][6]).abs().max_component());
                scale = space::max(scale, acc[i].abs().max_component());
            }

            Scalar convergence = diff / scale;

            DEBUG_MSG(true, "ck convergence: error= ", convergence);

            if (convergence >= last_convergence_)//begin to oscillate
            {
                resetLastConvergence();
                return true;
            } else {
                last_convergence_ = convergence;
                return convergence < convergent_limit_;
            }
        }

        /**
         *
         * @param BTab
         * @param acc
         * @return
         */
        Scalar calcuBError(const RadauTab &BTab, const VectorArray &acc) {
            size_t size = acc.size();
            Scalar diff = 0;
            Scalar scale = 0;
            for (size_t i = 0; i < size; ++i) {
                diff = space::max(diff, BTab[i][6].abs().max_component());
                scale = space::max(scale, acc[i].abs().max_component());
            }
            return diff / (scale * relativeError_);
        }

        /**
         *
         * @param error
         * @return
         */
        inline Scalar optimalStepCoef(Scalar error) {
            if (error == 0)
                return Radau::maxStepCof;
            else
                return space::max(Radau::minStepCof, space::min(pow(0.55 / error, 1.0 / 7), Radau::maxStepCof));
        }

        /**
         *
         */
        inline void resetLastConvergence() {
            last_convergence_ = Radau::maxFloat;
        }

    private:
        Integrator integrator_;

        Scalar relativeError_{1e-10};

        Scalar convergent_limit_{1e-16};

        Scalar last_convergence_{Radau::maxFloat};

        size_t particleNum_{ParticSys::arraySize};

        constexpr static size_t max_iter_ = 12;
    };

}
#endif
