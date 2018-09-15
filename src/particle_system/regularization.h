
#ifndef REGULARIZATION_H
#define REGULARIZATION_H

#include "../core_computation.h"
#include "../dev_tools.h"

namespace SpaceH {

    /**
     * @brief Abstract class of regularization, Never use it directly.
     * @tparam TypeClass
     */
    template<typename TypeClass>
    class Regularization {
    public:
        /* Typedef */
        using type        = TypeClass;
        using Scalar      = typename type::Scalar;
        using VectorArray = typename type::VectorArray;
        using ScalarArray = typename type::ScalarArray;
        /* Typedef */

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (NAME, TYPE, MEMBER). This macro will create an interface:
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *
         *  See macros definition in 'devTools.h'.
         */
        SPACEHUB_READ_INTERFACES_FOR_SCALAR(omega, Scalar, omega_);
        SPACEHUB_READ_INTERFACES_FOR_SCALAR(bindE, Scalar, bindE_);
        SPACEHUB_WRITE_INTERFACES_FOR_SCALAR(omega, Scalar, omega_);
        SPACEHUB_WRITE_INTERFACES_FOR_SCALAR(bindE, Scalar, bindE_);

        /**
         *
         * @param partc
         * @return
         */
        template <typename Particles>
        inline Scalar capitalOmega(const Particles &partc) {
            return -getPotentialEnergy(partc.mass(), partc.pos());
        }

        /**
         *
         * @param partc
         */
        template <typename Particles>
        void init(const Particles &partc) {
            omega_ = capitalOmega(partc);
            bindE_ = -getTotalEnergy(partc.mass(), partc.pos(), partc.vel());
        }

        /** @brief Advance the Omega.
         *  @param velIndepAcc Velocity independent acceleration array.
         *  @param vel         Velocity array.
         *  @param stepSize    Time stepSize.
         */
        inline void
        advanceOmega(const VectorArray &domega, const VectorArray &vel, const ScalarArray &mass, Scalar stepSize) {
            size_t size = mass.size();
            Scalar sum = 0;

            for (size_t i = 0; i < size; ++i)
                sum += dot(domega[i], vel[i]) * mass[i];

            omega_ += sum * stepSize;
        }

        /** @brief Advance the bindE.
         *  @param velDepAcc Velocity dependent acceleration array.
         *  @param vel       Velocity array.
         *  @param stepSize  Time stepSize.
         */
        inline void
        advanceBindE(const VectorArray &dbindE, const VectorArray &vel, const ScalarArray &mass, Scalar stepSize) {
            size_t particleNum = mass.size();
            Scalar sum = 0;

            for (size_t i = 0; i < particleNum; ++i)
                sum -= dot(dbindE[i], vel[i]) * mass[i];

            bindE_ += sum * stepSize;
        }

    private:
        Scalar omega_;
        Scalar bindE_;
    };

/**  @brief logH extention algorithmatic regularization interface
  *
  *  See detials in https://link.springer.com/article/10.1023%2FA%3A1008368322547 and
  *                 http://iopscience.iop.org/article/10.1086/301102/meta .
  */
    template<typename TypeClass>
    class LogH : public Regularization<TypeClass> {
    public:
        /* Typedef */
        using Base   = Regularization<TypeClass>;
        using type   = typename Base::type;
        using Scalar = typename type::Scalar;
        /* Typedef */
        using Base::bindE;

        static constexpr bool nonTrivial{true};

        /** @brief Calculate the physical time for position advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size. This could not be the physical time.
         */
        template <typename Particles>
        inline Scalar getPhysicalPosTime(const Particles &partc, Scalar stepSize) {
            return stepSize / (bindE() + getKineticEnergy(partc.mass(), partc.vel()));
        }

        /** @brief Calculate the physical time for velocity advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size. This could not be the physical time.
         */
        template <typename Particles>
        inline Scalar getPhysicalVelTime(const Particles &partc, Scalar stepSize) {
            return stepSize / -getPotentialEnergy(partc.mass(), partc.pos());
            //return stepSize/-getPotentialEnergy(partc.mass(), partc.pos(), partc.chain_pos(), partc.chain_index());
        }

        /*private:
            template<typename Partic_>
            inline typename std::enable_if<Partic_::dataStruct == SpaceH::DATASTRUCT::PLAIN, Scalar>::type
            getPotEng(const Partic_ &partc) {
                return -getPotentialEnergy(partc.mass(), partc.pos());
            }

            template<typename Partic_>
            inline typename std::enable_if<Partic_::dataStruct == SpaceH::DATASTRUCT::CHAIN, Scalar>::type
            getPotEng(const Partic_ &partc) {
                return -getPotentialEnergy(partc.mass(), partc.pos(), partc.chainPos(), partc.chainIndex());
            }*/
    };

/**  @brief Time Transform Leapfrog algorithmatic regularization interface
 *
 *  See detials in https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
 */
    template<typename TypeClass>
    class TTL : public Regularization<TypeClass> {
    public:
        /* Typedef */
        using Base   = Regularization<TypeClass>;
        using type   = typename Base::type;
        using Scalar = typename type::Scalar;
        /* Typedef */

        using Base::omega;
        using Base::capitalOmega;

        static constexpr bool nonTrivial{true};

        /** @brief Calculate the physical time for position advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size. This could not be the physical time. Look
         *                  references for details in class despriction.
         */
        template <typename Particles>
        inline Scalar getPhysicalPosTime(const Particles &partc, Scalar stepSize) {
            return stepSize / omega();
        }

        /** @brief Calculate the physical time for velocity advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size. This could not be the physical time. Look
         *                  references for details in class despriction.
         */
        template <typename Particles>
        inline Scalar getPhysicalVelTime(const Particles &partc, Scalar stepSize) {
            return stepSize / capitalOmega();
        }
    };

/**  @brief Ordinary algorithmatic regularization interface
 *
 *   No regularization.
 */
    template<typename TypeClass>
    class NoRegu : public Regularization<TypeClass> {
    public:
        /* Typedef */
        using Base        = Regularization<TypeClass>;
        using type        = typename Base::type;
        using Scalar      = typename type::Scalar;
        using VectorArray = typename type::VectorArray;
        using ScalarArray = typename type::ScalarArray;
        /* Typedef */

        static constexpr bool nonTrivial{false};

        /** @brief Calculate the physical time for position advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size.
         */
        template <typename Particles>
        inline Scalar getPhysicalPosTime(const Particles &partc, Scalar stepSize) {
            return stepSize;
        }

        /** @brief Calculate the physical time for velocity advance from integration step size
         *  @param partc    Dynamic system contains position, velocity and regularization variables.
         *                  See example class in particles.h.
         *  @param stepSize Integration step size.
         */
        template <typename Particles>
        inline Scalar getPhysicalVelTime(const Particles &partc, Scalar stepSize) {
            return stepSize;
        }

        /** @brief Trivial interface for non-regularized regularizator. Will be deleted after inlining.*/
        inline void
        advanceOmega(const VectorArray &domega, const VectorArray &vel, const ScalarArray &mass, Scalar stepSize) {}

        /** @brief Trivial interface for non-regularized regularizator. Will be deleted after inlining.*/
        inline void
        advanceBindE(const VectorArray &dbindE, const VectorArray &vel, const ScalarArray &mass, Scalar stepSize) {}
    };
}
#endif

