
#ifndef INTERACTION_H
#define INTERACTION_H

#include "../protoType.h"

namespace SpaceH {

    /**
     * Acceleration evaluator. This class is created on purpose that each new force manages its own acceleration array.
     * @tparam Forcefunc Callable object to provide the method of acceleration evaluation.
     * @tparam Dtype     Scalar type.
     * @tparam ArraySize Array size.
     */
    template<typename Forcefunc, typename Dtype, size_t ArraySize>
    struct AccEvaluator {
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        using Vector = typename type::Vector;
        using VectorArray = typename type::VectorArray;
        /* Typedef */

        void addAccTo(VectorArray &acc) {
            const size_t size = acc.size();

            for (size_t i = 0; i < size; ++i)
                acc[i] += this_acc_[i];
        }

        const VectorArray &acc() const { return this_acc_; }

        const Vector &acc(size_t i) const { return this_acc_[i]; }

        void resize(size_t size) {
            this_acc_.resize(size);
        }

        void reserve(size_t size) {
            this_acc_.reserve(size);
        }

        template<typename ParticleSystem>
        inline void evaluateAcc(const ParticleSystem &partc) {
            acc_evaluator_(partc, this_acc_);
        }

    private:
        VectorArray this_acc_;
        Forcefunc acc_evaluator_;
    };

    template<typename VelIndep, typename VelDep = void, typename ExtVelIndep = void, typename ExtVelDep = void>
    class Interactions {
    public:
        /* Typedef */
        using type        = typename VelIndep::type;
        using Scalar      = typename type::Scalar;
        using Vector      = typename type::Vector;
        using VectorArray = typename type::VectorArray;
        using ScalarArray = typename type::ScalarArray;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/
        constexpr static bool isVelDep{!std::is_void<VelDep>::value | !std::is_void<ExtVelDep>::value};

        inline const VectorArray &totalAcc() const { return acc_; }

        inline const VectorArray &velIndepAcc() const { return vel_indep_.acc(); }

        inline const VectorArray &velDepAcc() const { return vel_dep_.acc(); }

        inline const VectorArray &extVelIndepAcc() const { return ext_vel_indep_.acc(); }

        inline const VectorArray &extVelDepAcc() const { return ext_vel_dep_.acc(); }

        inline const Vector &totalAcc(size_t i) const { return acc_[i]; }

        inline const Vector &velIndepAcc(size_t i) const { return vel_indep_.acc(i); }

        inline const Vector &velDepAcc(size_t i) const { return vel_dep_.acc(i); }

        inline const Vector &extVelIndepAcc(size_t i) const { return ext_vel_indep_.acc(i); }

        inline const Vector &extVelDepAcc(size_t i) const { return ext_vel_dep_.acc(i); }

        inline void calcuVelIndepAcc(const Particles &partc) {
            vel_indep_.evaluateAcc(partc);
        }

        inline void calcuVelDepAcc(const Particles &partc) {
            if constexpr (!std::is_void<VelDep>::value) {
                vel_dep_.evaluateAcc(partc);
            }
        }

        inline void calcuExtVelIndepAcc(const Particles &partc) {
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.evaluateAcc(partc);
            }
        }

        inline void calcuExtVelDepAcc(const Particles &partc) {
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.evaluateAcc(partc);
            }
        }

        void zeroTotalAcc() {
            for (auto &a : acc_)
                a.setZero();
        }

        void sumTotalAcc() {
            acc_ = vel_indep_.acc();

            if constexpr (!std::is_void<VelDep>::value) {
                vel_dep_.addAccTo(acc_);
            }
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.addAccTo(acc_);
            }
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.addAccTo(acc_);
            }
        }

        void resize(size_t new_siz) {
            vel_indep_.resize(new_siz);

            if constexpr (!std::is_void<VelDep>::value) {
                vel_dep_.resize(new_siz);
            }
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.resize(new_siz);
            }
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.resize(new_siz);
            }
        }

        void reserve(size_t new_siz) {
            vel_indep_.reserve(new_siz);

            if constexpr (!std::is_void<VelDep>::value) {
                vel_dep_.reserve(new_siz);
            }
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.reserve(new_siz);
            }
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.reserve(new_siz);
            }
        }
    private:
        VectorArray acc_;

        AccEvaluator<VelIndep, Scalar, type::arraySize> vel_indep_;

        AccEvaluator<VelDep, Scalar, type::arraySize> vel_dep_;

        AccEvaluator<ExtVelIndep, Scalar, type::arraySize> ext_vel_indep_;

        AccEvaluator<ExtVelDep, Scalar, type::arraySize> ext_vel_dep_;
    };
}

#endif
