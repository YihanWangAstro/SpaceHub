
#ifndef INTERACTION_H
#define INTERACTION_H

#include "../protoType.h"
#include "../devTools.h"

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

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (TYPE, NAME, MEMBER). The first arg is the type of the array, the second
         *  arg is the name of the interface and the third arg is the private member that the interface connected. Each
         *  macros create five interfaces, they are :
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *  2. const typename TYPE::value_type & NAME (size_t i) const { return MEMBER[i];};
         *  3. void set_NAME (const TYPE &X) { MEMBER = X;};
         *  4. void set_NAME (size_t i, typename TYPE::value_type &X) { MEMBER[i] = X;};
         *  5. void swap_NAME (TYPE &X) { std::swap(X, MEMBER);};
         *
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_INTERFACES_FOR_ARRAY(VectorArray, acc, this_acc_);

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

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (TYPE, NAME, MEMBER). The first arg is the type of the array, the second
         *  arg is the name of the interface and the third arg is the private member that the interface connected. Each
         *  macros create five interfaces, they are :
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *  2. const typename TYPE::value_type & NAME (size_t i) const { return MEMBER[i];};
         *  3. void set_NAME (const TYPE &X) { MEMBER = X;};
         *  4. void set_NAME (size_t i, typename TYPE::value_type &X) { MEMBER[i] = X;};
         *  5. void swap_NAME (TYPE &X) { std::swap(X, MEMBER);};
         *
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_INTERFACES_FOR_ARRAY(VectorArray, totalAcc, acc_);

        /** @brief Interface adapter to inherit the interface of the data member
         *  The macros take four args (TYPE, MEMBER, NAME, NEWNAME). Each macros create five interfaces, they are:
         *
         *  1. const TYPE &NEWNAME () const { return MEMBER.NAME();};
         *  2. const typename TYPE::value_type & NEWNAME (size_t i) const { return MEMBER.NAME(i);};
         *  3. void set_NEWNAME (const TYPE &X) { MEMBER.set_NAME(X);};
         *  4. void set_NEWNAME (size_t i, typename TYPE::value_type &X) { MEMBER.set_NAME(i, X);};
         *  5. void swap_NEWNAME (TYPE &X) { MEMBER.swap_NAME(X)};
         *
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_INTERFACES_ADAPTER_FOR_ARRAY(VectorArray, vel_indep_,     acc, velIndepAcc);
        SPACEHUB_INTERFACES_ADAPTER_FOR_ARRAY(VectorArray, vel_dep_,       acc, velDepAcc);
        SPACEHUB_INTERFACES_ADAPTER_FOR_ARRAY(VectorArray, ext_vel_indep_, acc, extVelIndepAcc);
        SPACEHUB_INTERFACES_ADAPTER_FOR_ARRAY(VectorArray, ext_vel_dep_,   acc, extVelDepAcc);


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
