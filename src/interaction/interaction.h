
#ifndef INTERACTION_H
#define INTERACTION_H

#include "../type_class.h"
#include "../dev_tools.h"

namespace SpaceH {

    /**
     * Acceleration evaluator. This class is created on purpose that each new force manages its own acceleration array.
     * @tparam Forcefunc Callable object to provide the method of acceleration evaluation.
     * @tparam Dtype     Scalar type.
     * @tparam ArraySize Array size.
     */
    template<typename Forcefunc, typename TypeClass>
    struct AccEvaluator {
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        /* Typedef */

        void addAccTo(VectorArray &acc) {
            const size_t size = acc.size();

            for (size_t i = 0; i < size; ++i)
                acc[i] += this_acc_[i];
        }

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (NAME, TYPE, MEMBER). Each macros create two interfaces, they are :
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *  2. const typename TYPE::value_type & NAME (size_t i) const { return MEMBER[i];};
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(acc, VectorArray, this_acc_);

        void resize(size_t size) {
            if constexpr (Types::array_size == SpaceH::DYNAMICAL) {
                this_acc_.resize(size);
            } else {
                SPACEHUB_ABORT("Fixed particles number! Cannot be resized!");
            }
        }

        void reserve(size_t size) {
            if constexpr (Types::array_size == SpaceH::DYNAMICAL) {
                this_acc_.reserve(size);
            } else {
                SPACEHUB_ABORT("Fixed particles number! Cannot be reserved!");
            }
        }

        template<typename Particles>
        inline void evaluateAcc(const Particles &partc) {
            acc_callback_(partc, this_acc_);
        }

    private:
        VectorArray this_acc_;
        Forcefunc acc_callback_;
    };

    template<typename TypeClass>
    struct AccEvaluator<void, TypeClass> {
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        /* Typedef */

       /* const VectorArray& acc() const {
            //ERR_MSG("You are trying to access the acceleration of 'void' type force!");
        }*/

        void addAccTo(VectorArray &acc) {}

        void resize(size_t size) {}

        void reserve(size_t size) {}

        template<typename Particles>
        inline void evaluateAcc(const Particles &partc) {}
    };

    template<typename PairVelIndep, typename PairVelDep = void, typename ExtVelIndep = void, typename ExtVelDep = void>
    class Interactions {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(PairVelIndep);
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/
        constexpr static bool isVelDep{(!std::is_void<PairVelDep>::value) | (!std::is_void<ExtVelDep>::value)};

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (NAME, TYPE, MEMBER). Each macros create two interfaces, they are :
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *  2. const typename TYPE::value_type & NAME (size_t i) const { return MEMBER[i];};
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(acc, VectorArray, acc_);

        /** @brief Interface adapter to inherit the interface of the data member
         *  The macros take four args (NEWNAME, TYPE, MEMBER, NAME). Each macros create two interfaces, they are:
         *
         *  1. const TYPE &NEWNAME () const { return MEMBER.NAME();};
         *  2. const typename TYPE::value_type & NEWNAME (size_t i) const { return MEMBER.NAME(i);};
         *
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(pairVelIndepAcc, VectorArray, pair_vel_indep_, acc);
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(pairVelDepAcc,   VectorArray, pair_vel_dep_,   acc);
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(extVelIndepAcc,  VectorArray, ext_vel_indep_,  acc);
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(extVelDepAcc,    VectorArray, ext_vel_dep_,    acc);

        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuPairVelIndepAcc(const Particles &partc) {
            if constexpr (!std::is_void<PairVelIndep>::value) {
                pair_vel_indep_.evaluateAcc(partc);
            }
        }
        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuPairVelDepAcc(const Particles &partc) {
            if constexpr (!std::is_void<PairVelDep>::value) {
                pair_vel_dep_.evaluateAcc(partc);
            }
        }
        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuExtVelIndepAcc(const Particles &partc) {
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.evaluateAcc(partc);
            }
        }
        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuExtVelDepAcc(const Particles &partc) {
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.evaluateAcc(partc);
            }
        }
        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuVelIndepAcc(const Particles &partc) {
            calcuPairVelIndepAcc(partc);
            calcuExtVelIndepAcc(partc);
        }
        /**
         *
         * @tparam Particles
         * @param partc
         */
        template<typename Particles>
        inline void calcuVelDepAcc(const Particles &partc) {
            calcuPairVelDepAcc(partc);
            calcuExtVelDepAcc(partc);
        }

        inline void zeroTotalAcc() {
            for (auto &a : acc_)
                a.setZero();
        }

        void sumTotalAcc() {
            if constexpr  (!std::is_void<PairVelIndep>::value) {
                acc_ = pair_vel_indep_.acc();
            } else {
                zeroTotalAcc();
            }
            if constexpr (!std::is_void<PairVelDep>::value) {
                pair_vel_dep_.addAccTo(acc_);
            }
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.addAccTo(acc_);
            }
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.addAccTo(acc_);
            }
        }

        void resize(size_t new_siz) {
            if constexpr (!std::is_void<PairVelIndep>::value) {
                pair_vel_indep_.resize(new_siz);
            }
            if constexpr (!std::is_void<PairVelDep>::value) {
                pair_vel_dep_.resize(new_siz);
            }
            if constexpr (!std::is_void<ExtVelIndep>::value) {
                ext_vel_indep_.resize(new_siz);
            }
            if constexpr (!std::is_void<ExtVelDep>::value) {
                ext_vel_dep_.resize(new_siz);
            }
        }

        void reserve(size_t new_siz) {
            if constexpr (!std::is_void<PairVelIndep>::value) {
                pair_vel_indep_.reserve(new_siz);
            }
            if constexpr (!std::is_void<PairVelDep>::value) {
                pair_vel_dep_.reserve(new_siz);
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

        AccEvaluator<PairVelIndep, Types> pair_vel_indep_;

        AccEvaluator<PairVelDep,   Types> pair_vel_dep_;

        AccEvaluator<ExtVelIndep,  Types> ext_vel_indep_;

        AccEvaluator<ExtVelDep,    Types> ext_vel_dep_;
    };
}

#endif
