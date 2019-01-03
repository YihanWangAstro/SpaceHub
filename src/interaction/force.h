
#ifndef FORCES_H
#define FORCES_H

#include "../type_class.h"

namespace SpaceH {
    /**
     *
     * @tparam Forcefunc
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Forcefunc, typename Dtype, size_t ArraySize>
    struct Force {
        /* Typedef */
        using type = SpaceH::TypeSystem<Dtype, ArraySize>;
        using Vector = typename type::Vector;
        using VectorArray = typename type::VectorArray;
        /* Typedef */

        /**
         *
         * @param acc
         */
        void addTotal(VectorArray &acc) {
            const size_t size = acc.size();

            for (size_t i = 0; i < size; ++i)
                acc[i] += this_acc_[i];
        }

        /**
         *
         * @return
         */
        const VectorArray &acc() const { return this_acc_; }

        /**
         *
         * @param i
         * @return
         */
        const Vector &acc(size_t i) const { return this_acc_[i]; }

        /**
         *
         * @param size
         */
        void checkArraySize(size_t size) {
            if (size != this_acc_.size())
                this_acc_.resize(size);
        }

    protected:
        VectorArray this_acc_;
        Forcefunc force_;
    };

    /**
     *
     * @tparam PairForce
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct VelIndepForce : public Force<PairForce, Dtype, ArraySize> {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        using type        = typename Base::type;
        using Scalar      = typename type::Scalar;
        using Vector      = typename type::Vector;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        using IndexArray  = typename type::IndexArray;
        /* Typedef */
        using Base::checkArraySize;

        /**
         *
         * @param mass
         * @param pos
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 1; j < size; ++j)
                    Base::force_(Base::this_acc_[i], Base::this_acc_[j], mass[i], mass[j], pos[i], pos[j],
                                 pos[j] - pos[i]);
        }

        /**
         *
         * @param mass
         * @param pos
         * @param chainPos
         * @param chainInd
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &chainPos,
                      const IndexArray &chainInd) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size - 1; ++i)
                Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[i + 1]],
                             mass[chainInd[i]], mass[chainInd[i + 1]],
                             pos[chainInd[i]], pos[chainInd[i + 1]], chainPos[i]);

            for (size_t i = 0; i < size - 2; ++i)
                Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[i + 2]],
                             mass[chainInd[i]], mass[chainInd[i + 2]],
                             pos[chainInd[i]], pos[chainInd[i + 2]], chainPos[i] + chainPos[i + 1]);

            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 3; j < size; ++j)
                    Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[j]],
                                 mass[chainInd[i]], mass[chainInd[j]],
                                 pos[chainInd[i]], pos[chainInd[j]],
                                 pos[chainInd[j]] - pos[chainInd[i]]);
        }

    };

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct VelIndepForce<void, Dtype, ArraySize> {
        /* Typedef */
        using type        = SpaceH::TypeSystem<Dtype, ArraySize>;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        using IndexArray  = typename type::IndexArray;

        /* Typedef */
        void addTotal(VectorArray &acc) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &chainPos,
                      const IndexArray &chainInd) {}
    };

    /**
     *
     * @tparam ExtForce
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename ExtForce, typename Dtype, size_t ArraySize>
    struct ExtVelIndepForce : public Force<ExtForce, Dtype, ArraySize> {
        /* Typedef */
        using Base = Force<ExtForce, Dtype, ArraySize>;
        using type        = typename Base::type;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        /* Typedef */
        using Base::checkArraySize;

        /**
         *
         * @param mass
         * @param pos
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size; ++i)
                Base::force_(Base::this_acc_[i], mass[i], pos[i]);
        }

    };

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct ExtVelIndepForce<void, Dtype, ArraySize> {
        /* Typedef */
        using type        = SpaceH::TypeSystem<Dtype, ArraySize>;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;

        /* Typedef */
        void addTotal(VectorArray &acc) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos) {}
    };

    /**
     *
     * @tparam PairForce
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct VelDepForce : public Force<PairForce, Dtype, ArraySize> {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        using type        = typename Base::type;
        using Scalar      = typename type::Scalar;
        using Vector      = typename type::Vector;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        using IndexArray  = typename type::IndexArray;
        /* Typedef */
        using Base::checkArraySize;

        /**
         *
         * @param mass
         * @param pos
         * @param vel
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 1; j < size; ++j)
                    Base::force_(Base::this_acc_[i], Base::this_acc_[j], mass[i], mass[j],
                                 pos[i], pos[j], vel[i], vel[j],
                                 pos[j] - pos[i], vel[j] - vel[i]);
        }

        /**
         *
         * @param mass
         * @param pos
         * @param vel
         * @param chainPos
         * @param chainVel
         * @param chainInd
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel,
                      const VectorArray &chainPos, const VectorArray &chainVel, const IndexArray &chainInd) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size - 1; ++i)
                Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[i + 1]],
                             mass[chainInd[i]], mass[chainInd[i + 1]],
                             pos[chainInd[i]], pos[chainInd[i + 1]],
                             vel[chainInd[i]], vel[chainInd[i + 1]],
                             chainPos[i], chainVel[i]);


            for (size_t i = 0; i < size - 2; ++i)
                Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[i + 2]],
                             mass[chainInd[i]], mass[chainInd[i + 2]],
                             pos[chainInd[i]], pos[chainInd[i + 2]],
                             vel[chainInd[i]], vel[chainInd[i + 2]],
                             chainPos[i] + chainPos[i + 1], chainVel[i] + chainVel[i + 1]);

            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 3; j < size; ++j)
                    Base::force_(Base::this_acc_[chainInd[i]], Base::this_acc_[chainInd[j]],
                                 mass[chainInd[i]], mass[chainInd[j]],
                                 pos[chainInd[i]], pos[chainInd[j]],
                                 vel[chainInd[i]], vel[chainInd[j]],
                                 pos[chainInd[j]] - pos[chainInd[i]],
                                 vel[chainInd[j]] - vel[chainInd[i]]);
        }
    };

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct VelDepForce<void, Dtype, ArraySize> {
        /* Typedef */
        using type        = SpaceH::TypeSystem<Dtype, ArraySize>;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        using IndexArray  = typename type::IndexArray;

        /* Typedef */

        void addTotal(VectorArray &acc) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel,
                      const VectorArray &chainPos, const VectorArray &chainVel, const IndexArray &chainInd) {}
    };

    /**
     *
     * @tparam ExtForce
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename ExtForce, typename Dtype, size_t ArraySize>
    struct ExtVelDepForce : public Force<ExtForce, Dtype, ArraySize> {
        /* Typedef */
        using Base = Force<ExtForce, Dtype, ArraySize>;
        using type        = typename Base::type;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;
        /* Typedef */
        using Base::checkArraySize;

        /**
         *
         * @param mass
         * @param pos
         * @param vel
         */
        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel) {
            const size_t size = mass.size();
            checkArraySize(size);

            for (size_t i = 0; i < size; ++i)
                Base::this_acc_[i].setZero();

            for (size_t i = 0; i < size; ++i)
                Base::force_(Base::this_acc_[i], mass[i], pos[i], vel[i]);
        }

    };

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct ExtVelDepForce<void, Dtype, ArraySize> {
        /* Typedef */
        using type        = SpaceH::TypeSystem<Dtype, ArraySize>;
        using ScalarArray = typename type::ScalarArray;
        using VectorArray = typename type::VectorArray;

        /* Typedef */
        void addTotal(VectorArray &acc) {}

        void calcuAcc(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel) {}
    };

}

#endif
