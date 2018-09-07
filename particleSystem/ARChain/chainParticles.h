#ifndef DYNAMICCHAIN_H
#define DYNAMICCHAIN_H

#include "../reguSystem/regularParticles.h"
#include "chain.h"

namespace SpaceH {
    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     * @tparam IsVelDep
     */
    template<typename Dtype, size_t ArraySize, bool IsVelDep>
    class VelIndepChainParticles : public ReguParticles<Dtype, ArraySize, IsVelDep> {
    public:
        /* Typedef */
        using Base = ReguParticles<Dtype, ArraySize, IsVelDep>;
        using typename Base::type;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using IndexArray   = typename type::IndexArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        constexpr static SpaceH::DATASTRUCT dataStruct{SpaceH::DATASTRUCT::CHAIN};

        using Base::particleNumber;
        using Base::totalMass_;
        using Base::mass_;
        using Base::pos_;
        using Base::vel_;

        /** @brief Resize all containers if they are dynamical
         *  @param New size of container.
         */
        void resize(size_t new_siz) {
            static_cast<Base &>(*this).resize(new_siz);
            cartesian_pos_.resize(new_siz);
            cartesian_vel_.resize(new_siz);
            ch_index_.resize(new_siz);
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            static_cast<Base &>(*this).reserve(new_cap);
            cartesian_pos_.reserve(new_cap);
            cartesian_vel_.reserve(new_cap);
            ch_index_.reserve(new_cap);
        }

        /**  @brief Position array const interface. Reference to pos_*/
        inline const VectorArray &chainPos() const { return pos_; }

        /**  @brief Velocity array const interface. Reference to vel_*/
        inline const VectorArray &chainVel() const { return vel_; }

        /**  @brief Index array const interface. Reference to ch_index_*/
        inline const IndexArray &chainIndex() const { return ch_index_; }

        /**  @brief Position vector const interface. Reference to pos_[i]*/
        inline const Vector &chainPos(size_t i) const { return pos_[i]; }

        /**  @brief Velocity vecotr const interface. Reference to vel_[i]*/
        inline const Vector &chainVel(size_t i) const { return vel_[i]; }

        /**  @brief Index const interface. Reference to ch_index_[i]*/
        inline size_t chainIndex(size_t i) const { return ch_index_[i]; }

        /**  @brief Position array const interface. Reference to state.pos*/
        inline const VectorArray &pos() const { return cartesian_pos_; }

        /**  @brief Velocity array const interface. Reference to state.vel*/
        inline const VectorArray &vel() const { return cartesian_vel_; }

        /**  @brief Position vector const interface. Reference to state.pos[i]*/
        inline const Vector &pos(size_t i) const { return cartesian_pos_[i]; }

        /**  @brief Velocity vecotr const interface. Reference to state.vel[i]*/
        inline const Vector &vel(size_t i) const { return cartesian_vel_[i]; }

        /** @brief Advance the position array with internal velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(Scalar stepSize) {
            SpaceH::advanceVector(pos_, vel_, stepSize);
            SpaceH::chain::synCartesian(pos_, cartesian_pos_, ch_index_);

            Vector CMPos = SpaceH::calcuCMCoord(mass_, cartesian_pos_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_pos_, CMPos);
        }

        /** @brief Advance the position array with given velocity array.
         *  @param stepSize The advance step size.
         *  @param vel      The velocity array.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            SpaceH::advanceVector(cartesian_pos_, vel, stepSize);
            Vector CMPos = SpaceH::calcuCMCoord(mass_, cartesian_pos_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_pos_, CMPos);

            SpaceH::chain::synChain(cartesian_pos_, pos_, ch_index_);
        }

        /** @brief Advance the  velocity array with given acceleration array.
         *  @param stepSize The advance step size.
         *  @param cartesian_acc The acceleration array in Cartesian coordinates.
         */
        inline void advanceVel(const VectorArray &cartesian_acc, Scalar stepSize) {
            const size_t chain_num = this->particleNumber() - 1;

            VectorArray chain_acc;

            chain_acc.resize(chain_num + 1);

            for (size_t i = 0; i < chain_num; ++i)
                chain_acc[i] = cartesian_acc[ch_index_[i + 1]] - cartesian_acc[ch_index_[i]];

            chain_acc[chain_num].setZero();

            SpaceH::advanceVector(vel_, chain_acc, stepSize);
            SpaceH::chain::synCartesian(vel_, cartesian_vel_, ch_index_);

            Vector CMVel = SpaceH::calcuCMCoord(mass_, cartesian_vel_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_vel_, CMVel);
        }

        /* Update the way to connect the chain_ if needed*/
        void updateChain() {
            IndexArray new_index = ch_index_;
            SpaceH::chain::getChainIndex(cartesian_pos_, new_index);

            if (SpaceH::chain::isDiff(new_index, ch_index_)) {
                SpaceH::chain::updateChain(pos_, ch_index_, new_index);
                SpaceH::chain::updateChain(vel_, ch_index_, new_index);
                ch_index_ = new_index;
            }
        }

        /** @brief Input(Initialize) variables with istream.*/
        friend std::istream &operator>>(std::istream &is, VelIndepChainParticles &partc) {
            is >> static_cast<Base &>(partc);

            /*  @note Due to the downcast here, this class has to maintain the size of the additional Arraies by itself.*/
            if constexpr (type::arraySize == SpaceH::DYNAMICAL) {
                size_t num = partc.idn().size();
                partc.cartesian_pos_.resize(num);
                partc.cartesian_vel_.resize(num);
                partc.ch_index_.resize(num);
            }

            SpaceH::chain::getChainIndex(partc.pos_, partc.ch_index_);
            SpaceH::chain::synChain(partc.pos_, partc.cartesian_pos_, partc.ch_index_);
            SpaceH::chain::synChain(partc.vel_, partc.cartesian_vel_, partc.ch_index_);

            std::swap(partc.pos_, partc.cartesian_pos_);
            std::swap(partc.vel_, partc.cartesian_vel_);

            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const NbodyIO IO_flag = NbodyIO::STD) {
            size_t loc = static_cast<Base &>(*this).read(data, IO_flag);

            SpaceH::chain::synCartesian(pos_, cartesian_pos_, ch_index_);
            SpaceH::chain::synCartesian(vel_, cartesian_vel_, ch_index_);
            Vector CMPos = SpaceH::calcuCMCoord(mass_, cartesian_pos_, totalMass_);
            Vector CMVel = SpaceH::calcuCMCoord(mass_, cartesian_vel_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_pos_, CMPos);
            SpaceH::moveToCMCoord(cartesian_vel_, CMVel);

            return loc;
        }


    protected:
        /** @brief Chained position array of the particles. Element is 3D vector.*/
        VectorArray cartesian_pos_;

        /** @brief Chained velocity array of the particles. Element is 3D vector.*/
        VectorArray cartesian_vel_;

        /** @brief Index array from Cartesian coordinates to chain coordinates Element is 3D vector.*/
        IndexArray ch_index_;
    };

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     * @tparam IsVelDep
     */
    template<typename Dtype, size_t ArraySize, bool IsVelDep>
    class VelDepChainParticles : public VelIndepChainParticles<Dtype, ArraySize, IsVelDep> {
    public:
        /* Typedef */
        using Base = VelIndepChainParticles<Dtype, ArraySize, IsVelDep>;
        using typename Base::type;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        using Base::particleNumber;
        using Base::totalMass_;
        using Base::mass_;
        using Base::vel_;
        using Base::cartesian_vel_;
        using Base::ch_index_;
        using Base::auxi_vel_;

        /** @brief Resize all containers if they are dynamical
         *  @param New size of container.
         */
        void resize(size_t new_siz) {
            static_cast<Base &>(*this).resize(new_siz);
            cartesian_auxi_vel_.resize(new_siz);
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            static_cast<Base &>(*this).reserve(new_cap);
            cartesian_auxi_vel_.reserve(new_cap);
        }

        /**  @brief Velocity array const interface. Reference to vel_*/
        inline const VectorArray &chainAuxiVel() const { return auxi_vel_; }

        /**  @brief Velocity vecotr const interface. Reference to vel_[i]*/
        inline const Vector &chainAuxiVel(size_t i) const { return auxi_vel_[i]; }

        /**  @brief Auxiliary velocity array const interface. Reference to auxi_vel_*/
        inline const VectorArray &auxiVel() const { return cartesian_auxi_vel_; }

        /**  @brief Auxiliary velocity vecotr const interface. Reference to auxi_vel_[i] */
        inline const Vector &auxiVel(size_t i) const { return cartesian_auxi_vel_[i]; }

        /** @brief Advance the  velocity array with given acceleration array.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceAuxiVel(const VectorArray &cartesian_acc, Scalar stepSize) {
            const size_t chain_num = this->particleNumber() - 1;

            VectorArray chain_acc;
            chain_acc.resize(chain_num + 1);


            for (size_t i = 0; i < chain_num; ++i)
                chain_acc[i] = cartesian_acc[ch_index_[i + 1]] - cartesian_acc[ch_index_[i]];

            chain_acc[chain_num].setZero();

            SpaceH::advanceVector(auxi_vel_, chain_acc, stepSize);
            SpaceH::chain::synCartesian(auxi_vel_, cartesian_auxi_vel_, ch_index_);

            //SpaceH::advanceVector(auxi_vel_, acc, stepSize);

            Vector CMVel = SpaceH::calcuCMCoord(mass_, cartesian_auxi_vel_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_auxi_vel_, CMVel);

        }

        /**  @brief synchronize auxiVel_ with vel_ */
        inline void synAuxiVelwithVel() {
            auxi_vel_ = vel_;
            cartesian_auxi_vel_ = cartesian_vel_;
        }

        /** @brief Input(Initialize) variables with istream.*/
        friend std::istream &operator>>(std::istream &is, VelDepChainParticles &partc) {
            is >> static_cast<Base &>(partc);

            partc.auxi_vel_ = partc.vel_;//assign here automatically resize the chain_auxi_vel_
            partc.cartesian_auxi_vel_ = partc.cartesian_vel_;
            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const NbodyIO IO_flag = NbodyIO::STD) {
            size_t loc = static_cast<Base &>(*this).read(data, IO_flag);

            SpaceH::chain::synCartesian(auxi_vel_, cartesian_auxi_vel_, ch_index_);
            Vector CMVel = SpaceH::calcuCMCoord(mass_, cartesian_auxi_vel_, totalMass_);
            SpaceH::moveToCMCoord(cartesian_auxi_vel_, CMVel);
            return loc;
        }

    private:
        /** @brief Chained velocity array of the particles. Element is 3D vector.*/
        VectorArray cartesian_auxi_vel_;
    };

    template<typename Dtype, size_t ArraySize, bool IsVelDep>
    class ChainParticles : public VelIndepChainParticles<Dtype, ArraySize, IsVelDep> {
    };

    template<typename Dtype, size_t ArraySize>
    class ChainParticles<Dtype, ArraySize, true> : public VelDepChainParticles<Dtype, ArraySize, true> {
    };
}
#endif
