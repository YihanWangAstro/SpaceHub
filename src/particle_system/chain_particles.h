
#ifndef CHAINPARTICLES_H
#define CHAINPARTICLES_H

#include "../macros.h"
#include "../type_class.h"
#include "../dev_tools.h"
#include "../core_computation.h"
#include "chain.h"

namespace SpaceH {

    template<typename TypeClass>//not good design, try to optimize it
    class ChainCoord {
    public:
        template <typename T>
        friend class ChainParticles;
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        /* Typedef */
        void resize(size_t new_siz) {
            cartesian_.resize(new_siz);
            chain_.resize(new_siz);
        }

        void reserve(size_t new_cap) {
            cartesian_.reserve(new_cap);
            chain_.resize(new_cap);
        }

        const VectorArray &cartesian() {
            return cartesian_;
        }

        const VectorArray &raw() {
            return chain_;
        }

    private:
        void advanceCartesian(const VectorArray &cartesian_inc, Scalar stepSize) {
            const size_t size = cartesian_inc.size();
            VectorArray chain_inc;
            if constexpr (Types::array_size == SpaceH::DYNAMICAL)
                chain_inc.resize(size);
            SpaceH::chain::synChain(cartesian_inc, chain_inc, *index_);
            advanceChain(chain_inc, stepSize);
        }

        void advanceChain(const VectorArray &increament, Scalar stepSize) {
            SpaceH::advanceVector(chain_, increament, stepSize);
            synCartesian();
        }

        void synChain() {
            SpaceH::chain::synChain(cartesian_, chain_, *index_);
        }

        void synCartesian() {
            SpaceH::chain::synCartesian(chain_, cartesian_, *index_);
            //SpaceH::moveToCoM(*mass_, cartesian_);
        }

        void createChainIndex(IndexArray& index) {
            SpaceH::chain::getChainIndex(cartesian_, index);
        }

        void set_mass_index(const ScalarArray* mass, const IndexArray *index) {
            mass_  = mass;
            index_ = index;
        }
        VectorArray cartesian_;
        VectorArray chain_;
        const ScalarArray *mass_{nullptr};
        const IndexArray  *index_{nullptr};
    };

    /**
     * @brief Basic particles group with SoA memory laylout.
     * @tparam TypeClass Type class that define all the basic types.
     */
    template<typename TypeClass>
    class ChainParticles {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        using State        = ChainCoord<TypeClass>;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        ChainParticles() {
            pos_.set_mass_index(&mass_, &chain_index_);
            vel_.set_mass_index(&mass_, &chain_index_);
        }

        /** @brief Get the number of the particles.
         *  @return The particle number.
         */
        constexpr inline size_t particleNumber() const {
            return idn_.size();
        }

        /**  @brief Physical time scalar const interface. Reference to time_*/
        inline const Scalar &time() const {
            return time_;
        }

        /** @brief Advance the time.
         *  @param dt Time increament.
         */
        inline void advanceTime(Scalar dt) {
            time_ += dt;
        }

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (NAME, TYPE, MEMBER). Each macros create two interfaces, they are :
         *
         *  1. const TYPE &NAME () const { return MEMBER;};
         *  2. const typename TYPE::value_type & NAME (size_t i) const { return MEMBER[i];};
         *  See macros definition in 'devTools.h'.
         */
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(pos, VectorArray, pos_.cartesian_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(vel, VectorArray, vel_.cartesian_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(chain_pos, VectorArray, pos_.chain_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(chain_vel, VectorArray, vel_.chain_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(chain_index, IndexArray, chain_index_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(mass, ScalarArray, mass_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(radius, ScalarArray, radius_);
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(idn, IndexArray, idn_);

        inline const Vector& chain_ref_pos(size_t i) const {return pos_.cartesian_[chain_index_[i]];};
        inline const Vector& chain_ref_vel(size_t i) const {return vel_.cartesian_[chain_index_[i]];};
        inline Scalar chain_ref_mass(size_t i) const {return mass_[chain_index_[i]];};
        inline Scalar chain_ref_radius(size_t i) const {return radius_[chain_index_[i]];};
        inline size_t chain_ref_idn(size_t i) const {return idn_[chain_index_[i]];};

        /** Automaticlly create interfaces for data
        *  The macros takes three parameters (NAME, TYPE, MEMBER). This macros create one read interface :
        *  1. const TYPE &NAME () const { return MEMBER;};
        *  See macros definition in 'devTools.h'.
        */
        SPACEHUB_READ_INTERFACES_FOR_SCALAR(pos_state, State, pos_);
        SPACEHUB_READ_INTERFACES_FOR_SCALAR(vel_state, State, vel_);

        /** Automaticlly create interfaces for data
         *  The macros takes three parameters (NAME, TYPE, MEMBER). This macros create two interfaces :
         *
         *  1. void set_NAME (const TYPE& scalar) { MEMBER = scalar;};
         *  2. void swap_NAME (TYPE& scalar) { std::swap(MEMBER, scalar); }
         *  See macros definition in 'devTools.h'.
         */
        SPACEHUB_WRITE_INTERFACES_FOR_SCALAR(pos_state, State, pos_);
        SPACEHUB_WRITE_INTERFACES_FOR_SCALAR(vel_state, State, vel_);

        /** @brief Advance the position array with internal velocity array in real evolving coordinates.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(Scalar stepSize) {
            pos_.advanceChain(vel_.chain_, stepSize);
        }

        /** @brief Advance the position array with given velocity array in cartesian coordinates.
         *  @param vel The given velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            pos_.advanceCartesian(vel, stepSize);
        }

        /** @brief Advance the  velocity array with given acceleration array in cartesian coordinates.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceVel(const VectorArray &acc, Scalar stepSize) {
            vel_.advanceCartesian(acc, stepSize);
        }

        /** @brief Resize all containers if they are dynamical
         *  @param new_siz New size of container.
         */
        void resize(size_t new_siz) {
            if constexpr (Types::array_size == SpaceH::DYNAMICAL) {
                pos_.resize(new_siz);
                vel_.resize(new_siz);
                mass_.resize(new_siz);
                chain_index_.resize(new_siz);
                radius_.resize(new_siz);
                idn_.resize(new_siz);
            } else {
                SPACEHUB_ABORT("Fixed particles number! Cannot be resized!")
            }
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            if constexpr (Types::array_size == SpaceH::DYNAMICAL) {
                pos_.reserve(new_cap);
                vel_.reserve(new_cap);
                mass_.reserve(new_cap);
                chain_index_.reserve(new_cap);
                radius_.reserve(new_cap);
                idn_.reserve(new_cap);
            } else {
                SPACEHUB_ABORT("Fixed particles number! Cannot be reserved!")
            }
        }

        void moveToCoM() {
            SpaceH::moveToCoM(mass_, pos_.cartesian_);
            SpaceH::moveToCoM(mass_, vel_.cartesian_);
        }

        void moveToCoM(Scalar total_mass) {
            SpaceH::moveToCoM(mass_, pos_.cartesian_, total_mass);
            SpaceH::moveToCoM(mass_, vel_.cartesian_, total_mass);
        }

        void update() {
            IndexArray new_index = chain_index_;
            SpaceH::chain::getChainIndex(pos_.cartesian_, new_index);

            if (SpaceH::chain::isDiff(new_index, chain_index_)) {
                SpaceH::chain::updateChain(pos_.chain_, chain_index_, new_index);
                SpaceH::chain::updateChain(vel_.chain_, chain_index_, new_index);
                chain_index_ = new_index;
                pos_.synCartesian();
                vel_.synCartesian();
            }
        }

        /**
         * @brief Input variables from an IO stream
         * @param is The input IO stream.
         * @param fmt The units format of the input data, default to be all unity.
         * @param flag The IO flag for data. STD for all data, EVOLVED for evloved only data.
         * @return The partile number
         */
        size_t read(std::istream &is, const Unit::Fmt &fmt = Unit::UNITY_UNIT, IO_flag flag = IO_flag::STD) {
            is >> time_, time_ *= fmt.T;

            size_t particleNum = particleNumber();
            for (size_t loc = 0; loc < particleNum; loc++) {
                is >> pos_.cartesian_[loc], pos_.cartesian_[loc] *= fmt.L;
                is >> vel_.cartesian_[loc], vel_.cartesian_[loc] *= fmt.V;

                if (flag == IO_flag::STD) {
                    is >> mass_[loc], mass_[loc] *= fmt.M;
                    is >> radius_[loc], radius_[loc] *= fmt.L;
                    is >> idn_[loc];
                }
            }
            if (!is.good()) SPACEHUB_ABORT("Insufficent input data in initial file!");

            if (!chained) {
                moveToCoM();
                pos_.createChainIndex(chain_index_);
                chained = true;
            }
            pos_.synChain();
            vel_.synChain();

            return particleNum;
        }

        /**
         * @brief Write variables to an IO stream.
         * @param os The output IO stream.
         * @param fmt The units format of the input data, default to be all unity.
         * @param flag The IO flag for data. STD for all data, EVOLVED for evloved only data.
         * @return The partile number
         */
        size_t write(std::ostream &os, const Unit::Fmt &fmt = Unit::UNITY_UNIT, IO_flag flag = IO_flag::STD) const {
            os << ' ' << time_ / fmt.T << "\r\n";

            size_t particleNum = particleNumber();
            for (size_t i = 0; i < particleNum; ++i) {
                os << ' ' << pos_.cartesian_[i] / fmt.L
                   << ' ' << vel_.cartesian_[i] / fmt.V;
                if (flag == IO_flag::STD) {
                    os << ' ' << mass_[i] / fmt.M
                       << ' ' << radius_[i] / fmt.L
                       << ' ' << idn_[i];
                }
                os << "\r\n";
            }
            return particleNum;
        }

        /**
         * @brief Read variables from plain scalar buffer.
         * @param data
         * @param flag
         * @return
         */
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {
            size_t loc = 0;
            if (flag == IO_flag::EVOLVED) {
                if (!chained) SPACEHUB_ABORT(
                        "Chain Index hasn't been constructed by the positions in the Cartesian coordinates."
                        "Cannot read from raw chain data directly");
                //for locality, split into separate loops
                for (auto &p : pos_.chain_) {
                    p.x = data[loc++];
                    p.y = data[loc++];
                    p.z = data[loc++];
                }
                pos_.synCartesian();
                for (auto &v : vel_.chain_) {
                    v.x = data[loc++];
                    v.y = data[loc++];
                    v.z = data[loc++];
                }
                vel_.synCartesian();
            } else if (flag == IO_flag::STD) {
                for (auto &p : pos_.cartesian_) {
                    p.x = data[loc++];
                    p.y = data[loc++];
                    p.z = data[loc++];
                }
                if (!chained) {
                    pos_.createChainIndex(chain_index_);
                    chained = true;
                }
                for (auto &v : vel_.cartesian_) {
                    v.x = data[loc++];
                    v.y = data[loc++];
                    v.z = data[loc++];
                }
                for (auto &m : mass_)
                    m = data[loc++];
                for (auto &r : radius_)
                    r = data[loc++];
                for (auto &i : idn_)
                    i = static_cast<size_t >(data[loc++]);

                pos_.synChain();
                vel_.synChain();
            }
            time_ = data[loc++];
            return loc;
        }

        /**
         * @brief Write variables to plain scalar buffer.
         * @param data
         * @param flag
         * @return
         */
        size_t write(ScalarBuffer &data, const IO_flag flag = IO_flag::STD) const {
            size_t particleNum = particleNumber();

            data.clear();

            if (flag == IO_flag::EVOLVED) {
                data.reserve(particleNum * 6 + 1);
                //for locality, split into two loops
                for (const auto &p : pos_.chain_) {
                    data.emplace_back(p.x);
                    data.emplace_back(p.y);
                    data.emplace_back(p.z);
                }
                for (const auto &v : vel_.chain_) {
                    data.emplace_back(v.x);
                    data.emplace_back(v.y);
                    data.emplace_back(v.z);
                }
            } else if (flag == IO_flag::STD) {
                data.reserve(particleNum * 8 + 1);
                for (const auto &p : pos_.cartesian_) {
                    data.emplace_back(p.x);
                    data.emplace_back(p.y);
                    data.emplace_back(p.z);
                }
                for (const auto &v : vel_.cartesian_) {
                    data.emplace_back(v.x);
                    data.emplace_back(v.y);
                    data.emplace_back(v.z);
                }
                for (const auto &m : mass_)
                    data.emplace_back(m);
                for (const auto &r : radius_)
                    data.emplace_back(r);
                for (const auto &i : idn_)
                    data.emplace_back(static_cast<Scalar>(i));
            }

            data.emplace_back(time_);
            return data.size();
        }

    private:
        /** @brief Position array of the particles. Element is 3D vector.*/
        State pos_;

        /** @brief Velocity array of the particles. Element is 3D vector.*/
        State vel_;

        /** @brief Mass array of the particles. Element is Scalar.*/
        ScalarArray mass_;

        /** @brief Chain index array of the chain structure. Element is size_t.*/
        IndexArray chain_index_;

        /** @brief Radius array of the particles. Element is Scalar.*/
        ScalarArray radius_;

        /** @brief Id Array of the particles. Element is int.*/
        IndexArray idn_;

        /** @brief The physical time of the dynamic system*/
        Scalar time_;

        bool chained{false};
    };
}
#endif

