
#ifndef PARTICLES_H
#define PARTICLES_H

#include "macros.h"
#include "type_class.h"
#include "dev_tools.h"
#include "core_computation.h"

namespace SpaceH {

    template <typename TypeClass>
    class SimpleState : public TypeClass::VectorArray {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM(TypeClass);
        /* Typedef */
        const VectorArray & cartesian() {
            return *this;
        }
        const VectorArray & raw() {
            return *this;
        }
        VectorArray &data() {
            return *this;
        }
    };
 /**
  * @brief Basic particles group with SoA memory laylout.
  * @tparam TypeClass Type class that define all the basic types.
  */
    template<typename TypeClass>
    class Particles {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM(TypeClass);
        using State = SimpleState<TypeClass>;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

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
        SPACEHUB_READ_INTERFACES_FOR_ARRAY(pos, VectorArray, pos_);

        SPACEHUB_READ_INTERFACES_FOR_ARRAY(vel, VectorArray, vel_);

        SPACEHUB_READ_INTERFACES_FOR_ARRAY(mass, ScalarArray, mass_);

        SPACEHUB_READ_INTERFACES_FOR_ARRAY(radius, ScalarArray, radius_);

        SPACEHUB_READ_INTERFACES_FOR_ARRAY(idn, IndexArray, idn_);

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
            SpaceH::advanceVector(pos_.data(), vel_.data(), stepSize);
        }

        /** @brief Advance the position array with given velocity array in cartesian coordinates.
         *  @param vel The given velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            SpaceH::advanceVector(pos_.data(), vel, stepSize);
        }

        /** @brief Advance the  velocity array with given acceleration array in cartesian coordinates.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceVel(const VectorArray &acc, Scalar stepSize) {
            SpaceH::advanceVector(vel_.data(), acc, stepSize);
        }

        /** @brief Resize all containers if they are dynamical
         *  @param new_siz New size of container.
         */
        void resize(size_t new_siz) {
            if constexpr (Types::array_size == SpaceH::DYNAMICAL) {
                pos_.resize(new_siz);
                vel_.resize(new_siz);
                mass_.resize(new_siz);
                radius_.resize(new_siz);
                idn_.resize(new_siz);
            } else {
                SPACEHUB_ERR_MSG("Fixed particles number! Cannot be resized!")
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
                radius_.reserve(new_cap);
                idn_.reserve(new_cap);
            } else {
                SPACEHUB_ERR_MSG("Fixed particles number! Cannot be reserved!")
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
                is >> pos_[loc], pos_[loc] *= fmt.L;
                is >> vel_[loc], vel_[loc] *= fmt.V;

                if (flag == IO_flag::STD) {
                    is >> mass_[loc], mass_[loc] *= fmt.M;
                    is >> radius_[loc], radius_[loc] *= fmt.L;
                    is >> idn_[loc];
                }
            }
            if (!is.good())
                SPACEHUB_ERR_MSG("Insufficent input data in initial file!");
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
                os << ' ' << pos_[i] / fmt.L << ' ' << vel_[i] / fmt.V;
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
         * @note If the array size is dynamical, set the particle number with `resize()` before
         * reading data from the buffer.
         */
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {

            size_t loc = 0;
            //for locality, split into separate loops
            for (auto &p : pos_) {
                p.x = data[loc++];
                p.y = data[loc++];
                p.z = data[loc++];
            }
            for (auto &v : vel_) {
                v.x = data[loc++];
                v.y = data[loc++];
                v.z = data[loc++];
            }
            if (flag == IO_flag::STD) {
                for (auto &m : mass_)
                    m = data[loc++];
                for (auto &r : radius_)
                    r = data[loc++];
                for (auto &i : idn_)
                    i = static_cast<size_t>(data[loc++]);
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
            data.reserve(particleNum * 6 + 1);

            //for locality, split into two loops
            for (const auto &p : pos_) {
                data.emplace_back(p.x);
                data.emplace_back(p.y);
                data.emplace_back(p.z);
            }
            for (const auto &v : vel_) {
                data.emplace_back(v.x);
                data.emplace_back(v.y);
                data.emplace_back(v.z);
            }

            if (flag == IO_flag::STD) {
                data.reserve(particleNum * 8 + 1);
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
        void moveToCoM() {
            SpaceH::moveToCoM(mass_, pos_);
            SpaceH::moveToCoM(mass_, vel_);
        }
        void moveToCoM(Scalar total_mass) {
            SpaceH::moveToCoM(mass_, pos_, total_mass);
            SpaceH::moveToCoM(mass_, vel_, total_mass);
        }
        void update() {}
    private:
        /** @brief Position array of the particles. Element is 3D vector.*/
        State pos_;

        /** @brief Velocity array of the particles. Element is 3D vector.*/
        State vel_;

        /** @brief Mass array of the particles. Element is Scalar.*/
        ScalarArray mass_;

        /** @brief Radius array of the particles. Element is Scalar.*/
        ScalarArray radius_;

        /** @brief Id Array of the particles. Element is int.*/
        IndexArray idn_;

        /** @brief The physical time of the dynamic system*/
        Scalar time_;
    };
}
#endif

