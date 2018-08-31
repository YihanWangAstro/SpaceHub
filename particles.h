
#ifndef PARTICLES_H
#define PARTICLES_H

#include "macros.h"
#include "protoType.h"
#include "devTools.h"
#include "coreComputation.h"

namespace SpaceH {

/**
 *  @brief Basic particles group with SoA memory laylout.
 *  @tparam Dtype Type of scalar. e.g., float, double, kahanNumber...
 *  @tparam ArraySize The size of the arrays in whole system. SpaceH::DYNAMICAL for dynamical array.
 */
    template<typename Dtype, size_t ArraySize>
    class Particles {
    public:
        /* Typedef */
        using type         = SpaceH::ProtoType<Dtype, ArraySize>;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarArray  = typename type::ScalarArray;
        using IntArray     = typename type::IntArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        constexpr static size_t arraySize{type::arraySize};

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
        SPACEHUB_INTERFACES_FOR_ARRAY(VectorArray, pos,    pos_);
        SPACEHUB_INTERFACES_FOR_ARRAY(VectorArray, vel,    vel_);
        SPACEHUB_INTERFACES_FOR_ARRAY(ScalarArray, mass,   mass_);
        SPACEHUB_INTERFACES_FOR_ARRAY(ScalarArray, radius, radius_);
        SPACEHUB_INTERFACES_FOR_ARRAY(IntArray,    idn,    idn_);

        /** @brief Advance the position array with internal velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(Scalar stepSize) {
            SpaceH::advanceVector(pos_, vel_, stepSize);
        }

        /** @brief Advance the position array with given velocity array.
         *  @param vel The given velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            SpaceH::advanceVector(pos_, vel, stepSize);
        }

        /** @brief Advance the  velocity array with given acceleration array.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceVel(const VectorArray &acc, Scalar stepSize) {
            SpaceH::advanceVector(vel_, acc, stepSize);
        }

        /** @brief Resize all containers if they are dynamical
         *  @param new_siz New size of container.
         */
        void resize(size_t new_siz) {
            pos_.resize(new_siz);
            vel_.resize(new_siz);
            mass_.resize(new_siz);
            radius_.resize(new_siz);
            idn_.resize(new_siz);
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            pos_.reserve(new_cap);
            vel_.reserve(new_cap);
            mass_.reserve(new_cap);
            radius_.reserve(new_cap);
            idn_.reserve(new_cap);
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
                SpaceH::errMsg("Insufficent input data in initial file!", __FILE__, __LINE__);
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
        }

        /**
         * @brief Read variables from plain scalar buffer.
         * @param data
         * @param flag
         * @return
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
                    i = data[loc++];
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
                    data.emplace_back(i);
            }

            data.emplace_back(time_);
            return data.size();
        }

    private:
        /** @brief Position array of the particles. Element is 3D vector.*/
        VectorArray pos_;

        /** @brief Velocity array of the particles. Element is 3D vector.*/
        VectorArray vel_;

        /** @brief Mass array of the particles. Element is Scalar.*/
        ScalarArray mass_;

        /** @brief Radius array of the particles. Element is Scalar.*/
        ScalarArray radius_;

        /** @brief Id Array of the particles. Element is int.*/
        IntArray idn_;

        /** @brief The physical time of the dynamic system*/
        Scalar time_;
    };
}
#endif

