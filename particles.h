
#ifndef PARTICLES_H
#define PARTICLES_H

#include "macros.h"
#include "protoType.h"
#include "devTools.h"

namespace SpaceH {

/**
 *  @brief Basic particles group, wrapper on VelIndepParticles and VelDepParticles.
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
            return idn.size();
        }

        /** @brief Resize all containers if they are dynamical
         *  @param new_siz New size of container.
         */
        void resize(size_t new_siz) {
            pos.resize(new_siz);
            vel.resize(new_siz);
            mass.resize(new_siz);
            radius.resize(new_siz);
            idn.resize(new_siz);
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            pos.reserve(new_cap);
            vel.reserve(new_cap);
            mass.reserve(new_cap);
            radius.reserve(new_cap);
            idn.reserve(new_cap);
        }

        /**
         * @brief Input variables from a IO stream
         * @param is
         * @param fmt
         * @param flag
         * @return
         */
        size_t read(std::istream &is, const Unit::Fmt &fmt = Unit::UNITY_UNIT, IO_flag flag = IO_flag::STD) {
            is >> time, time *= fmt.T;

            size_t particleNum = particleNumber();
            for (size_t loc = 0; loc < particleNum ; loc++) {
                is >> pos[loc], pos[loc] *= fmt.L;
                is >> vel[loc], vel[loc] *= fmt.V;

                if (flag == IO_flag::STD) {
                    is >> mass[loc], mass[loc] *= fmt.M;
                    is >> radius[loc], radius[loc] *= fmt.L;
                    is >> idn[loc];
                }
            }
            if (!is.good())
                SpaceH::errMsg("Insufficent input data in initial file!", __FILE__, __LINE__);
            return particleNum;
        }

        /**
         * @brief Write variables to a IO stream.
         * @param os
         * @param fmt
         * @param flag
         * @return
         */
        size_t write(std::ostream &os, const Unit::Fmt &fmt = Unit::UNITY_UNIT, IO_flag flag = IO_flag::STD) const {
            os << ' ' << time / fmt.T << "\r\n";

            size_t particleNum = particleNumber();
            for (size_t i = 0; i < particleNum; ++i) {

                os << ' ' << pos[i] / fmt.L << ' ' << vel[i] / fmt.V;

                if (flag == IO_flag::STD) {
                    os << ' ' << mass[i] / fmt.M
                       << ' ' << radius[i] / fmt.L
                       << ' ' << idn[i];
                }
                os << "\r\n";
            }
        }

        /**
         * @brief Input variables from a plain scalar buffer.
         * @param data
         * @param IO_flag
         * @return
         */
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {

            size_t loc = 0;
            //for locality, split into separate loops
            for (auto &p : pos) {
                p.x = data[loc++];
                p.y = data[loc++];
                p.z = data[loc++];
            }
            for (auto &v : vel) {
                v.x = data[loc++];
                v.y = data[loc++];
                v.z = data[loc++];
            }

            if (flag == IO_flag::STD) {
                for (auto &m : mass)
                    m = data[loc++];
                for (auto &r : radius)
                    r = data[loc++];
                for (auto &i : idn)
                    i = data[loc++];
            }

            time = data[loc++];
            return loc;
        }

        /**
         * @brief Write variables to plain scalar buffer.
         * @param data
         * @param IO_flag
         * @return
         */
        size_t write(ScalarBuffer &data, const IO_flag flag = IO_flag::STD) const {
            size_t particleNum = particleNumber();

            data.clear();
            data.reserve(particleNum * 6 + 1);

            //for locality, split into two loops
            for (const auto &p : pos) {
                data.emplace_back(p.x);
                data.emplace_back(p.y);
                data.emplace_back(p.z);
            }
            for (const auto &v : vel) {
                data.emplace_back(v.x);
                data.emplace_back(v.y);
                data.emplace_back(v.z);
            }

            if (flag == IO_flag::STD) {
                data.reserve(particleNum * 8 + 1);
                for (const auto &m : mass)
                    data.emplace_back(m);
                for (const auto &r : radius)
                    data.emplace_back(r);
                for (const auto &i : idn)
                    data.emplace_back(i);
            }

            data.emplace_back(time);
            return data.size();
        }

    public:
        /** @brief Position array of the particles. Element is 3D vector.*/
        VectorArray pos;

        /** @brief Velocity array of the particles. Element is 3D vector.*/
        VectorArray vel;

        /** @brief Mass array of the particles. Element is Scalar.*/
        ScalarArray mass;

        /** @brief Radius array of the particles. Element is Scalar.*/
        ScalarArray radius;

        /** @brief Id Array of the particles. Element is int.*/
        IntArray idn;

        /** @brief The physical time of the dynamic system*/
        Scalar time;
    };
}
#endif

