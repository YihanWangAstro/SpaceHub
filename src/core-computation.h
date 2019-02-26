
#ifndef LIBS_H
#define LIBS_H

#include "own-math.h"

namespace SpaceH::calcu {


    template <typename Array>
    void array_set_zero(Array &arry){
        for(auto& a : arry){
            a = 0;
        }
    }

    template<typename ...Args>
    void set_arrays_zero(Args &...args){
        (..., (array_set_zero(args)));
    }

    template<typename Array, typename ...Args>
    auto array_dot(Array const &a, Array const &b, Args const &...args) {
        typename Array::value_type product{0};
        size_t size = a.size();
        for (size_t i = 0; i < size; ++i) {
            product += a[i] * b[i] * (args[i] * ...);
        }
    }

    template<typename Array, typename ...Args>
    inline void array_add(Array &dst, Array const &a, Array const &b, Args const &...args) {
        size_t size = dst.size();

        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] + b[i] + (args[i] + ...);
        }
    }

    template<typename Array, typename ...Args>
    inline void array_mul(Array &dst, Array const &a, Array const &b, Args const &...args) {
        size_t size = dst.size();

        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] * b[i] * (args[i] + ...);
        }
    }

    /** @brief Advance a vector with a given vector with a specific stepSize.
     *  @param var Vector need to be advanced.
     *  @param increase Given vector.
     *  @param stepSize Given stepSize.
     */
    template<typename Scalar, typename Array>
    inline void array_advance(Array &var, const Array &increase, Scalar stepSize) {
        size_t size = var.size();

        for (size_t i = 0; i < size; i++) {
            var[i] += increase[i] * stepSize;
        }
    }
    /**
     * @brief Return the summation of an array. The element of the array need to suppot one parameter initialization
     * with '0' and operator '+='.
     * @tparam Array
     * @param array
     * @return
     */
    template <typename Array>
    inline auto array_sum(Array &array) {
        typename  Array::value_type total = 0;
        for(auto & a : array) {
            total += a;
        }
        return total;
    }

    /** @brief Move variables to centr
     
     al mass coordinates
     *
     *  @param mass   Array of mass.
     *  @param phyVar Array of variables need to be moved.
     */
    template<typename ScalarArray, typename VectorArray>
    auto calcuCMCoord(const ScalarArray &mass, VectorArray &phyVar) {
        typename VectorArray::value_type centralMassVar(0.0, 0.0, 0.0);
        typename ScalarArray::value_type totalMass = 0;

        const size_t N = mass.size();

        for (size_t i = 0; i < N; ++i) {
            totalMass += mass[i];
            centralMassVar += phyVar[i] * mass[i];
        }

        centralMassVar /= totalMass;
        return centralMassVar;
    }

    /** @brief Move variables to central mass coordinates
     *
     *  @param mass   Array of mass.
     *  @param phyVar Array of variables need to be moved.
     *  @param totalMass The totalMass of the system
     */
    template<typename ScalarArray, typename VectorArray, typename Scalar>
    auto calcuCMCoord(const ScalarArray &mass, VectorArray &phyVar, const Scalar totalMass) {
        typename VectorArray::value_type centralMassVar(0.0, 0.0, 0.0);

        const size_t N = mass.size();

        for (size_t i = 0; i < N; ++i)
            centralMassVar += phyVar[i] * mass[i];

        centralMassVar /= totalMass;
        return centralMassVar;
    }

    /** @brief Move variables to central mass coordinates
     *
     *  @param mass   Array of mass.
     *  @param phyVar Array of variables need to be moved.
     */
    template<typename VectorArray>
    void moveToCMCoord(VectorArray &phyVar, const typename VectorArray::value_type &centralMassVar) {
        const size_t N = phyVar.size();

        for (auto& var : phyVar)
            var -= centralMassVar;
    }

    /** @brief Move variables to central mass coordinates
     *
     *  @param mass   Array of mass.
     *  @param phyVar Array of variables need to be moved.
     *  @param totalMass The totalMass of the system
     */
    template<typename ScalarArray, typename VectorArray, typename Scalar>
    void moveToCoM(const ScalarArray &mass, VectorArray &phyVar, const Scalar totalMass) {
        typename VectorArray::value_type centralMassVar(0.0, 0.0, 0.0);

        const size_t N = mass.size();

        for (size_t i = 0; i < N; ++i)
            centralMassVar += phyVar[i] * mass[i];

        centralMassVar /= totalMass;

        for (auto& var : phyVar)
            var -= centralMassVar;
    }

    /** @brief Move variables to central mass coordinates
     *
     *  @param mass   Array of mass.
     *  @param phyVar Array of variables need to be moved.
     *  @param totalMass The totalMass of the system
     */
    template<typename ScalarArray, typename VectorArray>
    void moveToCoM(const ScalarArray &mass, VectorArray &phyVar) {
        using Vector = typename VectorArray::value_type;
        using Scalar = typename Vector ::value_type;

        Vector centralMassVar(0.0, 0.0, 0.0);
        Scalar totalMass{0};

        const size_t N = mass.size();
        for (size_t i = 0; i < N; ++i) {
            centralMassVar += phyVar[i] * mass[i];
            totalMass += mass[i];
        }

        centralMassVar /= totalMass;

        for (auto& var : phyVar)
            var -= centralMassVar;
    }

    template<typename Particles>
    auto get_kinetic_energy(Particles const &partc) {
        typename Particles::Scalar kineticEnergy{0};

        size_t size = partc.number();

        for (size_t i = 0; i < size; ++i)
            kineticEnergy += 0.5 * partc.mass[i] * (partc.vx[i]*partc.vx[i] + partc.vy[i]*partc.vy[i] + partc.vz[i]*partc.vz[i]);

        return kineticEnergy;
    }

    template<typename Particles>
    auto get_potential_energy(Particles const&partc) {
        typename Particles::Scalar potentialEnergy{0};

        size_t size = partc.number();

        for (size_t i = 0; i < size; ++i)
            for (size_t j = i + 1; j < size; ++j){
                auto dx = partc.px[i] - partc.px[j];
                auto dy = partc.py[i] - partc.py[j];
                auto dz = partc.pz[i] - partc.pz[j];
                potentialEnergy -= partc.mass[i] * partc.mass[j] / sqrt(dx*dx + dy*dy + dz*dz);
            }
        return potentialEnergy;
    }

    template<typename Particles>
    auto get_total_energy(Particles const &partc) {
        typename Particles::Scalar potentialEnergy{0};
        typename Particles::Scalar kineticEnergy{0};
        size_t size = partc.number();

        for (size_t i = 0; i < size; ++i) {
            kineticEnergy += 0.5 * partc.mass[i] * (partc.vx[i]*partc.vx[i] + partc.vy[i]*partc.vy[i] + partc.vz[i]*partc.vz[i]);
            for (size_t j = i + 1; j < size; ++j){
                auto dx = partc.px[i] - partc.px[j];
                auto dy = partc.py[i] - partc.py[j];
                auto dz = partc.pz[i] - partc.pz[j];
                potentialEnergy -= partc.mass[i] * partc.mass[j] / sqrt(dx*dx + dy*dy + dz*dz);
            }
        }
        return potentialEnergy + kineticEnergy;
    }


    /** @brief Calculate the minimal fall free time of two particles
     *
     *  @param  mass mass array of particle.
     *  @param  pos  position array of particle.
     *  @return The minimal fall free time of the two particles
     */
    template<typename ScalarArray, typename VectorArray>
    inline auto minfallFreeTime(const ScalarArray &mass, const VectorArray &pos) {
        size_t size = mass.size();
        typename ScalarArray::value_type r = 0;
        typename ScalarArray::value_type min_fall_free = 0;
        typename ScalarArray::value_type fall_free = 0;

        for (size_t i = 0; i < size; i++) {
            for (size_t j = i + 1; j < size; j++) {
                r = (pos[i] - pos[j]).norm();
                fall_free = pow(r, 1.5) / (mass[i] + mass[j]);
                min_fall_free = min_fall_free == 0 ? fall_free : SpaceH::min(min_fall_free, fall_free);
            }
        }

        return 0.1 * min_fall_free;
    }

    /** @brief check if all elements in an vector array are zero vector
     *
     *  @param  array The array need to be checked.
     *  @return The boolean result.
     */
    template<typename VectorArray>
    bool isAllZero(const VectorArray &array) {
        for (auto& a : array) {
            if (a.norm() > 0)
                return false;
        }
        return true;
    }

    /** @brief Calculate the minimal fall free time of two particles
     *
     *  @param  mass mass array of particle.
     *  @param  pos  position array of particle.
     *  @param  vel  velocity array of particle.
     *  @return The minimal fall free time of the two particles
     */
    template<typename ScalarArray, typename VectorArray>
    inline auto minAccdot(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel) {
        using Vector = typename VectorArray::value_type;
        using Scalar = typename ScalarArray::value_type;
        size_t size = mass.size();
        Vector acc;
        Vector adot;
        Vector dr;
        Vector dv;

        Scalar r = 0;
        Scalar min_ds = 0;

        for (size_t i = 0; i < size; i++) {
            adot.setZero();
            acc.setZero();
            for (size_t j = 0; j < size; j++) {
                if (i != j) {
                    dr = pos[i] - pos[j];
                    dv = vel[i] - vel[j];
                    r = dr.norm();
                    acc += dr / (r * r * r) * mass[j];
                    adot += (dv / (r * r * r) + dr * (3 * dot(dv, dr) / (r * r * r * r * r))) * mass[j];
                }
            }
            Scalar ds = acc.norm() / adot.norm();
            min_ds = min_ds == 0 ? ds : SpaceH::min(min_ds, ds);
        }

        return min_ds;
    }




    /** @brief Calculate the potential energy of particles
     *
     *  @param  mass  Array of mass.
     *  @param  pos   Array of position.
     *  @param  chpos Array of chain position.
     *  @param  chind Array of chain index.
     *  @return The potential energy.
     */
    template<typename ScalarArray, typename VectorArray, typename IndexArray>
    auto getPotentialEnergy(const ScalarArray &mass, const VectorArray &pos, const VectorArray &chainPos,
                       const IndexArray &chainInd) {
        typename ScalarArray::value_type potentialEnergy = 0;

        size_t size = mass.size();

        for (size_t i = 0; i < size - 1; ++i)
            potentialEnergy -= mass[chainInd[i]] * mass[chainInd[i + 1]] / chainPos[i].norm();

        for (size_t i = 0; i < size - 2; ++i)
            potentialEnergy -= mass[chainInd[i]] * mass[chainInd[i + 2]] / (chainPos[i] + chainPos[i + 1]).norm();

        for (size_t i = 0; i < size; ++i)
            for (size_t j = i + 3; j < size; ++j)
                potentialEnergy -= mass[chainInd[i]] * mass[chainInd[j]] / distance(pos[chainInd[j]], pos[chainInd[i]]);

        return potentialEnergy;
    }



    /** @brief Calculate the energy error approximately for regularized system.
     *
     *  @param  mass  Array of mass.
     *  @param  pos   Array of position.
     *  @param  vel   Array of velocity.
     *  @param  bindE The binding energy of the regularized system.
     *  @return The approximated energy error.
     */
    template<typename ScalarArray, typename VectorArray, typename Scalar>
    inline Scalar
    getEnergyErr(const ScalarArray &mass, const VectorArray &pos, const VectorArray &vel, const Scalar bindE) {
        Scalar EK = get_kinetic_energy(mass, vel);
        Scalar EP = get_potential_energy(mass, pos);
        return log(abs((EK + bindE) / EP));
    }

}//end namespace SpaceH
#endif
