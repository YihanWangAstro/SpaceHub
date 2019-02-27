
#ifndef LIBS_H
#define LIBS_H

#include "own-math.h"

namespace SpaceH::calc {


    template <typename ...Args>
    auto add(Args && ...args){
        return (... + args);
    }

    template <typename ...Args>
    auto mul(Args && ...args){
        return (... * args);
    }

    template <typename ...Args>
    auto any(Args ...args){
        return (... || args);
    }

    template <typename ...Args>
    auto all(Args ...args){
        return (... && args);
    }

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
    void array_add(Array &dst, Array const &a, Array const &b, Args const &...args) {
        size_t size = dst.size();

        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] + b[i] + (args[i] + ...);
        }
    }

    template<typename Array, typename ...Args>
    void array_mul(Array &dst, Array const &a, Array const &b, Args const &...args) {
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
    void array_advance(Array &var, const Array &increase, Scalar stepSize) {
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
    auto array_sum(Array &array) {
        typename  Array::value_type total = 0;
        for(auto & a : array) {
            total += a;
        }
        return total;
    }


    template<typename Array1, typename Array2>
    auto calc_com(Array1 const & mass, Array2 const & var) {
        using T = typename Array2::value_type;
        using Scalar = typename Array1::value_type;

        T com_var{0};
        Scalar tot_mass{0};

        size_t const size = mass.size();

        for (size_t i = 0; i < size; ++i) {
            com_var += var[i] * mass[i];
            tot_mass += mass[i];
        }
        com_var /= tot_mass;
        return com_var;
    }

    template<typename Array>
    void move_to_com(Array &var, typename Array::value_type const &com_var) {
        for (auto& v : var)
            v -= com_var;
    }

    template<typename Array1, typename Array2>
    void move_to_com(Array1 const & mass, Array2 & var) {
        auto com_var = calc_com(mass, var);
        move_to_com(var, com_var);
    }

    template<typename Particles>
    auto calc_kinetic_energy(Particles const &ptc) {
        typename Particles::Scalar k_eng{0};

        size_t size = ptc.number();

        for (size_t i = 0; i < size; ++i)
            k_eng += 0.5 * ptc.mass(i) * (ptc.vx(i) * ptc.vx(i) + ptc.vy(i) * ptc.vy(i) + ptc.vz(i) * ptc.vz(i));

        return k_eng;
    }

    template<typename Particles>
    auto calc_potential_energy(Particles const &ptc) {
        typename Particles::Scalar p_eng{0};

        size_t size = ptc.number();

        for (size_t i = 0; i < size; ++i)
            for (size_t j = i + 1; j < size; ++j) {
                auto dx = ptc.px(i) - ptc.px(j);
                auto dy = ptc.py(i) - ptc.py(j);
                auto dz = ptc.pz(i) - ptc.pz(j);
                p_eng -= ptc.mass(i) * ptc.mass(j) / sqrt(dx * dx + dy * dy + dz * dz);
            }
        return p_eng;
    }

    template<typename Particles>
    auto calc_total_energy(Particles const &ptc) {
        typename Particles::Scalar p_eng{0};
        typename Particles::Scalar k_eng{0};
        size_t size = ptc.number();

        for (size_t i = 0; i < size; ++i) {
            k_eng += 0.5 * ptc.mass(i) * (ptc.vx(i) * ptc.vx(i) + ptc.vy(i) * ptc.vy(i) + ptc.vz(i) * ptc.vz(i));
            for (size_t j = i + 1; j < size; ++j) {
                auto dx = ptc.px(i) - ptc.px(j);
                auto dy = ptc.py(i) - ptc.py(j);
                auto dz = ptc.pz(i) - ptc.pz(j);
                p_eng -= ptc.mass(i) * ptc.mass(j) / sqrt(dx * dx + dy * dy + dz * dz);
            }
        }
        return p_eng + k_eng;
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
        Scalar EK = calc_kinetic_energy(mass, vel);
        Scalar EP = calc_potential_energy(mass, pos);
        return log(abs((EK + bindE) / EP));
    }

}//end namespace SpaceH
#endif
