

#ifndef KAHANNUMBER_h
#define KAHANNUMBER_h

#include "own-math.h"

namespace space {
/** @brief Kahan number
 *
 *  A way to reduce the round off error when adding a small number to a big one.
 *  See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
    template<typename T>
    struct kahan {
    public:
        using value_type = T;

        T real, err;

        kahan() : real(0), err(0) {};

        kahan(T r) : real(r), err(0) {};

        kahan(const kahan &k) : real(k.real), err(k.err) {};

        inline const kahan &operator=(const kahan &hs) {
            real = hs.real, err = 0;
            return *this;
        }

        inline operator T() {
            return real;
        }

        inline operator T() const {
            return real;
        }

        inline void zero_err() {
            err = 0;
        }

        friend inline kahan operator-(const kahan &hs) {
            return kahan(-hs.real);
        }

        friend inline const kahan &operator+=(kahan &lhs, const kahan &rhs) {
            T add = rhs.real - lhs.err;
            T sum = lhs.real + add;

            /*if (SpaceH::abs(add) < SpaceH::abs(lhs.real))
                lhs.err = (sum - lhs.real) - add;
            else
                lhs.err = (sum - add) - lhs.real;*/
            lhs.err  = (sum - lhs.real) - add;
            lhs.real = sum;
            return lhs;
        }

        friend inline const kahan &operator-=(kahan &lhs, const kahan &rhs) {
            T add = -rhs.real - lhs.err;
            T sum = lhs.real + add;

            /*if (SpaceH::abs(add) < SpaceH::abs(lhs.real))
                lhs.err = (sum - lhs.real) - add;
            else
                lhs.err = (sum - add) - lhs.real;*/
            lhs.err  = (sum - lhs.real) - add;
            lhs.real = sum;
            return lhs;
        }

        friend inline const kahan &operator/=(kahan &lhs, const kahan &rhs) {
            lhs.real /= rhs.real;
            return lhs;
        }

        friend inline const kahan &operator*=(kahan &lhs, const kahan &rhs) {
            lhs.real *= rhs.real;
            return lhs;
        }

        /** @brief Output to ostream */
        friend std::ostream &operator<<(std::ostream &output, const kahan &v) {
            output << v.real;
            return output;
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &input, kahan &v) {
            input >> v.real;
            v.err = 0;
            return input;
        }
    };

    using precise_d  = kahan<double>;
    using precise_f  = kahan<float>;
}//end namespace SpaceH
#endif /* kahanNumber_h */
