
#ifndef POSTNEWTONIAN_H
#define POSTNEWTONIAN_H

#include "../protoType.h"
#include "../macros.h"

namespace SpaceH {
    namespace Const {
        constexpr double INV_C = 1 / C;
        constexpr double INV_C2 = INV_C * INV_C;
        constexpr double INV_C3 = INV_C2 * INV_C;
        constexpr double INV_C4 = INV_C3 * INV_C;
        constexpr double INV_C5 = INV_C4 * INV_C;
    }

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct NewtonForce {
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;

        /* Typedef */

        inline void
        operator()(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2, const Vector &pos1, const Vector &pos2,
                   const Vector &dr) {
            Scalar inv_r = dr.reNorm();
            Scalar inv_r3 = inv_r * inv_r * inv_r;
            acc1 += dr * (inv_r3 * m2);
            acc2 -= dr * (inv_r3 * m1);
        }
    };

    /**
     *
     * @tparam T
     * @param x
     * @return
     */
    template<typename T>
    inline T KarmackFastInverseSquareRoot(T x) {
        return 1 / sqrt(x);
    }

    /**
     *
     * @param x
     * @return
     */
    template<>
    inline float KarmackFastInverseSquareRoot<float>(float x) {
        float xhalf = 0.5f * x;
        int i = *(int *) &x;
        //i = 0x5f3759df - (i >> 1);
        i = 0x5f375a86 - (i >> 1);
        x = *(float *) &i;
        x = x * (1.5f - xhalf * x * x);
        //x = x*(1.5f - xhalf*x*x);
        return x;
    }

    /**
     *
     * @param x
     * @return
     */
    template<>
    inline double KarmackFastInverseSquareRoot<double>(double x) {
        double xhalf = 0.5f * x;
        long long i = *(long long *) &x;
        i = 0x5fe6eb50c7aa19f9 - (i >> 1);
        x = *(double *) &i;
        x = x * (1.5f - xhalf * x * x);
        //x = x*(1.5f - xhalf*x*x);
        return x;
    }

    /**
     *
     * @tparam Dtype
     * @tparam ArraySize
     */
    template<typename Dtype, size_t ArraySize>
    struct KarmackNewtonian {
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;

        /* Typedef */

        inline void
        operator()(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2, const Vector &pos1, const Vector &pos2,
                   const Vector &dr) {
            Scalar r2 = dr * dr;
            Scalar inv_r = KarmackFastInverseSquareRoot<typename SpaceH::get_value_type<Scalar>::type>(r2);
            Scalar inv_r3 = inv_r * inv_r * inv_r;
            acc1 += dr * (inv_r3 * m2);
            acc2 -= dr * (inv_r3 * m1);
        }
    };


    /**
     * @brief Post newtonian pair interaction functor(c++ std11)*
     * @tparam Dtype
     * @tparam ArraySize
     * @tparam First
     * @tparam Second
     * @tparam Radiative
     */
    template<typename Dtype, size_t ArraySize, bool First, bool Second, bool Radiative>
    class PostNewtonianForce {
    public:
        /* Typedef */
        using type   = SpaceH::ProtoType<Dtype, ArraySize>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */

        /**
         *
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         */
        inline void operator()(Vector &acc1, Vector &acc2,
                               const Scalar m1, const Scalar m2,
                               const Vector &p1, const Vector &p2,
                               const Vector &v1, const Vector &v2,
                               const Vector &dr, const Vector &dv) {
            RadiativeOrder<Radiative>(acc1, acc2, m1, m2, p1, p2, v1, v2, dr, dv);
            SecondOrder<Second>(acc1, acc2, m1, m2, p1, p2, v1, v2, dr, dv);
            FirstOrder<First>(acc1, acc2, m1, m2, p1, p2, v1, v2, dr, dv);
        }

    private:
        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<!Swith>::type
        FirstOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                   const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                   const Vector &dr, const Vector &dv) {}

        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<Swith>::type
        FirstOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                   const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                   const Vector &dr, const Vector &dv) {
            Scalar inv_r = dr.reNorm();
            Scalar inv_r2 = inv_r * inv_r;
            Vector n12(dr * (-inv_r));
            Scalar nv1 = dot(n12, v1);
            Scalar nv2 = dot(n12, v2);
            Scalar A1 = 0.0, B1 = 0.0;

            A1 = ((5.0 * m1 + 4.0 * m2) * inv_r +
                  (1.5 * nv2 * nv2 - dot(v1, v1) + 4.0 * dot(v1, v2) - 2.0 * dot(v2, v2))) * m2 * inv_r2 *
                 Const::INV_C2;
            B1 = (4.0 * nv1 - 3.0 * nv2) * m2 * inv_r2 * Const::INV_C2;
            acc1 += n12 * A1 - dv * B1;
            A1 = ((5.0 * m2 + 4.0 * m1) * inv_r +
                  (1.5 * nv1 * nv1 - dot(v2, v2) + 4.0 * dot(v2, v1) - 2.0 * dot(v1, v1))) * m1 * inv_r2 *
                 Const::INV_C2;
            B1 = (-4.0 * nv2 + 3.0 * nv1) * m1 * inv_r2 * Const::INV_C2;
            acc2 -= n12 * A1 - dv * B1;
        }

        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<!Swith>::type
        SecondOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                    const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                    const Vector &dr, const Vector &dv) {}

        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<Swith>::type
        SecondOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                    const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                    const Vector &dr, const Vector &dv) {
            Scalar inv_r = dr.reNorm();
            Scalar inv_r2 = inv_r * inv_r;
            Vector n12(dr * (-inv_r));
            Scalar nv1 = dot(n12, v1);
            Scalar nv2 = dot(n12, v2);
            Scalar nv1s = nv1 * nv1;
            Scalar nv2s = nv2 * nv2;
            Scalar v1s = dot(v1, v1);
            Scalar v2s = dot(v2, v2);
            Scalar v12 = dot(v1, v2);
            Scalar C2 = inv_r2 * Const::INV_C4;
            Scalar A2 = 0.0, B2 = 0.0;

            A2 = (-2.0 * v2s * v2s + (4.0 * v2s - 2.0 * v12) * v12 +
                  (1.5 * v1s + 4.5 * v2s - 6.0 * v12 - 1.875 * nv2s) * nv2s
                  + m1 * inv_r * (-3.75 * v1s + 1.25 * v2s - 2.5 * v12 + 19.5 * nv1s - 39.0 * nv1 * nv2 + 8.5 * nv2s)
                  + m2 * inv_r * (4.0 * v2s - 8.0 * v12 + 2.0 * nv1s - 4.0 * nv1 * nv2 - 6.0 * nv2s)
                  + (-14.25 * m1 * m1 - 9.0 * m2 * m2 - 34.5 * m1 * m2) * inv_r2) * C2 * m2;
            B2 = ((4.0 * v2s - 4.0 * v12 - 6.0 * nv2s - 15.75 * m1 * inv_r - 2.0 * m2 * inv_r) * nv1
                  + (v1s - 5.0 * v2s + 4.0 * v12 + 4.5 * nv2s + 13.75 * m1 * inv_r - 2.0 * m2 * inv_r) * nv2) * C2 * m2;

            acc1 += n12 * A2 - dv * B2;

            A2 = (-2.0 * v1s * v1s + (4.0 * v1s - 2.0 * v12) * v12 +
                  (1.5 * v2s + 4.5 * v1s - 6.0 * v12 - 1.875 * nv1s) * nv1s
                  + m2 * inv_r * (-3.75 * v2s + 1.25 * v1s - 2.5 * v12 + 19.5 * nv2s - 39.0 * nv2 * nv1 + 8.5 * nv1s)
                  + m1 * inv_r * (4.0 * v1s - 8.0 * v12 + 2.0 * nv2s - 4.0 * nv2 * nv1 - 6.0 * nv1s)
                  + (-14.25 * m2 * m2 - 9.0 * m1 * m1 - 34.5 * m2 * m1) * inv_r2) * C2 * m1;
            B2 = -((4.0 * v1s - 4.0 * v12 - 6.0 * nv1s - 15.75 * m2 * inv_r - 2.0 * m1 * inv_r) * nv2
                   + (v2s - 5.0 * v1s + 4.0 * v12 + 4.5 * nv1s + 13.75 * m2 * inv_r - 2.0 * m1 * inv_r) * nv1) * C2 *
                 m1;

            acc2 -= n12 * A2 - dv * B2;
        }

        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<!Swith>::type
        RadiativeOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                       const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                       const Vector &dr, const Vector &dv) {}

        /**
         *
         * @tparam Swith
         * @param acc1
         * @param acc2
         * @param m1
         * @param m2
         * @param p1
         * @param p2
         * @param v1
         * @param v2
         * @param dr
         * @param dv
         * @return
         */
        template<bool Swith>
        inline typename std::enable_if<Swith>::type
        RadiativeOrder(Vector &acc1, Vector &acc2, const Scalar m1, const Scalar m2,
                       const Vector &p1, const Vector &p2, const Vector &v1, const Vector &v2,
                       const Vector &dr, const Vector &dv) {
            Scalar inv_r = dr.reNorm();
            Vector n12(dr * (-inv_r));
            Scalar nv12 = dot(n12, (v1 - v2));
            Scalar dvs = dot(dv, dv);
            Scalar C3 = 0.8 * inv_r * inv_r * inv_r * Const::INV_C5 * m1 * m2;
            Scalar A25 = 0.0, B25 = 0.0;

            A25 = nv12 * (3.0 * dvs + (52.0 / 3.0 * m2 - 6.0 * m1) * inv_r) * C3;
            B25 = (-dvs + (2.0 * m1 - 8.0 * m2) * inv_r) * C3;

            acc1 += n12 * A25 - dv * B25;

            A25 = nv12 * (3.0 * dvs + (52.0 / 3.0 * m1 - 6.0 * m2) * inv_r) * C3;
            B25 = (-dvs + (2.0 * m2 - 8.0 * m1) * inv_r) * C3;

            acc2 -= n12 * A25 - dv * B25;
        }
    };
}

#endif
