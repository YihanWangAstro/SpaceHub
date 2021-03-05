
#pragma once

#include "../lazy-evaluation/lazy_expr.h"
#include "vector3.hpp"
namespace hub::lazy {

    template <typename T>
    class Vec3;

    template <typename T>
    struct LazyVec3 : public Expr<LazyVec3<T>> {
       public:
        using value_type = T;

        union {
            struct {
                T x;
                T y;
                T z;
            };
            T data_[3];
        };

        LazyVec3() : x{0}, y{0}, z{0} {}

        LazyVec3(value_type s) : x(s), y(s), z(s) {}

        LazyVec3(value_type vx, value_type vy, value_type vz) : x(vx), y(vy), z(vz) {}

        LazyVec3(LazyVec3 const &v) { x = v.x, y = v.y, z = v.z; }

        ~LazyVec3() {}

        template <typename U>
        LazyVec3(const Expr<U> &expr) {
            const U &src = expr.cast();
            data_[0] = src.eval(0);
            data_[1] = src.eval(1);
            data_[2] = src.eval(2);
        }

        template <typename U>
        LazyVec3(const Vec3<U> &v) {
            x = v.x, y = v.y, z = v.z;
        }

        T const &eval(size_t i) const { return data_[i]; }

        LazyVec3 &operator=(LazyVec3 const &src) noexcept {
            this->x = src.x;
            this->y = src.y;
            this->z = src.z;
            return *this;
        }

        template <typename U>
        LazyVec3 &operator=(Expr<U> const &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            data_[0] = rhs.eval(0);
            data_[1] = rhs.eval(1);
            data_[2] = rhs.eval(2);
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator+=(Expr<U> const &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            data_[0] += rhs.eval(0);
            data_[1] += rhs.eval(1);
            data_[2] += rhs.eval(2);
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator-=(Expr<U> const &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            data_[0] -= rhs.eval(0);
            data_[1] -= rhs.eval(1);
            data_[2] -= rhs.eval(2);
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator*=(Expr<U> const &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            data_[0] *= rhs.eval(0);
            data_[1] *= rhs.eval(1);
            data_[2] *= rhs.eval(2);
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator/=(Expr<U> const &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            data_[0] /= rhs.eval(0);
            data_[1] /= rhs.eval(1);
            data_[2] /= rhs.eval(2);
            return *this;
        }

        template <typename U>
        LazyVec3 &operator=(Vec3<U> const &v) {
            x = v.x, y = v.y, z = v.z;
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator+=(Vec3<U> const &v) {
            x += v.x, y += v.y, z += v.z;
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator-=(Vec3<U> const &v) {
            x -= v.x, y -= v.y, z -= v.z;
            return *this;
        }

        template <typename U>
        LazyVec3 &operator*=(Vec3<U> const &v) {
            x *= v.x, y *= v.y, z *= v.z;
            return *this;
        }

        template <typename U>
        inline LazyVec3 &operator/=(Vec3<U> const &v) {
            x /= v.x, y /= v.y, z /= v.z;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &output, const LazyVec3 &v) {
            output << v.x << ',' << v.y << ',' << v.z;
            return output;
        }

        friend std::istream &operator>>(std::istream &input, LazyVec3 &v) {
            input >> v.x >> v.y >> v.z;
            return input;
        }
    };

    template <typename T1, typename T2>
    inline T1 distance(const LazyVec3<T1> &v1, const LazyVec3<T2> &v2) {
        auto dx = v1.x - v2.x;
        auto dy = v1.y - v2.y;
        auto dz = v1.z - v2.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    /** Calculate the inner product of two vectors */
    template <typename T1, typename T2>
    inline T1 dot(const LazyVec3<T1> &v1, const LazyVec3<T2> &v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    /** Calculate the length of a vector*/
    template <typename T>
    inline T norm(const LazyVec3<T> &v) {
        return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    /**  Calculate the inverse length of a vector*/
    template <typename T>
    inline T re_norm(const LazyVec3<T> &v) {
        return 1.0 / sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    /** Calculate the length square of a vector*/
    template <typename T>
    inline T norm2(const LazyVec3<T> &v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    /** Calculate the cross product of two vectors */
    template <typename T1, typename T2>
    inline LazyVec3<T1> cross(const LazyVec3<T1> &v1, const LazyVec3<T2> &v2) {
        return LazyVec3<T1>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

    template <typename T>
    inline LazyVec3<T> vec_abs(LazyVec3<T> const &v) {
        return LazyVec3<T>(fabs(v.x), fabs(v.y), fabs(v.z));
    }

    template <typename T>
    inline LazyVec3<T> vec_max(LazyVec3<T> const &v1, LazyVec3<T> const &v2) {
        return LazyVec3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    template <typename T>
    inline LazyVec3<T> vec_min(LazyVec3<T> const &v1, LazyVec3<T> const &v2) {
        return LazyVec3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    template <typename T>
    inline T max_abs(LazyVec3<T> const &v) {
        auto max = std::max(fabs(v.x), fabs(v.y));
        return std::max(max, fabs(v.z));
    }

}  // namespace hub::lazy
