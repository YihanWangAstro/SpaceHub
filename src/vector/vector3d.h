
#ifndef GENVECT3D_H
#define GENVECT3D_H

#ifdef __AVX__
#pragma message("Using AVX on vector3d")

#include <x86intrin.h>
namespace space
{

/** @brief Specilization of vector3d */
template<>
struct Vec3<double>
{
public:
    /* Typedef */
    using value_type = double;
    /* Typedef */

    union
    {
        struct
        {
            double x, y, z;
        };
        __m256d  mmvalue;
    } __attribute__((aligned(32)));

    Vec3(): mmvalue(_mm256_setzero_pd()) {};
    Vec3(double vx, double vy, double vz) : mmvalue(_mm256_set_pd(0.0, vz, vy, vx)) {};
    Vec3(double scalar) : mmvalue(_mm256_set1_pd(scalar)) {};
    Vec3(__m256d v) : mmvalue(v) {};
    Vec3(const Vec3& v) : mmvalue(v.mmvalue) {};
    /** @brief Addition by wise */
    inline Vec3 operator+(const Vec3& v) const
    {
        return Vec3(_mm256_add_pd(mmvalue, v.mmvalue));
    }
    /** @brief Subtraction by wise */
    inline Vec3 operator-(const Vec3& v) const
    {
        return Vec3(_mm256_sub_pd(mmvalue, v.mmvalue));
    }
    /** @brief Product by wise */
    inline Vec3 operator*(const Vec3& v) const
    {
        return Vec3(_mm256_mul_pd(mmvalue, v.mmvalue));
    }
    /** @brief Divition by wise */
    inline Vec3 operator/(const Vec3& v) const
    {
        return Vec3(_mm256_div_pd(mmvalue, v.mmvalue));
    }
    /** @brief Opposite vector */
    inline Vec3 operator-() const
    {
        //return Vec3(-x, -y, -z);
        return Vec3(-mmvalue);
    }
    /** @brief Absolute value by wise */
    inline Vec3 abs() const
    {
        //return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z);
        return Vec3(_mm256_max_pd(mmvalue,-mmvalue));
    }
    inline const Vec3& operator+=(const Vec3& v)
    {
        mmvalue = _mm256_add_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const Vec3& operator-=(const Vec3& v)
    {
        mmvalue = _mm256_sub_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const Vec3& operator*=(const Vec3& v)
    {
        mmvalue = _mm256_mul_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const Vec3& operator/=(const Vec3& v)
    {
        mmvalue = _mm256_div_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const Vec3& operator=(const Vec3& v)
    {
        mmvalue = v.mmvalue;
        return *this;
    }
    /** @brief Calculate the norm */
    inline double norm() const
    {
        Vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
        return sqrt(product.x + product.y + product.z);
    }
    /** @brief Calculate the norm */
    inline double norm2() const
    {
        Vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
        return product.x + product.y + product.z;
    }

    inline double max_component()
    {
        double max = (x > y ? x : y);
        return max > z ? max : z;
    }
    /** @brief Calculate the inverse of the norm */
    inline double reNorm() const
    {
        Vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
        return 1.0/sqrt(product.x + product.y + product.z);
    }
    inline void setZero()
    {
        mmvalue = _mm256_setzero_pd();
    }
    /** @brief Output to ostream */
    friend std::ostream& operator<<(std::ostream& output, const Vec3& v)
    {
        output << v.x << " " << v.y << " " << v.z;
        return output;
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, Vec3& v)
    {
        input >> v.x >> v.y >> v.z;
        return input;
    }
};

/** @brief Calculate the Euclid distance of two vectors */
template<>
inline double distance<double>(const Vec3<double>& v1, const Vec3<double>& v2)
{
    __m256d sub = _mm256_sub_pd(v1.mmvalue, v2.mmvalue);
    Vec3<double> product = _mm256_mul_pd(sub, sub);
    return sqrt(product.x + product.y + product.z);
}
/** @brief Calculate the inner product of two vectors */
template<>
inline double dot<double>(const Vec3<double>& v1, const Vec3<double>& v2)
{
    Vec3<double>  product = _mm256_mul_pd(v1.mmvalue, v2.mmvalue);
    return product.x + product.y + product.z;
}
/** @brief Calculate the cross product of two vectors */
template<>
inline Vec3<double> cross<double>(const Vec3<double>& v1, const Vec3<double>& v2)
{
    return Vec3<double>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

}
#endif
#endif

