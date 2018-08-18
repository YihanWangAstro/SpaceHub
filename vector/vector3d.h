
#ifndef GENVECT3D_H
#define GENVECT3D_H

#ifdef __AVX__
#pragma message("Using AVX on vector3d")
//#pragma G++ optimize("O3")
#include <x86intrin.h>
namespace SpaceH
{

/** @brief Specilization of vector3d */
template<>
struct vec3<double>
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

    vec3(): mmvalue(_mm256_setzero_pd()) {};
    vec3(double vx, double vy, double vz) : mmvalue(_mm256_set_pd(0.0, vz, vy, vx)) {};
    vec3(__m256d v) : mmvalue(v) {};
    vec3(const vec3& v) : mmvalue(v.mmvalue) {};
    /** @brief Addition by wise */
    inline vec3 operator+(const vec3& v) const
    {
        return vec3(_mm256_add_pd(mmvalue, v.mmvalue));
    }
    /** @brief Subtraction by wise */
    inline vec3 operator-(const vec3& v) const
    {
        return vec3(_mm256_sub_pd(mmvalue, v.mmvalue));
    }
    /** @brief Product by wise */
    inline vec3 operator*(const vec3& v) const
    {
        return vec3(_mm256_mul_pd(mmvalue, v.mmvalue));
    }
    /** @brief Divition by wise */
    inline vec3 operator/(const vec3& v) const
    {
        return vec3(_mm256_div_pd(mmvalue, v.mmvalue));
    }
    /** @brief Add scalar by wise */
    inline vec3 operator+(const double c) const
    {
        return vec3(_mm256_add_pd(mmvalue, _mm256_set1_pd(c)));
    }
    /** @brief Subtract scalar by wise */
    inline vec3 operator-(const double c) const
    {
        return vec3(_mm256_sub_pd(mmvalue, _mm256_set1_pd(c)));
    }
    /** @brief Multiply scalar by wise */
    inline vec3 operator*(const double c) const
    {
        return vec3(_mm256_mul_pd(mmvalue, _mm256_set1_pd(c)));
    }
    /** @brief Divide scalar by wise */
    inline vec3 operator/(const double c) const
    {
        return vec3(_mm256_div_pd(mmvalue, _mm256_set1_pd(c)));
    }
    /** @brief Opposite vector */
    inline vec3 operator-() const
    {
        //return vec3(-x, -y, -z);
        return vec3(-mmvalue);
    }
    /** @brief Absolute value by wise */
    inline vec3 abs() const
    {
        //return vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z);
        return vec3(_mm256_max_pd(mmvalue,-mmvalue));
    }
    inline const vec3& operator+=(const vec3& v)
    {
        mmvalue = _mm256_add_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator-=(const vec3& v)
    {
        mmvalue = _mm256_sub_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator*=(const vec3& v)
    {
        mmvalue = _mm256_mul_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator/=(const vec3& v)
    {
        mmvalue = _mm256_div_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator+=(const double c)
    {
        mmvalue = _mm256_add_pd(mmvalue, _mm256_set1_pd(c));
        return *this;
    }
    inline const vec3& operator-=(const double c)
    {
        mmvalue = _mm256_sub_pd(mmvalue, _mm256_set1_pd(c));
        return *this;
    }
    inline const vec3& operator*=(const double c)
    {
        mmvalue = _mm256_mul_pd(mmvalue, _mm256_set1_pd(c));
        return *this;
    }
    inline const vec3& operator/=(const double c)
    {
        mmvalue = _mm256_div_pd(mmvalue, _mm256_set1_pd(c));
        return *this;
    }
    inline const vec3& operator=(const vec3& v)
    {
        mmvalue = v.mmvalue;
        return *this;
    }
    /** @brief Calculate the norm */
    inline double norm() const
    {
        vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
        return sqrt(product.x + product.y + product.z);
    }
    /** @brief Calculate the norm */
    inline double norm2() const
    {
        vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
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
        vec3  product = _mm256_mul_pd(mmvalue, mmvalue);
        return 1.0/sqrt(product.x + product.y + product.z);
    }
    inline void setZero()
    {
        mmvalue = _mm256_setzero_pd();
    }
    friend vec3 operator+(const double c, const vec3& v)
    {
        return vec3(_mm256_add_pd(_mm256_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator-(const double c, const vec3& v)
    {
        return vec3(_mm256_sub_pd(_mm256_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator*(const double c, const vec3& v)
    {
        return vec3(_mm256_mul_pd(_mm256_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator/(const double c, const vec3& v)
    {
        return vec3(_mm256_div_pd(_mm256_set1_pd(c), v.mmvalue));
    }
    /** @brief Output to ostream */
    friend std::ostream& operator<<(std::ostream& output, const vec3& v)
    {
        output << v.x << " " << v.y << " " << v.z;
        return output;
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, vec3& v)
    {
        input >> v.x >> v.y >> v.z;
        return input;
    }
};

/** @brief Calculate the Euclid distance of two vectors */
template<>
inline double distance<double>(const vec3<double>& v1, const vec3<double>& v2)
{
    __m256d sub = _mm256_sub_pd(v1.mmvalue, v2.mmvalue);
    vec3<double> product = _mm256_mul_pd(sub, sub);
    return sqrt(product.x + product.y + product.z);
}
/** @brief Calculate the inner product of two vectors */
template<>
inline double dot<double>(const vec3<double>& v1, const vec3<double>& v2)
{
    vec3<double>  product = _mm256_mul_pd(v1.mmvalue, v2.mmvalue);
    return product.x + product.y + product.z;
}
/** @brief Calculate the cross product of two vectors */
template<>
inline vec3<double> cross<double>(const vec3<double>& v1, const vec3<double>& v2)
{
    return vec3<double>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

}
#endif
#endif

