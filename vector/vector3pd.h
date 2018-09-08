
#ifndef GENVECT3PD_H
#define GENVECT3PD_H

#ifdef __AVX512F__
#pragma message("Using AVX512 on vectorpd")
#include <x86intrin.h>
#include "../kahanNumber.h"

namespace SpaceH
{
/** @brief Specilization of vector3d */
template<>
struct vec3<precise_d>
{
public:
    /* Typedef */
    using value_type = precise_d;
    /* Typedef */

    
    vec3<double> __attribute__((aligned(32))) real;
    vec3<double> __attribute__((aligned(32))) err;

    vec3() : {};
    vec3(double vx, double vy, double vz) : mmvalue(_mm512_set_pd(vx, 0, vy, 0, vz, 0, 0, 0.0)) {};
    vec3(__m256d v) : mmvalue(v) {};
    vec3(const vec3& v) : mmvalue(v.mmvalue) {};
    /** @brief Addition by wise */
    inline vec3 operator+(const vec3& v) const
    {
        return vec3(_mm512_add_pd(mmvalue, v.mmvalue));
    }
    /** @brief Subtraction by wise */
    inline vec3 operator-(const vec3& v) const
    {
        return vec3(_mm512_sub_pd(mmvalue, v.mmvalue));
    }
    /** @brief Product by wise */
    inline vec3 operator*(const vec3& v) const
    {
        return vec3(_mm512_mul_pd(mmvalue, v.mmvalue));
    }
    /** @brief Divition by wise */
    inline vec3 operator/(const vec3& v) const
    {
        return vec3(_mm512_div_pd(mmvalue, v.mmvalue));
    }
    /** @brief Add scalar by wise */
    inline vec3 operator+(const precise_d c) const
    {
        return vec3(_mm512_add_pd(mmvalue, _mm512_set1_pd(c)));
    }
    /** @brief Subtract scalar by wise */
    inline vec3 operator-(const precise_d c) const
    {
        return vec3(_mm512_sub_pd(mmvalue, _mm512_set1_pd(c)));
    }
    /** @brief Multiply scalar by wise */
    inline vec3 operator*(const precise_d c) const
    {
        return vec3(_mm512_mul_pd(mmvalue, _mm512_set1_pd(c)));
    }
    /** @brief Divide scalar by wise */
    inline vec3 operator/(const precise_d c) const
    {
        return vec3(_mm512_div_pd(mmvalue, _mm512_set1_pd(c)));
    }
    /** @brief Opposite vector */
    inline vec3 operator-() const
    {
        return vec3(-x, -y, -z);
    }
    /** @brief Absolute value by wise */
    inline vec3 abs() const
    {
        return vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z);
    }
    inline const vec3& operator+=(const vec3& v)
    {
        mmvalue = _mm512_add_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator-=(const vec3& v)
    {
        mmvalue = _mm512_sub_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator*=(const vec3& v)
    {
        mmvalue = _mm512_mul_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator/=(const vec3& v)
    {
        mmvalue = _mm512_div_pd(mmvalue, v.mmvalue);
        return *this;
    }
    inline const vec3& operator+=(const precise_d c)
    {
        mmvalue = _mm512_add_pd(mmvalue, _mm512_set1_pd(c));
        return *this;
    }
    inline const vec3& operator-=(const precise_d c)
    {
        mmvalue = _mm512_sub_pd(mmvalue, _mm512_set1_pd(c));
        return *this;
    }
    inline const vec3& operator*=(const precise_d c)
    {
        mmvalue = _mm512_mul_pd(mmvalue, _mm512_set1_pd(c));
        return *this;
    }
    inline const vec3& operator/=(const precise_d c)
    {
        mmvalue = _mm512_div_pd(mmvalue, _mm512_set1_pd(c));
        return *this;
    }
    inline const vec3& operator=(const vec3& v)
    {
        mmvalue = v.mmvalue;
        return *this;
    }
    /** @brief Calculate the norm */
    inline precise_d norm() const
    {
        vec3  product = _mm512_mul_pd(mmvalue, mmvalue);
        return sqrt(product.x + product.y + product.z);
    }
    /** @brief Calculate the norm */
    inline precise_d norm2() const
    {
        vec3  product = _mm512_mul_pd(mmvalue, mmvalue);
        return product.x + product.y + product.z;
    }

    inline precise_d max_component()
    {
        precise_d max = (x > y ? x : y);
        return max > z ? max : z;
    }

    /** @brief Calculate the inverse of the norm */
    inline precise_d reNorm() const
    {
        vec3  product = _mm512_mul_pd(mmvalue, mmvalue);
        return 1.0/sqrt(product.x + product.y + product.z);
    }
    inline void setZero()
    {
        mmvalue = _mm512_setzero_pd();
    }
    friend vec3 operator+(const precise_d c, const vec3& v)
    {
        return vec3(_mm512_add_pd(_mm512_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator-(const precise_d c, const vec3& v)
    {
        return vec3(_mm512_sub_pd(_mm512_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator*(const precise_d c, const vec3& v)
    {
        return vec3(_mm512_mul_pd(_mm512_set1_pd(c), v.mmvalue));
    }
    friend vec3 operator/(const precise_d c, const vec3& v)
    {
        return vec3(_mm512_div_pd(_mm512_set1_pd(c), v.mmvalue));
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
inline precise_d distance<precise_d>(const vec3<precise_d>& v1, const vec3<precise_d>& v2)
{
    __m512d sub = _mm512_sub_pd(v1.mmvalue, v2.mmvalue);
    vec3<precise_d> product = _mm512_mul_pd(sub, sub);
    return sqrt(product.x + product.y + product.z);
}
/** @brief Calculate the inner product of two vectors */
template<>
inline precise_d dot<precise_d>(const vec3<precise_d>& v1, const vec3<precise_d>& v2)
{
    vec3<precise_d>  product = _mm512_mul_pd(v1.mmvalue, v2.mmvalue);
    return product.x + product.y + product.z;
}
/** @brief Calculate the cross product of two vectors */
template<>
inline vec3<precise_d> cross<precise_d>(const vec3<precise_d>& v1, const vec3<precise_d>& v2)
{
    return vec3<precise_d>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

}
#endif
#endif

