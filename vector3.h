
#ifndef GENVECT_H
#define GENVECT_H
#include <math.h>
#include <iostream>

/** @brief Self 3D vector class */
template<typename T>
struct vec3
{
public:
    /* Typedef */
    using value_type = T;
    /* Typedef */
    
    T  x;
    T  y;
    T  z;
    
    vec3() : x(0), y(0), z(0) {};
    vec3(T vx, T vy, T vz) : x(vx), y(vy), z(vz) {};
    vec3(const vec3& v) : x(v.x), y(v.y), z(v.z) {};
    /** @brief Addition by wise */
    inline vec3 operator+(const vec3& v) const
    {
        return vec3(x + v.x, y + v.y, z + v.z);
    }
    /** @brief Subtraction by wise */
    inline vec3 operator-(const vec3& v) const
    {
        return vec3(x - v.x, y - v.y, z - v.z);
    }
    /** @brief Divition by wise */
    inline vec3 operator/(const vec3& v) const
    {
        return vec3(x / v.x, y / v.y, z / v.z);
    }
    /** @brief Add scalar by wise */
    inline vec3 operator+(const T c) const
    {
        return vec3(x + c, y + c, z + c);
    }
    /** @brief Subtract scalar by wise */
    inline vec3 operator-(const T c) const
    {
        return vec3(x - c, y - c, z - c);
    }
    /** @brief Multiply scalar by wise */
    inline vec3 operator*(const T c) const
    {
        return vec3(x * c, y * c, z * c);
    }
    /** @brief Divide scalar by wise */
    inline vec3 operator/(const T c) const
    {
        return vec3(x / c, y / c, z / c);
    }
    /** @brief Opposite vector */
    inline vec3 operator-() const
    {
        return vec3(-x, -y, -z);
    }
    /** @brief Cross product */
    inline vec3 operator^(const vec3& v) const
    {
        return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    /** @brief Absolute value by wise */
    inline vec3 abs() const
    {
        return vec3(x > -x ? x : -x, y > -y ? y : -y, z > -z ? z : -z);
    }
    inline const vec3& operator+=(const vec3& v)
    {
        x += v.x, y += v.y, z += v.z;
        return *this;
    }
    inline const vec3& operator-=(const vec3& v)
    {
        x -= v.x, y -= v.y, z -= v.z;
        return *this;
    }
    inline const vec3& operator/=(const vec3& v)
    {
        x /= v.x, y /= v.y, z /= v.z;
        return *this;
    }
    inline const vec3& operator+=(const T c)
    {
        x += c, y += c, z += c;
        return *this;
    }
    inline const vec3& operator-=(const T c)
    {
        x -= c, y -= c, z -= c;
        return *this;
    }
    inline const vec3& operator*=(const T c)
    {
        x *= c, y *= c, z *= c;
        return *this;
    }
    inline const vec3& operator/=(const T c)
    {
        x /= c, y /= c, z /= c;
        return *this;
    }
    inline const vec3& operator=(const vec3& v)
    {
        x = v.x, y = v.y, z = v.z;
        return *this;
    }
    /** @brief Inner product */
    inline T operator*(const vec3& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    /** @brief Calculate the norm */
    inline T norm() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    /** @brief Calculate the inverse of the norm */
    inline T reNorm() const
    {
        return 1.0 / sqrt(x * x + y * y + z * z);
    }
    inline void setZero()
    {
        x = y = z = 0;
    }
    friend vec3 operator+(const T c, const vec3& v)
    {
        return vec3(v.x + c, v.y + c, v.z + c);
    }
    friend vec3 operator-(const T c, const vec3& v)
    {
        return vec3(c - v.x, c - v.y, c - v.z);
    }
    friend vec3 operator*(const T c, const vec3& v)
    {
        return vec3(v.x * c, v.y * c, v.z * c);
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
template<typename T>
inline T distance(const vec3<T>& v1, const vec3<T>& v2)
{
    return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z));
}
typedef vec3<double> vec3d;
typedef vec3<float>  vec3f;
typedef vec3<int>    vec3i;
typedef vec3<char>   vec3c;
#endif

