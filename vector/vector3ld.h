
#ifndef GENVECT_H
#define GENVECT_H
#include <math.h>
#include <iostream>
namespace SpaceH
{
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
        /** @brief Product by wise */
        inline vec3 operator*(const vec3& v) const
        {
            return vec3(x * v.x, y * v.y, z * v.z);
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
        /** @brief Calculate the norm */
        inline T norm() const
        {
            return sqrt(x * x + y * y + z * z);
        }
        /** @brief Calculate the norm */
        inline T norm2() const
        {
            return (x * x + y * y + z * z);
        }
        
        inline T max_component()
        {
            T max = (x > y ? x : y);
            return max > z ? max : z;
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
    /** @brief Calculate the inner product of two vectors */
    template<typename T>
    inline T dot(const vec3<T>& v1, const vec3<T>& v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
    /** @brief Calculate the cross product of two vectors */
    template<typename T>
    inline T cross(const vec3<T>& v1, const vec3<T>& v2)
    {
        return vec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }
    
    using vec3ld = vec3<long double>;
    using vec3d  = vec3<double>;
    using vec3f  = vec3<float>;
    using vec3i  = vec3<int>;
    using vec3c  = vec3<char>;
    using vec3b  = vec3<bool>;
}
#endif

