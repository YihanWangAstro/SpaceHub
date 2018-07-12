//
//  kahanNumber.h
//  SpaceHub
//
//  Created by 王艺涵 on 7/11/18.
//  Copyright © 2018 YihanWang. All rights reserved.
//

#ifndef KAHANNUMBER_h
#define KAHANNUMBER_h
namespace SpaceH
{
/** @brief Kahan number */
template<typename T>
struct kahan
{
public:
    using value_type = T;
    T real, err;
    kahan() {};
    kahan(T r) : real(r) {};
    kahan(const kahan& k) : real(k.real), err(k.err){};
    
    inline const kahan& operator=( const kahan& hs)
    {
        real = hs.real, err = 0;
        return *this;
    }
    
    inline operator T()
    {
        return real;
    }
    
    inline void zeroErr()
    {
        err = 0;
    }
    
    /** @brief Addition by wise */
    friend inline kahan operator+(const kahan& lhs, const kahan& rhs)
    {
        return kahan(lhs.real + rhs.real);
    }
    /** @brief Subtraction by wise */
    friend inline kahan operator-(const kahan& lhs, const kahan& rhs)
    {
        return kahan(lhs.real - rhs.real);
    }
    /** @brief Multiply by wise */
    friend inline kahan operator*(const kahan& lhs, const kahan& rhs)
    {
        return kahan(lhs.real * rhs.real);
    }
    /** @brief Divition by wise */
    friend inline kahan operator/(const kahan& lhs, const kahan& rhs)
    {
        return kahan(lhs.real / rhs.real);
    }
    /** @brief Opposite vector */
    friend inline kahan operator-(const kahan& hs)
    {
        return kahan(-hs.real);
    }
    
    friend inline const kahan& operator+=(kahan& lhs, const kahan& rhs)
    {
        T add = rhs.real - lhs.err;
        T sum = lhs.real + add;
        lhs.err = (sum - lhs.real) - add;
        lhs.real = sum;
        return lhs;
    }
    friend inline const kahan& operator-=(kahan& lhs, const kahan& rhs)
    {
        T add = -rhs.real - lhs.err;
        T sum = lhs.real + add;
        lhs.err = (sum - lhs.real) - add;
        lhs.real = sum;
        return lhs;
    }
    friend inline const kahan& operator/=(kahan& lhs, const kahan& rhs)
    {
        lhs.real /= rhs.real;
        return lhs;
    }
    
    friend inline const kahan& operator*=(kahan& lhs, const kahan& rhs)
    {
        lhs.real *= rhs.real;
        return lhs;
    }
    
    friend inline bool operator==(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real == rhs.real;
    }
    
    friend inline bool operator!=(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real != rhs.real;
    }
    
    friend inline bool operator>=(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real >= rhs.real;
    }
    
    friend inline bool operator<=(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real <= rhs.real;
    }
    
    friend inline bool operator>(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real > rhs.real;
    }
    
    friend inline bool operator<(const kahan& lhs, const kahan& rhs)
    {
        return lhs.real < rhs.real;
    }
    
    /** @brief Output to ostream */
    friend std::ostream& operator<<(std::ostream& output, const kahan& v)
    {
        output << v.real;
        return output;
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, kahan& v)
    {
        input >> v.real;
        return input;
    }
};

}
#endif /* kahanNumber_h */
