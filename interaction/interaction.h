////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:size_tercation.h                                                                                              //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef POSTNEWTONIAN_H
#define POSTNEWTONIAN_H
#include "../protoType.h"
namespace interact
{
    
constexpr double INV_C  = 1 / C;
constexpr double INV_C2 = INV_C * INV_C;
constexpr double INV_C3 = INV_C2 * INV_C;
constexpr double INV_C4 = INV_C3 * INV_C;
constexpr double INV_C5 = INV_C4 * INV_C;
    
template<typename Dtype, size_t ArraySize>
struct Empty : public type::ProtoType<Dtype, ArraySize>
{
    using Base = type::ProtoType<Dtype, ArraySize>;
    
    /* Typedef */
    template<typename T, size_t S>
    using Container   = typename Base::template Container<T, S>;
    
    using Scalar      = typename Base::Scalar;
    
    using Vector      = typename Base::Vector;
    
    using VectorArray = typename Base::VectorArray;
    
    using ScalarArray = typename Base::ScalarArray;
    
    using IntArray    = typename Base::IntArray;
    
    using SizeArray   = typename Base::SizeArray;
    
    constexpr static size_t arraySize{Base::arraySize};
    /* Typedef */
    
public:
    void operator()(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel) {}
    
    void addTotal(VectorArray& acc) {}
    
    const VectorArray& acc() {std::cout << "unaccessible data!" << "\r\n"; }
    
    const Vector& acc(size_t i) {std::cout << "unaccessible data!" << "\r\n"; }
    
    constexpr static bool isVelDep{false};
};

    
template<typename VelIndep, typename VelDep, typename ExtVelIndep, typename ExtVelDep>
class Interaction
{
public:
    /* Typedef */
    template<typename T, size_t S>
    using Container   = typename VelIndep::template Container<T, S>;
    
    using Scalar      = typename VelIndep::Scalar;
    
    using Vector      = typename VelIndep::Vector;
    
    using VectorArray = typename VelIndep::VectorArray;
    
    using ScalarArray = typename VelIndep::ScalarArray;
    
    using IntArray    = typename VelIndep::IntArray;
    
    using SizeArray   = typename VelIndep::SizeArray;
    
    constexpr static size_t arraySize{VelIndep::arraySize};
    /* Typedef */
    
    constexpr static bool isVelDep{ VelDep::isVelDep || ExtVelDep::isVelDep };
    
    void calcuVelIndepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        vel_indep_(mass, pos, vel);
    }
    
    void calcuVelDepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        vel_dep_(mass, pos, vel);
    }
    
    void calcuExtVelIndepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        ext_vel_indep_(mass, pos, vel);
    }
    
    void calcuExtVelDepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel )
    {
        ext_vel_dep_(mass, pos, vel);
    }
    
    const VectorArray& totalAcc()
    {
        return acc_;
    }
    
    const Vector& totalAcc(size_t i)
    {
        return acc_[i];
    }
    
    const VectorArray& velIndepAcc()
    {
        return vel_indep_.acc();
    }
    
    const VectorArray& velDepAcc()
    {
        return vel_dep_.acc();
    }
    
    const VectorArray& extVelIndepAcc()
    {
        return ext_vel_indep_.acc();
    }
    
    const VectorArray& extVelDepAcc()
    {
        return ext_vel_dep_.acc();
    }
    
    const Vector& velIndepAcc(size_t i)
    {
        return vel_indep_.acc(i);
    }
    
    const Vector& velDepAcc(size_t i)
    {
        return vel_dep_.acc(i);
    }
    
    const Vector& extVelIndepAcc(size_t i)
    {
        return ext_vel_indep_.acc(i);
    }
    
    const Vector& extVelDepAcc(size_t i)
    {
        return ext_vel_dep_.acc(i);
    }
    
    void zeroTotalAcc()
    {
        for(size_t i = 0 ; i < acc_.size() ; ++i)
            acc_[i].setZero();
    }
    
    void calcuTotalAcc()
    {
        acc_ = vel_indep_.acc();
        
        vel_dep_.addTotal(acc_);
        
        ext_vel_indep_.addTotal(acc_);
        
        ext_vel_dep_.addTotal(acc_);
    }
private:
    VectorArray acc_;
    
    VelIndep    vel_indep_;
    
    VelDep      vel_dep_;
    
    ExtVelIndep ext_vel_indep_;
    
    ExtVelDep   ext_vel_dep_;
};

template<typename Dtype, size_t ArraySize>
class NewtonGrav : public type::ProtoType<Dtype, ArraySize>
{
public:
    /* Typedef */
    using Base = type::ProtoType<Dtype, ArraySize>;
    
    template<typename T, size_t S>
    using Container   = typename Base::template Container<T, S>;
    
    using Scalar      = typename Base::Scalar;
    
    using Vector      = typename Base::Vector;
    
    using VectorArray = typename Base::VectorArray;
    
    using ScalarArray = typename Base::ScalarArray;
    
    using IntArray    = typename Base::IntArray;
    
    using SizeArray   = typename Base::SizeArray;
    
    constexpr static size_t arraySize{Base::arraySize};
    /* Typedef */
    
    constexpr static bool isVelDep{false};
    
    void operator()(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        Vector dr(0.0, 0.0, 0.0);
        Scalar inv_r  = 1;
        Scalar inv_r3 = 1;
        size_t size = mass.size();
        
        memset(&(newton_acc_[0]), 0, sizeof(Vector)*size);
        
        for(size_t i = 0 ; i < size ; ++i)
        {
            for(size_t j = i + 1 ; j < size ; ++j)
            {
                dr     = pos[j] - pos[i];
                inv_r  = dr.reNorm();
                inv_r3 = inv_r * inv_r * inv_r;
                newton_acc_[i] += dr * (inv_r3 * mass[j]);
                newton_acc_[j] -= dr * (inv_r3 * mass[i]);
            }
        }
    }
        
    void addTotal(VectorArray& totAcc)
    {
        for(size_t i = 0 ; i < totAcc.size() ;++i)
        {
            totAcc[i] += newton_acc_[i];
        }
    }
        
    const VectorArray& acc() { return newton_acc_; }
    
    const Vector& acc(size_t i) { return newton_acc_[i]; }
private:
    VectorArray newton_acc_;
};
    
/** @brief Post newtonian pair interaction functor(c++ std11)*/
template<typename Scalar>
class PostNewtonian
{
private:
    typedef vec3<Scalar> Vector;
public:
    /** @brief Update the velocity dependent acceleration of particle 1 and 2.
     *  @param m1   Mass of particle 1.
     *  @param m2   Mass of particle 2.
     *  @param dr   Relative position pos1 - pos2.
     *  @param dv   Relative velocity vel1 - vel2.
     *  @param v1   Velocity of particle 1.
     *  @param v2   Velocity of particle 2.
     *  @param acc1 Velocity dependent acceleration of particle 1 as return value.
     *  @param acc2 Velocity dependent acceleration of particle 1 as return value.
     */
    void operator()(Scalar m1, Scalar m2, Vector& dr, Vector& dv, Vector& v1, Vector& v2, Vector& acc1, Vector& acc2)
    {
        Scalar inv_r  = dr.reNorm();
        Scalar inv_r2 = inv_r * inv_r;
        Vector n12(dr * inv_r);
        Scalar nv1    = n12 * v1;
        Scalar nv2    = n12 * v2;
        Scalar A1 = 0.0, B1 = 0.0;
        A1 = ( (5.0 * m1 + 4.0 * m2) * inv_r + (1.5 * nv2 * nv2 - v1 * v1 + 4.0 * (v1 * v2) - 2.0 *
                                                (v2 * v2))) * m2 * inv_r2 * INV_C2;
        B1 = (4.0 * nv1 - 3.0 * nv2) * m2 * inv_r2 * INV_C2;
        acc1 += n12 * A1 + dv * B1;
        A1 = ( (5.0 * m2 + 4.0 * m1) * inv_r + (1.5 * nv1 * nv1 - v2 * v2 + 4.0 * (v2 * v1) - 2.0 *
                                                (v1 * v1))) * m1 * inv_r2 * INV_C2;
        B1 = (-4.0 * nv2 + 3.0 * nv1) * m1 * inv_r2 * INV_C2;
        acc2 -= n12 * A1 + dv * B1;
    }
};
    
}

#endif
