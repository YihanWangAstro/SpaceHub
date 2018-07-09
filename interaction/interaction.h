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
};

    
template<typename InVelIndep, typename InVelDep, typename OutVelIndep, typename OutVelDep>
class Interaction
{
public:
    /* Typedef */
    template<typename T, size_t S>
    using Container   = typename InVelIndep::template Container<T, S>;
    
    using Scalar      = typename InVelIndep::Scalar;
    
    using Vector      = typename InVelIndep::Vector;
    
    using VectorArray = typename InVelIndep::VectorArray;
    
    using ScalarArray = typename InVelIndep::ScalarArray;
    
    using IntArray    = typename InVelIndep::IntArray;
    
    using SizeArray   = typename InVelIndep::SizeArray;
    
    constexpr static size_t arraySize{InVelIndep::arraySize};
    /* Typedef */
    
    void calInnerVelIndepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        iVelIndepAcc(mass, pos, vel);
        iVelIndepAcc.addTotal(acc);
    }
    
    void calInnerVelDepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        iVelDepAcc(mass, pos, vel);
        iVelDepAcc.addTotal(acc);
    }
    
    void calOuterVelIndepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        oVelIndepAcc(mass, pos, vel);
        oVelIndepAcc.addTotal(acc);
    }
    
    void calOuterVelDepAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel )
    {
        oVelDepAcc(mass, pos, vel);
        oVelDepAcc.addTotal(acc);
    }
    
    const VectorArray& totalAcc()
    {
        return acc;
    }
    
    const Vector& totalAcc(size_t i)
    {
        return acc[i];
    }
    
    const VectorArray& innerVelIndepAcc()
    {
        return iVelIndepAcc.acc();
    }
    
    const VectorArray& innerVelDepAcc()
    {
        return iVelDepAcc.acc();
    }
    
    const VectorArray& outerVelIndepAcc()
    {
        return oVelIndepAcc.acc();
    }
    
    const VectorArray& outerVelDepAcc()
    {
        return oVelDepAcc.acc();
    }
    
    const Vector& innerVelIndepAcc(size_t i)
    {
        return iVelIndepAcc.acc(i);
    }
    
    const Vector& innerVelDepAcc(size_t i)
    {
        return iVelDepAcc.acc(i);
    }
    
    const Vector& outerVelIndepAcc(size_t i)
    {
        return oVelIndepAcc.acc(i);
    }
    
    const Vector& outerVelDepAcc(size_t i)
    {
        return oVelDepAcc.acc(i);
    }
    
    void zeroTotalAcc()
    {
        for(size_t i = 0 ; i < acc.size() ; ++i)
            acc[i].setZero();
    }
    
private:
    VectorArray acc;
    
    InVelIndep  iVelIndepAcc;
    
    InVelDep    iVelDepAcc;
    
    OutVelIndep oVelIndepAcc;
    
    OutVelDep   oVelDepAcc;
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
    
    void operator()(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
    {
        Vector dr(0.0, 0.0, 0.0);
        Scalar inv_r  = 1;
        Scalar inv_r3 = 1;
        size_t size = mass.size();
        
        memset(&(newtonAcc[0]), 0, sizeof(Vector)*size);
        
        for(size_t i = 0 ; i < size ; ++i)
        {
            for(size_t j = i + 1 ; j < size ; ++j)
            {
                dr     = pos[j] - pos[i];
                inv_r  = dr.reNorm();
                inv_r3 = inv_r * inv_r * inv_r;
                newtonAcc[i] += dr * (inv_r3 * mass[j]);
                newtonAcc[j] -= dr * (inv_r3 * mass[i]);
            }
        }
    }
        
    void addTotal(VectorArray& totAcc)
    {
        for(size_t i = 0 ; i < totAcc.size() ;++i)
        {
            totAcc[i] += newtonAcc[i];
        }
    }
        
    const VectorArray& acc()
    {
        return newtonAcc;
    }
    
    const Vector& acc(size_t i)
    {
        return newtonAcc[i];
    }
private:
    VectorArray newtonAcc;
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
