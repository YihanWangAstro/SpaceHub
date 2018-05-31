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
#include "../macros.h"
#include "../vector3.h"
namespace interact
{
    
constexpr double INV_C  = 1 / C;
constexpr double INV_C2 = INV_C * INV_C;
constexpr double INV_C3 = INV_C2 * INV_C;
constexpr double INV_C4 = INV_C3 * INV_C;
constexpr double INV_C5 = INV_C4 * INV_C;
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

/** @brief Marker of None velocity dependent force functor(c++ std11)*/
template<typename Scalar>
class Newtonian
{
public:
    void operator()() {}
};
    
}
#endif
