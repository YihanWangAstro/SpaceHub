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

constexpr double INV_C  = 1/C;
constexpr double INV_C2 = INV_C*INV_C;
constexpr double INV_C3 = INV_C2*INV_C;
constexpr double INV_C4 = INV_C3*INV_C;
constexpr double INV_C5 = INV_C4*INV_C;

template<typename Scalar>
class PN1th
{
private:
    typedef vec3<Scalar> Vector;
public:
    void operator()(Scalar m1, Scalar m2, Vector& p1, Vector& p2, Vector& v1, Vector& v2, Vector& acc1, Vector& acc2)
    {
        Vector dr(p1 - p2);
        Vector dv(v1 - v2);
        Scalar inv_r  = dr.reNorm();
        Scalar inv_r2 = inv_r*inv_r;
        Vector n12(dr*inv_r);
        Scalar nv1    = n12*v1;
        Scalar nv2    = n12*v2;
        Scalar A1 = 0.0, B1 = 0.0;
        A1 = ( (5.0*m1 + 4.0*m2)*inv_r + (1.5*nv2*nv2 - v1*v1 + 4.0*(v1*v2) - 2.0*(v2*v2)))*m2*inv_r2*INV_C2;
        B1 = (4.0*nv1 - 3.0*nv2)*m2*inv_r2*INV_C2;
        
        acc1 += n12*A1 + dv*B1;
        
        A1 = ( (5.0*m2 + 4.0*m1)*inv_r + (1.5*nv1*nv1 - v2*v2 + 4.0*(v2*v1) - 2.0*(v1*v1)))*m1*inv_r2*INV_C2;
        B1 = (-4.0*nv2 + 3.0*nv1)*m1*inv_r2*INV_C2;
        
        acc2 -= n12*A1 + dv*B1;
    }
};

template<typename Scalar>
class Newtonian
{
public:
    void operator()(){}
};
#endif
