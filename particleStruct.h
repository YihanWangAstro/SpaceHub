
#ifndef PARTICLESTRUCT_H
#define PARTICLESTRUCT_H

#include "protoType.h"

namespace particle
{
    /** @brief State of particles(evolved variables during evolution) */
    template<typename Dtype, size_t ArraySize>
    struct States : public type::ProtoType<Dtype, ArraySize>
    {
        using Base = type::ProtoType<Dtype, ArraySize>;
        
        /** @brief Position array of the particles. Element is 3D vector.*/
        typename Base::VectorArray pos;
        
        /** @brief Velocity array of the particles. Element is 3D vector.*/
        typename Base::VectorArray vel;
        
        /** @brief The physical time of the dynamic system*/
        typename Base::Scalar time;
        
        /** @brief The size of the particle number*/
        inline size_t particleNumber()
        {
            return pos.size();
        }
        
    };
    
    /** @brief Attributes of particles(const variables during evolution) */
    template<typename Dtype, size_t ArraySize>
    struct Attributes : public type::ProtoType<Dtype, ArraySize>
    {
        using Base = type::ProtoType<Dtype, ArraySize>;
        
        /** @brief Mass array of the particles. Element is Scalar.*/
        typename Base::ScalarArray mass;
        
        /** @brief Radius array of the particles. Element is Scalar.*/
        typename Base::ScalarArray radius;
        
        /** @brief Type Array of the particles. Element is int.*/
        typename Base::IntArray type;
        
        /** @brief Id Array of the particles. Element is int.*/
        typename Base::IntArray idn;
        
        /** @brief The total mass of the system*/
        typename Base::Scalar totalMass;
        
        /** @brief The size of the particle number*/
        inline size_t particleNumber()
        {
            return mass.size();
        }
        
    };
    
    /** @brief State of particles(evolved variables during evolution) */
    template<typename Dtype, size_t ArraySize>
    struct ReguState : public States<Dtype, ArraySize>
    {
        using Base = States<Dtype, ArraySize>;
        
        /** @brief The regularization scalar omega*/
        typename Base::Scalar omega;
        
        /** @brief The regularization scalar binding Energy*/
        typename Base::Scalar bindE;
        
        /** @brief The number of scalar*/
        constexpr static size_t scalarNumber{ArraySize*6+3};
    };
    
    /** @brief State of particles(evolved variables during evolution) */
    template<typename Dtype, size_t ArraySize>
    struct GAR : public States<Dtype, ArraySize>
    {
        using Base = States<Dtype, ArraySize>;
        
        /** @brief Auxiliary velocity array of the particles. Element is 3D vector.*/
        typename Base::VectorArray auxiVel;
    };
    
    /** @brief State of particles(evolved variables during evolution) */
    template<typename Dtype, size_t ArraySize>
    struct ReguGAR : public ReguState<Dtype, ArraySize>
    {
        using Base = ReguState<Dtype, ArraySize>;
        
        /** @brief Auxiliary velocity array of the particles. Element is 3D vector.*/
        typename Base::VectorArray auxiVel;
    };
}

#endif

