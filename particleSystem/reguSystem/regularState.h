#ifndef REGULARSTATE_H
#define REGULARSTATE_H
#include "../../particles.h"

/**
 *  @brief Class of dynamical system with regularization variables.
 *
 *  A simple extension of class dynamics in dynamicState.h. Used for regularization system. See detail
 *  in https://academic.oup.com/mnras/article/372/1/219/974304 .
 */
template<typename Dtype, size_t ArraySize>
class ReguParticles : public particles<Dtype, ArraySize>
{
public:
    /* Typedef */
    using Base = particles<Dtype, ArraySize>;
    
    using typename Base::type;
    
    template<typename T, size_t S>
    using Container = typename type::template Container<T, S>;
    
    using Scalar = typename type::Scalar;
    
    using Vector = typename type::Vector;
    
    using VectorArray = typename type::VectorArray;
    /* Typedef */

    constexpr static size_t activeScalar{6*type::arraySize + 3};
    
    using ActiveScalarArray = Container<Scalar, activeScalar>;
    

    /**  @brief Omega scalar const interface. Reference to state.time*/
    inline const Scalar& omega() const { return omega_; }
    
    /**  @brief BindE scalar const interface. Reference to state.time*/
    inline const Scalar& bindE() const { return bindE_; }
    
    /** @brief Advance the Omega.
     *  @param velIndepAcc Velocity independent acceleration array.
     *  @param vel         Velocity array.
     *  @param stepSize    Time stepSize.
     */
    inline void advanceOmega(const VectorArray& velIndepAcc, const VectorArray& vel, Scalar stepSize)
    {
        size_t size = this->particleNumber();
        Scalar dOmega = 0;
        
        for(size_t i = 0 ; i < size ; ++i)
            dOmega += (velIndepAcc[i] * vel[i]) * (this->mass_[i]);
        
        advanceScalar(omega_, dOmega*stepSize);
    }
    
    /** @brief Advance the bindE.
     *  @param velDepAcc Velocity dependent acceleration array.
     *  @param vel       Velocity array.
     *  @param stepSize  Time stepSize.
     */
    inline void advanceBindE(const VectorArray& velDepAcc, const VectorArray& vel, Scalar stepSize)
    {
        size_t particleNum = this->partc.particleNumber();
        Scalar dBindE = 0;
        
        for(size_t i = 0 ; i < particleNum ; ++i)
            dBindE -= (velDepAcc[i] * vel[i]) * (this->mass_[i]);
        
        advanceScalar(bindE_, dBindE*stepSize);
    }
    
    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, ReguParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        is >> partc.time_;
        
        partc.totalMass_ = 0;
        
        for(size_t i = 0 ; i < particleNum ; ++i)
        {
            is >> partc.idn_[i]
            >> partc.type_[i]
            >> partc.mass_[i]
            >> partc.radius_[i]
            >> partc.pos_[i]
            >> partc.vel_[i];
            
            partc.totalMass_ += partc.mass_[i];
        }
        
        partc.omega_ = partc.getCapitalOmega();
        partc.bindE_ = -getTotalEnergy(partc.mass(), partc.pos(), partc.vel());
        return is;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend ActiveScalarArray& operator>>(ActiveScalarArray& data, ReguParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            partc.pos_[i].x = data[d_loc++];
            partc.pos_[i].y = data[d_loc++];
            partc.pos_[i].z = data[d_loc++];
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            partc.vel_[i].x = data[d_loc++];
            partc.vel_[i].y = data[d_loc++];
            partc.vel_[i].z = data[d_loc++];
        }
        
        partc.time_  = data[d_loc++];
        partc.omega_ = data[d_loc++];
        partc.bindE_ = data[d_loc];
        
        return data;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend ActiveScalarArray& operator<<(ActiveScalarArray& data, const ReguParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            data[d_loc++] = partc.pos_[i].x;
            data[d_loc++] = partc.pos_[i].y;
            data[d_loc++] = partc.pos_[i].z;
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            data[d_loc++] = partc.vel_[i].x;
            data[d_loc++] = partc.vel_[i].y;
            data[d_loc++] = partc.vel_[i].z;
        }
        
        data[d_loc++] = partc.time_;
        data[d_loc++] = partc.omega_;
        data[d_loc]   = partc.bindE_;
        
        return data;
    }
    
protected:
    Scalar omega_;
    
    Scalar bindE_;
    
    /** @brief Calculate the regularized variable Omega.
     *  @return The value of capital omega.
     */
    Scalar getCapitalOmega()
    {
        return -getPotentialEnergy(this->mass(), this->pos());
    }
};
#endif

