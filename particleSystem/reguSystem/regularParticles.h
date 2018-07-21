#ifndef REGULARSTATE_H
#define REGULARSTATE_H
#include "../../particles.h"

namespace SpaceH
{
/**
 *  @brief Class of dynamical system with regularization variables.
 *
 *  A simple extension of class dynamics in dynamicState.h. Used for regularization system. See detail
 *  in https://academic.oup.com/mnras/article/372/1/219/974304 .
 */
template<typename Dtype, size_t ArraySize, bool IsVelDep>
class ReguParticles : public Particles<Dtype, ArraySize, IsVelDep>
{
public:
    /* Typedef */
    using Base = Particles<Dtype, ArraySize, IsVelDep>;
    using typename Base::type;
    using Scalar       = typename type::Scalar;
    using Vector       = typename type::Vector;
    using VectorArray  = typename type::VectorArray;
    using ScalarBuffer = typename type::ScalarBuffer;
    /* Typedef */
    
    /*Template parameter check*/
    CHECK_POD(Dtype)
    /*Template parameter check*/
    
    /**  @brief Omega scalar const interface. Reference to state.time*/
    inline const Scalar& omega() const { return omega_; }
    
    /**  @brief BindE scalar const interface. Reference to state.time*/
    inline const Scalar& bindE() const { return bindE_; }
    
    /** @brief Calculate the regularized variable Omega.
     *  @return The value of capital omega.
     */
    inline Scalar getCapitalOmega() const
    {
        return -getPotentialEnergy(this->mass(), this->pos());
    }
    
    /** @brief Advance the Omega.
     *  @param velIndepAcc Velocity independent acceleration array.
     *  @param vel         Velocity array.
     *  @param stepSize    Time stepSize.
     */
    inline void advanceOmega(const VectorArray& velIndepAcc, const VectorArray& velocity, Scalar stepSize)
    {
        size_t size = this->particleNumber();
        Scalar dOmega = 0;
        
        for(size_t i = 0 ; i < size ; ++i)
            dOmega += (velIndepAcc[i] * velocity[i]) * (this->mass_[i]);
        
        SpaceH::advanceScalar(omega_, dOmega*stepSize);
    }
    
    /** @brief Advance the bindE.
     *  @param velDepAcc Velocity dependent acceleration array.
     *  @param vel       Velocity array.
     *  @param stepSize  Time stepSize.
     */
    inline void advanceBindE(const VectorArray& velDepAcc, const VectorArray& velocity, Scalar stepSize)
    {
        size_t particleNum = this->particleNumber();
        Scalar dBindE = 0;
        
        for(size_t i = 0 ; i < particleNum ; ++i)
            dBindE -= (velDepAcc[i] * velocity[i]) * (this->mass_[i]);
        
        SpaceH::advanceScalar(bindE_, dBindE*stepSize);
    }
    
    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, ReguParticles& partc)
    {
        is >> static_cast<Base&>(partc);
    
        partc.omega_ = partc.getCapitalOmega();
        partc.bindE_ = -getTotalEnergy(partc.mass(), partc.pos(), partc.vel());
        
        return is;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend size_t operator>>(const ScalarBuffer& data, ReguParticles& partc)
    {
        size_t loc = data >> static_cast<Base&>(partc);
        
        partc.omega_ = data[loc++];
        partc.bindE_ = data[loc++];
        
        return loc;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend size_t operator<<(ScalarBuffer& data, const ReguParticles& partc)
    {
        size_t loc = data << static_cast<const Base&>(partc);
        
        data.reserve(loc + 2);
        
        data.emplace_back(partc.omega_);
        data.emplace_back(partc.bindE_);
        
        return data.size();
    }
    
protected:
    Scalar omega_;
    
    Scalar bindE_;
};
}
#endif

