#ifndef REGULARSTATE_H
#define REGULARSTATE_H
#include "../particles.h"

/**
 *  @brief Class of dynamical system with regularization variables.
 *
 *  A simple extension of class dynamics in dynamicState.h. Used for regularization system. See detail
 *  in https://academic.oup.com/mnras/article/372/1/219/974304 .
 */
template <typename States, typename Attributes>
class reguDynamics : public particles<States, Attributes>
{
    using Base = dynamics<States, Attributes>;
    
public:
    /**  @brief Omega scalar const interface. Reference to state.time*/
    inline const typename Base::Scalar& omega() const { return this->stat.omega; }
    
    /**  @brief BindE scalar const interface. Reference to state.time*/
    inline const typename Base::Scalar& bindE() const { return this->stat.bindE; }
    
    /** @brief Advance the omega.
     *  @param dOmega Increament of omega.
     */
    inline void advanceOmega(typename Base::Scalar dOmega, typename Base::Scalar stepSize)
    {
#ifdef KAHAN_SUMMATION
        KahanAdvance(this->stat.omega, dOmega*stepSize, this->err.omega);
#else
        this->stat.omega += dOmega;
#endif
    }
    
    /** @brief Advance the bindE.
     *  @param dBindE Increament of bindE.
     */
    inline void advanceBindE(typename Base::Scalar dBindE, typename Base::Scalar stepSize)
    {
#ifdef KAHAN_SUMMATION
        KahanAdvance(this->stat.bindE, dBindE*stepSize, this->err.bindE);
#else
        this->stat.bindE += dBindE;
#endif
    }
    
    /** @brief Initialize variables with istream.
     *  @param input istream.
     */
    void initialize(std::istream& input)
    {
        Base::initialize(input);
        this->stat.omega = getCapitalOmega();
        this->stat.bindE = getTotalEnergy(this->attr.mass, this->stat.pos, this->stat.vel);
        
    }
    
    /** @brief Calculate the regularized variable Omega.
     *  @return The value of capital omega.
     */
    typename Base::Scalar getCapitalOmega()
    {
        return -getPotentialEnergy(this->attr.mass, this->stat.pos);
    }
};
#endif

