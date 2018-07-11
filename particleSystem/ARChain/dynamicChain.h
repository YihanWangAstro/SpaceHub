
#ifndef DYNAMICCHAIN_H
#define DYNAMICCHAIN_H
#include "dynamics.h"
/**
 *  @brief Class of dynamical variable.
 *
 */
template <typename States, typename Attributes>
class dynamicChain : public dynamics< States, Attributes >
{
public:
    using Base = dynamics< States, Attributes >;
    
    /**  @brief Position array const interface. Reference to cartesian.pos*/
    inline const typename Base::VectorArray& pos() const { return cartesian.pos; }
    
    /**  @brief Velocity array const interface. Reference to cartesian.vel*/
    inline const typename Base::VectorArray& vel() const { return cartesian.vel; }
    
    /**  @brief Position vector const interface. Reference to cartesian.pos[i]*/
    inline const typename Base::Vector& pos(size_t i) const { return cartesian.pos[i]; }
    
    /**  @brief Velocity vecotr const interface. Reference to cartesian.vel[i]*/
    inline const typename Base::Vector& vel(size_t i) const { return cartesian.vel[i]; }
    
    
    /** @brief Advance the position array with internal velocity array.
     *  @param stepSize The advance step size.
     */
    inline void advancePos(typename Base::Scalar stepSize)
    {
        Base::advancePos(stepSize);
    }
    
    /** @brief Advance the  velocity array with given acceleration array.
     *  @param stepSize The advance step size.
     *  @param acc      The acceleration array.
     */
    inline void advanceVel(typename Base::VectorArray& acc, typename Base::Scalar stepSize)
    {
        Base::advanceVel(acc, stepSize);
    }

    /** @brief Initialize variables with istream.
     *  @param input istream.
     */
    void initialize(std::istream& input)
    {
        Base::initialize(input);
    };
    

private:
    /** @brief Attributes of particles(const variables during evolution) */
    States cartesian;
};
#endif

