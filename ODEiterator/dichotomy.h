////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:constIterator.h                                                                                            //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DICHOTOMY_H
#define DICHOTOMY_H
/** @brief Dichotomy iterator
 *
 *  Dichotomy iterator, use dichotomy to adjust the step size.
 */
template <typename ParticSys, typename Integrator>
class dichoIterator
{
public:
    typedef typename ParticSys::Scalar     Scalar;
    typedef typename ParticSys::PlainArray PlainArray;
    /** @brief interface to iterate particle system for one step
     *  @param particles  Particle system needs evolution.
     *  @param integrator Integrator to integrate the particle system.
     *  @param stepLength Macro step length for iteration(Here, the step length of the integrator).
     *  @return step length for next iteration.
     */
    Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength);
    
    /** @brief Set the local relative error*/
    void setRelativeError(Scalar relError)
    {
        relativeError = relError;
    }
    
    /** @brief Set the local absolute error*/
    void setAbsoluteError(Scalar absError)
    {
        absoluteError = absError;
    }
    
private:
    /** @brief The local partical system used to iterate.*/
    ParticSys localSystem1;
    
    /** @brief The local partical system used to iterate.*/
    ParticSys localSystem2;
    
    /** @brief The local partical system used to iterate.*/
    ParticSys localSystem3;
    
    /** @brief Local absolute error*/
    Scalar absoluteError{1e-15};
    
    /** @brief Local relative error*/
    Scalar relativeError{1e-15};
    
    /** @brief Calculate the error of two integration results*/
    Scalar getError(PlainArray& array1, PlainArray& array2) const;
};

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar dichoIterator<ParticSys, Integrator>::iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
{
    localSystem1 = particles;
    localSystem2 = particles;
    
    integrator.integrate(localSystem1, stepLength);
    stepLength *= 0.5;
    integrator.integrate(localSystem2, stepLength);
    localSystem3 = localSystem2;
    integrator.integrate(localSystem2, stepLength);
    
    Scalar err = getError(localSystem1.array(), localSystem2.array());
    
    if(err >= 1)
    {
        for(; getError(localSystem1.array(), localSystem2.array()) > 1 ; )
        {
            localSystem1 = localSystem3;
            localSystem2 = particles;
            stepLength *= 0.5
            integrator.integrate(localSystem2, stepLength);
            localSystem3 = localSystem2;
            integrator.integrate(localSystem2, stepLength);
        }
    }
    else
    {
        for(; getError(localSystem1.array(), localSystem2.array()) > 1 ; )
        {
            localSystem3 = localSystem2 = localSystem1;
            integrator.integrate(localSystem2, stepLength);
            localSystem1 = particles;
            stepLength  *= 2;
            integrator.integrate(localSystem1, stepLength);
        }
    }
    
    
    return stepLength;
}

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar dichoIterator<ParticSys, Integrator>::getError(PlainArray& array1, PlainArray& array2)
{
    size_t size = array1.size();
    Scalar maxError = 0;
    Scalar error    = 0;
    
    for(size_t i = 0 ; i < size; ++i)
    {
        error = abs(array1[i] - array2[i]) / (min(abs(array1[i]), array2[i]) ) * this->relativeError + this->absoluteError);
        maxError = max(maxError, error);
    }
    
    return maxError;
}
#endif
