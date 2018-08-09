
#ifndef PREDICTORCORRECTOR_H
#define PREDICTORCORRECTOR_H
/** @brief Basic predictor corrector iterator
 *
 */
template <typename ParticSys, typename Integrator>
class PCIterator
{
public:
    /* Typedef */
    typedef typename ParticSys::Scalar     Scalar;
    typedef typename ParticSys::PlainArray PlainArray;
    /* Typedef */
    
    /*Template parameter check*/
    CHECK_TYPE(ParticSys, Integrator)
    /*Template parameter check*/
    
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
    
    /** @brief Local absolute error*/
    Scalar absoluteError{1e-15};
    
    /** @brief Local relative error*/
    Scalar relativeError{1e-15};
};


#endif
