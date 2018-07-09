
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:constIterator.h                                                                                            //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef BSITERATOR_H
#define BSITERATOR_H
#include "../libs.h"

/** @brief Bulirsch-Stoer extrapolation algorithm*/
template <typename ParticSys, typename Integrator>
class BSIterator
{
    
public:
    using ActiveScalarArray = typename ParticSys::ActiveScalarArray;
    using Scalar            = typename ParticSys::Scalar;
    
    template<typename T, size_t S>
    using Container = typename ParticSys::template Container<T, S>;

    /** @brief Constructor for initializing cost, nSteps, fmin and CC*/
    BSIterator();
    
    /** @brief Interface of ODE iterator
     *  @note  BSiterator will force use the internal mid-point integrator as the basic integrator.
     */
    Scalar iterate(ParticSys& particles, Scalar stepLength);
    
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
    
public:
    /** @brief The Maximum iteration depth*/
    static const size_t MaxDepth{8};
    
    /** @brief The local partical system used to iterate.*/
    ParticSys localSystem;
    
    /** @brief Extrapolation table.*/
    Container < ActiveScalarArray, (MaxDepth + 1)* (MaxDepth + 2) / 2 > extrapTab;
    
    /** @brief Extrapolation coefficient.*/
    Container < Scalar, (MaxDepth + 1)* (MaxDepth + 2) / 2 > CC;
    
    /** @brief Macro step length for different iteration depth.*/
    Container < Scalar, MaxDepth + 1 > macroStepLength;
    
    /** @brief The work(calculation quantities) per integration length of each iteration depth.*/
    Container < Scalar, MaxDepth + 1 > work;
    
    /** @brief The minimal coeffient of integration step estimation.*/
    Container < Scalar, MaxDepth + 1 > fmin;
    
    /** @brief The work(calculation quantities) of each iteration depth.*/
    Container < size_t, MaxDepth + 1 > cost;
    
    /** @brief Steps of integration of each iteration depth.*/
    Container < size_t, MaxDepth + 1 > nSteps;
    
    /** @brief Local absolute error*/
    Scalar absoluteError{1e-15};
    
    /** @brief Local relative error*/
    Scalar relativeError{1e-15};
    
    /** @brief Coefficient of new step length estimation. See "Numerical recipes" on page 926.*/
    Scalar s1{0.94};
    
    /** @brief Coefficient of new step length estimation*/
    Scalar s2{0.95};
    
    /** @brief Coefficient of new step length estimation*/
    Scalar s3{0.02};
    
    /** @brief Coefficient of new step length estimation*/
    Scalar s4{4.0};
    
    /** @brief Current iteraation depth.*/
    size_t iterDepth{7};
    
private:
    /** @brief Internal integrator */
    Integrator integrator;
    
    /** @brief Check the rejection criteria for current iteration.*/
    bool checkRejection(Scalar error, size_t k) const;
    
    /** @brief Extrapolate the kth row to the right end.*/
    void extrapolate(size_t k);
    
    /** @brief Calculate the error of the k row of the extrapolation table.*/
    Scalar getError(size_t k) const;
    
    /** @brief Calculate the new iteration integration step coefficient*/
    Scalar getTimeStepCoef(Scalar error, size_t k);
    
    /** @brief Calculate the new iteration integration step and new iteration depth for next iteration.*/
    Scalar prepareForNewIteration(size_t k, bool lastRejection);
};

/** @brief Constructor for initializing cost, nSteps, fmin and CC*/
template <typename ParticSys, typename Integrator>
BSIterator<ParticSys, Integrator>::BSIterator()
{
    fmin[0]      = s3;
    Scalar ratio = 1;
    nSteps[0]    = 1;
    cost[0]      = 1;
    CC[0]        = 1;

    for(size_t i = 1 ; i < MaxDepth + 1; ++i)
    {
        nSteps[i] = 2 * i;
        cost[i]   = cost[i - 1] + nSteps[i];
        
        fmin[i] = pow(s3, 1.0 / (2*i + 1));
        
        for(size_t j = 0 ; j < i; ++j)
        {
            ratio = (Scalar)nSteps[i] / (Scalar)nSteps[i - j - 1];
            CC[i * (i + 1) / 2 + j] = 1.0 / (ratio * ratio - 1);
        }
    }
}

/** @brief Interface of ODE iterator
 *  @param particles Particle system need iteration.
 *  @param integrator Basica integrator used to evolve, but here BS iterator will force use internal mid-point integrator.
 *  @param stepLength Macro integration step length.
 *  @return The next macro integration step length.
 */
template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::iterate(ParticSys& particles,
        Scalar stepLength)
{
    Scalar error  = 0;
    Scalar H      = stepLength;
    Scalar h      = stepLength;
    bool reject   = false;

    for(;;)
    {
        for(size_t k = 0 ; k <= iterDepth + 1 ; ++k)
        {
            localSystem = particles;
            h = H / nSteps[k];
            
            /*for(size_t i = 0 ; i < nSteps[k]; i++)
                integrator.integrate(localSystem,h);*/
            
            localSystem.drift(0.5 * h);
            for(size_t i = 1 ; i < nSteps[k]; i++)
            {
                localSystem.kick(h);
                localSystem.drift(h);
            }
            localSystem.kick(h);
            localSystem.drift(0.5 * h);
            
            
            extrapTab[k * (k + 1) / 2] << localSystem;//copyDataToExtrapTab;
            
            extrapolate(k);
            error = getError(k);
            macroStepLength[k] = H * getTimeStepCoef(error, k);
            work[k] = cost[k] / macroStepLength[k];
            //if(iterDepth!=7)
            //std::cout << cost[k] << " "<< k << " " << iterDepth << " " << error << " " << H << "\r\n";
            if(k == iterDepth - 1 || k == iterDepth || k == iterDepth + 1)
            {
                if(error <= 1)
                {
                    reject = false;
                    H = prepareForNewIteration(k, reject);
                    extrapTab[k * (k + 1) / 2 + k] >> localSystem;
                    particles = localSystem;
                    return H;
                }
                else if(checkRejection(error, k))
                {
                    reject = true;
                    H = prepareForNewIteration(k, reject);
                    break;
                }
            }
        }
    }
}

/** @brief Check the rejection criteria for current iteration.
 *  @param k The kth iteration depth.
 *  @return If current iteration is rejected.
 */
template <typename ParticSys, typename Integrator>
bool BSIterator<ParticSys, Integrator>::checkRejection(Scalar error, size_t k) const
{
    Scalar stepLengthRatio = 1;

    if(k == iterDepth - 1)
        stepLengthRatio = nSteps[k + 1] * nSteps[k + 2] / (Scalar)(nSteps[0] * nSteps[0]);
    else if(k == iterDepth)
        stepLengthRatio = nSteps[k + 1] / (Scalar)nSteps[0];
    else
        return true;//k == iterDepth+1 and error >1 reject directly

    return error > stepLengthRatio * stepLengthRatio;
}

/** @brief Extrapolate the k-th row to the right end.
 *  @param k The kth row of extrapolation table.
 */
template <typename ParticSys, typename Integrator>
void BSIterator<ParticSys, Integrator>::extrapolate(size_t k)
{
    if(k > 0)
    {
        size_t n    = k * (k + 1) / 2;
        size_t pn   = (k - 1) * k / 2;
        size_t size = extrapTab[n].size();

        for(size_t j = 0 ; j < k; ++j)
        {
            size_t now  = n + j;
            size_t next = now + 1;
            size_t last = pn + j;

            for(size_t i = 0 ; i < size ; ++i)
                extrapTab[next][i] = extrapTab[now][i] + (extrapTab[now][i] - extrapTab[last][i]) * CC[now];
        }
    }
    else
    {
        return;
    }
}
/** @brief Calculate the error of the k row of the extrapolation table.
 *  @param k The kth row of extrapolation table.
 */
template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::getError(size_t k) const
{
    if(k != 0)
    {
        size_t topk     = k * (k + 1) / 2 + k;
        size_t subk     = topk - 1;
        size_t size     = extrapTab[topk].size();
        Scalar maxError = 0;
        Scalar error    = 0;

        for(size_t i = 0 ; i < size; ++i)
        {
            error = abs(extrapTab[topk][i] - extrapTab[subk][i])
                    / (min(abs(extrapTab[topk][i]), abs(extrapTab[subk][i]) ) * this->relativeError + this->absoluteError);
            maxError = max(maxError, error);
        }

        return maxError;
    }
    else
    {
        return 1;
    }
}

/** @brief Calculate the new iteration integration step coefficient
 *  @param error The error of the k-th row of extrapolation table.
 *  @param order The order of the error in k-th row of extrapolation table.
 *  @return The new iteration integration step length coefficient.
 */
template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::getTimeStepCoef(Scalar error, size_t k)
{
    if(error != 0)
        return max(fmin[k] / s4, min( s1 * pow(s2 / error, 1.0 / (2 * k + 1)), 1.0 / fmin[k]));
    else
        return 1.0 / fmin[k];
}

/** @brief Calculate the new iteration integration step and new iteration depth for next iteration.
 *  @param k The k-th iteration depth.
 *  @param lastRejection The rejection status of last trail of iteration.
 *  @return The new iteration step length.
 */
template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::prepareForNewIteration(size_t k, bool lastRejection)
{
    Scalar newH = 0;

    if(k == iterDepth - 1)
    {
        iterDepth = min(MaxDepth - 1, max(2, k) );
        newH      = macroStepLength[iterDepth];

        if(work[k - 1] < 0.8 * work[k] && k >= 2)
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 1) );
            newH      = macroStepLength[iterDepth];
        }
        else if(work[k] < 0.9 * work[k - 1] || iterDepth <= 2)
        {
            iterDepth = min(MaxDepth - 1, max(2, k + 1) );
            newH      = macroStepLength[iterDepth - 1] * cost[iterDepth] / static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }
    else if(k == iterDepth)
    {
        iterDepth = min(MaxDepth - 1, max(2, k) );
        newH      = macroStepLength[iterDepth];

        if(work[k - 1] < 0.8 * work[k])
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 1) );
            newH      = macroStepLength[iterDepth];
        }
        else if(!lastRejection && work[k] < 0.9 * work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k + 1) );
            newH      = macroStepLength[iterDepth - 1] * cost[iterDepth] / static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }
    else if(k == iterDepth + 1)
    {
        iterDepth = min(MaxDepth - 1, max(2, k - 1) );
        newH      = macroStepLength[iterDepth];

        if(work[k - 2] < 0.8 * work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 2) );
            newH      = macroStepLength[iterDepth];
        }
        else if(!lastRejection && work[k] < 0.9 * work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k) );
            newH      = macroStepLength[iterDepth - 1] * cost[iterDepth] / static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }

    return newH;
}
#endif
