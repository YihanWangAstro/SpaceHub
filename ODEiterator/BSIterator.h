
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


template <typename ParticSys, typename Integrator>
class BSIterator
{
    //////////////////////////////////Interface/////////////////////////////////////
public:
    typedef typename ParticSys::Scalar Scalar;
    template <typename Scalar, size_t N> using scalarArray = std::array<Scalar,N>;
    
    BSIterator();
    Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength);
    void setRelativeError(Scalar relError){ relativeError = relError; }
    void setAbsoluteError(Scalar absError){ absoluteError = absError; }
    ///////////////////////////////Member variables/////////////////////////////////
private:
    static const size_t MaxDepth{8};
    ParticSys        localSystem;
    scalarArray<scalarArray<Scalar, ParticSys::volume()>, (MaxDepth + 1)*(MaxDepth + 2)/2> extrapTab;
    scalarArray<Scalar, (MaxDepth + 1)*(MaxDepth + 2)/2> CC;
    scalarArray<Scalar, MaxDepth + 1> macroStepLength;
    scalarArray<Scalar, MaxDepth + 1> work;
    scalarArray<Scalar, MaxDepth + 1> fmin;
    scalarArray<size_t, MaxDepth + 1> cost;
    scalarArray<size_t, MaxDepth + 1> nSteps;
    Scalar        absoluteError{1e-15};
    Scalar        relativeError{1e-15};
    Scalar        s1{0.94};
    Scalar        s2{0.95};
    Scalar        s3{0.02};
    Scalar        s4{4.0};
    size_t        iterDepth{7};
    //////////////////////////////Private Function//////////////////////////////////
private:
      void copyDataToExtrapTab(size_t k);
      bool checkRejection(Scalar error, size_t k) const;
      void extrapolate(size_t k);
    Scalar getError(size_t k) const;
    Scalar getTimeStepCoef(Scalar error, size_t order);
    Scalar prepareForNewIteration(size_t k, bool lastRejection);
};
    /////////////////////////////size_terface Implement////////////////////////////////
template <typename ParticSys, typename Integrator>
BSIterator<ParticSys, Integrator>::BSIterator()
{
    Scalar ratio = 1;
    nSteps[0]    = 1;
    cost[0]      = 1;
    CC[0]        = 1;
    
    for(size_t i = 1 ; i < MaxDepth + 1; ++i)
    {
        nSteps[i] = 2 * i;
        cost[i]   = cost[i - 1] + nSteps[i];
        
        for(size_t j = 0 ; j < i; ++j)
        {
            ratio               = (Scalar)nSteps[i]/(Scalar)nSteps[i - j - 1];
            CC[i*(i + 1)/2 + j] = 1.0/(ratio*ratio - 1);
        }
    }
    
    for(size_t i = 0 ; i < MaxDepth*2 + 2; ++i)
    {
        fmin[i] = pow(s3, 1.0/i);
    }
}

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
{
    Scalar error  = 0;
    Scalar H      = stepLength;
    Scalar h      = stepLength;
    bool   reject = false;
    
    for(;;)
    {
        for(size_t k = 0 ; k <= iterDepth + 1 ; ++k)
        {
            localSystem = particles;

            h = H/nSteps[k];
            
            localSystem.advancePos(0.5*h);
            for(size_t i = 1 ; i < nSteps[k]; i++)
            {
                localSystem.advanceVel(h);
                localSystem.advancePos(h);
            }
            localSystem.advanceVel(h);
            localSystem.advancePos(0.5*h);
            
            copyDataToExtrapTab(k);
            extrapolate(k);
            error = getError(k);
            macroStepLength[k] = H*getTimeStepCoef(error, 2*k+1);
            work[k] = cost[k] / macroStepLength[k];
            //std::cout << localSystem.time << " "<< k << " " << iterDepth << " " << error << " " << H << "\r\n";
            if(k == iterDepth - 1 || k == iterDepth || k == iterDepth + 1)
            {
                if(error <= 1)
                {
                    reject = false;
                    H = prepareForNewIteration(k, reject);
                    particles.load(extrapTab[k*(k + 1)/2 + k]);
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

//////////////////////////////Private Implement/////////////////////////////////
template <typename ParticSys, typename Integrator>
void BSIterator<ParticSys, Integrator>::copyDataToExtrapTab(size_t k)
{
    size_t n = k*(k + 1)/2;
    extrapTab[n] = localSystem.array();
}

template <typename ParticSys, typename Integrator>
bool BSIterator<ParticSys, Integrator>::checkRejection(Scalar error, size_t k) const
{
    Scalar stepLengthRatio = 1;
    
    if(k == iterDepth - 1)
        stepLengthRatio = nSteps[k + 1]*nSteps[k + 2]/(Scalar)(nSteps[0]*nSteps[0]);
    else if(k == iterDepth)
        stepLengthRatio = nSteps[k + 1]/(Scalar)nSteps[0];
    else
        return true;//k == iterDepth+1 and error >1 reject directly
    
    return error > stepLengthRatio*stepLengthRatio;
}

template <typename ParticSys, typename Integrator>
void BSIterator<ParticSys, Integrator>::extrapolate(size_t k)
{
    if(k > 0)
    {
        size_t n    = k*(k + 1)/2;
        size_t pn   = (k - 1)*k/2;
        size_t size = extrapTab[n].size();
        
        for(size_t j = 0 ; j < k; ++j)
        {
            size_t now  = n + j;
            size_t next = n + j + 1;
            size_t last = pn + j;
            for(size_t i = 0 ; i < size ; ++i)
                extrapTab[next][i] = extrapTab[now][i] + (extrapTab[now][i] - extrapTab[last][i])*CC[now];
        }
    }
    else
    {
        return;
    }
}

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::getError(size_t k) const
{
    if(k != 0)
    {
        size_t topk     = k*(k+1)/2 + k;
        size_t subk     = topk - 1;
        size_t size     = extrapTab[topk].size();
        Scalar maxError = 0;
        Scalar error    = 0;
        for(size_t i = 0 ; i < size; ++i)
        {
            error = abs(extrapTab[topk][i] - extrapTab[subk][i])
            /(min(abs(extrapTab[topk][i]), abs(extrapTab[subk][i]) )*this->relativeError + this->absoluteError);
            maxError = max(maxError, error);
            
        }
        return maxError;
    }
    else
    {
        return 1;
    }
}

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::getTimeStepCoef(Scalar error, size_t order)
{
    if(error != 0)
        return max(fmin[order]/s4, min( s1*pow(s2/error, 1.0/order), 1/fmin[order]));
    else
        return 1/fmin[order];
}

template <typename ParticSys, typename Integrator>
typename ParticSys::Scalar BSIterator<ParticSys, Integrator>::prepareForNewIteration(size_t k, bool lastRejection)
{
    Scalar newH = 0;
    if(k == iterDepth-1)
    {
        iterDepth = min(MaxDepth - 1, max(2, k) );
        newH      = macroStepLength[iterDepth];
        
        if(work[k - 1] < 0.8*work[k] && k >= 2)
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 1) );
            newH      = macroStepLength[iterDepth];
        }
        else if(work[k] < 0.9*work[k - 1] || iterDepth <= 2)
        {
            iterDepth = min(MaxDepth - 1, max(2, k + 1) );
            newH      = macroStepLength[iterDepth - 1]*cost[iterDepth]/static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }
    else if(k == iterDepth)
    {
        iterDepth = min(MaxDepth - 1, max(2, k) );
        newH      = macroStepLength[iterDepth];
        
        if(work[k - 1] < 0.8*work[k])
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 1) );
            newH      = macroStepLength[iterDepth];
        }
        else if(!lastRejection && work[k] < 0.9*work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k + 1) );
            newH      = macroStepLength[iterDepth - 1]*cost[iterDepth]/static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }
    else if(k == iterDepth + 1)
    {
        iterDepth = min(MaxDepth - 1, max(2, k - 1) );
        newH      = macroStepLength[iterDepth];
        
        if(work[k-2] < 0.8*work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k - 2) );
            newH      = macroStepLength[iterDepth];
        }
        else if(!lastRejection && work[k] < 0.9*work[k - 1] )
        {
            iterDepth = min(MaxDepth - 1, max(2, k) );
            newH      = macroStepLength[iterDepth - 1]*cost[iterDepth]/static_cast<Scalar>(cost[iterDepth - 1]);
        }
    }
    return newH;
}
#endif
