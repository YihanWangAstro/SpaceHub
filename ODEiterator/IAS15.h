
#ifndef IAS15ITERATOR_H
#define IAS15ITERATOR_H

namespace SpaceH
{
   
    namespace Radau
    {
        constexpr double maxStepCof = 1.74867862159014;// = pow(0.02,1/7)
    }
    
    /** @brief Gauss RADAU iterator
     *
     */
    template <typename ParticSys, typename Integrator>
    class IAS15
    {
    public:
        /* Typedef */
        using type         = typename ParticSys::type
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        using RadauArray   = typename Integrator::RadauArray;
        using RadauTab     = typename Integrator::RadauTab;
        /* Typedef */
        /* Typedef */
        
        /*Template parameter check*/
        CHECK_TYPE(ParticSys, Integrator)
        /*Template parameter check*/
        
        /** @brief interface to iterate particle system for one step
         *  @param particles  Particle system needs evolution.
         *  @param integrator Integrator to integrate the particle system.
         *  @param stepLength Macro step length for iteration(Here, the step length of the integrator).
         *  @return step length for next iteration.
         */
        Scalar iterate(ParticSys& particles, Scalar stepLength)
        {
            checkTabVolume(particles.particleNumber());
            Scalar iterH = stepLength;
            RadauTab iterBTab = BTab0 = integrator_.getBTab();
            for(size_t k = 0; k < max_iter_; ++k)
            {
                localSystem_ = particles;
                //integrator.setBTab(iterBTab);
                integrator_.integrate(localSystem_, iterH);//advance system with Gauss-Radau method and update all coresponding B values in interpolation.
                
                if(isConvergent(iterBTab, integrator_.getBTab(), localSystem_.acc()))
                {
                    Scaclar error = calcuBError(integrator_.getBTab(), localSystem_.acc());
                    iterH *= optimalStepCoef(error);
                    if(error < 1)
                    {
                        particles = localSystem_;
                        return iterH;
                    }
                    else//current stepSize is too large, restart the iteration with smaller iterH that been determined by current error.
                    {
                        integrator.setBTab(BTab0);
                        k = 0;
                    }
                }
                iterBTab = integrator_.getBTab();
            }
            SpaceH::errMsg("IAS15: iteration exceed the max iteration depth!", __FILE__, __LINE__);
        }
        
        void setRelativeError(Scalar rel_tol)
        {
            relativeError_ = rel_tol;
        }
        
        void SetConvergentLimit(Scalar limit)
        {
            convergent_limit_ = limit;
        }
        
    private:
        
        void checkTabVolume(size_t particleNum)
        {
            if(particleNum_ < particleNum)
            {
                BTab_.resize(particleNum);
                GTab_.resize(particleNum);
                particleNum_ = particleNum;
            }
        }
        
        bool isConvergent(RadauTab& BTab, RadauTab& newBTab, VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
                msr_err += ((BTab[i][6] - newBTab[i][6])/acc[i]).norm2();
            
            return sqrt(msr_err/size/3) < relativeError_;
        }
        
        Scalar calcuBError(RadauTab& BTab, VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
                msr_err  += (BTab[i][6]/acc[i]).norm2();
            
            return sqrt(msr_err/size/3)/convergent_limit_;
        }
        
        inline Scalar optimalStepCoef(Scalar error)
        {
            if(error == 0)
                return Radau::maxStepCof;
            else
                return pow(1/error, 1.0/7);
        }
        
    private:
        ParticSys localSystem_;
        
        Integrator integrator_;
        
        Scalar relativeError_{1e-16};
        
        Scalar convergent_limit_{1e-9};
        
        size_t particleNum_{ParticSys::arraySize};
        
        constexpr static size_t max_iter_ = 12；
    };

}
#endif
