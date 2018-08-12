
#ifndef IAS15ITERATOR_H
#define IAS15ITERATOR_H

namespace SpaceH
{
   
    namespace Radau
    {
        constexpr double maxStepCof = 1.74867862159014;// = pow(0.02,1/7)
    }
    
    /** @brief IAS15 iterator see details in https://arxiv.org/abs/1409.4779 .
     *
     */
    template <typename ParticSys, typename Integrator>
    class IAS15Iterator
    {
    public:
        /* Typedef */
        using type         = typename ParticSys::type;
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
            Scalar iterH = stepLength;
            RadauTab iterBTab = integrator_.getBTab();
            RadauTab BTab0 = iterBTab;
            
            integrator_.template checkTabVolume<ParticSys::arraySize>(particles.particleNumber());
            for(size_t k = 0; k < max_iter_; ++k)
            {
                integrator_.calcuBTab(particles, iterH);
                if(isConvergent(iterBTab, integrator_.getBTab(), localSystem_.acc()))
                {
                    Scalar error = calcuBError(integrator_.getBTab(), localSystem_.acc());
                    iterH *= optimalStepCoef(error);
                    if(error < 1)
                    {
                        integrator_.evaluateSystemAt(particles, iterH, Integrator::finalPoint);
                        return iterH;
                    }
                    else//current stepSize is too large, restart the iteration with smaller iterH that has been determined by current error.
                    {
                        integrator_.setBTab(BTab0);
                        k = 0;
                    }
                }
                iterBTab = integrator_.getBTab();
            }
            SpaceH::errMsg("IAS15: iteration exceed the max iteration depth!", __FILE__, __LINE__);
            return 0;
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
        bool isConvergent(const RadauTab& BTab, const RadauTab& newBTab, const VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
                msr_err += ((BTab[i][6] - newBTab[i][6])/acc[i]).norm2();
            
            return sqrt(msr_err/size/3) < relativeError_;
        }
        
        Scalar calcuBError(const RadauTab& BTab, const VectorArray& acc)
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
        
        constexpr static size_t max_iter_ = 12;
    };

}
#endif
