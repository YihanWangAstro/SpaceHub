
#ifndef IAS15ITERATOR_H
#define IAS15ITERATOR_H

namespace SpaceH
{
   
    namespace Radau
    {
        constexpr double maxStepCof = 1.74867862159014;// = 1/pow(0.02,1/7)
        constexpr double minStepCof = 0.5718603679678214;// = pow(0.02,1/7)
        constexpr double maxFloat   = 1e100;
    }
    
    /** @brief IAS15 iterator see details in https://arxiv.org/abs/1409.4779 .
     *
     */
    template <typename ParticSys, typename Integrator>
    class IAS15
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
            integrator_.template checkTabVolume<ParticSys::arraySize>(particles.particleNumber());
            
            Scalar iterH      = stepLength;
            RadauTab iterBTab = integrator_.getBTab();//get the b value of the last step.

            //std::cout << '\n';
            for(size_t k = 0; k < max_iter_; ++k)
            {
                integrator_.calcuBTab(particles, iterH);
        
                if(isConvergent(iterBTab, integrator_.getBTab(), particles.acc(), k))
                {
                    Scalar error = calcuBError(integrator_.getBTab(), particles.acc());
                    Scalar stepQ = optimalStepCoef(error);
                    
                    if(error < 1)
                    {
                        integrator_.evaluateSystemAt(particles, iterH, Integrator::finalPoint);
                        iterH *= stepQ;
                        integrator_.predictNewB(stepQ);
                        return iterH;
                    }
                    else//current stepSize is too large, restart the iteration with smaller iterH that has been determined by current error.
                    {
                        iterH *= stepQ;
                        integrator_.predictNewB(stepQ);
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
        bool isConvergent(const RadauTab& BTab, const RadauTab& newBTab, const VectorArray& acc, size_t iter)
        {
            size_t size = acc.size();
            Scalar maxdb = 0;
            Scalar maxa  = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                maxdb = SpaceH::max(maxdb, (BTab[i][6] - newBTab[i][6]).abs().max_component());
                maxa  = SpaceH::max(maxa,  acc[i].abs().max_component());
            }
            if(maxdb >= last_delta_b)//begin to oscillate
            {
                resetLastDeltaB();
                return true;
            }
            else
            {
                last_delta_b = maxdb;
                return maxdb/maxa < relativeError_;
            }
            /*Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
                msr_err += ((BTab[i][6] - newBTab[i][6])/(acc[i]+absoluteError_)).norm2();

            return sqrt(msr_err/size/3) < relativeError_;*/
        }
        
        Scalar calcuBError(const RadauTab& BTab, const VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar maxb = 0;
            Scalar maxa = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                maxb = SpaceH::max(maxb, BTab[i][6].abs().max_component());
                maxa = SpaceH::max(maxa, acc[i].abs().max_component());
            }
            return maxb/maxa/convergent_limit_;
            /*Scalar msr_err = 0;
             for(size_t i = 0 ; i < size ; ++i)
             msr_err  += (BTab[i][6]/(acc[i]+absoluteError_)).norm2();
             
             return sqrt(msr_err/size/3)/convergent_limit_;*/
        }
        
        inline Scalar optimalStepCoef(Scalar error)
        {
            if(error == 0)
                return Radau::maxStepCof;
            else
                return SpaceH::min(pow(0.95/error, 1.0/7), Radau::maxStepCof);
        }
        
        inline void resetLastDeltaB()
        {
            last_delta_b = Radau::maxFloat;
        }
    private:
        Integrator integrator_;
        
        Scalar relativeError_{1e-16};
        
        Scalar convergent_limit_{1e-9};
        
        Scalar last_delta_b{Radau::maxFloat};
        
        size_t particleNum_{ParticSys::arraySize};
        
        constexpr static size_t max_iter_ = 20;
    };

}
#endif
