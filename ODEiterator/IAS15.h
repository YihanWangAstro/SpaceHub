
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

            resetLastConvergence();
            for(size_t k = 0; k < max_iter_; ++k)
            {
                integrator_.calcuBTab(particles, iterH);
        
                if(isConvergent(iterBTab, integrator_.getBTab(), integrator_.localAcc()))
                {
                    Scalar error = calcuBError(integrator_.getBTab(), integrator_.localAcc());
                    Scalar stepQ = optimalStepCoef(error);
                    
                    if(error < 1)
                    {
                        //std::cout << "acept:" << error << ' ' << k << ' ' << iterH << '\n';
                        integrator_.evaluateSystemAt(particles, iterH, Integrator::finalPoint);
                        integrator_.predictNewB(stepQ);
                        iterH *= stepQ;
                        return iterH;
                    }
                    else//current stepSize is too large, restart the iteration with smaller iterH that has been determined by current error.
                    {
                        //std::cout << "reject:" << error << ' ' << k << ' ' << iterH << '\n';
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
        bool isConvergent(const RadauTab& BTab, const RadauTab& newBTab, const VectorArray& acc)
        {
            size_t size = acc.size();
            
            Scalar diff = 0;
            Scalar scale  = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                diff  = SpaceH::max(diff, (BTab[i][6] - newBTab[i][6]).abs().max_component());
                scale = SpaceH::max(scale,  acc[i].abs().max_component());
            }
            
            Scalar convergence = diff/scale;
            
            if(convergence  >= last_convergence_)//begin to oscillate
            {
                resetLastConvergence();
                return true;
            }
            else
            {
                last_convergence_ = convergence;
                return convergence < convergent_limit_;
            }
        }
        
        Scalar calcuBError(const RadauTab& BTab, const VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar diff  = 0;
            Scalar scale = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                //if(vel[i].norm()*dt < pos[i].norm()*1e-8)
                  //  continue;
                
                diff = SpaceH::max(diff, BTab[i][6].abs().max_component());
                scale = SpaceH::max(scale, acc[i].abs().max_component());
            }
            return diff/(scale*relativeError_);
        }
        
        inline Scalar optimalStepCoef(Scalar error)
        {
            if(error == 0)
                return Radau::maxStepCof;
            else
                return SpaceH::max(Radau::minStepCof ,SpaceH::min(pow(0.55/error, 1.0/7), Radau::maxStepCof));
        }
        
        inline void resetLastConvergence()
        {
            last_convergence_ = Radau::maxFloat;
        }
    private:
        Integrator integrator_;
        
        Scalar relativeError_{1e-9};
        
        Scalar convergent_limit_{1e-16};
        
        Scalar last_convergence_{Radau::maxFloat};
        
        size_t particleNum_{ParticSys::arraySize};
        
        constexpr static size_t max_iter_ = 12;
    };

}
#endif
