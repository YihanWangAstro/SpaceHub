#ifndef BSITERATOR_H
#define BSITERATOR_H
#include "../libs.h"
#include "../devTools.h"
namespace SpaceH
{
    
    /** @brief Bulirsch-Stoer extrapolation algorithm*/
    template <typename ParticSys, typename Integrator>
    class BSIterator
    {
    public:
        /* Typedef */
        using type   = typename ParticSys::type;
        using Scalar = typename type::Scalar;
        using ScalarBuffer = typename type::ScalarBuffer;
        
        template<typename T, size_t S>
        using Container = typename type::template Container<T, S>;
        /* Typedef */
        
        /*Template parameter check*/
        CHECK_TYPE(ParticSys, Integrator)
        /*Template parameter check*/
        
        /** @brief Constructor for initializing work_, nSteps, fmin and CC*/
        BSIterator();
        
        /** @brief Interface of ODE iterator
         *  @note  BSiterator will force use the internal mid-point integrator as the basic integrator.
         */
        Scalar iterate(ParticSys& particles, Scalar stepLength);
        
        /** @brief Set the local relative error*/
        void setRelativeError(Scalar relError)
        {
            relativeError_ = relError;
        }
        
        /** @brief Set the local absolute error*/
        void setAbsoluteError(Scalar absError)
        {
            absoluteError_ = absError;
        }
        
        Scalar getRejRate()
        {
            return static_cast<float> (rej_num_)/iter_num_;
        }
    private:
        /** @brief The Maximum iteration depth*/
        static const size_t MaxDepth{8};
        
        /** @brief The local partical system used to iterate.*/
        ParticSys localSystem;
        
        /** @brief Extrapolation table.*/
        Container < ScalarBuffer, (MaxDepth + 1)* (MaxDepth + 2) / 2 > extrapTab;
        
        /** @brief Extrapolation coefficient.*/
        Container < Scalar, (MaxDepth + 1)* (MaxDepth + 2) / 2 > extrap_coef_;
        
        /** @brief The optimal step size array.*/
        Container < Scalar, MaxDepth + 1 > optimal_H_;
        
        /** @brief The work(computation resource) needed to converge at column k.*/
        Container < Scalar, MaxDepth + 1 > work_per_len_;
        
        /** @brief The exponent of error estimate at column k.*/
        Container < Scalar, MaxDepth + 1 > err_expon_;
        
        /** @brief The minimal coeffient of integration step estimation.*/
        Container < Scalar, MaxDepth + 1 > stepsize_limiter;
        
        /** @brief The work(computation resource) per step size of each iteration depth.*/
        Container < size_t, MaxDepth + 1 > work_;
        
        /** @brief Steps of integration of each iteration depth.*/
        Container < size_t, MaxDepth + 1 > sub_steps_;
        
        /** @brief Local absolute error*/
        Scalar absoluteError_{1e-15};
        
        /** @brief Local relative error*/
        Scalar relativeError_{1e-15};
        
        /** @brief Current iteraation depth.*/
        size_t optimal_iter_{7};
        
        /** @brief Total volume of extrapolation table(in scalar).*/
        size_t extrapTabVolume_{0};
        
        /** @brief Rejectin number*/
        size_t rej_num_{0};
        
        /** @brief Total iteration number*/
        size_t iter_num_{0};
        
    private:
        /** @brief Internal integrator */
        Integrator integrator;
        
        /** @brief Check the rejection criteria for current iteration.*/
        bool divergedInOrderWindow(Scalar error, size_t k) const;
        
        /** @brief Extrapolate the kth row to the right end.*/
        void extrapolate(size_t k);
        
        /** @brief Check/resize the extrapolation table volume.*/
        void checkExtrapVolume();
        
        /** @brief Calculate the error of the k row of the extrapolation table.*/
        Scalar calcuError(size_t k) const;
        
        /** @brief Calculate the new iteration integration step coefficient*/
        Scalar calcuIdealStepCoef(Scalar error, size_t k);
        
        /** @brief Prepare the next iteration.*/
        Scalar prepareNextIteration(size_t iter, bool last_rejected);
        
        /** @brief Evolve the local system*/
        void evolveLocalSys(size_t iter, Scalar step);
        
        /** @brief step size reductio after rejection.*/
        inline Scalar stepSizeReduction(size_t iter)
        {
            if(iter == optimal_iter_ - 1)
                return optimal_H_[iter];
            else
                return optimal_H_[optimal_iter_];
        }
        
        /** @brief check if current iteration can converge in current order window*/
        inline bool isInConvergenceWindow(size_t iter)
        {
            return iter == optimal_iter_ - 1 || iter == optimal_iter_ || iter == optimal_iter_ + 1;
        }
        
        /** @brief return convergence column in allowed range*/
        inline size_t allowedRange(size_t i)
        {
            return SpaceH::min(MaxDepth - 1, SpaceH::max(2, i));
        }
        
        /** @brief Step size limiter*/
        inline Scalar stepSizeLimiter(Scalar H, size_t k)
        {
            return SpaceH::max(stepsize_limiter[k] / 4, SpaceH::min(1.0*H, 1.0/stepsize_limiter[k]) );
        }
    };
    
    /** @brief Constructor for initializing cost, nSteps, fmin and CC*/
    template <typename ParticSys, typename Integrator>
    BSIterator<ParticSys, Integrator>::BSIterator()
    {
        Scalar ratio        = 1;
        stepsize_limiter[0] = 0.02;
        sub_steps_[0]       = 1;
        work_[0]            = 1;
        extrap_coef_[0]     = 1;
        err_expon_[0]       = 1;
        extrap_coef_[0]     = 1;
        for(size_t i = 1 ; i < MaxDepth + 1; ++i)
        {
            sub_steps_[i] = 2*i;
            work_[i]  = work_[i - 1] + sub_steps_[i];
            err_expon_[i] = static_cast<Scalar>(1.0 / (2*i - 1));
            stepsize_limiter[i] = pow(0.02, err_expon_[i]);
            
            for(size_t j = 1 ; j <= i; ++j)
            {
                ratio = static_cast<Scalar>(sub_steps_[i]) / static_cast<Scalar>(sub_steps_[i - j]);
                extrap_coef_[i * (i + 1) / 2 + j] = 1.0 / (ratio * ratio - 1);
            }
        }
        
        absoluteError_ = 5*std::numeric_limits<typename SpaceH::get_value_type<Scalar>::type>::epsilon();
        relativeError_ = 5*std::numeric_limits<typename SpaceH::get_value_type<Scalar>::type>::epsilon();
        
    }

    template <typename ParticSys, typename Integrator>
    void BSIterator<ParticSys, Integrator>::evolveLocalSys(size_t iter, Scalar macroStepSize)
    {
        size_t steps = sub_steps_[iter];
        Scalar h = macroStepSize / steps;
        
        localSystem.drift(0.5 * h);
        for(size_t i = 1 ; i < steps; i++)
        {
            localSystem.kick(h);
            localSystem.drift(h);
        }
        localSystem.kick(h);
        localSystem.drift(0.5 * h);
    }
    
    /** @brief Interface of ODE iterator
     *  @param particles Particle system need iteration.
     *  @param integrator Basica integrator used to evolve, but here BS iterator will force use internal mid-point integrator.
     *  @param stepLength Macro integration step length.
     *  @return The next macro integration step length.
     */
    template <typename ParticSys, typename Integrator>
    typename ParticSys::type::Scalar BSIterator<ParticSys, Integrator>::iterate(ParticSys& particles, Scalar stepLength)
    {
        Scalar error     = 0;
        Scalar iter_H    = stepLength;
        bool last_rejected = false;
        
        for(;;)
        {
            //std::cout << '\n' << "H ="<<iter_H <<"\r\n";
            iter_num_ ++;
            
            localSystem = particles;
            evolveLocalSys(0, iter_H);
            localSystem.write(extrapTab[0], NbodyIO::ACTIVE);//copyDataToExtrapTab;
            
            checkExtrapVolume();
            
            for(size_t iter = 1 ; iter <= optimal_iter_ + 1 ; ++iter)
            {
                localSystem = particles;
                evolveLocalSys(iter, iter_H);
                
                localSystem.write(extrapTab[iter * (iter + 1) / 2], NbodyIO::ACTIVE);//copyDataToExtrapTab;
                
                extrapolate(iter);
                
                error = calcuError(iter);
                
                Scalar step_cof = calcuIdealStepCoef(error, iter);
                
                optimal_H_[iter] = iter_H * step_cof;
                
                work_per_len_[iter] = work_[iter] / step_cof;
                
                 //std::cout << iter << ' ' << optimal_iter_ << " err: " << error << " H: " << optimal_H_[iter] << " W:" << work_per_len_[iter] << " T cof:" << calcuIdealStepCoef(error, iter) << "\n";
                
                    
                if(isInConvergenceWindow(iter))
                {
                    if(error < 1.0)
                    {
                        last_rejected = false;
                        
                        iter_H = prepareNextIteration(iter, last_rejected);
                        localSystem.read(extrapTab[iter * (iter + 1) / 2 + iter], NbodyIO::ACTIVE);
                        particles = localSystem;
                        return iter_H;
                    }
                    else if(divergedInOrderWindow(error, iter))
                    {
                        rej_num_ ++;
                
                        //iter_H = prepareNextIteration(iter, true);
                        iter_H = stepSizeReduction(iter);
                        last_rejected = true;
                        break;
                    }
                }
            }
        }
    }
    
    /** @brief Check the last_rejected criteria for current iteration.
     *  @param k The kth iteration depth.
     *  @return If current iteration is rejected.
     */
    template <typename ParticSys, typename Integrator>
    bool BSIterator<ParticSys, Integrator>::divergedInOrderWindow(Scalar error, size_t iter) const
    {
        Scalar r = 1;
        
        if(iter == optimal_iter_ - 1)
            r = static_cast<Scalar>(sub_steps_[iter + 1] * sub_steps_[iter + 2]) / static_cast<Scalar>(sub_steps_[1] * sub_steps_[1]);
        else if(iter == optimal_iter_)
            r = static_cast<Scalar>(sub_steps_[iter + 1]) / static_cast<Scalar>(sub_steps_[1]);
        else
            return  error > 1.0;;//k == iterDepth+1 and error >1 reject directly
        
        return error > r * r;
    }
    
    /** @brief Extrapolate the k-th row to the right end.
     *  @param k The kth row of extrapolation table.
     */
    template <typename ParticSys, typename Integrator>
    void BSIterator<ParticSys, Integrator>::extrapolate(size_t k)
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
                extrapTab[next][i] = extrapTab[now][i] + (extrapTab[now][i] - extrapTab[last][i]) * extrap_coef_[next];
        }
    }
    
    /** @brief Calculate the error of the k row of the extrapolation table.
     *  @param k The kth row of extrapolation table.
     */
    template <typename ParticSys, typename Integrator>
    typename ParticSys::type::Scalar BSIterator<ParticSys, Integrator>::calcuError(size_t k) const
    {
        if(k != 0)
        {
            size_t topk     = k * (k + 1) / 2 + k;
            size_t subk     = topk - 1;
            size_t size     = extrapTab[topk].size();
            Scalar max_err  = 0;
            
            for(size_t i = 0 ; i < size; ++i)
            {
                Scalar scale = (SpaceH::max(SpaceH::abs(extrapTab[topk][i]), SpaceH::abs(extrapTab[subk][i]) ) * relativeError_ + absoluteError_);
                Scalar d = extrapTab[topk][i] - extrapTab[subk][i];
                max_err = SpaceH::max(1.0*max_err, SpaceH::abs(d/scale));
            }
            return max_err;
        }
        else
            return 1e10;
    }
    
    /** @brief Calculate the new iteration integration step coefficient
     *  @param error The error of the k-th row of extrapolation table.
     *  @param order The order of the error in k-th row of extrapolation table.
     *  @return The new iteration integration step length coefficient.
     */
    template <typename ParticSys, typename Integrator>
    typename ParticSys::type::Scalar BSIterator<ParticSys, Integrator>::calcuIdealStepCoef(Scalar error, size_t k)
    {
        if(error != 0)
            return SpaceH::max(stepsize_limiter[k] / 4, SpaceH::min( 0.94*pow(0.95 / error, err_expon_[k]), 1.0/stepsize_limiter[k] ) );
        else
            return  1.0 / stepsize_limiter[k];
    }
    
    /** @brief Calculate the new iteration depth for next iteration.
     *  @param iter The current iteration index depth.
     *  @param lastRejection The rejection status of last trail of iteration.
     *  @return The new iteration step length.
     */
    template <typename ParticSys, typename Integrator>
    typename ParticSys::type::Scalar BSIterator<ParticSys, Integrator>::prepareNextIteration(size_t iter, bool last_rejected)
    {
        switch (static_cast<int>(iter - optimal_iter_))
        {
            case -1:
                
                if(work_per_len_[iter] < 0.9 * work_per_len_[iter - 1] || optimal_iter_ <= 2)//optimal_iter_ <= 2 here to avoid none calculated work_per_len[1-1=0]
                {
                    return optimal_H_[iter] * static_cast<Scalar>(work_[iter+1]) / static_cast<Scalar>(work_[iter]);
                }
                else
                {
                    optimal_iter_ = allowedRange(optimal_iter_ - 1);
                    return optimal_H_[optimal_iter_];//reduce order
                }
                
            case 0:
                
                if(work_per_len_[iter - 1] < 0.8 * work_per_len_[iter])
                {
                    optimal_iter_ = allowedRange(optimal_iter_ - 1);
                    return optimal_H_[optimal_iter_];//order decrease
                }
                else if(!last_rejected && work_per_len_[iter] < 0.9 * work_per_len_[iter - 1])
                {
                    optimal_iter_ = allowedRange(optimal_iter_ + 1);
                    return optimal_H_[iter] * static_cast<Scalar>(work_[iter+1]) / static_cast<Scalar>(work_[iter]);//order increase
                }
                else
                {
                    return optimal_H_[optimal_iter_];//keep order
                }
                
            case 1:
                
                if(work_per_len_[iter - 2] < 0.8 * work_per_len_[iter - 1] )
                {
                    optimal_iter_ = allowedRange(optimal_iter_ - 1);
                }
                if(!last_rejected && work_per_len_[iter] < 0.9 * work_per_len_[optimal_iter_] )
                {
                    optimal_iter_ = allowedRange(optimal_iter_ + 1);
                }
                return optimal_H_[optimal_iter_];
                
            default:
                SpaceH::errMsg("unexpected iteration index!",__FILE__, __LINE__);
                exit(0);
        }
    }
    
    /** @brief Check/resize the extrapolation table.*/
    template <typename ParticSys, typename Integrator>
    void BSIterator<ParticSys, Integrator>::checkExtrapVolume()
    {
        size_t tab_size = extrapTab.size();
        size_t array_size = extrapTab[0].size();
        
        if( extrapTabVolume_ < tab_size * array_size )
        {
            for(size_t i = 1 ; i < tab_size ; i++)
            {
                extrapTab[i].resize(array_size);
            }
            extrapTabVolume_ = tab_size * array_size;
        }
    }
}

#endif
