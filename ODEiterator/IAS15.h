
#ifndef IAS15ITERATOR_H
#define IAS15ITERATOR_H

namespace SpaceH
{
   
    namespace Radau
    {
        /** Gauss radau stepping. see details in
         https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
         */
        
        // all const computation here need to be rewrite with template metaprogramming
        
        constexpr double H[8] = {0, 0.05626256053692215, 0.18024069173689236, 0.35262471711316964, 0.54715362633055538, 0.73421017721541053, 0.88532094683909577, 0.97752061356128750};
        
        constexpr double c[21] =
        {
            -0.05626256053692215,  0.01014080283006363, -0.00357589772925162,  0.00195656540994722, -0.00143653023637089,  0.00127179030902687,
            -0.23650325227381452,  0.09353769525946207, -0.05475538688906869,  0.04215852772126871, -0.03876035791590677,
            -0.58912796938698420,  0.41588120008230689, -0.36009959650205681,  0.36096224345284600,
            -1.13628159571753962,  1.25015071184069093, -1.46688420840042699,
            -1.87049177293294999,  2.90613625930842900,
            -2.75581271977204567
        };
        
        constexpr double d[21] =
        {
            0.05626256053692215, 0.00316547571817083, 0.00017809776922174, 0.00001002023652233, 0.00000056376416393, 0.00000003171881540,
            0.23650325227381452, 0.04579298550602792, 0.00843185715352570, 0.00152978400250047, 0.00027629309098265,
            0.58912796938698420, 0.25353406905456927, 0.09783423653244401, 0.03602855398373646,
            1.13628159571753962, 0.87525466468409108, 0.57673300027707874,
            1.87049177293294999, 2.24858876076915948,
            2.75581271977204567
        };
        
        constexpr double r[28] =
        {
            17.77380891407799979, 5.54813671853721679, 2.83587607864443880, 1.82764026751759778, 1.36200781606246957, 1.12953387533678984, 1.02299632982348676,
            8.06593864838188779, 3.37424997696263551, 2.03711183535858442, 1.47504021756041159, 1.20618766605844563, 1.08547219393864247,
            5.80100155926406202, 2.72544221180822621, 1.80515358014025140, 1.41827826373473909, 1.25426462228187785,
            5.14062410581093232, 2.62064492638703506, 1.87724249618681016, 1.60026654949081615,
            5.34597689987110947, 2.95711601729045581, 2.32359830021969449,
            6.61766201370242158, 4.10997577834455807,
            10.84602619023684689
        };
        
        template<typename ScalarArray>
        inline void transB2G(ScalarArray& B, ScalarArray& G)
        {
            B[0] = G[0] + c[0]*G[1] + c[1]*G[2] + c[2] *G[3] + c[3] *G[4] + c[4] *G[5] + c[5] *G[6];
            B[1] =             G[1] + c[6]*G[2] + c[7] *G[3] + c[8] *G[4] + c[9] *G[5] + c[10]*G[6];
            B[2] =                         G[2] + c[11]*G[3] + c[12]*G[4] + c[13]*G[5] + c[14]*G[6];
            B[3] =                                      G[3] + c[15]*G[4] + c[16]*G[5] + c[17]*G[6];
            B[4] =                                                   G[4] + c[18]*G[5] + c[19]*G[6];
            B[5] =                                                                G[5] + c[20]*G[6];
            B[6] =                                                                             G[6];
        }
        
        template<typename ScalarArray>
        inline void transG2B(ScalarArray& G, ScalarArray& B)
        {
            G[0] = B[0] + c[0]*B[1] + c[1]*B[2] + c[2] *B[3] + c[3] *B[4] + c[4] *B[5] + c[5] *B[6];
            G[1] =             B[1] + c[6]*B[2] + c[7] *B[3] + c[8] *B[4] + c[9] *B[5] + c[10]*B[6];
            G[2] =                         B[2] + c[11]*B[3] + c[12]*B[4] + c[13]*B[5] + c[14]*B[6];
            G[3] =                                      B[3] + c[15]*B[4] + c[16]*B[5] + c[17]*B[6];
            G[4] =                                                   B[4] + c[18]*B[5] + c[19]*B[6];
            G[5] =                                                                B[5] + c[20]*B[6];
            G[6] =                                                                             B[6];
        }
        
        
    }
    
    /** @brief Gauss RADAU iterator
     *
     */
    template <typename ParticSys, typename Integrator>
    class IAS15
    {
    public:
        /* Typedef */
        using type = typename ParticSys::type
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        using VectorArray = typename type::VectorArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        
        
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
        Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
        {
            checkTabVolume(particles.particleNumber());
            Scalar err = 1;
            for(size_t k = 0 ; k < max_iter_; ++k)
            {
                for(size_t i = 0 ; i < 7; ++i)
                {
                    localSystem_ = particles;
                    Scalar h = stepLength*Radau::H[i];
                    integrator_.mountBTab(Btab_);
                    integrator_.integrate(localSystem_, h);
                    calcuGTab(particles.acc(), localSystem_.acc(), i);
                }
                Radau::transG2B(Gtab_, newBtab_);
                err = checkErr(newBtab_, Btab_);
                Btab_ = newBtab_;
                if(err < 1)
                {
                    return nextStepLength(stepLength, err);
                }
            }
            SpaceH::errMsg("IAS15: iteration exceed the max iteration depth!", __FILE__, __LINE__);
        }
        
        
        IAS15()
        {
            //init tab
        }
        
    private:
        
        void checkTabVolume(size_t particleNum)
        {
            if(particleNum_ < particleNum)
            {
                cofBTab_.resize(particleNum);
                cofGTab_.resize(particleNum);
                particleNum_ = particleNum;
            }
        }
        
        void calcuGTab(VectorArray& a0, VectorArray& a, size_t step)
        {
            
        }
        
        Scalar checkErr(RadauTab& Btab, RadauTab& newBtab)
        {
            
        }
        
        Scalar nextStepLength(Scalar currentStepLength, Sclar error)
        {
            return currentStepLength;
        }
        
    private:
        Integrator integrator_;
        ParticSys localSystem_;
        RadauTab    Btab_;
        RadauTab    newBtab_;
        RadauTab    Gtab_;
        size_t      particleNum_{0};
        constexpr static size_t max_iter_ = 12；
    };

}
#endif
