
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
        
        constexpr double interplt[8] = {0.05626256053692215, 0.18024069173689236, 0.35262471711316964, 0.54715362633055538, 0.73421017721541053, 0.88532094683909577, 0.97752061356128750};
        
        constexpr double c[21] =
        {
            -5.62625605369221488e-02, 1.01408028300636298e-02, -3.57589772925161718e-03, 1.95656540994722109e-03, -1.43653023637089149e-03, 1.27179030902686775e-03,
            -2.36503252273814524e-01, 9.35376952594620670e-02, -5.47553868890686898e-02, 4.21585277212687057e-02, -3.87603579159067693e-02,
            -5.89127969386984196e-01, 4.15881200082306890e-01, -3.60099596502056807e-01, 3.60962243452845999e-01,
            -1.13628159571753962e+00, 1.25015071184069093e+00, -1.46688420840042699e+00,
            -1.87049177293294999e+00, 2.90613625930842900e+00,
            -2.75581271977204567e+00
        };
        
        constexpr double d[21] =
        {
            5.62625605369221488e-02, 3.16547571817082972e-03, 1.78097769221743430e-04, 1.00202365223291294e-05, 5.63764163931820866e-07, 3.17188154017613769e-08,
            2.36503252273814524e-01, 4.57929855060279223e-02, 8.43185715352570177e-03, 1.52978400250046594e-03, 2.76293090982647678e-04,
            5.89127969386984196e-01, 2.53534069054569267e-01, 9.78342365324440105e-02, 3.60285539837364582e-02,
            1.13628159571753962e+00, 8.75254664684091077e-01, 5.76733000277078744e-01,
            1.87049177293294999e+00, 2.24858876076915948e+00,
            2.75581271977204567e+00
        };
        
        /*constexpr double r[28] =
        {
            1.77738089140779998e+01,
            5.54813671853721679e+00, 8.06593864838188779e+00,
            2.83587607864443880e+00, 3.37424997696263551e+00, 5.80100155926406202e+00,
            1.82764026751759778e+00, 2.03711183535858442e+00, 2.72544221180822621e+00, 5.14062410581093232e+00,
            1.36200781606246957e+00, 1.47504021756041159e+00, 1.80515358014025140e+00, 2.62064492638703506e+00, 5.34597689987110947e+00,
            1.12953387533678984e+00, 1.20618766605844563e+00, 1.41827826373473909e+00, 1.87724249618681016e+00, 2.95711601729045581e+00, 6.61766201370242158e+00,
            1.02299632982348676e+00, 1.08547219393864247e+00, 1.25426462228187785e+00, 1.60026654949081615e+00, 2.32359830021969449e+00, 4.10997577834455807e+00, 1.08460261902368469e+01
        };*/
        
        /*coef for computing g(derived from r)*/
        constexpr double gg[28] =
        {
            1.77738089140779998e+01,
            4.47509303845560006e+01, 8.06593864838188779e+00,
            5.55095216749226980e+01, 1.95740293777069740e+01, 5.80100155926406202e+00,
            5.21625022561530060e+01, 2.85409022679298947e+01, 1.40104739330160325e+01, 5.14062410581093232e+00,
            5.08080910907447725e+01, 3.73038175637124385e+01, 2.52900342103279819e+01, 1.40099072392295135e+01, 5.34597689987110947e+00,
            7.09853803416487306e+01, 6.28448441358017932e+01, 5.21020450666394531e+01, 3.67361232269326457e+01, 1.95691943377340427e+01, 6.61766201370242158e+00,
            2.30858165231426691e+02, 2.25668615322657305e+02, 2.07899029180855715e+02, 1.65753721732680277e+02, 1.03578820531755125e+02, 4.45769049331641467e+01, 1.08460261902368469e+01
        };
        
        template<typename RadauTab>
        inline void transB2G(RadauTab& B, RadauTab& G)
        {
            size_t size = B.size();
            for(size_t i = 0 ; i < size; ++i)
            {
                B[i][0] = G[i][0] + c[0]*G[i][1] + c[1]*G[i][2] + c[2] *G[i][3] + c[3] *G[i][4] + c[4] *G[i][5] + c[5] *G[i][6];
                B[i][1] =                G[i][1] + c[6]*G[i][2] + c[7] *G[i][3] + c[8] *G[i][4] + c[9] *G[i][5] + c[10]*G[i][6];
                B[i][2] =                               G[i][2] + c[11]*G[i][3] + c[12]*G[i][4] + c[13]*G[i][5] + c[14]*G[i][6];
                B[i][3] =                                               G[i][3] + c[15]*G[i][4] + c[16]*G[i][5] + c[17]*G[i][6];
                B[i][4] =                                                               G[i][4] + c[18]*G[i][5] + c[19]*G[i][6];
                B[i][5] =                                                                               G[i][5] + c[20]*G[i][6];
                B[i][6] =                                                                                               G[i][6];
            }
        }
        
        template<typename RadauTab>
        inline void transG2B(RadauTab& G, RadauTab& B)
        {
            size_t size = B.size();
            for(size_t i = 0 ; i < size; ++i)
            {
                G[i][0] = B[i][0] + d[0]*B[i][1] + d[1]*B[i][2] + d[2] *B[i][3] + d[3] *B[i][4] + d[4] *B[i][5] + d[5] *B[i][6];
                G[i][1] =                B[i][1] + d[6]*B[i][2] + d[7] *B[i][3] + d[8] *B[i][4] + d[9] *B[i][5] + d[10]*B[i][6];
                G[i][2] =                               B[i][2] + d[11]*B[i][3] + d[12]*B[i][4] + d[13]*B[i][5] + d[14]*B[i][6];
                G[i][3] =                                               B[i][3] + d[15]*B[i][4] + d[16]*B[i][5] + d[17]*B[i][6];
                G[i][4] =                                                               B[i][4] + d[18]*B[i][5] + d[19]*B[i][6];
                G[i][5] =                                                                               B[i][5] + d[20]*B[i][6];
                G[i][6] =                                                                                               B[i][6];
            }
        }
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
            RadauTab iterBTab = BTab0 = BTab_;
            for(size_t k = 0; k < max_iter_; ++k)
            {
                for(size_t i = 0 ; i < 7; ++i)
                {
                    localSystem_ = particles;
                    Scalar h = iterH*Radau::interplt[i];
                    integrator.mountBTab(iterBTab);
                    integrator.integrate(localSystem_, h);
                    calcuGTab(particles.acc(), localSystem_.acc(), i);
                }
                
                Radau::transG2B(Gtab_, iterBtab);
                
                if(isConvergent(iterBTab, particles.acc()))
                {
                    Scaclar error = calcuBError(iterBTab, particles.acc());
                    iterH *= optimalStepCoef(error);
                    if(error < 1)
                    {
                        BTab_ = iterBTab;
                        particles = localSystem;
                        return iterH;
                    }
                    else//current stepSize is too large, restart the iteration with smaller iterH that been determined by current error.
                    {
                        iterBTab = BTab0;
                        k = 0;
                    }
                }
                BTab_ = iterBTab;
            }
            SpaceH::errMsg("IAS15: iteration exceed the max iteration depth!", __FILE__, __LINE__);
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
        
        void calcuGTab(VectorArray& a0, VectorArray& a, size_t iter)
        {
            using namespace Radau;
            size_t size  = a0.size();
            size_t offset = iter*(iter+1)/2 + 1;
            for(size_t i = 0 ; i < size; ++i)
            {
                GTab_[i][iter] = (a[i] - a0[i])*gg[offset - 1];
                
                for(size_t j = 0 ; j < iter; ++j)
                {
                    GTab_[i][iter] -= GTab_[i][j]*gg[j + offset];
                }
            }
        }
        
        bool isConvergent(RadauTab& iterBTab, VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                msr_err += ((BTab_[i][6] - iterBTab[i][6])/acc[i]).norm2();
            }
            return sqrt(msr_err/size/3) < relativeError_;
        }
        
        Scalar calcuBError(RadauTab& iterBTab, VectorArray& acc)
        {
            size_t size = acc.size();
            Scalar msr_err = 0;
            for(size_t i = 0 ; i < size ; ++i)
            {
                msr_err  += (iterBTab[i][6]/acc[i]).norm2();
            }
            return sqrt(msr_err/size/3)/convergent_limit_;
        }
        
        inline Scalar optimalStepCoef(Scalar error)
        {
            if(error == 0)
                return Radau::maxStepCof;
            else
                return pow(1/error, 1.0/7);
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
        ParticSys localSystem_;
        
        RadauTab BTab_;
        
        RadauTab GTab_;
        
        /** @brief Local relative error*/
        Scalar relativeError_{1e-16};
        
        Scalar convergent_limit_{1e-9};
        
        size_t particleNum_{0};
        
        constexpr static size_t max_iter_ = 12；
    };

}
#endif
